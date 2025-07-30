# pathway_parameters.R

pathwayParametersUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    box(
      title = "Analysis Parameters", width = 12, status = "primary",
      selectInput(ns("pathway_omics"), "Omics Layer",
                  choices = c("Transcriptomics", "Proteomics", "Metabolomics"),
                  selected = "Transcriptomics"),
      
      uiOutput(ns("id_type_ui")),
      
      selectInput(ns("pathway_strata"), "Stratify By", 
                  choices = c("None")),
      selectInput(ns("pathway_species"), "Species",
                  choices = c("Human" = "hsa", "Mouse" = "mmu", "Rat" = "rno"),
                  selected = "hsa"),
      numericInput(ns("pathway_pval"), "P-value Threshold", 
                   value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput(ns("pathway_qval"), "Q-value Threshold", 
                   value = 0.1, min = 0, max = 1, step = 0.01),
      numericInput(ns("pathway_fc"), "Minimum Fold Change", 
                   value = 1.5, min = 0, step = 0.1),
      numericInput(ns("top_pathways"), "Number of Top Pathways to Show",
                   value = 10, min = 1, max = 50, step = 1),
      actionButton(ns("run_pathway"), "Run Pathway Analysis",
                   class = "btn-primary")
    )
  )
}

pathwayParametersServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update stratification dropdown based on metadata
    observe({
      req(values$metadata)
      meta_vars <- colnames(values$metadata)[-1] # exclude SampleID
      updateSelectInput(session, "pathway_strata", choices = c("None", meta_vars))
    })
    
    # Dynamic UI for ID type selection
    output$id_type_ui <- renderUI({
      if (input$pathway_omics == "Transcriptomics") {
        selectInput(ns("pathway_idtype"), "Gene ID Type",
                    choices = c("Gene Symbol" = "SYMBOL", 
                                "Ensembl Gene" = "ENSEMBL",
                                "Entrez ID" = "ENTREZID"),
                    selected = "SYMBOL")
      } else if (input$pathway_omics == "Proteomics") {
        selectInput(ns("pathway_idtype"), "Protein ID Type",
                    choices = c("UniProt" = "UNIPROT",
                                "Gene Symbol" = "SYMBOL",
                                "Entrez ID" = "ENTREZID"),
                    selected = "UNIPROT")
      } else if (input$pathway_omics == "Metabolomics") {
        selectInput(ns("pathway_idtype"), "Metabolite ID Type",
                    choices = c("HMDB" = "HMDB",
                                "KEGG" = "KEGG",
                                "ChEBI" = "CHEBI"),
                    selected = "HMDB")
      }
    })
    
    # Helper function to extract accession numbers
    extract_accession <- function(ids) {
      sapply(ids, function(x) {
        if (grepl("\\(", x)) {
          sub("\\(.*", "", x)
        } else {
          x
        }
      })
    }
    
    # Enhanced compound ID conversion function
    convert_compound_ids <- function(cpd_list, id_type, species) {
      if (is.null(cpd_list)) return(NULL)
      
      # If IDs are already in KEGG format (C#####)
      if (all(grepl("^C\\d+$", names(cpd_list)))) {
        return(cpd_list)
      }
      
      # Convert HMDB to KEGG IDs
      if (id_type == "HMDB") {
        tryCatch({
          # Get the first part of HMDB IDs (HMDB00001 -> 00001)
          clean_ids <- sub("^HMDB", "", names(cpd_list))
          
          # Convert using pathview
          kegg_ids <- pathview::cpdidmap(
            in.ids = clean_ids,
            in.type = "hmdb",
            out.type = "kegg"
          )
          
          # Keep only successfully mapped compounds
          mapped_idx <- match(clean_ids, kegg_ids[,1])
          mapped_cpds <- cpd_list[!is.na(mapped_idx)]
          names(mapped_cpds) <- kegg_ids[na.omit(mapped_idx), 2]
          
          # Filter for valid KEGG IDs
          valid_cpds <- mapped_cpds[grepl("^C\\d+$", names(mapped_cpds))]
          
          if (length(valid_cpds) == 0) {
            showNotification("No valid KEGG compound IDs found after conversion", 
                             type = "warning")
            return(NULL)
          }
          
          message(paste(
            "Converted", length(valid_cpds), "of", length(cpd_list), "compounds to KEGG IDs",
            "\nExample:", names(cpd_list)[1], "->", names(valid_cpds)[1]
          ))
          
          return(valid_cpds)
        }, error = function(e) {
          showNotification(paste("Compound ID conversion failed:", e$message), 
                           type = "error")
          return(NULL)
        })
      }
      return(cpd_list) # Return unchanged if not HMDB
    }
    
    # Helper function to prepare omics data
    prepare_omics_data <- function(omics_data, strata_var, groups) {
      group1_samples <- values$metadata$SampleID[strata_var == groups[1]]
      group2_samples <- values$metadata$SampleID[strata_var == groups[2]]
      
      # Calculate fold changes with pseudocount
      group1_means <- rowMeans(omics_data[, colnames(omics_data) %in% group1_samples, drop = FALSE], na.rm = TRUE)
      group2_means <- rowMeans(omics_data[, colnames(omics_data) %in% group2_samples, drop = FALSE], na.rm = TRUE)
      fc <- log2((group2_means + 1)/(group1_means + 1))
      names(fc) <- extract_accession(omics_data[,1])
      
      # Calculate p-values
      pvals <- apply(omics_data[,-1, drop = FALSE], 1, function(x) {
        tryCatch(
          t.test(x[colnames(omics_data)[-1] %in% group1_samples],
                 x[colnames(omics_data)[-1] %in% group2_samples])$p.value,
          error = function(e) NA
        )
      })
      
      # Filter significant results
      qvals <- p.adjust(pvals, method = "BH")
      sig_idx <- which(pvals < input$pathway_pval & 
                         qvals < input$pathway_qval & 
                         abs(fc) > log2(input$pathway_fc))
      fc[sig_idx]
    }
    
    # Helper function to convert gene IDs to ENTREZID
    convert_gene_ids <- function(gene_list, id_type) {
      if (id_type == "ENTREZID") return(gene_list)
      
      gene_ids <- names(gene_list)
      id_map <- tryCatch({
        clusterProfiler::bitr(gene_ids,
                              fromType = id_type,
                              toType = "ENTREZID",
                              OrgDb = "org.Hs.eg.db")
      }, error = function(e) {
        showNotification(paste("ID conversion error:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(id_map)) return(NULL)
      
      # Keep only genes that were successfully mapped
      mapped_idx <- match(gene_ids, id_map[[id_type]])
      mapped_genes <- gene_list[!is.na(mapped_idx)]
      names(mapped_genes) <- id_map$ENTREZID[mapped_idx[!is.na(mapped_idx)]]
      
      # Remove duplicates by averaging values
      if (any(duplicated(names(mapped_genes)))) {
        mapped_genes <- tapply(mapped_genes, names(mapped_genes), mean)
      }
      
      mapped_genes
    }
    
    # Run pathway analysis
    pathway_results <- eventReactive(input$run_pathway, {
      req(!is.null(values$transcript_data) || !is.null(values$protein_data) || !is.null(values$metabo_data))
      
      tryCatch({
        showNotification("Running pathway analysis...", type = "message", duration = NULL)
        
        # Initialize data lists
        gene_list <- NULL
        cpd_list <- NULL
        
        # Handle stratified vs non-stratified analysis
        if (input$pathway_strata != "None") {
          strata_var <- values$metadata[[input$pathway_strata]]
          groups <- unique(na.omit(strata_var))
          
          if (length(groups) != 2) {
            showNotification("Need exactly 2 groups for comparison", type = "error")
            return(NULL)
          }
          
          # Process selected omics layer
          if (input$pathway_omics == "Transcriptomics" && !is.null(values$transcript_data)) {
            gene_list <- prepare_omics_data(values$transcript_data, strata_var, groups)
            if (!is.null(gene_list)) {
              gene_list <- convert_gene_ids(gene_list, input$pathway_idtype)
              if (is.null(gene_list)) return(NULL)
            }
            
            # Add metabolites if available
            if (!is.null(values$metabo_data)) {
              cpd_list <- prepare_omics_data(values$metabo_data, strata_var, groups)
              if (!is.null(cpd_list)) {
                cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
              }
            }
          } 
          else if (input$pathway_omics == "Proteomics" && !is.null(values$protein_data)) {
            gene_list <- prepare_omics_data(values$protein_data, strata_var, groups)
            if (!is.null(gene_list)) {
              gene_list <- convert_gene_ids(gene_list, input$pathway_idtype)
              if (is.null(gene_list)) return(NULL)
            }
            
            # Add metabolites if available
            if (!is.null(values$metabo_data)) {
              cpd_list <- prepare_omics_data(values$metabo_data, strata_var, groups)
              if (!is.null(cpd_list)) {
                cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
              }
            }
          }
          else if (input$pathway_omics == "Metabolomics" && !is.null(values$metabo_data)) {
            cpd_list <- prepare_omics_data(values$metabo_data, strata_var, groups)
            if (!is.null(cpd_list)) {
              cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
            }
          }
        } else {
          # Non-stratified analysis - just use means
          if (input$pathway_omics == "Transcriptomics" && !is.null(values$transcript_data)) {
            gene_list <- rowMeans(values$transcript_data[,-1, drop = FALSE], na.rm = TRUE)
            names(gene_list) <- extract_accession(values$transcript_data[,1])
            gene_list <- convert_gene_ids(gene_list, input$pathway_idtype)
            if (is.null(gene_list)) return(NULL)
            
            # Add metabolites if available
            if (!is.null(values$metabo_data)) {
              cpd_list <- rowMeans(values$metabo_data[,-1, drop = FALSE], na.rm = TRUE)
              names(cpd_list) <- extract_accession(values$metabo_data[,1])
              cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
            }
          } 
          else if (input$pathway_omics == "Proteomics" && !is.null(values$protein_data)) {
            gene_list <- rowMeans(values$protein_data[,-1, drop = FALSE], na.rm = TRUE)
            names(gene_list) <- extract_accession(values$protein_data[,1])
            gene_list <- convert_gene_ids(gene_list, input$pathway_idtype)
            if (is.null(gene_list)) return(NULL)
            
            # Add metabolites if available
            if (!is.null(values$metabo_data)) {
              cpd_list <- rowMeans(values$metabo_data[,-1, drop = FALSE], na.rm = TRUE)
              names(cpd_list) <- extract_accession(values$metabo_data[,1])
              cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
            }
          }
          else if (input$pathway_omics == "Metabolomics" && !is.null(values$metabo_data)) {
            cpd_list <- rowMeans(values$metabo_data[,-1, drop = FALSE], na.rm = TRUE)
            names(cpd_list) <- extract_accession(values$metabo_data[,1])
            cpd_list <- convert_compound_ids(cpd_list, input$pathway_idtype, input$pathway_species)
          }
        }
        
        # Run KEGG enrichment if we have gene data
        kegg_res <- NULL
        if (!is.null(gene_list) && length(gene_list) > 0) {
          kegg_res <- clusterProfiler::enrichKEGG(
            gene = names(gene_list),
            organism = input$pathway_species,
            pvalueCutoff = input$pathway_pval,
            qvalueCutoff = input$pathway_qval
          )
          
          # Sort by p-value and select top pathways
          if (!is.null(kegg_res)) {
            kegg_res <- kegg_res[order(kegg_res$p.adjust),]
            kegg_res <- kegg_res[1:min(input$top_pathways, nrow(kegg_res)),]
          }
        }
        
        showNotification("Pathway analysis completed", type = "message")
        
        # Return results with converted compound IDs
        list(
          kegg_res = kegg_res,
          gene_list = gene_list,
          cpd_list = cpd_list,  # This now contains KEGG-formatted IDs
          species = input$pathway_species,
          omics = input$pathway_omics,
          metadata = values$metadata
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        return(NULL)
      })
    })
    
    return(pathway_results)
  })
}