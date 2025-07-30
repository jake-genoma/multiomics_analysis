# pathway_selection.R

pathwaySelectionUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    box(
      title = "Pathway Selection", width = 12, status = "info",
      DTOutput(ns("pathway_table")),
      fluidRow(
        column(6, uiOutput(ns("download_selected_ui"))),
        column(6, uiOutput(ns("download_all_ui")))
      )
    )
  )
}

pathwaySelectionServer <- function(id, pathway_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Track generated files for cleanup
    generated_files <- reactiveVal(character())
    
    # Cleanup function for pathview files
    cleanup_pathview_files <- function() {
      if (length(generated_files()) > 0) {
        unlink(generated_files())
        generated_files(character())
      }
      
      # Clean any remaining pathview files
      remaining_files <- list.files(
        pattern = paste0("^", pathway_results()$species, "\\d+.*\\.(png|xml)$"),
        full.names = TRUE
      )
      if (length(remaining_files) > 0) {
        unlink(remaining_files)
      }
    }
    
    # Clean up when module is destroyed
    onStop(function() {
      cleanup_pathview_files()
    })
    
    # Display pathway table
    output$pathway_table <- renderDT({
      req(pathway_results())
      
      if (is.null(pathway_results()$kegg_res)) {
        return(datatable(data.frame(Note = "No pathway enrichment results. Try different thresholds.")))
      }
      
      pathway_df <- as.data.frame(pathway_results()$kegg_res)
      pathway_df <- pathway_df[, c("ID", "Description", "pvalue", "p.adjust", "geneID")]
      colnames(pathway_df) <- c("Pathway ID", "Pathway", "P-value", "Adj. P-value", "Genes")
      
      pathway_df$Genes <- vapply(pathway_df$Genes, function(x) {
        genes <- strsplit(x, "/")[[1]]
        if (length(genes) > 5) paste0(paste(genes[1:5], collapse = ", "), ", +", length(genes)-5, " more")
        else paste(genes, collapse = ", ")
      }, character(1))
      
      datatable(pathway_df, 
                options = list(
                  pageLength = 10,
                  scrollX = TRUE,
                  dom = 'tip'
                ),
                rownames = FALSE,
                selection = 'single') %>%
        formatSignif(columns = c("P-value", "Adj. P-value"), digits = 3)
    })
    
    # Reactive for selected pathway index
    selected_pathway_index <- reactive({
      req(input$pathway_table_rows_selected)
      input$pathway_table_rows_selected
    })
    
    # Download button UIs
    output$download_selected_ui <- renderUI({
      req(pathway_results(), selected_pathway_index())
      downloadButton(ns("download_selected"), "Download Selected", 
                     class = "btn-primary", style = "width: 100%;")
    })
    
    output$download_all_ui <- renderUI({
      req(pathway_results())
      downloadButton(ns("download_all"), "Download All", 
                     class = "btn-success", style = "width: 100%;")
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
    
    # Safe pathview wrapper with proper compound handling
    safe_pathview <- function(pathway_id, gene_list, cpd_list, species, suffix) {
      # Clean up previous files first
      cleanup_pathview_files()
      
      tryCatch({
        pathway_id_clean <- sub(paste0("^", species), "", pathway_id)
        
        # Verify compound IDs are in KEGG format
        if (!is.null(cpd_list)) {
          if (any(grepl("^HMDB", names(cpd_list)))) {
            showNotification("Compound IDs not properly converted to KEGG format", 
                             type = "error", duration = 10)
            return(NULL)
          }
          
          # Filter out any invalid entries
          valid_cpds <- cpd_list[grepl("^C\\d+$", names(cpd_list))]
          if (length(valid_cpds) == 0) {
            message("No valid KEGG compound IDs found")
            cpd_list <- NULL
          } else {
            cpd_list <- valid_cpds
          }
        }
        
        # Debug output
        message("\n=== Pathway Visualization Parameters ===")
        message("Pathway: ", pathway_id)
        message("Species: ", species)
        message("Genes: ", length(gene_list))
        if (!is.null(cpd_list)) {
          message("Compounds: ", length(cpd_list))
          message("Example compound: ", names(cpd_list)[1], " = ", cpd_list[1])
        }
        
        # Generate visualization
        suppressMessages(
          pathview::pathview(
            gene.data = gene_list,
            cpd.data = cpd_list,
            pathway.id = pathway_id_clean,
            species = species,
            kegg.native = TRUE,
            gene.idtype = "entrez",
            cpd.idtype = "kegg",
            out.suffix = suffix,
            low = list(gene = "green", cpd = "blue"),
            mid = list(gene = "gray", cpd = "gray"),
            high = list(gene = "red", cpd = "yellow"),
            na.col = "white"
          )
        )
        
        # Find and track generated files
        img_pattern <- paste0(species, pathway_id_clean, ".*", suffix, "\\.png$")
        generated_img <- list.files(pattern = img_pattern, full.names = TRUE)
        
        if (length(generated_img) > 0) {
          generated_files(c(generated_files(), generated_img))
          message("Successfully generated visualization: ", generated_img[1])
          return(generated_img[1])
        } else {
          message("No image file was generated")
          return(NULL)
        }
      }, error = function(e) {
        message("Pathview error: ", e$message)
        showNotification(paste("Visualization failed:", e$message), 
                         type = "error", duration = 10)
        return(NULL)
      })
    }
    
    # Download handler for selected pathway
    output$download_selected <- downloadHandler(
      filename = function() {
        selected_pathway <- pathway_results()$kegg_res$ID[selected_pathway_index()]
        paste0("pathway_", selected_pathway, ".zip")
      },
      content = function(file) {
        req(pathway_results(), selected_pathway_index())
        
        # Create temp directory for packaging
        temp_dir <- tempfile()
        dir.create(temp_dir)
        on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
        
        selected_pathway <- pathway_results()$kegg_res$ID[selected_pathway_index()]
        pathway_name <- pathway_results()$kegg_res$Description[selected_pathway_index()]
        
        # Generate visualization
        img_path <- safe_pathview(
          pathway_id = selected_pathway,
          gene_list = pathway_results()$gene_list,
          cpd_list = pathway_results()$cpd_list,
          species = pathway_results()$species,
          suffix = "selected"
        )
        
        files_to_zip <- c()
        
        if (!is.null(img_path)) {
          new_img_path <- file.path(temp_dir, paste0("pathway_", selected_pathway, ".png"))
          if (file.copy(img_path, new_img_path)) {
            files_to_zip <- c(files_to_zip, new_img_path)
          }
        }
        
        # Create CSV with pathway data
        pathway_df <- as.data.frame(pathway_results()$kegg_res[selected_pathway_index(), ])
        
        # Extract mapped gene data
        gene_ids_in_pathway <- unlist(strsplit(pathway_df$geneID, "/"))
        gene_data <- data.frame(
          GeneID = gene_ids_in_pathway,
          Log2FC = NA,
          Mapped = "No"
        )
        
        if (!is.null(pathway_results()$gene_list)) {
          matched_genes <- gene_ids_in_pathway[gene_ids_in_pathway %in% names(pathway_results()$gene_list)]
          if (length(matched_genes) > 0) {
            gene_data$Log2FC[match(matched_genes, gene_data$GeneID)] <- 
              pathway_results()$gene_list[matched_genes]
            gene_data$Mapped[match(matched_genes, gene_data$GeneID)] <- "Yes"
          }
        }
        
        # Extract mapped compound data if available
        cpd_data <- NULL
        if (!is.null(pathway_results()$cpd_list)) {
          cpd_data <- data.frame(
            OriginalID = if (any(grepl("^HMDB", names(pathway_results()$cpd_list)))) {
              names(pathway_results()$cpd_list)
            } else {
              NA
            },
            KEGG_ID = names(pathway_results()$cpd_list),
            Log2FC = pathway_results()$cpd_list,
            Mapped = ifelse(grepl("^C\\d+$", names(pathway_results()$cpd_list)), "Yes", "No")
          )
        }
        
        # Save data files
        csv_file1 <- file.path(temp_dir, paste0("pathway_", selected_pathway, "_summary.csv"))
        write.csv(pathway_df, csv_file1, row.names = FALSE)
        files_to_zip <- c(files_to_zip, csv_file1)
        
        csv_file2 <- file.path(temp_dir, paste0("pathway_", selected_pathway, "_gene_data.csv"))
        write.csv(gene_data, csv_file2, row.names = FALSE)
        files_to_zip <- c(files_to_zip, csv_file2)
        
        if (!is.null(cpd_data)) {
          csv_file3 <- file.path(temp_dir, paste0("pathway_", selected_pathway, "_compound_data.csv"))
          write.csv(cpd_data, csv_file3, row.names = FALSE)
          files_to_zip <- c(files_to_zip, csv_file3)
        }
        
        # Create README file
        readme_content <- paste(
          "Pathway Analysis Results",
          "=======================",
          paste("Pathway ID:", selected_pathway),
          paste("Pathway Name:", pathway_name),
          paste("Omics Layer:", pathway_results()$omics),
          paste("Species:", pathway_results()$species),
          "",
          "Files included:",
          ifelse(!is.null(img_path), "1. Pathway visualization (PNG)", "1. [No visualization generated]"),
          "2. Pathway summary (CSV)",
          "3. Mapped gene data (CSV)",
          ifelse(!is.null(cpd_data), "4. Mapped compound data (CSV)", ""),
          "",
          "Generated on:", Sys.time(),
          sep = "\n"
        )
        
        readme_file <- file.path(temp_dir, "README.txt")
        writeLines(readme_content, readme_file)
        files_to_zip <- c(files_to_zip, readme_file)
        
        # Create ZIP file
        if (length(files_to_zip) > 0) {
          old_wd <- setwd(temp_dir)
          on.exit(setwd(old_wd))
          zip::zip(file, files = basename(files_to_zip))
        } else {
          showNotification("No files were generated for download", type = "warning")
        }
      }
    )
    
    # Download handler for all pathways
    output$download_all <- downloadHandler(
      filename = function() {
        "all_pathways.zip"
      },
      content = function(file) {
        req(pathway_results())
        
        # Create temp directory
        temp_dir <- tempfile()
        dir.create(temp_dir)
        on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
        
        # Create subdirectories
        pathways_dir <- file.path(temp_dir, "pathways")
        data_dir <- file.path(temp_dir, "data")
        dir.create(pathways_dir)
        dir.create(data_dir)
        
        # Create progress bar
        progress <- shiny::Progress$new()
        progress$set(message = "Generating pathway package", value = 0)
        on.exit(progress$close(), add = TRUE)
        
        pathway_ids <- pathway_results()$kegg_res$ID
        n_pathways <- min(length(pathway_ids), nrow(pathway_results()$kegg_res))
        
        # Track successfully processed pathways
        successful_pathways <- 0
        
        # Process each pathway
        for (i in 1:n_pathways) {
          progress$inc(1/n_pathways, detail = paste("Processing pathway", i, "of", n_pathways))
          
          pathway_id <- pathway_ids[i]
          pathway_name <- pathway_results()$kegg_res$Description[i]
          
          # Generate visualization
          img_path <- safe_pathview(
            pathway_id = pathway_id,
            gene_list = pathway_results()$gene_list,
            cpd_list = pathway_results()$cpd_list,
            species = pathway_results()$species,
            suffix = paste0("all_", i)
          )
          
          # Only proceed if visualization was successful
          if (!is.null(img_path) && file.exists(img_path)) {
            new_img_path <- file.path(pathways_dir, paste0("pathway_", pathway_id, ".png"))
            
            if (file.copy(img_path, new_img_path)) {
              successful_pathways <- successful_pathways + 1
            }
          }
          
          # Create data files regardless of visualization success
          pathway_df <- as.data.frame(pathway_results()$kegg_res[i, ])
          
          # Gene data
          gene_ids_in_pathway <- unlist(strsplit(pathway_df$geneID, "/"))
          gene_data <- data.frame(
            GeneID = gene_ids_in_pathway,
            Log2FC = NA,
            Mapped = "No"
          )
          
          if (!is.null(pathway_results()$gene_list)) {
            matched_genes <- gene_ids_in_pathway[gene_ids_in_pathway %in% names(pathway_results()$gene_list)]
            if (length(matched_genes) > 0) {
              gene_data$Log2FC[match(matched_genes, gene_data$GeneID)] <- 
                pathway_results()$gene_list[matched_genes]
              gene_data$Mapped[match(matched_genes, gene_data$GeneID)] <- "Yes"
            }
          }
          
          # Compound data if available
          cpd_data <- NULL
          if (!is.null(pathway_results()$cpd_list)) {
            cpd_data <- data.frame(
              OriginalID = if (any(grepl("^HMDB", names(pathway_results()$cpd_list)))) {
                names(pathway_results()$cpd_list)
              } else {
                NA
              },
              KEGG_ID = names(pathway_results()$cpd_list),
              Log2FC = pathway_results()$cpd_list,
              Mapped = ifelse(grepl("^C\\d+$", names(pathway_results()$cpd_list)), "Yes", "No")
            )
          }
          
          # Save data files
          data_subdir <- file.path(data_dir, pathway_id)
          dir.create(data_subdir, recursive = TRUE, showWarnings = FALSE)
          
          # Only write files if directory was created successfully
          if (dir.exists(data_subdir)) {
            write.csv(pathway_df, file.path(data_subdir, "pathway_summary.csv"), row.names = FALSE)
            write.csv(gene_data, file.path(data_subdir, "gene_data.csv"), row.names = FALSE)
            if (!is.null(cpd_data)) {
              write.csv(cpd_data, file.path(data_subdir, "compound_data.csv"), row.names = FALSE)
            }
          }
        }
        
        # Create master CSV with all pathway results
        if (n_pathways > 0) {
          all_pathways_df <- as.data.frame(pathway_results()$kegg_res)
          write.csv(all_pathways_df, file.path(temp_dir, "all_pathways_summary.csv"), row.names = FALSE)
        }
        
        # Create README file
        readme_content <- paste(
          "Pathway Analysis Results - Complete Set",
          "=====================================",
          paste("Number of Pathways Processed:", successful_pathways, "of", n_pathways),
          paste("Omics Layer:", pathway_results()$omics),
          paste("Species:", pathway_results()$species),
          "",
          "Folder Structure:",
          "1. pathways/ - Contains pathway visualizations (PNG)",
          "2. data/ - Contains subfolders for each pathway with data files",
          "   - Each pathway has its own folder with:",
          "     - pathway_summary.csv",
          "     - gene_data.csv",
          "     - compound_data.csv (if available)",
          "3. all_pathways_summary.csv - Summary of all pathways",
          "",
          "Generated on:", Sys.time(),
          sep = "\n"
        )
        
        writeLines(readme_content, file.path(temp_dir, "README.txt"))
        
        # Create ZIP file
        tryCatch({
          # Change to temp directory to create proper relative paths
          old_wd <- setwd(temp_dir)
          on.exit(setwd(old_wd), add = TRUE)
          
          # Get all files recursively
          all_files <- list.files(recursive = TRUE, full.names = FALSE)
          
          if (length(all_files) > 0) {
            zip::zip(file, files = all_files)
          } else {
            showNotification("No valid files were generated for download", type = "warning")
          }
        }, error = function(e) {
          showNotification(paste("Error creating ZIP file:", e$message), type = "error")
        })
      }
    )
    
    # Return the selected pathway index
    return(selected_pathway_index)
  })
}