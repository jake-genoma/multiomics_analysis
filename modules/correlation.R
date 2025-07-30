# modules/correlation.R

correlationUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Stratified Correlation Analysis"),
    fluidRow(
      box(
        title = "Analysis Parameters", width = 4, status = "primary",
        selectInput(ns("corr_method"), "Correlation Method",
                    choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                    selected = "pearson"),
        selectInput(ns("corr_strata"), "Stratify By (Optional)", 
                    choices = c("None")),
        uiOutput(ns("omics1_ui")),
        uiOutput(ns("omics2_ui")),
        numericInput(ns("corr_threshold"), "Significance Threshold (FDR)", 
                     value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("min_abs_corr"), "Minimum Absolute Correlation", 
                     value = 0.5, min = 0, max = 1, step = 0.1),
        checkboxInput(ns("filter_significant"), "Filter Significant Correlations Only", value = TRUE),
        actionButton(ns("run_correlation"), "Run Correlation Analysis", class = "btn-primary"),
        tags$hr(),
        downloadButton(ns("download_corr"), "Download Results")
      ),
      
      box(
        title = "Results", width = 8, status = "info",
        tabsetPanel(
          tabPanel("Correlation Matrix", 
                   plotlyOutput(ns("corr_heatmap")) %>% withSpinner(),
                   uiOutput(ns("corr_strata_header"))),
          tabPanel("Top Correlations", 
                   DTOutput(ns("top_corr_table")) %>% withSpinner()),
          tabPanel("Stratified Scatter", 
                   fluidRow(
                     column(6, uiOutput(ns("scatter_x_ui"))),
                     column(6, uiOutput(ns("scatter_y_ui")))
                   ),
                   plotlyOutput(ns("stratified_scatter_plot")) %>% withSpinner(),
                   verbatimTextOutput(ns("scatter_stats"))),
          tabPanel("Biological Interpretation",
                   h4("Pathway Analysis of Correlated Features"),
                   selectInput(ns("pathway_db"), "Pathway Database",
                               choices = c("KEGG", "Reactome", "GO Biological Process"),
                               selected = "KEGG"),
                   actionButton(ns("run_pathway"), "Analyze Pathways"),
                   plotlyOutput(ns("pathway_plot")) %>% withSpinner(),
                   DTOutput(ns("pathway_table")) %>% withSpinner())
        )
      )
    )
  )
}

correlationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive function to get available omics data types
    available_omics <- reactive({
      omics_types <- c()
      if (!is.null(values$transcript_data)) omics_types <- c(omics_types, "Transcriptomics")
      if (!is.null(values$protein_data)) omics_types <- c(omics_types, "Proteomics")
      if (!is.null(values$metabo_data)) omics_types <- c(omics_types, "Metabolomics")
      
      if (length(omics_types) < 2) {
        showNotification("At least two omics datasets are required for correlation analysis", 
                         type = "warning")
      }
      omics_types
    })
    
    # Dynamic UI for omics selection
    output$omics1_ui <- renderUI({
      omics <- available_omics()
      if (length(omics) >= 2) {
        selectInput(ns("corr_omics1"), "First Omics Layer",
                    choices = omics,
                    selected = omics[1])
      } else {
        disabled(selectInput(ns("corr_omics1"), "First Omics Layer (need 2+ datasets)",
                             choices = omics))
      }
    })
    
    output$omics2_ui <- renderUI({
      omics <- available_omics()
      if (length(omics) >= 2) {
        selectInput(ns("corr_omics2"), "Second Omics Layer",
                    choices = omics,
                    selected = omics[2])
      } else {
        disabled(selectInput(ns("corr_omics2"), "Second Omics Layer (need 2+ datasets)",
                             choices = omics))
      }
    })
    
    # Update stratification dropdown when metadata changes
    observe({
      req(values$metadata)
      updateSelectInput(session, "corr_strata", 
                        choices = c("None", setdiff(colnames(values$metadata), "SampleID")))
    })
    
    # Safe data extraction function
    safe_extract_matrix <- function(omics_data) {
      tryCatch({
        mat <- as.matrix(omics_data[,-1])
        rownames(mat) <- omics_data[,1]
        t(mat)
      }, error = function(e) {
        showNotification(paste("Error processing omics data:", e$message), type = "error")
        NULL
      })
    }
    
    # Reactive function to calculate correlation with p-values
    calculate_correlation <- function(mat1, mat2, method) {
      tryCatch({
        # Ensure matrices have the same samples in the same order
        common_samples <- intersect(rownames(mat1), rownames(mat2))
        
        if (length(common_samples) < 3) {
          stop("Insufficient common samples (need at least 3) for correlation analysis")
        }
        
        mat1 <- mat1[common_samples, , drop = FALSE]
        mat2 <- mat2[common_samples, , drop = FALSE]
        
        # Calculate correlation matrix
        cor_mat <- cor(mat1, mat2, method = method, use = "pairwise.complete.obs")
        
        # Calculate p-values
        n <- nrow(mat1)
        cor_p <- matrix(NA, nrow = ncol(mat1), ncol = ncol(mat2))
        rownames(cor_p) <- colnames(mat1)
        colnames(cor_p) <- colnames(mat2)
        
        for (i in 1:ncol(mat1)) {
          for (j in 1:ncol(mat2)) {
            test <- try(cor.test(mat1[,i], mat2[,j], method = method), silent = TRUE)
            if (!inherits(test, "try-error")) {
              cor_p[i,j] <- test$p.value
            }
          }
        }
        
        # Adjust p-values for multiple testing
        cor_fdr <- matrix(p.adjust(cor_p, method = "fdr"), 
                          nrow = nrow(cor_p), ncol = ncol(cor_p))
        rownames(cor_fdr) <- rownames(cor_p)
        colnames(cor_fdr) <- colnames(cor_p)
        
        return(list(cor = cor_mat, p = cor_p, fdr = cor_fdr))
      }, error = function(e) {
        showNotification(paste("Error calculating correlations:", e$message), type = "error")
        NULL
      })
    }
    
    # Main correlation analysis
    observeEvent(input$run_correlation, {
      req(input$corr_omics1, input$corr_omics2)
      
      # Validate that we have different omics types selected
      if (input$corr_omics1 == input$corr_omics2) {
        showNotification("Please select two different omics layers for comparison", 
                         type = "warning")
        return()
      }
      
      tryCatch({
        showNotification("Running correlation analysis...", type = "message")
        
        # Get the two omics datasets to correlate
        omics1_data <- switch(input$corr_omics1,
                              "Transcriptomics" = values$transcript_data,
                              "Proteomics" = values$protein_data,
                              "Metabolomics" = values$metabo_data)
        
        omics2_data <- switch(input$corr_omics2,
                              "Transcriptomics" = values$transcript_data,
                              "Proteomics" = values$protein_data,
                              "Metabolomics" = values$metabo_data)
        
        if (is.null(omics1_data) || is.null(omics2_data)) {
          stop("Selected omics data not available")
        }
        
        # Prepare data (samples x features)
        mat1 <- safe_extract_matrix(omics1_data)
        mat2 <- safe_extract_matrix(omics2_data)
        
        if (is.null(mat1) || is.null(mat2)) return()
        
        # Calculate correlation with p-values
        cor_method <- input$corr_method
        cor_res <- calculate_correlation(mat1, mat2, cor_method)
        if (is.null(cor_res)) return()
        
        # Filter results based on thresholds
        sig_mask <- if (input$filter_significant) {
          cor_res$fdr < input$corr_threshold & abs(cor_res$cor) >= input$min_abs_corr
        } else {
          abs(cor_res$cor) >= input$min_abs_corr
        }
        
        filtered_cor <- cor_res$cor * sig_mask
        filtered_cor[filtered_cor == 0] <- NA
        
        # Prepare top correlations table
        cor_df <- tryCatch({
          melted <- reshape2::melt(filtered_cor, na.rm = TRUE)
          colnames(melted) <- c("Feature1", "Feature2", "Correlation")
          
          # Add p-values and FDR
          melted$PValue <- cor_res$p[as.matrix(melted[,1:2])]
          melted$FDR <- cor_res$fdr[as.matrix(melted[,1:2])]
          
          # Add feature annotations
          melted$Feature1_Type <- input$corr_omics1
          melted$Feature2_Type <- input$corr_omics2
          
          # Order by absolute correlation
          melted[order(abs(melted$Correlation), decreasing = TRUE),]
        }, error = function(e) {
          showNotification(paste("Error processing correlation results:", e$message), type = "error")
          NULL
        })
        
        if (is.null(cor_df)) return()
        
        # Store results
        values$corr_results <- list(
          matrix = filtered_cor,
          full_matrix = cor_res$cor,
          method = cor_method,
          omics1 = input$corr_omics1,
          omics2 = input$corr_omics2,
          all_correlations = cor_df,
          top_correlations = head(cor_df, 100)
        )
        
        showNotification("Correlation analysis completed", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in correlation analysis:", e$message), type = "error")
      })
    })
    
    # Dynamic UI for scatter plot feature selection
    output$scatter_x_ui <- renderUI({
      req(values$corr_results)
      selectInput(ns("scatter_x"), "X-axis Feature", 
                  choices = rownames(values$corr_results$matrix))
    })
    
    output$scatter_y_ui <- renderUI({
      req(values$corr_results)
      selectInput(ns("scatter_y"), "Y-axis Feature", 
                  choices = colnames(values$corr_results$matrix))
    })
    
    # Heatmap visualization
    output$corr_heatmap <- renderPlotly({
      req(values$corr_results)
      
      tryCatch({
        cor_mat <- values$corr_results$matrix
        
        if (nrow(cor_mat) == 0 || ncol(cor_mat) == 0) {
          stop("No significant correlations found with current thresholds")
        }
        
        # Cluster the matrix for better visualization
        row_order <- try(hclust(dist(cor_mat, method = "euclidean"), method = "complete")$order)
        col_order <- try(hclust(dist(t(cor_mat), method = "euclidean"), method = "complete")$order)
        
        if (!inherits(row_order, "try-error")) {
          cor_mat <- cor_mat[row_order, , drop = FALSE]
        }
        if (!inherits(col_order, "try-error")) {
          cor_mat <- cor_mat[, col_order, drop = FALSE]
        }
        
        # Create custom hover text
        hover_text <- matrix(paste0(
          "Feature 1: ", rownames(cor_mat), "<br>",
          "Feature 2: ", rep(colnames(cor_mat), each = nrow(cor_mat)), "<br>",
          "Correlation: ", round(cor_mat, 3)),
          nrow = nrow(cor_mat))
        
        plot_ly(
          x = colnames(cor_mat),
          y = rownames(cor_mat),
          z = cor_mat,
          type = "heatmap",
          colors = colorRamp(c("blue", "white", "red")),
          colorbar = list(title = paste(toupper(values$corr_results$method), "Correlation")),
          hoverinfo = "text",
          text = hover_text
        ) %>% layout(
          title = paste(values$corr_results$omics1, "vs", values$corr_results$omics2, 
                        "Correlation (|r| >=", input$min_abs_corr, ")"),
          xaxis = list(title = values$corr_results$omics2, tickangle = 45),
          yaxis = list(title = values$corr_results$omics1),
          margin = list(l = 150, b = 150)
        )
      }, error = function(e) {
        showNotification(paste("Error generating heatmap:", e$message), type = "error")
        NULL
      })
    })
    
    # Top correlations table
    output$top_corr_table <- renderDT({
      req(values$corr_results)
      
      tryCatch({
        datatable(values$corr_results$top_correlations, 
                  options = list(
                    scrollX = TRUE,
                    pageLength = 10,
                    dom = 'Bfrtip',
                    buttons = c('csv', 'excel')
                  ),
                  extensions = 'Buttons',
                  rownames = FALSE,
                  filter = 'top') %>%
          formatRound(columns = c("Correlation", "PValue", "FDR"), digits = 4)
      }, error = function(e) {
        showNotification(paste("Error rendering correlation table:", e$message), type = "error")
        NULL
      })
    })
    
    # Stratified scatter plot
    output$stratified_scatter_plot <- renderPlotly({
      req(values$corr_results, input$scatter_x, input$scatter_y)
      
      tryCatch({
        # Get expression values for selected features
        x_feature <- input$scatter_x
        y_feature <- input$scatter_y
        
        # Get the appropriate data for each feature
        get_feature_values <- function(feature, omics_type) {
          data <- switch(omics_type,
                         "Transcriptomics" = values$transcript_data,
                         "Proteomics" = values$protein_data,
                         "Metabolomics" = values$metabo_data)
          
          if (is.null(data)) {
            stop(paste(omics_type, "data not available"))
          }
          
          if (!feature %in% data[,1]) {
            stop(paste("Feature", feature, "not found in", omics_type, "data"))
          }
          
          data[data[,1] == feature, -1]
        }
        
        x_vals <- get_feature_values(x_feature, values$corr_results$omics1)
        y_vals <- get_feature_values(y_feature, values$corr_results$omics2)
        
        plot_data <- data.frame(
          x = as.numeric(x_vals),
          y = as.numeric(y_vals),
          Sample = colnames(values$transcript_data)[-1]
        )
        
        # Merge with metadata if stratification is requested
        if (input$corr_strata != "None" && !is.null(values$metadata)) {
          strata_var <- values$metadata[[input$corr_strata]]
          plot_data$Strata <- strata_var[match(plot_data$Sample, values$metadata$SampleID)]
        }
        
        # Calculate correlation stats for display
        cor_test <- cor.test(plot_data$x, plot_data$y, method = values$corr_results$method)
        cor_stats <- paste0(
          "Correlation: ", round(cor_test$estimate, 3), "\n",
          "P-value: ", format.pval(cor_test$p.value, digits = 3), "\n",
          "N: ", nrow(plot_data)
        )
        
        output$scatter_stats <- renderText(cor_stats)
        
        # Create plot
        p <- ggplot(plot_data, aes(x = x, y = y, text = paste0(
          "Sample: ", Sample, "\n",
          x_feature, ": ", round(x, 2), "\n",
          y_feature, ": ", round(y, 2)
        ))) +
          geom_point(size = 3, alpha = 0.7) +
          geom_smooth(method = "lm", se = FALSE, color = "darkred") +
          theme_minimal() +
          labs(x = paste0(x_feature, " (", values$corr_results$omics1, ")"),
               y = paste0(y_feature, " (", values$corr_results$omics2, ")"))
        
        if (input$corr_strata != "None" && !is.null(values$metadata)) {
          p <- p + aes(color = Strata) +
            scale_color_brewer(palette = "Set1") +
            theme(legend.position = "bottom")
        }
        
        ggplotly(p, tooltip = "text") %>%
          layout(legend = list(orientation = "h", x = 0.5, y = -0.2))
      }, error = function(e) {
        showNotification(paste("Error generating scatter plot:", e$message), type = "error")
        NULL
      })
    })
    
    # Pathway analysis of correlated features
    observeEvent(input$run_pathway, {
      req(values$corr_results)
      
      tryCatch({
        showNotification("Running pathway analysis...", type = "message")
        
        # Get top correlated features
        top_features <- unique(c(
          values$corr_results$top_correlations$Feature1,
          values$corr_results$top_correlations$Feature2
        ))
        
        # Convert to gene symbols (for transcriptomics/proteomics)
        if (values$corr_results$omics1 %in% c("Transcriptomics", "Proteomics") || 
            values$corr_results$omics2 %in% c("Transcriptomics", "Proteomics")) {
          
          # Extract gene symbols from feature names
          gene_symbols <- tryCatch({
            unique(unlist(lapply(top_features, function(f) {
              # Handle different naming conventions
              if (grepl("\\(.*\\)", f)) {
                # For protein names like "P38398(BRCA1)"
                gsub(".*\\((.*)\\).*", "\\1", f)
              } else if (grepl("^ENSG", f)) {
                # For ENSEMBL IDs
                f
              } else {
                # Assume it's already a gene symbol
                f
              }
            })))
          }, error = function(e) {
            showNotification(paste("Error extracting gene symbols:", e$message), type = "error")
            NULL
          })
          
          if (is.null(gene_symbols)) return()
          
          # Run pathway analysis
          pathway_res <- tryCatch({
            switch(input$pathway_db,
                   "KEGG" = enrichKEGG(
                     gene = gene_symbols,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1
                   ),
                   "Reactome" = enrichPathway(
                     gene = gene_symbols,
                     organism = "human",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1,
                     readable = TRUE
                   ),
                   "GO Biological Process" = enrichGO(
                     gene = gene_symbols,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1,
                     readable = TRUE
                   )
            )
          }, error = function(e) {
            showNotification(paste("Error in pathway analysis:", e$message), type = "error")
            NULL
          })
          
          if (!is.null(pathway_res) && nrow(pathway_res) > 0) {
            pathway_df <- as.data.frame(pathway_res)
            values$pathway_results <- pathway_df
            
            # Create pathway plot
            output$pathway_plot <- renderPlotly({
              plot_df <- head(pathway_df, 20)
              plot_df$Description <- factor(plot_df$Description, 
                                            levels = rev(plot_df$Description))
              
              p <- ggplot(plot_df, aes(x = -log10(p.adjust), y = Description, 
                                       text = paste0("Description: ", Description, "\n",
                                                     "Genes: ", gsub("/", ", ", geneID), "\n",
                                                     "P.adjust: ", format.pval(p.adjust)))) +
                geom_col(fill = "steelblue") +
                theme_minimal() +
                labs(x = "-log10(Adjusted P-value)", y = "")
              
              ggplotly(p, tooltip = "text")
            })
            
            output$pathway_table <- renderDT({
              datatable(pathway_df[, c("Description", "GeneRatio", "pvalue", "p.adjust", "geneID")],
                        options = list(scrollX = TRUE, pageLength = 5))
            })
          } else {
            showNotification("No significant pathways found", type = "warning")
          }
        } else {
          showNotification("Pathway analysis is only available for transcriptomics/proteomics data", 
                           type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Error in pathway analysis:", e$message), type = "error")
      })
    })
    
    # Download handler
    output$download_corr <- downloadHandler(
      filename = function() {
        paste0("correlation_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(values$corr_results)
        tryCatch({
          write.csv(values$corr_results$all_correlations, file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Error downloading results:", e$message), type = "error")
        })
      }
    )
  })
}