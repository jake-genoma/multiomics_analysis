# modules/differential.R

differentialUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Differential Analysis by Metadata"),
    fluidRow(
      box(
        title = "Analysis Setup", width = 4, status = "primary",
        selectInput(ns("diff_omics"), "Omics Layer",
                    choices = c("Transcriptomics", "Proteomics", "Metabolomics"),
                    selected = "Transcriptomics"),
        selectInput(ns("diff_group"), "Comparison Group", choices = NULL),
        selectInput(ns("diff_ref"), "Reference Level", choices = NULL),
        selectInput(ns("diff_covariates"), "Adjust for Covariates", 
                    choices = NULL, multiple = TRUE),
        numericInput(ns("diff_pval"), "P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("diff_fc"), "Min Fold Change", value = 1.5, min = 1, step = 0.1),
        radioButtons(ns("diff_method"), "Method",
                     choices = c("limma (microarray-like)" = "limma",
                                 "DESeq2 (RNA-seq)" = "deseq2",
                                 "EdgeR (RNA-seq)" = "edger"),
                     selected = "limma"),
        actionButton(ns("run_diff"), "Run Differential Analysis", class = "btn-primary"),
        downloadButton(ns("download_results"), "Download Results", class = "btn-success")
      ),
      
      box(
        title = "Results", width = 8, status = "info",
        tabsetPanel(
          tabPanel("Volcano Plot", plotlyOutput(ns("volcano_plot"))),
          tabPanel("Top Features", 
                   DTOutput(ns("diff_results_table")),
                   downloadButton(ns("download_table"), "Download Table")),
          tabPanel("Feature Plot", 
                   selectInput(ns("diff_feature"), "Select Feature", choices = NULL),
                   plotlyOutput(ns("feature_boxplot"))),
          tabPanel("QC Plots",
                   plotOutput(ns("pval_histogram")),
                   plotOutput(ns("ma_plot")))
        )
      )
    )
  )
}

differentialServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update dropdowns when metadata changes
    observeEvent(values$metadata, {
      req(values$metadata)
      meta_vars <- colnames(values$metadata)[-1] # exclude SampleID
      updateSelectInput(session, "diff_group", choices = meta_vars)
      updateSelectInput(session, "diff_covariates", choices = meta_vars)
    })
    
    # Update reference level when group changes
    observeEvent(input$diff_group, {
      req(values$metadata, input$diff_group)
      group_levels <- unique(values$metadata[[input$diff_group]])
      updateSelectInput(session, "diff_ref", choices = group_levels)
    })
    
    # Run differential analysis
    observeEvent(input$run_diff, {
      req(values$metadata, input$diff_group, input$diff_ref)
      
      tryCatch({
        showModal(modalDialog("Running differential analysis...", footer = NULL))
        
        # Get the appropriate omics data
        omics_data <- switch(input$diff_omics,
                             "Transcriptomics" = values$transcript_data,
                             "Proteomics" = values$protein_data,
                             "Metabolomics" = values$metabo_data)
        
        if (is.null(omics_data)) {
          stop(paste(input$diff_omics, "data not loaded"))
        }
        
        # Prepare metadata - ensure SampleID column exists
        if (!"SampleID" %in% colnames(values$metadata)) {
          stop("Metadata must contain 'SampleID' column")
        }
        meta <- values$metadata
        rownames(meta) <- meta$SampleID
        
        # Ensure samples match between omics and metadata
        samples <- colnames(omics_data)[-1]  # Assuming first column is feature IDs
        common_samples <- intersect(samples, meta$SampleID)
        
        if (length(common_samples) == 0) {
          stop("No matching samples between omics data and metadata")
        }
        
        # Subset data to common samples
        meta <- meta[common_samples, , drop = FALSE]
        omics_data <- omics_data[, c(1, which(colnames(omics_data) %in% common_samples))]
        
        # Extract expression matrix and ensure numeric
        expr_mat <- as.matrix(omics_data[, -1, drop = FALSE])
        rownames(expr_mat) <- omics_data[, 1]
        expr_mat <- apply(expr_mat, 2, as.numeric)  # Ensure all values are numeric
        rownames(expr_mat) <- omics_data[, 1]  # Restore rownames
        
        # Create group factor - handle NAs
        group_var <- meta[[input$diff_group]]
        if (any(is.na(group_var))) {
          showNotification("NA values found in group variable - these samples will be removed", 
                           type = "warning")
          keep <- !is.na(group_var)
          group_var <- group_var[keep]
          expr_mat <- expr_mat[, keep, drop = FALSE]
          meta <- meta[keep, , drop = FALSE]
        }
        
        # Ensure we have at least 2 samples per group after NA removal
        group_var <- factor(group_var)
        if (length(levels(group_var)) < 2) {
          stop("Comparison group must have at least 2 levels after NA removal")
        }
        
        # Relevel with reference
        group_var <- relevel(group_var, ref = input$diff_ref)
        
        # Handle covariates if specified
        if (!is.null(input$diff_covariates)) {
          # Check covariates exist and have no NAs
          covars <- meta[, input$diff_covariates, drop = FALSE]
          if (any(is.na(covars))) {
            stop("NA values found in covariates - please handle these first")
          }
          covar_formula <- paste(input$diff_covariates, collapse = " + ")
          design_formula <- as.formula(paste("~", covar_formula, "+ group_var"))
        } else {
          design_formula <- ~ group_var
        }
        
        # Create design matrix
        design <- model.matrix(design_formula, data = meta)
        
        # Check for rank deficiency
        if (qr(design)$rank < ncol(design)) {
          stop("Design matrix is rank deficient - check for colinear variables")
        }
        
        # Perform limma analysis
        if (input$diff_method == "limma") {
          fit <- lmFit(expr_mat, design)
          fit <- eBayes(fit)
          
          # Find the coefficient for the group comparison
          coef_name <- paste0("group_var", levels(group_var)[2])
          if (!coef_name %in% colnames(design)) {
            coef_name <- "group_var"  # Fallback if simple design
          }
          
          results <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
          
        } else if (input$diff_method %in% c("deseq2", "edger")) {
          stop("Selected method requires count data. Please use limma for normalized data.")
        }
        
        # Add feature names and format results
        results$Feature <- rownames(results)
        results$OmicsLayer <- input$diff_omics
        results$ComparisonGroup <- input$diff_group
        results$ReferenceLevel <- input$diff_ref
        
        # Calculate -log10 p-value and significance
        results$negLogP <- -log10(results$P.Value)
        results$Significant <- ifelse(results$adj.P.Val < input$diff_pval & 
                                        abs(results$logFC) > log2(input$diff_fc), 
                                      "Yes", "No")
        
        # Store results
        values$diff_results <- results
        
        # Update feature selection dropdown
        sig_features <- results$Feature[results$Significant == "Yes"]
        if (length(sig_features) > 0) {
          updateSelectInput(session, "diff_feature", 
                            choices = sig_features[1:min(50, length(sig_features))])
        } else {
          updateSelectInput(session, "diff_feature", 
                            choices = results$Feature[1:min(50, nrow(results))])
          showNotification("No significant features found at current thresholds", 
                           type = "warning")
        }
        
        removeModal()
        showNotification("Differential analysis completed", type = "message")
        
      }, error = function(e) {
        removeModal()
        showNotification(paste("Error in differential analysis:", e$message), type = "error")
      })
    })
    
    # Volcano plot
    output$volcano_plot <- renderPlotly({
      req(values$diff_results)
      
      tryCatch({
        results <- values$diff_results
        
        p <- ggplot(results, aes(x = logFC, y = negLogP, 
                                 color = Significant, 
                                 text = paste("Feature:", Feature, "<br>",
                                              "logFC:", round(logFC, 2), "<br>",
                                              "P-value:", format.pval(P.Value), "<br>",
                                              "Adj.P-value:", format.pval(adj.P.Val)))) +
          geom_point(alpha = 0.6, size = 2) +
          scale_color_manual(values = c("Yes" = "red", "No" = "gray60")) +
          geom_hline(yintercept = -log10(input$diff_pval), linetype = "dashed", color = "black") +
          geom_vline(xintercept = c(-log2(input$diff_fc), log2(input$diff_fc)), 
                     linetype = "dashed", color = "black") +
          theme_minimal() +
          labs(x = "log2 Fold Change", 
               y = "-log10 p-value",
               title = paste("Volcano Plot:", input$diff_group, "vs", input$diff_ref),
               color = "Significant")
        
        ggplotly(p, tooltip = "text") %>%
          layout(hoverlabel = list(bgcolor = "white"))
        
      }, error = function(e) {
        plotly_empty() %>%
          layout(title = paste("Error creating volcano plot:", e$message))
      })
    })
    
    # Results table
    output$diff_results_table <- renderDT({
      req(values$diff_results)
      
      # Ensure dplyr is available
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("dplyr package is required for this functionality")
      }
      
      # Extract results and format
      results <- values$diff_results
      results <- results[, c("Feature", "logFC", "P.Value", "adj.P.Val", "Significant", 
                             "OmicsLayer", "ComparisonGroup", "ReferenceLevel")]
      
      # Round numeric columns
      numeric_cols <- sapply(results, is.numeric)
      results[numeric_cols] <- lapply(results[numeric_cols], function(x) round(x, 4))
      
      # Create the datatable
      DT::datatable(
        results,
        filter = 'top',
        extensions = 'Buttons',
        options = list(
          scrollX = TRUE,
          pageLength = 10,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        )
      ) %>%
        DT::formatStyle(
          'Significant',
          backgroundColor = DT::styleEqual(c("Yes", "No"), c('#FF9999', '#CCCCCC'))
        ) %>%
        DT::formatStyle(
          c('P.Value', 'adj.P.Val'),
          color = DT::styleInterval(c(0.05), c('red', 'black'))
        )
    })
    
    # Feature boxplot
    output$feature_boxplot <- renderPlotly({
      req(values$metadata, input$diff_feature, input$diff_group)
      
      # Get the appropriate omics data
      omics_data <- switch(input$diff_omics,
                           "Transcriptomics" = values$transcript_data,
                           "Proteomics" = values$protein_data,
                           "Metabolomics" = values$metabo_data)
      
      feature_row <- which(omics_data[,1] == input$diff_feature)
      if (length(feature_row) == 0) return(NULL)
      
      plot_data <- data.frame(
        Expression = as.numeric(omics_data[feature_row, -1]),
        Group = values$metadata[[input$diff_group]][match(colnames(omics_data)[-1], values$metadata$SampleID)],
        Sample = colnames(omics_data)[-1]
      )
      
      # Remove NA values
      plot_data <- na.omit(plot_data)
      
      p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group, text = Sample)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        theme_minimal() +
        labs(title = input$diff_feature, 
             y = ifelse(input$diff_omics == "Transcriptomics", "log2 Expression", "Intensity"),
             x = input$diff_group) +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = "text")
    })
    
    # P-value histogram
    output$pval_histogram <- renderPlot({
      req(values$diff_results)
      
      ggplot(values$diff_results, aes(x = P.Value)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white") +
        theme_minimal() +
        labs(title = "P-value Distribution",
             x = "P-value",
             y = "Count")
    })
    
    # MA plot
    output$ma_plot <- renderPlot({
      req(values$diff_results)
      
      tryCatch({
        # Get the appropriate omics data
        omics_data <- switch(input$diff_omics,
                             "Transcriptomics" = values$transcript_data,
                             "Proteomics" = values$protein_data,
                             "Metabolomics" = values$metabo_data)
        
        if (is.null(omics_data)) {
          stop("Omics data not available for MA plot")
        }
        
        # Prepare the expression matrix
        expr_mat <- as.matrix(omics_data[, -1])
        rownames(expr_mat) <- omics_data[, 1]
        
        # Calculate average expression
        results <- values$diff_results
        results$AveExpr <- rowMeans(expr_mat, na.rm = TRUE)
        
        # Create the plot
        ggplot(results, aes(x = AveExpr, y = logFC, color = Significant)) +
          geom_point(alpha = 0.6) +
          scale_color_manual(values = c("Yes" = "red", "No" = "gray60")) +
          geom_hline(yintercept = 0, color = "black") +
          theme_minimal() +
          labs(title = "MA Plot",
               x = "Average Expression",
               y = "log2 Fold Change")
        
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = paste("Error creating MA plot:", e$message)) +
          theme_void()
      })
    })
    
    # Download handlers
    output$download_results <- downloadHandler(
      filename = function() {
        paste("differential_results_", input$diff_omics, "_", input$diff_group, "_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(values$diff_results, file, row.names = FALSE)
      }
    )
    
    output$download_table <- downloadHandler(
      filename = function() {
        paste("differential_table_", input$diff_omics, "_", input$diff_group, "_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(values$diff_results, file, row.names = FALSE)
      }
    )
  })
}