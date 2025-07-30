
# modules/machine_learning.R

machineLearningUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Metadata-Informed Machine Learning"),
    fluidRow(
      box(
        title = "Model Configuration", width = 4, status = "primary",
        selectInput(ns("ml_task"), "Prediction Task",
                    choices = c("Classification", "Regression"),
                    selected = "Classification"),
        selectInput(ns("ml_outcome"), "Outcome Variable", choices = NULL),
        selectInput(ns("ml_algorithm"), "Algorithm",
                    choices = c("Random Forest", "SVM", "Elastic Net", "XGBoost"),
                    selected = "Random Forest"),
        selectInput(ns("ml_omics"), "Omics Data Source",
                    choices = c("Transcriptomics", "Proteomics", "Metabolomics", "Integrated"),
                    selected = "Integrated"),
        selectInput(ns("ml_covariates"), "Include Covariates", 
                    choices = NULL, multiple = TRUE),
        numericInput(ns("ml_cv"), "Cross-Validation Folds", value = 5),
        actionButton(ns("run_ml"), "Run Analysis")
      ),
      
      box(
        title = "Results", width = 8, status = "info",
        tabsetPanel(
          tabPanel("Model Performance", 
                   verbatimTextOutput(ns("ml_summary")),
                   plotlyOutput(ns("ml_perf_plot"))),
          tabPanel("Feature Importance", 
                   plotlyOutput(ns("feature_importance_plot")),
                   DTOutput(ns("feature_importance_table"))),
          tabPanel("Covariate Effects", 
                   plotlyOutput(ns("covariate_effect_plot")))
        )
      )
    )
  )
}

machineLearningServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_ml, {
      req(values$transcript_data, values$metadata, input$ml_outcome)
      
      tryCatch({
        # Prepare data
        outcome_var <- values$metadata[[input$ml_outcome]]
        features <- t(values$transcript_data[,-1])
        colnames(features) <- values$transcript_data[,1]
        
        # Add covariates if selected
        if (!is.null(input$ml_covariates)) {
          covariates <- values$metadata[, input$ml_covariates, drop = FALSE]
          features <- cbind(features, covariates)
        }
        
        # Remove samples with missing outcome
        complete_cases <- !is.na(outcome_var)
        features <- features[complete_cases,]
        outcome_var <- outcome_var[complete_cases]
        
        # Create train control
        train_control <- trainControl(
          method = "cv",
          number = input$ml_cv,
          savePredictions = "final",
          classProbs = input$ml_task == "Classification"
        )
        
        # Train model based on selected algorithm
        if (input$ml_algorithm == "Random Forest") {
          model <- caret::train(
            x = features,
            y = outcome_var,
            method = "rf",
            trControl = train_control,
            importance = TRUE
          )
        } else if (input$ml_algorithm == "SVM") {
          model <- caret::train(
            x = features,
            y = outcome_var,
            method = "svmRadial",
            trControl = train_control,
            preProcess = c("center", "scale")
          )
        } else if (input$ml_algorithm == "Elastic Net") {
          model <- caret::train(
            x = features,
            y = outcome_var,
            method = "glmnet",
            trControl = train_control,
            tuneLength = 5
          )
        } else if (input$ml_algorithm == "XGBoost") {
          model <- caret::train(
            x = features,
            y = outcome_var,
            method = "xgbTree",
            trControl = train_control,
            tuneLength = 3
          )
        }
        
        values$ml_model <- model
        showNotification("Machine learning analysis completed", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in machine learning:", e$message), type = "error")
      })
    })
    
    output$ml_summary <- renderPrint({
      req(values$ml_model)
      print(values$ml_model)
    })
    
    output$ml_perf_plot <- renderPlotly({
      req(values$ml_model)
      ggplotly(plot(values$ml_model))
    })
    
    output$feature_importance_plot <- renderPlotly({
      req(values$ml_model)
      
      imp <- varImp(values$ml_model)$importance
      imp$Feature <- rownames(imp)
      imp <- imp[order(imp$Overall, decreasing = TRUE),]
      imp <- head(imp, 20)
      
      p <- ggplot(imp, aes(x = Overall, y = reorder(Feature, Overall))) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal() +
        labs(x = "Importance", y = "Feature")
      
      ggplotly(p)
    })
    
    output$feature_importance_table <- renderDT({
      req(values$ml_model)
      
      imp <- varImp(values$ml_model)$importance
      imp$Feature <- rownames(imp)
      imp <- imp[order(imp$Overall, decreasing = TRUE),]
      
      datatable(imp, options = list(scrollX = TRUE))
    })
  })
}