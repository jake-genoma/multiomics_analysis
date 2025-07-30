# modules/dataupload_proteomics.R

proteomicsUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "Proteomics Data", width = 4, status = "warning",
    selectInput(ns("species"), "Species",
                choices = c("Human" = "human", "Mouse" = "mouse", "Rat" = "rat"),
                selected = "human"),
    fileInput(ns("file"), "Choose CSV/TSV File", accept = c(".csv", ".tsv")),
    checkboxInput(ns("header"), "Header", TRUE),
    selectInput(ns("sep"), "Separator",
                choices = c(Comma = ",", Semicolon = ";", Tab = "\t", Space = " "), 
                selected = ","),
    selectInput(ns("id_type"), "Identifier Type",
                choices = c("UniProt ID" = "uniprot",
                            "Gene Symbol" = "symbol",
                            "Ensembl Protein ID" = "ensembl_protein")),
    actionButton(ns("load"), "Load Proteomics", 
                 class = "btn-primary"),
    actionButton(ns("view_as_symbol"), "View as Symbol", 
                 class = "btn-info"),
    uiOutput(ns("status")),
    conditionalPanel(
      condition = paste0("output['", ns("conversion_needed"), "']"),
      ns = ns,
      checkboxInput(ns("aggregate"), "Aggregate multiple matches (mean)"),
      verbatimTextOutput(ns("conversion_stats"))
    )
  )
}

proteomicsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for tracking status
    status <- reactiveValues(
      loaded = FALSE,
      converted = FALSE
    )
    
    # Load data
    observeEvent(input$load, {
      df <- load_omics_data(input$file, 
                            input$header, 
                            input$sep,
                            "proteomics",
                            input$id_type,
                            input$species)
      
      if (!is.null(df)) {
        values$protein_data <- df
        status$loaded <- TRUE
        status$converted <- FALSE
      }
    })
    
    # Convert to symbols
    observeEvent(input$view_as_symbol, {
      req(values$protein_data, status$loaded)
      
      tryCatch({
        result <- convert_to_symbols(values$protein_data, input$aggregate)
        
        values$protein_data <- result$data
        status$converted <- TRUE
        
        output$conversion_stats <- renderPrint({
          cat("Conversion Statistics:\n")
          cat("Original IDs:", result$stats$original_ids, "\n")
          cat("Mapped to symbols:", result$stats$mapped, "\n")
          cat("Unique symbols:", result$stats$unique_symbols, "\n")
          cat("Multi-mapped IDs:", result$stats$multi_mapped, "\n")
        })
        
        showNotification("Successfully converted to gene symbols", type = "message")
      }, error = function(e) {
        showNotification(e$message, type = "error")
      })
    })
    
    # Status indicator
    output$status <- renderUI({
      if (status$loaded) {
        if (status$converted) {
          tags$div(icon("check-circle"), "Data loaded and converted to symbols", style = "color: green;")
        } else {
          tags$div(icon("check-circle"), "Data loaded", style = "color: green;")
        }
      } else {
        tags$div(icon("info-circle"), "Data not loaded", style = "color: gray;")
      }
    })
    
    # Determine if conversion is needed
    output$conversion_needed <- reactive({
      status$loaded && !status$converted
    })
    outputOptions(output, "conversion_needed", suspendWhenHidden = FALSE)
  })
}