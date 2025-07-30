# modules/data_upload_metabolomics.R

metabolomicsUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "Metabolomics Data", width = 4, status = "success",
    fileInput(ns("file"), "Choose CSV/TSV File", accept = c(".csv", ".tsv")),
    checkboxInput(ns("header"), "Header", TRUE),
    selectInput(ns("sep"), "Separator",
                choices = c(Comma = ",", Semicolon = ";", Tab = "\t", Space = " "), 
                selected = ","),
    selectInput(ns("id_type"), "Identifier Type",
                choices = c("HMDB ID" = "hmdb",
                            "KEGG ID" = "kegg",
                            "PubChem CID" = "pubchem")),
    actionButton(ns("load"), "Load Metabolomics", 
                 class = "btn-primary"),
    uiOutput(ns("status"))
  )
}

metabolomicsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for tracking status
    status <- reactiveValues(loaded = FALSE)
    
    # Load data
    observeEvent(input$load, {
      df <- load_omics_data(input$file, 
                            input$header, 
                            input$sep,
                            "metabolomics",
                            input$id_type)
      
      if (!is.null(df)) {
        values$metabo_data <- df
        status$loaded <- TRUE
      }
    })
    
    # Status indicator
    output$status <- renderUI({
      if (status$loaded) {
        tags$div(icon("check-circle"), "Data loaded", style = "color: green;")
      } else {
        tags$div(icon("info-circle"), "Data not loaded", style = "color: gray;")
      }
    })
  })
}