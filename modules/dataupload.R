# modules/dataupload.R

source("modules/dataupload_helpers.R")
source("modules/dataupload_transcriptomics.R")
source("modules/dataupload_proteomics.R")
source("modules/dataupload_metabolomics.R")

dataUploadUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Upload and Process Omics Data"),
    fluidRow(
      transcriptomicsUI(ns("transcriptomics")),
      proteomicsUI(ns("proteomics")),
      metabolomicsUI(ns("metabolomics"))
    ),
    
    fluidRow(
      box(
        title = "Data Summary", width = 12, status = "info",
        tabsetPanel(
          id = ns("summary_tabs"),
          tabPanel("Overview",
                   h4("Loaded Datasets"),
                   tableOutput(ns("dataset_summary"))),
          tabPanel("Transcriptomics", 
                   DTOutput(ns("transcript_preview_table"))),
          tabPanel("Proteomics", 
                   DTOutput(ns("protein_preview_table"))),
          tabPanel("Metabolomics", 
                   DTOutput(ns("metabo_preview_table")))
        )
      )
    )
  )
}

dataUploadServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Call submodule servers
    transcriptomicsServer("transcriptomics", values)
    proteomicsServer("proteomics", values)
    metabolomicsServer("metabolomics", values)
    
    # Dataset summary
    output$dataset_summary <- renderTable({
      datasets <- c("Transcriptomics", "Proteomics", "Metabolomics")
      loaded <- c(!is.null(values$transcript_data),
                  !is.null(values$protein_data),
                  !is.null(values$metabo_data))
      features <- c(ifelse(loaded[1], nrow(values$transcript_data), NA),
                    ifelse(loaded[2], nrow(values$protein_data), NA),
                    ifelse(loaded[3], nrow(values$metabo_data), NA))
      samples <- c(ifelse(loaded[1], ncol(values$transcript_data)-1, NA),
                   ifelse(loaded[2], ncol(values$protein_data)-1, NA),
                   ifelse(loaded[3], ncol(values$metabo_data)-1, NA))
      
      data.frame(Dataset = datasets,
                 Loaded = ifelse(loaded, "Yes", "No"),
                 Features = features,
                 Samples = samples)
    })
    
    # Data previews
    output$transcript_preview_table <- renderDT({
      req(values$transcript_data)
      datatable(values$transcript_data, options = list(scrollX = TRUE, pageLength = 5))
    })
    
    output$protein_preview_table <- renderDT({
      req(values$protein_data)
      datatable(values$protein_data, options = list(scrollX = TRUE, pageLength = 5))
    })
    
    output$metabo_preview_table <- renderDT({
      req(values$metabo_data)
      datatable(values$metabo_data, options = list(scrollX = TRUE, pageLength = 5))
    })
  })
}