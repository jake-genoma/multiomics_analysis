# modules/metadata.R
metadataUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Sample Metadata Management"),
    fluidRow(
      box(
        title = "Metadata Upload", width = 4, status = "info",
        fileInput(ns("metadata_file"), "Upload Metadata File", accept = ".csv"),
        checkboxInput(ns("metadata_header"), "Header", TRUE),
        selectInput(ns("metadata_sep"), "Separator",
                    choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), 
                    selected = ","),
        actionButton(ns("metadata_preview"), "Preview Metadata"),
        hr(),
        h4("Data Status"),
        uiOutput(ns("sample_consistency")),
        uiOutput(ns("sample_mismatch"))
      ),
      
      box(
        title = "Metadata Summary", width = 8, status = "info",
        verbatimTextOutput(ns("metadata_summary")),
        plotlyOutput(ns("meta_distribution_plot")),
        DTOutput(ns("metadata_preview_table"))
      )
    ),
    
    fluidRow(
      box(
        title = "Metadata Visualization", width = 12, status = "info",
        selectInput(ns("meta_var_x"), "X-axis Variable", choices = NULL),
        selectInput(ns("meta_var_y"), "Y-axis Variable", choices = NULL),
        selectInput(ns("meta_var_color"), "Color By", choices = NULL),
        plotlyOutput(ns("meta_scatter_plot"))
      )
    )
  )
}

metadataServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$metadata_preview, {
      values$metadata <- read_and_validate(
        input$metadata_file, 
        input$metadata_header, 
        input$metadata_sep,
        "metadata"
      )
      
      # Update all dropdowns that depend on metadata
      update_metadata_dropdowns(session, values$metadata)
    })
    
    output$metadata_preview_table <- renderDT({
      req(values$metadata)
      datatable(values$metadata, options = list(scrollX = TRUE))
    })
    
    # Sample consistency check
    output$sample_consistency <- renderUI({
      req(values$transcript_data, values$metadata)
      
      transcript_samples <- colnames(values$transcript_data)[-1]
      metadata_samples <- values$metadata$SampleID
      
      if (setequal(transcript_samples, metadata_samples)) {
        tags$div(
          style = "color: green;",
          icon("check-circle"), 
          "All samples match between transcriptomics data and metadata"
        )
      } else {
        tags$div(
          style = "color: red;",
          icon("exclamation-triangle"), 
          "Sample mismatch between transcriptomics data and metadata"
        )
      }
    })
    
    # Metadata summary and visualization
    output$metadata_summary <- renderPrint({
      req(values$metadata)
      summary(values$metadata)
    })
    
    output$meta_distribution_plot <- renderPlotly({
      req(values$metadata, input$meta_var_x)
      
      if (is.numeric(values$metadata[[input$meta_var_x]])) {
        p <- ggplot(values$metadata, aes_string(x = input$meta_var_x)) +
          geom_histogram(fill = "steelblue", bins = 30) +
          theme_minimal()
      } else {
        p <- ggplot(values$metadata, aes_string(x = input$meta_var_x)) +
          geom_bar(fill = "steelblue") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      
      ggplotly(p)
    })
    
    output$meta_scatter_plot <- renderPlotly({
      req(values$metadata, input$meta_var_x, input$meta_var_y)
      
      p <- ggplot(values$metadata, aes_string(x = input$meta_var_x, y = input$meta_var_y)) +
        theme_minimal()
      
      if (!is.null(input$meta_var_color)) {
        p <- p + geom_point(aes_string(color = input$meta_var_color), size = 3)
      } else {
        p <- p + geom_point(size = 3)
      }
      
      ggplotly(p)
    })
  })
}