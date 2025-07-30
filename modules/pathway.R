# pathway.R

source("modules/pathway_parameters.R")
source("modules/pathway_selection.R")
source("modules/pathway_visualization.R")

pathwayUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Pathway Analysis and Visualization"),
    fluidRow(
      column(4, pathwayParametersUI(ns("params"))),
      column(4, pathwaySelectionUI(ns("selection"))),
      column(4, pathwayVisualizationUI(ns("visualization")))
    )
  )
}

pathwayServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Run pathway analysis and get results
    pathway_results <- pathwayParametersServer("params", values)
    
    # Display pathway selection and get selected index
    selected_pathway_index <- pathwaySelectionServer("selection", pathway_results)
    
    # Visualize selected pathway
    pathwayVisualizationServer("visualization", pathway_results, selected_pathway_index)
  })
}