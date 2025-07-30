# modules/network.R

networkUI <- function(id) {
  ns <- NS(id)
  
  tabItem(
    tabName = id,
    h2("Metadata-Annotated Networks"),
    fluidRow(
      box(
        title = "Network Parameters", width = 4, status = "primary",
        selectInput(ns("network_type"), "Network Type",
                    choices = c("Co-expression", "Multi-omics", "PPI"),
                    selected = "Multi-omics"),
        selectInput(ns("network_strata"), "Color Nodes By", choices = NULL),
        numericInput(ns("network_threshold"), "Edge Threshold", 
                     value = 0.7, min = 0, max = 1, step = 0.05),
        actionButton(ns("generate_network"), "Generate Network")
      ),
      
      box(
        title = "Network Visualization", width = 8, status = "info",
        visNetworkOutput(ns("network_plot"), height = "600px"),
        downloadButton(ns("download_network"), "Download Network")
      )
    ),
    
    fluidRow(
      box(
        title = "Network Statistics by Metadata", width = 12, status = "info",
        plotlyOutput(ns("network_stats_plot")),
        DTOutput(ns("network_stats_table")))
    )
  )
}

networkServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$generate_network, {
      req(values$corr_results, values$metadata)
      
      tryCatch({
        # Create a correlation network
        cor_mat <- values$corr_results$matrix
        
        # Apply threshold
        cor_mat[abs(cor_mat) < input$network_threshold] <- 0
        
        # Create nodes and edges
        nodes <- data.frame(
          id = c(rownames(cor_mat), colnames(cor_mat)),
          label = c(rownames(cor_mat), colnames(cor_mat)),
          group = c(rep(values$corr_results$omics1, nrow(cor_mat)),
                    rep(values$corr_results$omics2, ncol(cor_mat)))
        )
        
        edges <- reshape2::melt(cor_mat)
        edges <- edges[edges$value != 0, ]
        colnames(edges) <- c("from", "to", "value")
        edges$width <- abs(edges$value) * 5
        edges$color <- ifelse(edges$value > 0, "green", "red")
        
        # Add metadata if stratification is selected
        if (input$network_strata != "None") {
          # This is a simplified version - in practice you'd need to map features to samples
          # and then to metadata, which can be complex
          nodes$title <- input$network_strata
        }
        
        values$network <- list(nodes = nodes, edges = edges)
        showNotification("Network generated", type = "message")
      }, error = function(e) {
        showNotification(paste("Error generating network:", e$message), type = "error")
      })
    })
    
    output$network_plot <- renderVisNetwork({
      req(values$network)
      
      visNetwork(values$network$nodes, values$network$edges) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(solver = "forceAtlas2Based", 
                   forceAtlas2Based = list(gravitationalConstant = -50))
    })
    
    output$download_network <- downloadHandler(
      filename = function() {
        paste("network-", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
        net <- visNetwork(values$network$nodes, values$network$edges) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visPhysics(solver = "forceAtlas2Based", 
                     forceAtlas2Based = list(gravitationalConstant = -50))
        
        visSave(net, file)
      }
    )
  }) 
}