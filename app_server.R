# app_server.R
source("globals.R")
source("helpers.R")
source("modules/dataupload.R")
source("modules/metadata.R")
source("modules/correlation.R")
source("modules/differential.R")
source("modules/pathway.R")
source("modules/machine_learning.R")
source("modules/network.R")
source("modules/ai_insights.R")

server <- function(input, output, session) {
  # Create reactive values to store data
  values <- reactiveValues(
    transcript_data = NULL,
    protein_data = NULL,
    metabo_data = NULL,
    metadata = NULL,
    merged_data = NULL,
    diff_results = NULL,
    corr_results = NULL,
    pathway_results = NULL,
    ml_model = NULL,
    network = NULL,
    ai_insights = NULL
  )
  
  # Call module servers
  dataUploadServer("upload", values)
  metadataServer("metadata", values)
  correlationServer("correlation", values)
  differentialServer("differential", values)
  pathwayServer("pathway", values)
  machineLearningServer("ml", values)
  networkServer("network", values)
  aiInsightsServer("ai_insights", values)
}