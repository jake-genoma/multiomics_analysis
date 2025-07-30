# app_ui.R
source("globals.R")
source("modules/dataupload.R")
source("modules/metadata.R")
source("modules/correlation.R")
source("modules/differential.R")
source("modules/pathway.R")
source("modules/machine_learning.R")
source("modules/network.R")
source("modules/ai_insights.R")

ui <- dashboardPage(
  dashboardHeader(title = "Multi-Omics Analyzer AI"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Metadata", tabName = "metadata", icon = icon("table")),
      menuItem("Correlation Analysis", tabName = "correlation", icon = icon("project-diagram")),
      menuItem("Differential Analysis", tabName = "differential", icon = icon("chart-bar")),
      menuItem("Pathway Analysis", tabName = "pathway", icon = icon("sitemap")),
      menuItem("Machine Learning", tabName = "ml", icon = icon("brain")),
      menuItem("Network Visualization", tabName = "network", icon = icon("network-wired")),
      menuItem("AI Insights", tabName = "ai_insights", icon = icon("robot"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    useShinyFeedback(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css/styles.css"),
      tags$script(src = "js/custom.js")
    ),
    
    tabItems(
      dataUploadUI("upload"),
      metadataUI("metadata"),
      correlationUI("correlation"),
      differentialUI("differential"),
      pathwayUI("pathway"),
      machineLearningUI("ml"),
      networkUI("network"),
      aiInsightsUI("ai_insights")
    )
  )
)