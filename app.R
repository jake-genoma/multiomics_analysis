# app.R
source("app_ui.R")
source("app_server.R")

# Run the application
shinyApp(ui = ui, server = server)