library(shiny)

source("preprocessing_tab.R")
source("quantitative_analysis_tab.R")
source("network_analysis_tab.R")        
source("functional_analysis_tab.R")

ui <- fluidPage(
  titlePanel("PiProteline Shiny App"),
  tabsetPanel(
    preprocessing_ui(),
    quantitative_analysis_ui(),
    network_analysis_ui(),        
    functional_analysis_ui()      
  )
)

server <- function(input, output, session) {
  shared_data <- preprocessing_server(input, output, session)
  
  quantitative_analysis_server(input, output, session, shared_data)
  network_analysis_server(input, output, session, shared_data)     
  functional_analysis_server(input, output, session, shared_data)  
}

shinyApp(ui = ui, server = server)
