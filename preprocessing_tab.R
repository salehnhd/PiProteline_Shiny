# preprocessing_tab.R
library(shiny)
library(DT)
library(openxlsx)
library(shinyWidgets)
library(PiProteline)

preprocessing_ui <- function() {
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        fileInput("dataset", "Upload dataset (.xlsx, .xls, .csv)",
                  accept = c(".xlsx", ".xls", ".csv")),
        uiOutput("gene_column_ui"),   
        
        textInput("group_names", "Group names (comma-separated)",
                  value = "PSM_C,PSM_PD,PSM_PG"),
        
        selectInput("norm_type", "Normalization type",
                    choices = c("ln","Znorm","MinMax","Robust","UnitVector",
                                "TotSigNorm","MaxSigNorm","RowSigmaNorm"),
                    selected = "Znorm"),
        
        actionBttn("preprocess", "Preprocess data", style = "fill", color = "primary"),
        
        tags$hr(),
        actionButton("show_data_unique", "Show data_unique"),
        actionButton("show_data_norm",   "Show data_norm"),
        tags$hr(),
        actionBttn("desc_stats", "Descriptive statistics", style = "fill", color = "success"),
        actionButton("show_desc_stats_col", "Show DS_col"),
        actionButton("show_desc_stats_row", "Show DS_row")
      ),
      mainPanel(
        verbatimTextOutput("debug_output"),
        
        conditionalPanel("input.show_data_unique > 0",
                         h4("data_unique"), DTOutput("data_unique_table"),
                         downloadButton("download_data_unique","Download")
        ),
        conditionalPanel("input.show_data_norm > 0",
                         h4("data_norm"), DTOutput("data_norm_table"),
                         downloadButton("download_data_norm","Download")
        ),
        conditionalPanel("input.show_desc_stats_col > 0",
                         h4("Descriptive Stats – Columns"), DTOutput("desc_stats_col_table"),
                         downloadButton("download_desc_stats_col","Download")
        ),
        conditionalPanel("input.show_desc_stats_row > 0",
                         h4("Descriptive Stats – Rows"), DTOutput("desc_stats_row_table"),
                         downloadButton("download_desc_stats_row","Download")
        )
      )
    )
  )
}

preprocessing_server <- function(input, output, session) {
  
  shared_data <- reactiveValues(
    dataset      = NULL,
    group_names  = NULL,    
    gene_column  = NULL,    
    preproc_data = NULL,
    desc_stats   = NULL
  )
  
  output$debug_output <- renderText({
    paste("Temp dir:", tempdir())
  })
  
  # Load dataset
  observeEvent(input$dataset, {
    req(input$dataset)
    ext <- tolower(tools::file_ext(input$dataset$name))
    fp  <- input$dataset$datapath
    shared_data$dataset <- switch(
      ext,
      xlsx = openxlsx::read.xlsx(fp),
      xls  = openxlsx::read.xlsx(fp),
      csv  = read.csv(fp, check.names = FALSE),
      { showNotification("Unsupported file type", type = "error"); NULL }
    )
  })
  
  
  output$gene_column_ui <- renderUI({
    req(shared_data$dataset)
    cols <- colnames(shared_data$dataset)
    default <- if ("GeneName" %in% cols) "GeneName" else cols[1]
    selectInput("gene_column", "Select gene column", choices = cols, selected = default)
  })
  
 
  observe({ req(input$group_names); shared_data$group_names <- trimws(strsplit(input$group_names, ",")[[1]]) })
  observe({ shared_data$gene_column <- input$gene_column })
  
 
  observeEvent(input$preprocess, {
    req(shared_data$dataset, shared_data$group_names, shared_data$gene_column, input$norm_type)
    
    withProgress(message = "Preprocessing", value = 0, {
      tryCatch({
        res <- PiProteline::preprocessing_data(
          dataset         = shared_data$dataset,
          names_of_groups = shared_data$group_names,
          gene_column     = shared_data$gene_column,
          norm_type       = input$norm_type
        )
        shared_data$preproc_data <- res
        
       
        shared_data$gene_column   <- colnames(res$data_norm)[1]
        
        showNotification("Preprocessing done", type = "message")
      }, error = function(e) {
        showNotification(paste("Preprocessing error:", e$message), type = "error")
        shared_data$preproc_data <- NULL
      })
    })
  })
  
  
  observeEvent(input$desc_stats, {
    req(shared_data$preproc_data)
    tryCatch({
      result <- PiProteline::descriptive_statistics(
        shared_data$preproc_data$data_unique,
        shared_data$preproc_data$data_grouped_full
      )
      shared_data$desc_stats <- result
      showNotification("Descriptive statistics computed", type = "message")
    }, error = function(e) {
      showNotification(paste("Desc. stats error:", e$message), type = "error")
      shared_data$desc_stats <- NULL
    })
  })
  
  
  output$data_unique_table <- renderDT({ req(shared_data$preproc_data, input$show_data_unique); shared_data$preproc_data$data_unique })
  output$data_norm_table   <- renderDT({ req(shared_data$preproc_data, input$show_data_norm);   shared_data$preproc_data$data_norm })
  output$desc_stats_col_table <- renderDT({ req(shared_data$desc_stats, input$show_desc_stats_col); shared_data$desc_stats$DS_col })
  output$desc_stats_row_table <- renderDT({ req(shared_data$desc_stats, input$show_desc_stats_row); shared_data$desc_stats$DS_row })
  
  output$download_data_unique <- downloadHandler("data_unique.csv", function(f) write.csv(shared_data$preproc_data$data_unique, f, row.names = FALSE))
  output$download_data_norm   <- downloadHandler("data_norm.csv",   function(f) write.csv(shared_data$preproc_data$data_norm,   f, row.names = FALSE))
  output$download_desc_stats_col <- downloadHandler("DS_col.csv", function(f) write.csv(shared_data$desc_stats$DS_col, f, row.names = FALSE))
  output$download_desc_stats_row <- downloadHandler("DS_row.csv", function(f) write.csv(shared_data$desc_stats$DS_row, f, row.names = FALSE))
  
  return(shared_data)
}
