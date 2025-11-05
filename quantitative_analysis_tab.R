# quantitative_analysis_tab.R
library(shiny)
library(DT)
library(shinyWidgets)
library(PiProteline)

quantitative_analysis_ui <- function() {
  tabPanel(
    "Quantitative Analysis",
    sidebarLayout(
      sidebarPanel(
        numericInput("significance_manova", "Significance threshold (MANOVA)",
                     value = 0.05, min = 0.001, max = 1, step = 0.001),
        actionBttn("run_analysis", "Run Quantitative Analysis", style = "fill", color = "primary"),
        
        tags$hr(),
        actionButton("show_manova", "Show MANOVA Results"),
        
        tags$hr(),
        actionButton("show_volcano_1", "Show Volcano: Group1 vs Group2"),
        actionButton("show_volcano_2", "Show Volcano: Group1 vs Group3"),
        actionButton("show_volcano_3", "Show Volcano: Group2 vs Group3"),
        
        tags$hr(),
        actionButton("show_mds",   "Show MDS Plot"),
        actionButton("show_mds_1", "Show Pairwise MDS: Group1 vs Group2"),
        actionButton("show_mds_2", "Show Pairwise MDS: Group1 vs Group3")
      ),
      mainPanel(
        verbatimTextOutput("qa_debug"),
        
        conditionalPanel("input.show_manova > 0",
                         DTOutput("manova_table"),
                         downloadButton("download_manova", "Download MANOVA Results")
        ),
        
        conditionalPanel("input.show_volcano_1 > 0", plotOutput("volcano1"), downloadButton("download_volcano1", "Download Volcano 1")),
        conditionalPanel("input.show_volcano_2 > 0", plotOutput("volcano2"), downloadButton("download_volcano2", "Download Volcano 2")),
        conditionalPanel("input.show_volcano_3 > 0", plotOutput("volcano3"), downloadButton("download_volcano3", "Download Volcano 3")),
        
        conditionalPanel("input.show_mds > 0",   plotOutput("mds_plot"), downloadButton("download_mds", "Download MDS Plot")),
        conditionalPanel("input.show_mds_1 > 0", plotOutput("mds1"),     downloadButton("download_mds1", "Download MDS Pairwise 1")),
        conditionalPanel("input.show_mds_2 > 0", plotOutput("mds2"),     downloadButton("download_mds2", "Download MDS Pairwise 2"))
      )
    )
  )
}

quantitative_analysis_server <- function(input, output, session, shared_data) {
  
  qa_result <- reactiveVal(NULL)
  
  observe({
    req(shared_data$group_names)
    g <- shared_data$group_names
    if (length(g) >= 3) {
      updateActionButton(session, "show_volcano_1", label = paste("Show Volcano:", g[1], "vs", g[2]))
      updateActionButton(session, "show_volcano_2", label = paste("Show Volcano:", g[1], "vs", g[3]))
      updateActionButton(session, "show_volcano_3", label = paste("Show Volcano:", g[2], "vs", g[3]))
      updateActionButton(session, "show_mds_1",     label = paste("Show Pairwise MDS:", g[1], "vs", g[2]))
      updateActionButton(session, "show_mds_2",     label = paste("Show Pairwise MDS:", g[1], "vs", g[3]))
    }
  })
  
  observeEvent(input$run_analysis, {
    req(shared_data$preproc_data, shared_data$group_names)
    
    data_norm   <- shared_data$preproc_data$data_norm
    gene_column <- shared_data$gene_column
    if (is.null(gene_column) || !(gene_column %in% colnames(data_norm))) gene_column <- 1
    
    withProgress(message = "Running Quantitative Analysis", value = 0, {
      setProgress(0.35, detail = "Computing MANOVA, volcano, MDS ...")
      tryCatch({
        result <- PiProteline::quantitative_analysis(
          dataset             = data_norm,                 # correct arg name
          names_of_groups     = shared_data$group_names,
          gene_column         = gene_column,               # name or 1
          significance_manova = input$significance_manova
        )
        
        if (is.null(result$manova_pairw_results)) {
          if (!is.null(result$manova_pairw))            result$manova_pairw_results <- result$manova_pairw
          if (!is.null(result$manova_pairwise_results)) result$manova_pairw_results <- result$manova_pairwise_results
        }
        
        qa_result(result)
        shared_data$quant_results <- result   
        
        setProgress(1)
        showNotification("Quantitative analysis completed.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in quantitative analysis:", e$message), type = "error")
        qa_result(NULL)
      })
    })
  })
  
  output$qa_debug <- renderPrint({
    if (is.null(qa_result())) return("No analysis yet.")
    list(
      result_elements = names(qa_result()),
      volcano_keys    = tryCatch(names(qa_result()$volcano_plots), error = function(...) NULL),
      pairwise_keys   = tryCatch(names(qa_result()$mds_plot_pairw), error = function(...) NULL)
    )
  })
  
  output$manova_table <- renderDT({ req(qa_result()); qa_result()$manova_results })
  
  output$volcano1 <- renderPlot({ req(qa_result()); qa_result()$volcano_plots[[1]] })
  output$volcano2 <- renderPlot({ req(qa_result()); qa_result()$volcano_plots[[2]] })
  output$volcano3 <- renderPlot({ req(qa_result()); qa_result()$volcano_plots[[3]] })
  
  output$mds_plot <- renderPlot({ req(qa_result()); qa_result()$mds_plot })
  output$mds1     <- renderPlot({ req(qa_result()); qa_result()$mds_plot_pairw[[1]] })
  output$mds2     <- renderPlot({ req(qa_result()); qa_result()$mds_plot_pairw[[2]] })
  
  # Downloads
  output$download_manova <- downloadHandler("manova_results.csv", function(file) write.csv(qa_result()$manova_results, file, row.names = FALSE))
  output$download_volcano1 <- downloadHandler("volcano1.png", function(file) { png(file); print(qa_result()$volcano_plots[[1]]); dev.off() })
  output$download_volcano2 <- downloadHandler("volcano2.png", function(file) { png(file); print(qa_result()$volcano_plots[[2]]); dev.off() })
  output$download_volcano3 <- downloadHandler("volcano3.png", function(file) { png(file); print(qa_result()$volcano_plots[[3]]); dev.off() })
  output$download_mds      <- downloadHandler("mds_plot.png", function(file) { png(file); print(qa_result()$mds_plot); dev.off() })
  output$download_mds1     <- downloadHandler("mds_pairwise1.png", function(file) { png(file); print(qa_result()$mds_plot_pairw[[1]]); dev.off() })
  output$download_mds2     <- downloadHandler("mds_pairwise2.png", function(file) { png(file); print(qa_result()$mds_plot_pairw[[2]]); dev.off() })
}
