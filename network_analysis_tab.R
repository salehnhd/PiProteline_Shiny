# =============================
# Network Analysis Tab (runs pipeline)
# =============================
library(shiny)
library(shinyWidgets)
library(CoPPIs)
library(igraph)
library(dplyr)
library(PiProteline)
library(DT)

network_analysis_ui <- function() {
  tabPanel(
    "Network Analysis",
    sidebarLayout(
      sidebarPanel(
        
        textInput("na_save_prefix", "Save results as (prefix)", value = "Parkinson_sampled_test"),
        actionBttn("na_run_pipeline", "Run pipeline", style = "fill", color = "primary"),
        helpText("Run once after preprocessing. Functional Analysis tab will use these results."),
        tags$hr(),
        
        
        h4("Centralities"),
        selectInput("na_cent_group", "Group", choices = c("PSM_C","PSM_PD","PSM_PG")),
        actionButton("na_show_unweighted_cent", "Unweighted_centralities"),
        actionButton("na_show_weighted_cent",   "Weighted_centralities"),
        tags$hr(),
        
        
        h4("Critical Nodes"),
        selectInput("na_cn_scope", "Scope", choices = c("NotSpecific","Specific","CentralitySpecific")),
        selectInput("na_cn_group", "Group", choices = c("PSM_C","PSM_PD","PSM_PG")),
        actionButton("na_show_unweighted_cn", "Unweighted_criticalNodes"),
        actionButton("na_show_weighted_cn",   "Weighted_criticalNodes"),
        tags$hr(),
        
        
        h4("PPI graphs"),
        selectInput("na_ppi_group", "Group", choices = c("PSM_C","PSM_PD","PSM_PG")),
        actionButton("na_plot_ppi_unweighted",   "Show PPI_unweighted"),
        actionButton("na_plot_ppi_correlations","Show PPI_correlations")
      ),
      mainPanel(
        verbatimTextOutput("na_debug"),
        
        h4(textOutput("na_table_title")),
        DTOutput("na_table"),
        downloadButton("na_table_dl", "Download CSV"),
        tags$hr(),
        
        h4(textOutput("na_plot_title")),
        plotOutput("na_plot", height = 500),
        downloadButton("na_plot_dl", "Download Plot")
      )
    )
  )
}

network_analysis_server <- function(input, output, session, shared_data) {
  
  
  get_path <- function(x, path) {
    out <- x
    for (nm in path) {
      if (is.null(out) || is.null(out[[nm]])) return(NULL)
      out <- out[[nm]]
    }
    out
  }
  
  
  observeEvent(input$na_run_pipeline, {
    req(shared_data$dataset, shared_data$group_names)
    
    withProgress(message = "Running pipeline", value = 0, {
      setProgress(0.25, detail = "Building interactome...")
      g_interactome <- CoPPIs::interactome.hs %>%
        CoPPIs::filter_interactome(scores_threshold = c(experimental = 150, database = 300)) %>%
        dplyr::select(3, 4) %>%
        igraph::graph_from_data_frame(directed = FALSE)
      
      setProgress(0.55, detail = "Calling PiProteline::pipeline()...")
      set.seed(123)
      res <- PiProteline::pipeline(
        dataset         = shared_data$dataset,         
        names_of_groups = shared_data$group_names,     
        g_interactome   = g_interactome,
        gene_column     = 1,
        save_results_as = input$na_save_prefix
      )
      
      
      shared_data$pipelineResults <- res
      
      setProgress(1, detail = "Done.")
      showNotification("Pipeline finished. Functional Analysis tab is now ready.", type = "message")
    })
  })
  
 
  current_table <- reactiveVal(NULL)
  current_table_name <- reactiveVal("")
  
  observeEvent(input$na_show_unweighted_cent, {
    req(shared_data$pipelineResults)
    grp <- req(input$na_cent_group)
    tab <- get_path(shared_data$pipelineResults, c("networkAnalysis","Unweighted_centralities", grp))
    validate(need(!is.null(tab), "Unweighted_centralities table not found for this group."))
    current_table(as.data.frame(tab))
    current_table_name(paste("Unweighted_centralities -", grp))
  })
  
  observeEvent(input$na_show_weighted_cent, {
    req(shared_data$pipelineResults)
    grp <- req(input$na_cent_group)
    tab <- get_path(shared_data$pipelineResults, c("networkAnalysis","Weighted_centralities", grp))
    validate(need(!is.null(tab), "Weighted_centralities table not found for this group."))
    current_table(as.data.frame(tab))
    current_table_name(paste("Weighted_centralities -", grp))
  })
  
  output$na_table_title <- renderText({ req(current_table_name()); current_table_name() })
  output$na_table <- renderDT({ req(current_table()); current_table() },
                              options = list(pageLength = 10, scrollX = TRUE))
  output$na_table_dl <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_]+","_", current_table_name()), ".csv"),
    content  = function(file) write.csv(current_table(), file, row.names = FALSE)
  )
  
  
  show_cn <- function(weighted = FALSE) {
    req(shared_data$pipelineResults)
    scope <- req(input$na_cn_scope)   
    grp   <- req(input$na_cn_group)   
    path <- if (!weighted) c("networkAnalysis","Unweighted_criticalNodes", scope, grp)
    else            c("networkAnalysis","Weighted_criticalNodes",   scope, grp)
    obj <- get_path(shared_data$pipelineResults, path)
    validate(need(!is.null(obj), "Critical nodes not found for this selection."))
    rn <- rownames(obj)
    validate(need(length(rn) > 0, "No critical nodes available."))
    
    df <- data.frame(Node = rn, stringsAsFactors = FALSE)
    current_table(df)
    current_table_name(paste(ifelse(weighted,"Weighted","Unweighted"), "criticalNodes -", scope, "-", grp))
  }
  
  observeEvent(input$na_show_unweighted_cn, { show_cn(FALSE) })
  observeEvent(input$na_show_weighted_cn,   { show_cn(TRUE)  })
  
  
  output$na_plot_title <- renderText("")
  output$na_plot <- renderPlot({})
  
  plot_ppi <- function(kind = c("PPI_unweighted","PPI_correlations")) {
    req(shared_data$pipelineResults)
    kind <- match.arg(kind)
    grp  <- req(input$na_ppi_group)
    gobj <- get_path(shared_data$pipelineResults, c("networkAnalysis", kind, grp))
    validate(need(!is.null(gobj), paste("Graph not found for", kind, "—", grp)))
    output$na_plot_title <- renderText(paste(kind, "-", grp))
    output$na_plot <- renderPlot({ plot(gobj) })
    output$na_plot_dl <- downloadHandler(
      filename = function() paste0(gsub("[^A-Za-z0-9_]+","_", paste(kind, grp)), ".png"),
      content  = function(file) { png(file, 1400, 900, res = 150); plot(gobj); dev.off() }
    )
  }
  
  observeEvent(input$na_plot_ppi_unweighted,   { plot_ppi("PPI_unweighted")   })
  observeEvent(input$na_plot_ppi_correlations, { plot_ppi("PPI_correlations") })
  
  
  output$na_debug <- renderPrint({
    list(
      have_pipeline = !is.null(shared_data$pipelineResults),
      cent_group = input$na_cent_group %||% NA,
      cn_scope   = input$na_cn_scope %||% NA,
      cn_group   = input$na_cn_group %||% NA,
      ppi_group  = input$na_ppi_group %||% NA
    )
  })
}
