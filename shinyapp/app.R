# Load required libraries
library(shiny)
library(shinythemes)
library(globpso)
library(DiagrammeR)

# Load external files
source("util.R")  # Load utility functions
source("kStageP2A_Objective.R")  # Load objective functions for clinical trial design

# Define the UI
ui <- fluidPage(
  theme = shinytheme("united"),  # Use United theme
  tags$head(
    tags$style(HTML("
    .btn-orange {
      background-color: orange;
      color: white;
      border: none;
    }
    .btn-orange:hover {
      background-color: #cc8400;
      color: white;
    }
    .table { font-size: 20px !important; }
    pre { font-size: 14px !important; }
    .summary-table {
      font-size: 18px;
      margin-top: 20px;
      border-collapse: collapse;
      width: 100%;
    }
    .summary-table th, .summary-table td {
      border: 1px solid #ddd;
      text-align: left;
      padding: 8px;
    }
    .summary-table th {
      background-color: #f2f2f2;
      font-weight: bold;
    }
    .cpu-time-box {
      background-color: #f9f9f9;
      border: 1px solid #ddd;
      padding: 10px;
      margin-top: 10px;
      font-size: 20px;
      font-weight: bold;
      color: #333;
      border-radius: 5px;
    }
  "))
  ),
  
  titlePanel("PSO for K-Stage Clinical Trial Design"),
  
  sidebarLayout(
    sidebarPanel(
      # Group basic parameters into a well panel
      wellPanel(
        h4("Basic Parameters:"),
        numericInput("p0", "Baseline Response Rate (p0):", value = 0.1, min = 0, max = 1, step = 0.05),
        helpText("The baseline response rate (p0) is the expected response rate under the null hypothesis."),
        numericInput("p1", "Target Response Rate (p1):", value = 0.3, min = 0, max = 1, step = 0.05),
        helpText("The target response rate (p1) is the expected response rate under the alternative hypothesis."),
        numericInput("alpha", "Significance Level (α):", value = 0.1, min = 0, max = 1, step = 0.01),
        helpText("The significance level (α) represents the probability of a Type I error (false positive)."),
        numericInput("beta", "Type II Error (β):", value = 0.1, min = 0, max = 1, step = 0.01),
        helpText("The Type II error (β) represents the probability of failing to reject the null hypothesis when it is false.")
      ),
      
      # Use dropdown menus for interim analyses and PSO type selections
      selectInput("nStage", "Number of Interim Analyses:",
                  choices = list("2 Stages" = 2, "3 Stages" = 3, "4 Stages" = 4, "5 Stages" = 5),
                  selected = 3),
      helpText("Select the number of interim analyses to be conducted during the trial."),
      
      selectInput("psoType", "PSO Type:",
                  choices = list("Basic" = "basic", "Quantum" = "quantum", "LCRI" = "lcri", "Comp" = "comp", "Dexp" = "dexp"),
                  selected = "basic"),
      helpText("Select the type of Particle Swarm Optimization (PSO) to use for the trial design optimization."),
      
      # Sliders for swarm size and iterations
      sliderInput("nSwarm", "Swarm Size for Optimization:", min = 100, max = 3000, value = 300, step = 50),
      helpText("Set the number of particles in the swarm for the optimization process."),
      sliderInput("maxIter", "Number of Optimization Iterations:", min = 100, max = 3000, value = 300, step = 50),
      helpText("Set the maximum number of iterations for the PSO optimization process."),
      
      # Button to run PSO optimization
      actionButton("run", "Run PSO", class = "btn-orange")
    ),
    
    mainPanel(
      tabsetPanel(
        # Display optimization results table
        tabPanel("Optimization Results Table", 
                 tableOutput("resultsTable"),
                 hr(),
                 h4("Summary:"),
                 tableOutput("summaryTable"),
                 hr(),
                 h4("CPU Time:"),
                 tags$div(class = "cpu-time-box", textOutput("cpuTime"))
        ),
        # Display detailed optimization metrics
        tabPanel("Detailed Optimization Metrics", verbatimTextOutput("details")),
        # Display clinical trial flowchart
        tabPanel("Clinical Trial Flowchart", grVizOutput("flowchart", width = "90%", height = "600px"))
      )
    )
  )
)

# Define the server
server <- function(input, output, session) {
  # Reactive values for stage and PSO type
  nStage <- reactive({ input$nStage })
  psoType <- reactive({ input$psoType })
  
  # Run PSO optimization
  observeEvent(input$run, {
    cliRequirement <- list(
      p0 = input$p0,
      p1 = input$p1,
      alpha = input$alpha,
      beta = input$beta
    )
    
    # Get current stage and settings
    nStageValue <- as.numeric(nStage())
    nSwarm <- input$nSwarm
    maxIter <- input$maxIter
    n1Min <- 1
    nrMin <- 1
    nMaxRange <- c(15, 70)
    nMinEachInterim <- c(n1Min, rep(nrMin, nStageValue - 1))
    upper <- c(nMaxRange[2], rep(0.5 * pi, nStageValue - 1), rep(1, nStageValue))
    lower <- c(nMaxRange[1], rep(0.0 * pi, nStageValue - 1), rep(0, nStageValue))
    
    # PSO settings
    algSetting <- getPSOInfo(
      nSwarm = nSwarm,
      maxIter = maxIter,
      psoType = psoType()
    )
    
    # Run global PSO optimization
    optimRes <- globpso(
      objFunc = kStageOptimObj, 
      PSO_INFO = algSetting, 
      lower = lower, 
      upper = upper, 
      seed = NULL, 
      verbose = TRUE, 
      nMin = nMinEachInterim, 
      cliRequirement = cliRequirement
    )
    
    # Get optimal design
    optimDesign <- kStageFreqCrit(
      nPolarized = optimRes$par[2:nStageValue],  
      rProportion = optimRes$par[(nStageValue + 1):length(optimRes$par)], 
      nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement)
    
    # Optimization results table
    results <- data.frame(
      Stage = 1:nStageValue,
      SampleSizes = optimDesign$nseq,
      StoppingCutoffs = optimDesign$rseq
    )
    
    output$resultsTable <- renderTable({
      results
    })
    
    # Summary information table
    early_termination_labels <- paste0("Early Termination (Stage ", 1:(nStageValue - 1), ")")
    summary_data <- data.frame(
      Metric = c("Expected Total Sample Size", "Significance Level", "Statistical Power (1 - β)", early_termination_labels),
      Value = c(round(optimDesign$en, 2),
                format(optimDesign$t1e, scientific = FALSE), 
                format(1 - optimDesign$t2e, scientific = FALSE),
                round(optimDesign$pet_seq, 2))
    )
    
    output$summaryTable <- renderTable({
      summary_data
    }, bordered = TRUE, striped = TRUE, hover = TRUE)
    
    # Display CPU time
    output$cpuTime <- renderText({
      paste0(round(optimRes$cputime["elapsed"], 2), "s")
    })
    
    # Display detailed optimization metrics
    output$details <- renderPrint({
      list(
        TypeIError = optimDesign$t1e,
        TypeIIError = optimDesign$t2e,
        ExpectedSampleSize = round(optimDesign$en, 2),
        ObjectiveValue = optimRes$val,
        CPUTime = optimRes$cputime,
        EarlyTerminationProbabilities = optimDesign$pet_seq,
        SampleSizesByStage = optimDesign$nseq,
        StoppingCutoffsByStage = optimDesign$rseq
      )
    })
    
    # Render clinical trial flowchart
    output$flowchart <- renderGrViz({
      nodes <- c()
      edges <- c()
      
      for (i in 1:nStageValue) {
        # Add nodes for each stage
        nodes <- c(nodes, 
                   paste0("Stage_", i, 
                          " [label=<<table border='0' cellborder='0' cellspacing='0'>",
                          "<tr><td><b>Stage ", i, "</b></td></tr>",
                          "<tr><td><font point-size='12'>Sample: ", optimDesign$nseq[i], "</font></td></tr>",
                          "<tr><td><font point-size='12'>Cutoff: ", optimDesign$rseq[i], "</font></td></tr>",
                          "</table>>, shape=box, style=filled, fillcolor=white, color=orange, fontsize=14]"))
        
        if (i < nStageValue) {
          # Add intermediate stage decision branches
          edges <- c(edges, 
                     paste0("Stage_", i, " -> Fail_", i, 
                            " [label='≤ ", optimDesign$rseq[i], "', color=black, fontsize=12, labeldistance=2, labelangle=90]"),
                     paste0("Stage_", i, " -> Stage_", i + 1, 
                            " [label='> ", optimDesign$rseq[i], "', color=black, fontsize=12, labeldistance=2, labelangle=-90]"),
                     paste0("Fail_", i, 
                            " [label='Fail', shape=ellipse, style=filled, fillcolor=white, color=gray, fontsize=12]"))
        } else {
          # Final stage decision branches
          edges <- c(edges, 
                     paste0("Stage_", i, " -> Success", 
                            " [label='> ", optimDesign$rseq[i], "', color=black, fontsize=12]"),
                     paste0("Stage_", i, " -> Fail_Final", 
                            " [label='≤ ", optimDesign$rseq[i], "', color=black, fontsize=12]"),
                     paste0("Success [label='Effective', shape=ellipse, style=filled, fillcolor=white, color=lightblue, fontsize=12]"),
                     paste0("Fail_Final [label='Fail', shape=ellipse, style=filled, fillcolor=white, color=red, fontsize=12]"))
        }
      }
      
      # Generate flowchart
      grViz(
        paste0(
          "digraph G {",
          "graph [rankdir=TB, nodesep=2, ranksep=1, splines=line]",  # Increase node spacing
          paste(nodes, collapse = "\n"),
          "\n",
          paste0(edges, collapse = "\n"),
          "edge [arrowhead=vee, arrowtail=dot, penwidth=1.5, color=black, arrowsize=1.2]",  # Adjust arrow length and width
          "}"
        )
      )
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
