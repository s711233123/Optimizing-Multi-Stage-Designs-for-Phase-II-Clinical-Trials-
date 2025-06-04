# === Setup ===
library(shiny)
library(shinyWidgets)
library(shinyBS)
library(DT)
library(globpso)
library(GA)
library(DEoptim)
library(ABCoptim)
library(shinythemes)

source("util.R")
source("kStageP2A_Objective.R")
source("run_pso_stage.R")
source("run_ga_stage.R")
source("run_de_stage.R")
source("run_abc_stage.R")

# === UI ===
ui <- fluidPage(
  theme = shinytheme("united"),
  withMathJax(),
  tags$head(
    tags$style(HTML("
    
    
    .shiny-notification {
    width: 290px !important;  /* 改這裡 */
    font-size: 16px !important;
    font-weight: bold;
    padding: 15px;
  }
    .shiny-notification .progress-bar {
      height: 10px !important;
    }
    
    body { font-family: 'Segoe UI', 'Noto Sans', sans-serif; }
    .shiny-input-container { margin-bottom: 10px; }
    .btn { font-size: 15px; }
    .tab-content { padding-top: 10px; }
    .well { background-color: #fffef9; border-color: #ffc107; }
    .dataTables_wrapper .dataTables_paginate .paginate_button { padding: 2px 10px; }

    /* 主進度條（黃色） */
    .irs-bar {
      background-color: #ffc107 !important;
    }

    /* 左邊填色邊緣（同樣黃色） */
    .irs-bar-edge {
      background-color: #ffc107 !important;
    }

    /* 滑塊（圓圈）樣式 */
    .irs-handle {
      border: 3px solid #000000 !important;
      background-color: white !important;
    }

    /* 數值顯示框：單值與 range 左右值 */
    .irs-single, .irs-from, .irs-to {
      background-color: #000000 !important;
      color: white !important;
      font-weight: bold;
      border: none;
    }
    
    /* ✅ Checkbox 勾勾顏色（藍改黃） */
  input[type='checkbox']:checked {
    accent-color: #ffc107 !important;
  }

  /* ✅ Checkbox label 字體樣式 */
  .checkbox label {
    font-weight: bold;
    font-size: 15px;
    color: #333333;
  "))
  )
  ,
  tags$div(
    style = "background-color: #fff8e1; padding: 20px; text-align: center;",
    tags$h2("Efficient Swarm Intelligence Algorithm for Optimizing Multi-Stage Designs for Phase II Clinical Trials"),
    tags$p("Yu-Hung Chou"),
    tags$p("Department of Statistics, National Taipei University, Taiwan")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("⚙️ Algorithm & Model Settings", style = "margin-top:5px; margin-bottom:8px;"),
      pickerInput(
        inputId = "algo",
        label = "Choose Algorithm:",
        choices = c("PSO" = "pso", "GA" = "ga", "DE" = "de", "ABC" = "abc"),
        selected = "pso",
        options = list(
          style = "btn-warning"
        )
      ),
      uiOutput("algo_desc"),
      tags$hr(),
      checkboxGroupInput(
        inputId = "design_type",
        label = "Select Designs to Search:",
        choices = c("Optimal", "Minimax", "Admissible"),
        selected = c("Optimal", "Minimax")
      ),
      conditionalPanel(
        condition = "input.design_type.includes('Admissible')",
        sliderInput("q_adm", "Admissible Design Weight (q2):", min = 0, max = 1, value = 0.3, step = 0.05)
      ),
      tags$hr(),
      radioGroupButtons(
        inputId = "nStage",
        label = "Number of Stages (K):",
        choices = c("2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6),
        selected = 3,
        status = "warning",
        size = "sm",
        justified = TRUE
      ),
      sliderInput("nMaxRange", "Total Sample Size Range (nMax):", min = 10, max = 100, value = c(10, 70)),
      tags$hr(),
      h4("🧪 Design Parameters"),
      sliderInput("p0range", "Range of p₀:", min = 0.05, max = 0.95, value = c(0.05, 0.75), step = 0.05),
      fluidRow(
        column( 6, numericInput("pstep", "Step size for p₀:", value = 0.05, min = 0.01, step = 0.01)),
        column( 6, numericInput("p_diff", "Difference (p₁ - p₀):", value = 0.2, min = 0.05, step = 0.01))
      ),
      
      fluidRow(
        column(6, numericInput("alpha", "Type I Error (α):", value = 0.1, min = 0, max = 1, step = 0.05)),
        column(6, numericInput("beta", "Type II Error (β):", value = 0.1, min = 0, max = 1, step = 0.05))
      ),
      
      tags$hr(),
      h4("🔁 Optimization Settings"),
      sliderInput("nParticles", "Particles / Individuals:", value = 256, min = 32, max = 1024, step = 32),
      sliderInput("maxIter1", "Initial Iterations:", min = 50, max = 1000, value = 100, step = 50),
      sliderInput("maxIter2", "Refined Iterations:", min = 100, max = 10000, value = 500, step = 100),
      tags$hr(),
      
      actionButton("runBtn", "🚀 Run Simulation", class = "btn btn-warning btn-lg", width = "100%")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Optimal Design",
                 withMathJax(
                   div(
                     style = "background-color: #fff8f0; border: 1px solid #ffc107; border-radius: 10px; padding: 15px; margin-bottom: 15px;
             box-shadow: 2px 2px 6px rgba(0,0,0,0.1); font-size:16px; word-wrap:break-word; white-space:normal;",
                     HTML("
      <b>🔍 Optimal Design Objective:</b><br/>
      The goal is to minimize the expected sample size under the null hypothesis \\( H_0 \\).<br/>
      Here, \\( PET_s(p_0) \\) denotes the probability of early termination at stage \\( s \\) under \\( H_0 \\).<br/>
      <div style='overflow-x:auto;'>
        $$E(N \\mid H_0) = n_1 + \\sum_{k=2}^{K} n_k \\left(1 - \\sum_{s=1}^{k-1} PET_s(p_0)\\right)$$
      </div>
      The objective function is:<br/>
      <div style='overflow-x:auto;'>
        $$f_{opt}(\\xi) = \\min\\ E(N \\mid H_0) + M \\cdot \\max \\{ I(\\text{Type I} > \\alpha), I(\\text{Type II} > \\beta) \\}$$
      </div>
      A large penalty \\( M \\), e.g., \\(10^8\\), is used to exclude invalid designs that violate error constraints.
    ")
                   )
                 ),
                 bsCollapse(multiple = TRUE,
                            bsCollapsePanel("📋 Table", dataTableOutput("opt_table"), style = "info"),
                            bsCollapsePanel("📊 Plot", plotOutput("opt_plot", height = "300px"), style = "warning"),
                            bsCollapsePanel("⏱️ CPU Time", uiOutput("opt_time"), style = "danger")
                 )
        ),
        tabPanel("Minimax Design",
                 withMathJax(
                   div(
                     style = "background-color: #fff8f0; border: 1px solid #ffc107; border-radius: 10px; padding: 15px; margin-bottom: 15px;
             box-shadow: 2px 2px 6px rgba(0,0,0,0.1); font-size:16px; word-wrap:break-word; white-space:normal;",
                     HTML("
      <b>🔍 Minimax Design Objective:</b><br/>
      Traditional minimax design minimizes only the maximum sample size \\( N_{\\text{max}} \\), potentially at the cost of higher average size.<br/>
      This study introduces a revised objective combining both criteria.<br/>
      Let \\( PET_s(p_0) \\) be the early termination probability at stage \\( s \\) under the null.<br/>
      <div style='overflow-x:auto;'>
        $$E(N \\mid H_0) = n_1 + \\sum_{k=2}^{K} n_k \\left(1 - \\sum_{s=1}^{k-1} PET_s(p_0)\\right)$$
      </div>
      The revised objective is:<br/>
      <div style='overflow-x:auto;'>
        $$f_{\\text{minimax}}(\\xi) = \\min \\left[N_{\\text{max}} + \\frac{E(N \\mid H_0)}{N_{\\text{max}}}\\right] + M \\cdot \\max \\{ I(\\text{Type I} > \\alpha), I(\\text{Type II} > \\beta) \\}$$
      </div>
      The ratio term penalizes designs with inefficient expected sample size, promoting both compactness and feasibility.
    ")
                   )
                 ),
                 bsCollapse(multiple = TRUE,
                            bsCollapsePanel("📋 Table", dataTableOutput("min_table"), style = "info"),
                            bsCollapsePanel("📊 Plot", plotOutput("min_plot", height = "300px"), style = "warning"),
                            bsCollapsePanel("⏱️ CPU Time", uiOutput("min_time"), style = "danger")
                 )
        ),
        tabPanel("Admissible Design",
                 withMathJax(
                   div(
                     style = "background-color: #fff8f0; border: 1px solid #ffc107; border-radius: 10px; padding: 15px; margin-bottom: 15px;
             box-shadow: 2px 2px 6px rgba(0,0,0,0.1); font-size:16px; word-wrap:break-word; white-space:normal;",
                     HTML("
      <b>🔍 Admissible Design Objective:</b><br/>
      The first stage adopts the Minimax design (i.e., \\(q_1 = 1\\)) to determine the maximum sample size;<br/>
      The second stage uses a user-specified weight \\(q_2\\) to balance between Minimax and Optimal designs.<br/><br/>
      Let \\( PET_s(p_0) \\) denote the early stopping probability at stage \\( s \\) under the null hypothesis, and define:<br/>
      <div style='overflow-x:auto;'>
        $$E(N \\mid H_0) = n_1 + \\sum_{k=2}^{K} n_k \\left(1 - \\sum_{s=1}^{k-1} PET_s(p_0)\\right)$$
      </div>
      Then, the combined objective function is:<br/>
      <div style='overflow-x:auto;'>
        $$f_{\\text{adm}}(\\xi; q) = \\min \\left[
        q \\cdot \\left( \\frac{E(N \\mid H_0)}{N_{\\text{max}}} + N_{\\text{max}} \\right) +
        (1 - q) \\cdot E(N \\mid H_0) +
        M \\cdot \\max \\{ I(t1e > \\alpha), I(t2e > \\beta) \\}
        \\right]$$
      </div>
      where \\(q \\in [0, 1]\\) reflects the user's preference trade-off.<br/>
      This design is suitable for scenarios requiring a balance between resource constraints and average efficiency.
    ")
                   )
                 )
                 
                 ,
                 bsCollapse(multiple = TRUE,
                            bsCollapsePanel("📋 Table", dataTableOutput("adm_table"), style = "info"),
                            bsCollapsePanel("📊 Plot", plotOutput("adm_plot", height = "300px"), style = "warning"),
                            bsCollapsePanel("⏱️ CPU Time", uiOutput("adm_time"), style = "danger")
                 )
        )
      )
    )
  ),
  tags$div(
    style = "background-color: #fff8e1; padding: 15px 20px; text-align: right; font-size: 14px;",
    HTML('<span style="font-size: 10px;">&#9711;</span> Maintainer: Yu-Hung Chou (<a href="mailto:jayden000819@gmail.com" style="color:#e74c3c;">jayden000819@gmail.com</a>)')
  )
)


# === SERVER ===
server <- function(input, output, session) {
  
  output$algo_desc <- renderUI({
    desc <- switch(input$algo,
                   "pso" = withMathJax(HTML('
  <div style="font-size:18px; word-wrap:break-word; white-space:normal;">
    <b>🐦 PSO (Particle Swarm Optimization)</b><br/>
    Each particle updates its velocity and position based on cognitive and social components:<br/>
    <div style="overflow-x:auto;">
      $$v_i(t+1) = \\omega v_i(t) + c_1 r_1 (p_i - x_i) + c_2 r_2 (g - x_i)$$
    </div>
    Suitable for optimization in continuous search spaces.
  </div>
')),
                   "ga" = withMathJax(HTML('
      <div style="font-size:18px; word-wrap:break-word; white-space:normal;">
        <b>🧬 GA (Genetic Algorithm)</b><br/>
        New individuals are generated through crossover and mutation:<br/>
        <div style="overflow-x:auto;">
        $$\\text{offspring} = x^{(1)}_{1:k} \\cup x^{(2)}_{k+1:n}$$
        </div>
        Capable of handling complex structures or discrete variables.
      </div>
    ')),
                   "de" = withMathJax(HTML('
      <div style="font-size:18px; word-wrap:break-word; white-space:normal;">
        <b>🌱 DE (Differential Evolution)</b><br/>
        Generates variation by perturbing solutions with vector differences:<br/>
        <div style="overflow-x:auto;">
        $$v_i = x_{r_1} + F \\cdot (x_{r_2} - x_{r_3})$$
        </div>
        Balances exploration and exploitation in real-valued optimization.
      </div>
    ')),
                   "abc" = withMathJax(HTML('
      <div style="font-size:18px; word-wrap:break-word; white-space:normal;">
        <b>🐝 ABC (Artificial Bee Colony)</b><br/>
        Mimics the foraging behavior of honey bees:<br/>
        <div style="overflow-x:auto;">
        $$v_{i,j} = x_{i,j} + \\phi(x_{i,j} - x_{k,j})$$
        </div>
        Utilizes employed, onlooker, and scout bees to search efficiently in complex spaces.
      </div>
    '))
    )
  })
  
  
  
  
  
  
  observeEvent(input$runBtn, {
    nStage <- as.numeric(input$nStage)
    nMaxRange <- input$nMaxRange
    seed <- 1
    q2 <- input$q_adm
    maxIter1 <- input$maxIter1
    maxIter2 <- input$maxIter2
    nParticles <- input$nParticles
    alpha <- input$alpha
    beta <- input$beta
    
    
    
    # === p0, p1 組合自定義 ===
    p0_seq <- seq(input$p0range[1], input$p0range[2], by = input$pstep)
    p1_seq <- p0_seq + input$p_diff
    p_combinations <- data.frame(p0 = p0_seq, p1 = p1_seq)
    
    # 過濾不合法組合 (避免 p1 > 1)
    p_combinations <- subset(p_combinations, p1 <= 1)
    
    
    
    results <- list()
    
    withProgress(message = "模擬中...", value = 0, {
      total_progress_steps <- nrow(p_combinations) * length(input$design_type)
      progress_count <- 0
      
      for (i in 1:nrow(p_combinations)) {
        p0 <- p_combinations$p0[i]
        p1 <- p_combinations$p1[i]
        cliRequirement <- list(p0 = p0, p1 = p1, alpha = alpha, beta = beta)
        
        if (input$algo == "pso") {
          if ("Optimal" %in% input$design_type) {
            start_opt <- Sys.time()
            opt_design <- run_pso_two_stage(seed, 1, 0, cliRequirement, nStage, maxIter1, maxIter2, nSwarm = nParticles, nMaxRange, FALSE)
            end_opt <- Sys.time()
            opt_cputime <- round(as.numeric(difftime(end_opt, start_opt, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("PSO - Optimal: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            opt_design <- NULL
            opt_cputime <- NA
          }
          
          if ("Minimax" %in% input$design_type) {
            start_min <- Sys.time()
            min_design <- run_pso_two_stage(seed, 1, 1, cliRequirement, nStage, maxIter1, maxIter2, nSwarm = nParticles, nMaxRange, TRUE)
            end_min <- Sys.time()
            min_cputime <- round(as.numeric(difftime(end_min, start_min, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("PSO - Minimax: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            min_design <- NULL
            min_cputime <- NA
          }
          
          if ("Admissible" %in% input$design_type) {
            start_adm <- Sys.time()
            adm_design <- run_pso_two_stage(seed, 1, q2, cliRequirement, nStage, maxIter1, maxIter2, nSwarm = nParticles, nMaxRange, FALSE)
            end_adm <- Sys.time()
            adm_cputime <- round(as.numeric(difftime(end_adm, start_adm, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("PSO - Admissible: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            adm_design <- NULL
            adm_cputime <- NA
          }
        }
        
        if (input$algo == "ga") {
          if ("Optimal" %in% input$design_type) {
            start_opt <- Sys.time()
            opt_design <- run_ga_two_stage(seed, 1, 0, cliRequirement, nStage, maxIter1, maxIter2, popSize = nParticles, nMaxRange, FALSE)
            end_opt <- Sys.time()
            opt_cputime <- round(as.numeric(difftime(end_opt, start_opt, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("GA - Optimal: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            opt_design <- NULL
            opt_cputime <- NA
          }
          
          if ("Minimax" %in% input$design_type) {
            start_min <- Sys.time()
            min_design <- run_ga_two_stage(seed, 1, 1, cliRequirement, nStage, maxIter1, maxIter2, popSize = nParticles, nMaxRange, TRUE)
            end_min <- Sys.time()
            min_cputime <- round(as.numeric(difftime(end_min, start_min, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("GA - Minimax: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            min_design <- NULL
            min_cputime <- NA
          }
          
          if ("Admissible" %in% input$design_type) {
            start_adm <- Sys.time()
            adm_design <- run_ga_two_stage(seed, 1, q2, cliRequirement, nStage, maxIter1, maxIter2, popSize = nParticles, nMaxRange, FALSE)
            end_adm <- Sys.time()
            adm_cputime <- round(as.numeric(difftime(end_adm, start_adm, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("GA - Admissible: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            adm_design <- NULL
            adm_cputime <- NA
          }
        }
        
        if (input$algo == "de") {
          if ("Optimal" %in% input$design_type) {
            start_opt <- Sys.time()
            opt_design <- run_de_two_stage(seed, 1, 0, cliRequirement, nStage, maxIter1, maxIter2, NP = nParticles, nMaxRange, FALSE)
            end_opt <- Sys.time()
            opt_cputime <- round(as.numeric(difftime(end_opt, start_opt, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("DE - Optimal: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            opt_design <- NULL
            opt_cputime <- NA
          }
          
          if ("Minimax" %in% input$design_type) {
            start_min <- Sys.time()
            min_design <- run_de_two_stage(seed, 1, 1, cliRequirement, nStage, maxIter1, maxIter2, NP = nParticles, nMaxRange, TRUE)
            end_min <- Sys.time()
            min_cputime <- round(as.numeric(difftime(end_min, start_min, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("DE - Minimax: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            min_design <- NULL
            min_cputime <- NA
          }
          
          if ("Admissible" %in% input$design_type) {
            start_adm <- Sys.time()
            adm_design <- run_de_two_stage(seed, 1, q2, cliRequirement, nStage, maxIter1, maxIter2, NP = nParticles, nMaxRange, FALSE)
            end_adm <- Sys.time()
            adm_cputime <- round(as.numeric(difftime(end_adm, start_adm, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("DE - Admissible: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            adm_design <- NULL
            adm_cputime <- NA
          }
        }
        
        if (input$algo == "abc") {
          if ("Optimal" %in% input$design_type) {
            start_opt <- Sys.time()
            opt_design <- run_abc_two_stage(seed, 1, 0, cliRequirement, nStage, maxIter1, maxIter2, foodNumber = nParticles, nMaxRange, FALSE)
            end_opt <- Sys.time()
            opt_cputime <- round(as.numeric(difftime(end_opt, start_opt, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("ABC - Optimal: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            opt_design <- NULL
            opt_cputime <- NA
          }
          
          if ("Minimax" %in% input$design_type) {
            start_min <- Sys.time()
            min_design <- run_abc_two_stage(seed, 1, 1, cliRequirement, nStage, maxIter1, maxIter2, foodNumber = nParticles, nMaxRange, TRUE)
            end_min <- Sys.time()
            min_cputime <- round(as.numeric(difftime(end_min, start_min, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("ABC - Minimax: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            min_design <- NULL
            min_cputime <- NA
          }
          
          if ("Admissible" %in% input$design_type) {
            start_adm <- Sys.time()
            adm_design <- run_abc_two_stage(seed, 1, q2, cliRequirement, nStage, maxIter1, maxIter2, foodNumber = nParticles, nMaxRange, FALSE)
            end_adm <- Sys.time()
            adm_cputime <- round(as.numeric(difftime(end_adm, start_adm, units = "secs")), 3)
            progress_count <- progress_count + 1
            setProgress(value = progress_count / total_progress_steps,
                        message = sprintf("ABC - Admissible: p₀=%.2f, p₁=%.2f", p0, p1))
          } else {
            adm_design <- NULL
            adm_cputime <- NA
          }
        }
        
        results[[i]] <- list(
          p0 = p0, p1 = p1,
          optimal = opt_design, minimax = min_design, admissible = adm_design,
          time_opt = opt_cputime, time_min = min_cputime, time_adm = adm_cputime
        )
      }
    })
    
    opt_list <- lapply(results, function(res) {
      if (is.null(res) || is.null(res$optimal)) return(NULL)
      row <- data.frame(p0 = round(res$p0, 2), p1 = round(res$p1, 2))
      for (k in 1:nStage) {
        row[[paste0("stage", k)]] <- paste0(res$optimal$rseq[k], "/", res$optimal$nseq[k])
      }
      row$EN_opt <- formatC(res$optimal$en, format = "f", digits = 2)
      for (k in 1:(nStage - 1)) {
        row[[paste0("pet", k, "_opt")]] <- formatC(res$optimal$pet_seq[k], format = "f", digits = 2)
      }
      return(row)
    })
    opt_df <- do.call(rbind, Filter(Negate(is.null), opt_list))
    
    
    min_df <- do.call(rbind, lapply(Filter(function(res) {
      !is.null(res) && !is.null(res$minimax)
    }, results), function(res) {
      row <- data.frame(p0 = round(res$p0, 2), p1 = round(res$p1, 2))
      for (k in 1:nStage) {
        row[[paste0("stage", k)]] <- paste0(res$minimax$rseq[k], "/", res$minimax$nseq[k])
      }
      row$EN_minmax <- formatC(res$minimax$en, format = "f", digits = 2)
      for (k in 1:(nStage - 1)) {
        row[[paste0("pet", k, "_min")]] <- formatC(res$minimax$pet_seq[k], format = "f", digits = 2)
      }
      row
    }))
    
    adm_df <- do.call(rbind, lapply(Filter(function(res) {
      !is.null(res) && !is.null(res$admissible)
    }, results), function(res) {
      row <- data.frame(p0 = round(res$p0, 2), p1 = round(res$p1, 2))
      for (k in 1:nStage) {
        row[[paste0("stage", k)]] <- paste0(res$admissible$rseq[k], "/", res$admissible$nseq[k])
      }
      row$EN_adm <- formatC(res$admissible$en, format = "f", digits = 2)
      for (k in 1:(nStage - 1)) {
        row[[paste0("pet", k, "_adm")]] <- formatC(res$admissible$pet_seq[k], format = "f", digits = 2)
      }
      row
    }))
    
    output$opt_table <- renderDataTable({
      datatable(opt_df, options = list(
        pageLength = 5,
        lengthMenu = c(3, 5, 10, 20),
        scrollX = TRUE   # 👉 這行讓它可以橫向捲動
      ))
    })
    
    output$opt_time <- renderUI({
      total_time <- sum(sapply(results, function(r) r$time_opt), na.rm = TRUE)
      div(style = "margin-bottom: 10px; font-weight: bold;",
          sprintf("⏱️ Total CPU Time (Optimal): %.2f sec", total_time))
    })
    
    output$opt_plot <- renderPlot({
      if (length(results) == 0) return()
      labels <- sapply(results, function(r) paste0("(", r$p0, ",", r$p1, ")"))
      cputimes <- sapply(results, function(r) r$time_opt)
      
      
      barplot(
        height = cputimes,
        names.arg = labels,
        col = "#f4a261", # 橘色
        border = NA,
        las = 1,
        ylab = "Time (sec)",
        main = "CPU Time for Optimal Design",
        ylim = c(0, max(cputimes, na.rm = TRUE) * 1.2),
        cex.names = 0.8
      )
      
    })
    
    output$min_table <- renderDataTable({
      datatable(min_df, options = list(
        pageLength = 5,
        lengthMenu = c(3, 5, 10, 20),
        scrollX = TRUE   # 👉 這行讓它可以橫向捲動
      ))
    })
    
    output$min_time <- renderUI({
      total_time <- sum(sapply(results, function(r) r$time_min), na.rm = TRUE)
      div(style = "margin-bottom: 10px; font-weight: bold;",
          sprintf("⏱️ Total CPU Time (Minimax): %.2f sec", total_time))
    })
    
    output$min_plot <- renderPlot({
      if (length(results) == 0) return()
      labels <- sapply(results, function(r) paste0("(", r$p0, ",", r$p1, ")"))
      cputimes <- sapply(results, function(r) r$time_min)
      
      barplot(
        height = cputimes,
        names.arg = labels,
        col = "#f4a261", # 橘色
        border = NA,
        las = 1,
        ylab = "Time (sec)",
        main = "CPU Time for Minimax Design",
        ylim = c(0, max(cputimes, na.rm = TRUE) * 1.2),
        cex.names = 0.8
      )
    })
    
    output$adm_table <- renderDataTable({
      datatable(adm_df, options = list(
        pageLength = 5,
        lengthMenu = c(3, 5, 10, 20),
        scrollX = TRUE   # 👉 這行讓它可以橫向捲動
      ))
    })
    
    output$adm_time <- renderUI({
      total_time <- sum(sapply(results, function(r) r$time_adm), na.rm = TRUE)
      div(style = "margin-bottom: 10px; font-weight: bold;",
          sprintf("⏱️ Total CPU Time (Admissible): %.2f sec", total_time))
    })
    
    
    output$adm_plot <- renderPlot({
      if (length(results) == 0) return()
      labels <- sapply(results, function(r) paste0("(", r$p0, ",", r$p1, ")"))
      cputimes <- sapply(results, function(r) r$time_adm)
      barplot(
        height = cputimes,
        names.arg = labels,
        col = "#f4a261", # 橘色
        border = NA,
        las = 1,
        ylab = "Time (sec)",
        main = "CPU Time for Admissible Design",
        ylim = c(0, max(cputimes, na.rm = TRUE) * 1.2),
        cex.names = 0.8
      )
    })
    
    total_time_opt <- sum(sapply(results, function(r) r$time_opt), na.rm = TRUE)
    total_time_min <- sum(sapply(results, function(r) r$time_min), na.rm = TRUE)
    total_time_adm <- sum(sapply(results, function(r) r$time_adm), na.rm = TRUE)
    
    showNotification("Simulation completed!", type = "message")
    
  })
}

# 啟動 APP
shinyApp(ui = ui, server = server)


