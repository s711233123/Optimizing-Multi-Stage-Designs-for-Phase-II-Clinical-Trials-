
# Optimizing Multi-Stage Designs for Phase II Clinical Trials

This repository implements a flexible framework in R to generate multi-stage Phase II clinical trial designs using metaheuristic algorithms. It supports design generation based on **PSO**, **GA**, **DE**, and **ABC** algorithms, and produces **optimal**, **minimax**, and **admissible** trial configurations across different $(p_0, p_1)$ response rate scenarios.

---

## 📌 Key Features

- Supports **2 to 5 stages** for trial design.
- Implements **three design types**: Optimal, Minimax, Admissible.
- Modular structure using `util.R`, `kStageP2A_Objective.R`, and algorithm wrappers.
- Compatible with PSO (globpso), GA (GA), DE (DEoptim), and ABC (ABCoptim) packages.
- Postprocessing tools to extract EN, PET, rejection rules, and sample sizes.

---

## 📁 Structure

```
.
├── run/
│   ├── util.R
│   ├── kStageP2A_Objective.R
│   ├── run_pso_stage.R
│   ├── run_ga_stage.R
│   ├── run_de_stage.R
│   └── run_abc_stage.R
├── 1227.R            # Demonstration script
└── README.md
```

---

## ⚙️ Setup Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/s711233123/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-.git
cd Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-
```

### 2. Install Required R Packages

```r
install.packages("globpso")
install.packages("GA")
install.packages("DEoptim")
install.packages("ABCoptim")
```

### 3. Set Working Directory

```r
# For Windows:
setwd("D:/NTPU/PSO/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")

# For macOS/Linux:
setwd("~/your_project/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
```

---

## 🚀 Run the Demo Script

Edit `algo` to switch algorithms:

```r
algo <- "pso"  # options: "pso", "ga", "de", "abc"
```

Then run the main file `1227.R`. It will generate designs for multiple $(p_0, p_1)$ settings and return:

- Expected Sample Size (`en`)
- PET values (`pet_seq`)
- Stage-wise cutoffs (`rseq`) and sample sizes (`nseq`)

---

## 📊 Output Tables

Summary tables `opt_df`, `min_df`, `adm_df` are created with columns:

- `p0`, `p1`, `stage1`, `stage2`, ..., `EN`, `PET` values per stage

They can be exported as:

```r
write.csv(opt_df, "opt_results.csv", row.names = FALSE)
saveRDS(results, "all_designs.rds")
```

---

## 📚 Documentation

Full explanation is available in the Appendix of the associated thesis (see `\chapter{R Implementation for Generating Multi-Stage Trial Designs}`), covering:

- Global parameter settings (e.g., `alpha`, `beta`, `q_adm`, `nStage`)
- Modular script responsibilities
- Main loop workflow and optimization logic
- Result postprocessing and export

---

## 🧠 Remarks

- You can adjust `nStage` to generate 3-stage, 4-stage, or 5-stage designs.
- The trade-off parameter `q_adm ∈ [0,1]` balances minimax and optimal designs for admissible cases.
- Designed for reproducibility and extensibility in future adaptive trial frameworks.

---

## 📬 Contact

Maintainer: 周昱宏  
Institution: National Taipei University – Graduate Institute of Statistics  
GitHub: [@s711233123](https://github.com/s711233123)



---

## 🌐 Interactive App (Shiny)

We provide a live demonstration Shiny App that allows users to interactively explore multi-stage design outputs generated by this framework:

🔗 **[Launch Shiny App](https://yhchou.shinyapps.io/multi-stage-design-via-metaheuristic/)**

This app allows:

- Selection of (p₀, p₁) scenarios
- Real-time display of Optimal, Minimax, and Admissible designs
- Visualization of stage-wise cutoffs, sample sizes, and PETs
- Downloadable result tables for local analysis

The app is fully built using R + Shiny, and integrates outputs from the underlying optimization scripts included in this repository.

