setwd("~/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
shiny::runApp("app.R", launch.browser = TRUE)

#架設環境
library(rsconnect)
rsconnect::setAccountInfo(name='yhchou',
                          token='59F419316A590FBA496D2B409879889F',
                          secret='+9sqQMYSNfdlmiKXmM1JYxaWQWt834ACd6zbJALi')
rsconnect::deployApp(
  appDir = "/home/yhchou/kstagemain/Admissible-main/app",
  appName = "multi-stage-design-via-metaheuristic"
)
