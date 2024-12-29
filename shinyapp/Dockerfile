# 使用官方 R Shiny 基礎映像
FROM rocker/shiny:latest

# 系統更新與必要依賴安裝（globpso 和其他 R 包可能需要）
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk-dev

# 安裝必要的 R 包
RUN R -e "install.packages(c('shiny', 'shinythemes', 'DiagrammeR', 'devtools'), repos='https://cloud.r-project.org/')"

# 安裝 globpso 包
# 方法 1: 如果 globpso 存在於 GitHub 上
RUN R -e "devtools::install_github('PingYangChen/globpso')"

# 或 方法 2: 使用本地的 tar.gz 文件安裝 globpso
# COPY globpso.tar.gz /tmp/globpso.tar.gz
# RUN R -e "install.packages('/tmp/globpso.tar.gz', repos=NULL, type='source')"

# 複製應用程序文件到容器
COPY . /srv/shiny-server/

# 設定文件的權限（重要）
RUN chmod -R 755 /srv/shiny-server

# 確保應用程序入口文件（如 app.R）存在於正確路徑
WORKDIR /srv/shiny-server

# 暴露 3838 端口（Shiny 預設端口）
EXPOSE 3838

# 啟動 Shiny Server
CMD ["/usr/bin/shiny-server"]
