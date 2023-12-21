# Use the official R base image
FROM r-base:latest

# Maintainer label
LABEL maintainer="nikipavel14@gmail.com"

# Update and install system libraries
RUN apt-get update && \
    apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    procps \
    g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Bioconductor and the desired packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('CAGEr', 'BSgenome'))"