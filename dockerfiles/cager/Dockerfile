# Use the official R base image
FROM r-base:latest

# Maintainer label
LABEL maintainer="your.email@example.com"

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

RUN chmod -R 777 /usr/local/lib/R/site-library/

# RUN sudo usermod -a -G staff your_user_name
# Command to run R by default when the container starts
# CMD ["R"]ls -