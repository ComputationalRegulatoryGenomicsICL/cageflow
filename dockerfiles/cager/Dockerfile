# Use the official R base image
FROM r-base:4.3.3

# Maintainer label
LABEL maintainer="nikipavel14@gmail.com"

# Set up folder structure
WORKDIR /opt/software

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

# Install R dependencies
COPY requirements.R /opt/software/requirements.R
RUN Rscript requirements.R
