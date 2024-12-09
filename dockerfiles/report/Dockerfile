# Use the official R base image
FROM r-base:4.4.1

# Maintainer label
LABEL maintainer="katalitf@uio.no"

# Set up folder structure
WORKDIR /opt/software

# Update and install system libraries
RUN apt-get update && \
    apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    littler \
    procps \
    g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install R dependencies
COPY requirements.R /opt/software/requirements.R
RUN Rscript requirements.R
CMD ["/bin/bash"]
