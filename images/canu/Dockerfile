# Use Miniforge as the base image
FROM condaforge/miniforge3:26.9.4-2@sha256:3a0d9b6b5d8b9059aca9b82ee7660b5f558b71fe764b94dea294549e88287bcf

# Set environment variables to avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    g++ \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    default-jdk \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set Canu version
ENV CANU_VERSION=2.2

# Download and install Canu
RUN wget https://github.com/marbl/canu/releases/download/v$CANU_VERSION/canu-$CANU_VERSION.Linux-amd64.tar.xz && \
    tar -xvf canu-$CANU_VERSION.Linux-amd64.tar.xz && \
    mv canu-$CANU_VERSION /usr/local/canu && \
    chmod +x /usr/local/canu/bin/canu && \
    ln -s /usr/local/canu/bin/canu /usr/local/bin/canu

# Add Canu to PATH
ENV PATH="/usr/local/canu/bin:$PATH"

# Copy the Conda environment YAML file into the container
COPY env.yaml .

# Create the Conda environment
RUN conda env create -f env.yaml

# Set the active Conda environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
