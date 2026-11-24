# Start from a base image that includes Miniforge
FROM condaforge/miniforge3:26.9.4-2@sha256:3a0d9b6b5d8b9059aca9b82ee7660b5f558b71fe764b94dea294549e88287bcf

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY env.yaml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f env.yaml

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
ENV CONDA_PREFIX=/opt/conda/envs/env
