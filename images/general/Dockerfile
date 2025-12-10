# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

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
