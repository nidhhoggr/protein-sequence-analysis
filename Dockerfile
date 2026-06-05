FROM ubuntu:22.04

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    vim \
    nano \
    perl \
    ghostscript \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda and accept license
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /root/miniconda3 && \
    rm miniconda.sh && \
    /root/miniconda3/bin/conda clean --all -y && \
    /root/miniconda3/bin/conda config --system --prepend channels conda-forge && \
    /root/miniconda3/bin/conda config --system --prepend channels bioconda && \
    /root/miniconda3/bin/conda config --system --set auto_update_conda false && \
    /root/miniconda3/bin/conda config --system --set solver libmamba

# Accept Anaconda ToS and configure channels properly
RUN /root/miniconda3/bin/conda config --system --remove-key channels && \
    /root/miniconda3/bin/conda config --system --add channels conda-forge && \
    /root/miniconda3/bin/conda config --system --add channels bioconda

# Create conda environment and install bioinformatics tools using conda directly
RUN /root/miniconda3/bin/conda create -n bio -y \
    -c conda-forge -c bioconda \
    --override-channels \
    python=3.11 \
    biopython \
    mafft \
    fasttree \
    iqtree \
    raxml \
    cd-hit \
    seqtk \
    seqmagick \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    dendropy && \
    /root/miniconda3/bin/conda clean --all -y

# Make conda environment available by default
RUN echo "conda activate bio" >> ~/.bashrc

# Set working directory
WORKDIR /data

# Create data directories
RUN mkdir -p /data/sequences /data/results /data/scripts /data/config

# Copy entrypoint script
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]
