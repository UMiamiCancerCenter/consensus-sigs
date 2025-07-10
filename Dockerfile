# Start with rocker using your R version
FROM rocker/r-ver:4.4.1


# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    cmake \
    wget \
    bzip2 \
    libssl-dev \
    libffi-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libpng-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -afy

# Put conda in PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Install Python 3.11 with conda
RUN conda install -y python=3.11 && \
    conda clean -afy

# Install renv for R
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')"

# Install Poetry
RUN pip install poetry==2.1.1 && poetry config virtualenvs.create false

ENV RETICULATE_PYTHON=$CONDA_DIR/bin/python

# Set working directory
WORKDIR /app

# Copy lockfiles
COPY renv.lock renv/activate.R ./
COPY pyproject.toml poetry.lock ./

# Restore R and Python envs
RUN R -e "renv::restore()"
RUN poetry install --no-root

# Copy scripts
COPY scripts/ /app/

# Entrypoint
ENTRYPOINT ["Rscript", "main.R"]
