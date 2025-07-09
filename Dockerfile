# Start with rocker using your R version
FROM rocker/r-ver:4.4.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Base build tools
    build-essential \
    curl \
    git \
    # Python & venv tools
    python3 python3-pip python3-venv \
    # SSL and crypto
    libssl-dev libffi-dev \
    # R system libraries
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
    # Clean up apt cache to reduce image size
    && rm -rf /var/lib/apt/lists/*

# Install env managers
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')"
RUN pip install poetry==2.1.1 && poetry config virtualenvs.create false

# Set the working directory inside the container
WORKDIR /app

# Copy environment files
COPY renv.lock renv/activate.R ./
COPY pyproject.toml poetry.lock ./

# Restore the environment exactly as in the lockfile
RUN R -e "renv::restore()"
RUN poetry install --no-root

# Copy scripts and other necessary files
COPY scripts/ /app/

ENV RETICULATE_PYTHON=/usr/bin/python3

# Default to interactive R session
ENTRYPOINT ["Rscript", "main.R"]
