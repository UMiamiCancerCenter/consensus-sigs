# Start with rocker using your R version
FROM rocker/r-ver:4.4.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install renv (no need to snapshot here)
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')"

# Create a working directory in the container
WORKDIR /app

# Copy renv.lock and renv/activate.R into the image
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R

# Restore the environment exactly as in the lockfile
RUN R -e "renv::restore()"

COPY scripts/ /app/

# Default to interactive R session
ENTRYPOINT ["Rscript", "deconstructsigs_assignment.R"]
