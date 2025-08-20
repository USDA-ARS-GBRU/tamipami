# Use a more specific base image for better caching
FROM condaforge/mambaforge:24.9.2-0 as base
LABEL maintainer="adam.rivers@usda.gov"

# Set environment variables for better caching and performance
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_BANNER=1 \
    CONDA_ALWAYS_YES=true

WORKDIR /app

# Install system dependencies in a single layer with cleanup
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        openjdk-8-jre-headless \
        wget \
        curl \
        git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install conda dependencies early for better caching
RUN mamba install -y -c bioconda bbmap=39.28 && \
    mamba clean -afy

# Copy requirements first for better Docker layer caching
COPY requirements.txt pyproject.toml setup.py setup.cfg /app/
COPY tamipami/__init__.py /app/tamipami/__init__.py

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy source code (this layer changes most frequently, so put it last)
COPY . /app

# Install the package
RUN pip install --no-cache-dir .

# Add labels for better container management
LABEL org.opencontainers.image.source="https://github.com/usda-ars-gbru/tamipami"
LABEL org.opencontainers.image.description="Tamipami application"
LABEL org.opencontainers.image.version="${VERSION_CLEAN:-latest}"

EXPOSE 8501

# Add a comprehensive healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# Use exec form for better signal handling
ENTRYPOINT ["streamlit", "run", "/app/tamipami/app.py", "--server.port=8501", "--server.address=0.0.0.0"]