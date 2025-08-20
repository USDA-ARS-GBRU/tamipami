FROM condaforge/mambaforge:24.9.2-0

LABEL maintainer="adam.rivers@usda.gov"

# Environment setup
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_BANNER=1 \
    CONDA_ALWAYS_YES=true

ARG SETUPTOOLS_SCM_PRETEND_VERSION
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${SETUPTOOLS_SCM_PRETEND_VERSION}

WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        openjdk-8-jre-headless \
        wget \
        curl \
        git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install conda packages
RUN mamba install -y -c conda-forge -c bioconda \
        bbmap=39.28 \
        biopython \
        scipy \
        pandas \
        matplotlib \
        scikit-bio \
        pyyaml \
        streamlit \
        setuptools \
        pip && \
    mamba clean -afy

# Install pip-only dependencies
RUN pip install --no-cache-dir \
        ckmeans \
        treelib \
        textdistance \
        logomaker \
        altair \
        tables

# Copy and install package
COPY . /app
RUN pip install --no-cache-dir .

# Verify installation
RUN python -c "from tamipami._version import __version__; print(f'Version: {__version__}')"

EXPOSE 8501
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1
ENTRYPOINT ["streamlit", "run", "/app/tamipami/app.py", "--server.port=8501", "--server.address=0.0.0.0"]