FROM condaforge/miniforge3:25.11.0-1
LABEL maintainer="adam.rivers@usda.gov"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_BANNER=1 \
    CONDA_ALWAYS_YES=true



WORKDIR /app

# Install system dependencies (curl is used by HEALTHCHECK; Java for BBMerge comes from
# the bioconda bbmap package and all pip requirements ship binary wheels)
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        curl && \
    rm -rf /var/lib/apt/lists/*

# Copy requirement files
COPY conda-requirements.txt .
COPY pip-requirements.txt .

# Install conda packages
RUN mamba install -y -c conda-forge -c bioconda --file conda-requirements.txt && \
    mamba clean -afy

# Install pip packages
RUN pip install --no-cache-dir -r pip-requirements.txt

# 4. VERSION ARGS
ARG SETUPTOOLS_SCM_PRETEND_VERSION
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${SETUPTOOLS_SCM_PRETEND_VERSION}

# Copy and install the package (application files stay root-owned and read-only to the app user;
# uploads and logs are written to the system temp directory)
COPY . /app
RUN pip install --no-cache-dir . && rm -rf build

# Unprivileged user to run the app
RUN useradd --create-home --system --shell /usr/sbin/nologin tamipami

# Verify installation
RUN python -c "from tamipami._version import __version__; print(f'Version: {__version__}')"

# Streamlit environment variables
ENV STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_MAX_UPLOAD_SIZE=200 \
    STREAMLIT_CLIENT_SHOW_ERROR_DETAILS=type \
    STREAMLIT_SERVER_FILE_WATCHER_TYPE=none \
    TAMIPAMI_SOLVER_WORKERS=4

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# Run the container as the non-root user instead of root
USER tamipami

ENTRYPOINT ["streamlit", "run", "/app/tamipami/app.py"]