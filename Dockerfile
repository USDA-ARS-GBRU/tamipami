FROM condaforge/miniforge3:25.3.1-0

LABEL maintainer="adam.rivers@usda.gov"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_BANNER=1 \
    CONDA_ALWAYS_YES=true

ARG SETUPTOOLS_SCM_PRETEND_VERSION
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${SETUPTOOLS_SCM_PRETEND_VERSION}

WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        openjdk-8-jre-headless \
        wget \
        curl \
        git && \
    rm -rf /var/lib/apt/lists/*

# Copy requirement files
COPY conda-requirements.txt .
COPY pip-requirements.txt .

# Install conda packages removing packages after install is done to keep image in sync with Mac Arm installs where ortools requires its wn versions of these libs
RUN mamba install -y -c conda-forge -c bioconda --file conda-requirements.txt && \
    mamba remove -p /opt/conda --force -y libabseil glog gflags protobuf && \ 
    mamba clean -afy

# Install pip packages
RUN pip install --no-cache-dir -r pip-requirements.txt

# Copy and install your package
COPY . /app
RUN pip install --no-cache-dir .

# Verify installation
RUN python -c "from tamipami._version import __version__; print(f'Version: {__version__}')"

EXPOSE 8501
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1
ENTRYPOINT ["streamlit", "run", "/app/tamipami/app.py", "--server.port=8501", "--server.address=0.0.0.0"]