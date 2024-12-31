FROM condaforge/mambaforge:24.9.2-0

WORKDIR /app

# # Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential openjdk-8-jre-headless wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Install conda dependencies
RUN mamba install -y -c bioconda bbmap

# Copy the itsxpress package files and install dependencies
LABEL org.opencontainers.image.source="https://github.com/usda-ars-gbru/tamipami"
# COPY . /app
RUN cd /app && git clone https://github.com/USDA-ARS-GBRU/tamipami.git && cd tamipami 
RUN pip install --no-cache-dir -r requirements.txt && pip install .

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

# Download and expand bbtools
#RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.13.tar.gz && \
#    tar -xzf BBMap_39.13.tar.gz && \
#    rm BBMap_39.13.tar.gz

#ENV PATH="${PATH}:/bbmap" \
#LC_ALL=C

    # Set the default command to run itsxpress
ENTRYPOINT ["streamlit", "run", "/app/tamipami/app.py", "--server.port=8501", "--server.address=0.0.0.0"]