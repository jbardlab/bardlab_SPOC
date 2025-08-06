FROM condaforge/miniforge3:24.11.3-0

ENV CONDA_DIR=/opt/conda
# Use bash as a login shell so that conda is initialized
SHELL ["/bin/bash", "--login", "-c"]

# Copy environment file
COPY environment.yml /tmp/environment.yml

RUN apt-get -y update && apt-get install -y libgomp1 && \
    conda env create -f /tmp/environment.yml && \
    conda clean --tarballs --index-cache --packages --yes && \
    conda clean --force-pkgs-dirs --all --yes && \
    rm -rf /var/lib/apt/lists/*


# Copy the entrypoint script into the container
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]