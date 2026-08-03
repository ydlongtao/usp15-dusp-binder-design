FROM mambaorg/micromamba:2.3.0

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV CONDA_OVERRIDE_CUDA=12.9

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates tini \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge \
      "openmm=8.5.2" \
      "cuda-version=12.9" \
      "pdbfixer>=1.11" \
      "mdtraj>=1.10" \
      "mdanalysis>=2.8" \
      "parmed>=4.3" \
      "ambertools>=24" \
      "numpy<3" \
      "pandas>=2.2" \
      "scipy>=1.14" \
    && micromamba clean --all --yes

ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/_entrypoint.sh"]
CMD ["python", "-c", "import openmm; print(openmm.__version__); print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"]
