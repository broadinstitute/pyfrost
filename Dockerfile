FROM mambaorg/micromamba:0.25.1 AS builder
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/pyfrost/
USER root
RUN apt-get update && apt-get install -y \
    make && rm -rf /var/lib/{apt,dpkg,cache,log}

USER $MAMBA_USER
RUN micromamba install -y -f /tmp/pyfrost/env-build.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --upgrade build && cd /tmp/pyfrost/ && python3 -m build

FROM mambaorg/micromamba:0.25.1 AS install
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /tmp/pyfrost/ /tmp/pyfrost
RUN micromamba install -y -f /tmp/pyfrost/env.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install /tmp/pyfrost/dist/*.whl
