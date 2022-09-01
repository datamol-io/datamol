FROM mambaorg/micromamba:focal


# Install system-wide packages with `apt`
USER root
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update --fix-missing && \
    apt install --no-install-recommends -y wget git openssh-client && \
    apt clean
USER $MAMBA_USER

# Copy env.yml file
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

# NOTE(hadim): see `build.args` in `devcontainer.json` to get `***_VERSION` variables
# in order to customize the installation of some deps such as python and rdkit.

# Install the mamba env
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

# Activate mamba for other below commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install git prompt
# NOTE(hadim): likely a bit hacky
COPY --chown=$MAMBA_USER:$MAMBA_USER .devcontainer/bashrc /tmp/bashrc
RUN cat /tmp/bashrc >> ~/.bashrc
