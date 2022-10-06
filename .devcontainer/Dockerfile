FROM mambaorg/micromamba:focal

# Take root for the system-wide setup
USER root

# Install system-wide packages with `apt`
RUN export DEBIAN_FRONTEND=noninteractive && apt update --fix-missing && \
    apt install --no-install-recommends -y \
    ca-certificates \
    moreutils \
    bash-completion \
    less \
    watch \
    git \
    patch \
    git-lfs \
    nano \
    jq \
    gnupg2 \
    openssh-client \
    file \
    htop \
    zip \
    unzip \
    p7zip-full \
    curl \
    wget \
    lsof \
    iputils-ping \
    iproute2 \
    net-tools \
    dnsutils \
    socat \
    telnet \
    && apt clean && rm -rf /var/lib/apt/lists/*

# Install bash completions
RUN : \
    && mkdir -p /etc/bash_completion.d \
    && sh -c "curl -L https://raw.githubusercontent.com/docker/compose/1.29.2/contrib/completion/bash/docker-compose > /etc/bash_completion.d/docker-compose" \
    && sh -c "curl -L https://raw.githubusercontent.com/docker/cli/v20.10.13/contrib/completion/bash/docker > /etc/bash_completion.d/docker" \
    && sh -c "curl -L https://raw.githubusercontent.com/git/git/v2.35.1/contrib/completion/git-completion.bash > /etc/bash_completion.d/git" \
    ;

# Go back to regular user
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
