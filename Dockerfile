FROM ubuntu:16.04
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"

ARG DENTIST_VERSION
ARG PREFIX="/usr/local"
ENV BINDIR="$PREFIX/bin"
ENV PATH="$BINDIR:$PATH"
SHELL ["/bin/bash", "-c"]

# Configure tzdata (timezone)
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Fix error message about locales (at least sometimes)
ENV LANGUAGE=C \
    LANG=C \
    LC_ALL=C

# Install DENTIST using conda
RUN bail_out() { echo "ERROR" "$@"; exit 1; } >&2; \
    ( [[ -n "${DENTIST_VERSION:+set}" ]] || bail_out "missing required --build-arg DENTIST_VERSION=..." ) \
    && apt-get update && apt-get install -y git curl patch \
    && curl https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh > install.sh \
    && chmod +x install.sh \
    && printf "\nyes\n/opt/miniconda3\n\n" | ./install.sh \
    && rm install.sh \
    && ( \
        export PATH="/opt/miniconda3/bin:$PATH" \
        && eval "$(conda shell.bash hook)" \
        && conda create -y -p "$PREFIX" -c a_ludi -c bioconda dentist-core=="$DENTIST_VERSION" jq==1.6 \
    ) \
    && rm -rf /opt/miniconda3 \
    && apt-get remove -y curl && apt-get autoremove -y && apt-get clean
COPY snakemake/Snakefile ./snakemake/Snakefile
COPY tests/test-commands.sh ./tests/test-commands.sh
RUN dentist -d && ls -R snakemake tests && pwd && ./tests/test-commands.sh

# Provide alternative location for bash
RUN ln -s /bin/bash /usr/bin/bash
