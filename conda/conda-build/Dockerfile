FROM ubuntu:16.04
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"
SHELL ["/bin/dash", "-c"]

RUN apt-get update && apt-get install -y git curl patch \
    && curl https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh > install.sh \
    && chmod +x install.sh \
    && printf "\nyes\n/opt/miniconda3\n\n" | ./install.sh \
    && rm install.sh \
    && apt-get remove -y curl && apt-get autoremove -y && apt-get clean
ADD activate-conda.sh /app/activate-conda.sh
RUN . /app/activate-conda.sh && \
    conda update conda && \
    conda install conda-build git

VOLUME /recipe
ADD ./entrypoint.sh /app/entrypoint.sh
ENTRYPOINT ["/app/entrypoint.sh"]
