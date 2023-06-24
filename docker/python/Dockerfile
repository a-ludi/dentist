FROM ubuntu:16.04
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"

RUN apt-get update && apt-get install -y binutils python3 python3-pip
RUN pip3 install pipenv

VOLUME /sources
WORKDIR /sources
