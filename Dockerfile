FROM ubuntu:16.04
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"

ARG NCPUS=1
ARG PREFIX="/usr/local"
ENV BINDIR="$PREFIX/bin"
ENV PATH="$BINDIR:$PATH"

# Configure tzdata (timezone)
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Fix error message about locales (at least sometimes)
ENV LANGUAGE=C \
    LANG=C \
    LC_ALL=C

# Provide our convenient build script to reduce verbosity
COPY ./build-and-install.sh /opt/

# Build runtime dependencies

# libmaus2 & daccord
RUN BUILD_DEPS="build-essential autoconf pkg-config zlib1g-dev git libtool libgmp-dev" && \
    RUNTIME_DEPS="libgmp10 libgomp1" && \
    apt-get update && apt-get install -y $BUILD_DEPS $RUNTIME_DEPS && \
    REPO=https://gitlab.com/german.tischler/libmaus2.git \
    BRANCH=2.0.724-release-20200702192714 \
    PREBUILD='autoupdate && autoreconf -i -f && ./configure --prefix="$PREFIX" --with-gmp && make' \
    INSTALL_CMD='make install' \
    CLEANBUILD=0 \
    /opt/build-and-install.sh libmaus2 make && \
    \
    REPO=https://gitlab.com/german.tischler/daccord.git \
    BRANCH=0.0.18-release-20200702195851 \
    PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:}$PKG_CONFIG_PATH" \
    MAKEFLAGS='LDFLAGS=-static-libtool-libs' \
    PREBUILD='autoreconf -i -f && ./configure --prefix="$PREFIX"' \
    INSTALL_CMD='make install' \
    /opt/build-and-install.sh daccord make && \
    make -C /opt/libmaus2 uninstall && rm -rf /opt/libmaus2 && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DAZZ_DB
RUN BUILD_DEPS="build-essential git zlib1g-dev" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/DAZZ_DB.git \
    BRANCH=d22ae58d32a663d09325699f17373ccf8c6f93a0 \
    /opt/build-and-install.sh DAZZ_DB make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DASCRUBBER
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/DASCRUBBER.git \
    BRANCH=a53dbe879a716e7b08338f397de5a0403637641e \
    /opt/build-and-install.sh DASCRUBBER make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DAMASKER
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/DAMASKER.git \
    BRANCH=22139ff1c2b2c0ff2589fbc9cc948370be799827 \
    PREBUILD='sed -i -E '\''s/\bDB_CSS\b/DB_CCS/g'\'' *.c *.h' \
    /opt/build-and-install.sh DAMASKER make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DALIGNER
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/DALIGNER.git \
    BRANCH=c2b47da6b3c94ed248a6be395c5b96a4e63b3f63 \
    /opt/build-and-install.sh DALIGNER make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DAMAPPER
RUN BUILD_DEPS="build-essential git" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    REPO=https://github.com/thegenemyers/DAMAPPER.git \
    BRANCH=b2c9d7fd64bb4dd2dde7c69ff3cc8a04cbeeebbc \
    /opt/build-and-install.sh DAMAPPER make && \
    apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean

# DENTIST
ARG DMD_VERSION="2.096.0-0"
ARG DUB_VERSION="1.23.0"
ARG DENTIST_BUILD=release
COPY ./ /opt/dentist
RUN BUILD_DEPS="build-essential git jq curl libcurl3" && \
    apt-get update && apt-get install -y $BUILD_DEPS && \
    curl -o /tmp/dmd.deb -L "https://s3.us-west-2.amazonaws.com/downloads.dlang.org/releases/2021/dmd_${DMD_VERSION}_amd64.deb" && \
    dpkg -i /tmp/dmd.deb && \
    rm /tmp/dmd.deb && \
    BUILD_DEPS="$BUILD_DEPS dmd" && \
    curl -o /tmp/dub.tar.gz -L "https://github.com/dlang/dub/releases/download/v${DUB_VERSION}/dub-v${DUB_VERSION}-linux-x86_64.tar.gz" && \
    tar -C /usr/local/bin -xzf /tmp/dub.tar.gz && \
    rm /tmp/dub.tar.gz && \
    REPO=https://github.com/a-ludi/dentist.git \
    BRANCH='' \
    BUILD=unittest \
    CLEANBUILD=0 \
    INSTALL_CMD='install dentist "$BINDIR/dentist-unittest"' \
    /opt/build-and-install.sh dentist dub && \
    dentist-unittest && \
    rm -f "$BINDIR/dentist-unittest" && \
    REPO=https://github.com/a-ludi/dentist.git \
    BRANCH='' \
    BUILD="$DENTIST_BUILD" \
    CLEANBUILD=1 \
    /opt/build-and-install.sh dentist dub && \
    dentist -d && \
    rm -f /usr/local/bin/dub && rm -rf "$HOME/.dub" && apt-get remove -y $BUILD_DEPS && apt-get autoremove -y && apt-get clean
