FROM alpine:3
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"

# Install dependencies (build & runtime) via apk
RUN echo 'http://dl-cdn.alpinelinux.org/alpine/edge/main' >> /etc/apk/repositories && \
    echo 'http://dl-cdn.alpinelinux.org/alpine/edge/community' >> /etc/apk/repositories && \
    apk update && \
    # install build dependencies
    apk add \
        build-base abuild binutils dmd dub zlib-static zlib-dev avr-libc \
        git autoconf automake libtool gmp-dev && \
    # install runtime (and build) dependencies
    apk add bash jq
# Provide our convenient build script to reduce verbosity
COPY ./build-and-install.sh /opt/

# Build runtime dependencies
RUN REPO=https://gitlab.com/german.tischler/libmaus2.git \
    BRANCH=2.0.724-release-20200702192714 \
    PREBUILD='autoupdate && autoreconf -i -f && ./configure --with-gmp && echo --- BEGIN && cat Makefile && echo --- END && make -d' \
    INSTALL_CMD='make install' \
    /opt/build-and-install.sh libmaus2 make
RUN REPO=https://gitlab.com/german.tischler/daccord.git \
    BRANCH=0.0.18-release-20200702195851 \
    PREBUILD='autoreconf -i -f && ./configure --with-libmaus2=/usr/local' \
    INSTALL_CMD='make install' \
    /opt/build-and-install.sh daccord make
RUN REPO=https://github.com/thegenemyers/DAZZ_DB.git \
    BRANCH=d22ae58d32a663d09325699f17373ccf8c6f93a0 \
    /opt/build-and-install.sh DAZZ_DB make
RUN REPO=https://github.com/thegenemyers/DASCRUBBER.git \
    BRANCH=a53dbe879a716e7b08338f397de5a0403637641e \
    /opt/build-and-install.sh DASCRUBBER make
RUN REPO=https://github.com/thegenemyers/DAMASKER.git \
    BRANCH=22139ff1c2b2c0ff2589fbc9cc948370be799827 \
    PREBUILD="sed -i -E 's/\\bDB_CSS\\b/DB_CCS/g' *.c *.h" \
    /opt/build-and-install.sh DAMASKER make
RUN REPO=https://github.com/thegenemyers/DALIGNER.git \
    BRANCH=c2b47da6b3c94ed248a6be395c5b96a4e63b3f63 \
    /opt/build-and-install.sh DALIGNER make
RUN REPO=https://github.com/thegenemyers/DAMAPPER.git \
    BRANCH=b2c9d7fd64bb4dd2dde7c69ff3cc8a04cbeeebbc \
    /opt/build-and-install.sh DAMAPPER make
COPY ./ /opt/dentist
RUN REPO=https://github.com/a-ludi/dentist.git \
    BRANCH='' \
    BUILD=release \
    /opt/build-and-install.sh dentist dub
# Check if dependencies are correctly installed and remove build dependencies
# and artifacts
RUN rm -rf /opt/build-and-install.sh \
           "$HOME/.dub" && \
    # prevent required shared objects from removal
    apk add \
        llvm-libunwind libstdc++ libgcc gmp libgomp && \
    apk del \
        build-base gcc abuild binutils dmd dub zlib-static zlib-dev avr-libc \
        git autoconf automake libtool gmp-dev && \
    # let dentist check if the dependencies are present
    dentist -d
