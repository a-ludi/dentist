FROM ubuntu:16.04
LABEL maintainer="Arne Ludwig <ludwig@mpi-cbg.de>"
SHELL ["/bin/dash", "-c"]

ENV DLANG_ROOT=/opt/dlang

# ADD ./install-dmd.sh /opt/install-dmd.sh
# RUN /opt/install-dmd.sh


# install prerequisites
RUN apt-get update
RUN apt-get install -y git curl patch xz-utils

# install dlang installer
RUN mkdir -p "$DLANG_ROOT"
RUN curl https://dlang.org/install.sh > "$DLANG_ROOT/install.sh"
RUN chmod +x "$DLANG_ROOT/install.sh"

# NOTE: older GnuPG versions do not support latest trends in
#       signature algos. This, we sacrifice verification for
#       for a working installation.
RUN sed -i -E 's/^verify\(\)\s*\{/\0\n    return/' "$DLANG_ROOT/install.sh"

# download dmd and dub
ARG DMD_VERSION
RUN "$DLANG_ROOT/install.sh" -p "$DLANG_ROOT" install "dmd-$DMD_VERSION"

ENV LIBRARY_PATH="$DLANG_ROOT/dmd-$DMD_VERSION/linux/lib64${LIBRARY_PATH:+:}${LIBRARY_PATH:-}"
ENV LD_LIBRARY_PATH="$DLANG_ROOT/dmd-$DMD_VERSION/linux/lib64${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH:-}"
ENV PATH="$DLANG_ROOT/dmd-$DMD_VERSION/linux/bin64${PATH:+:}${PATH:-}"
ENV DMD=dmd
ENV DC=dmd

# test
RUN dmd --version
RUN dub --version
RUN git --version

# clean up
RUN apt-get remove -y curl patch xz-utils
RUN apt-get clean


VOLUME /sources
WORKDIR /sources
CMD ["dub", "build"]
