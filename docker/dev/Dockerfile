FROM ubuntu:focal

ENV DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-lc"]

ARG JULIA_VERSION=1.12.4
ARG TEMPO_REPO=https://github.com/AlexBatrakov/tempo-13.103_ddstg
ARG TEMPO_COMMIT=755a4502825fd49822fccba9bfded556a321b3b0
ARG TEMPO2_REPO=https://github.com/AlexBatrakov/tempo2
ARG TEMPO2_COMMIT=0dfe8e244949c1e8c6182cdbd2b6abac72d76012

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PSRHOME=/software \
    TEMPO=/software/tempo \
    TEMPO2=/software/tempo2/T2runtime \
    TEMPO_CMD=/software/tempo/bin/tempo \
    TEMPO2_CMD=/software/tempo2/T2runtime/bin/tempo2 \
    PGPLOT_DIR=/usr/lib/pgplot5 \
    PGPLOT_FONT=/usr/lib/pgplot5/grfont.dat \
    PGPLOT_INCLUDES=/usr/include \
    PGPLOT_BACKGROUND=white \
    PGPLOT_FOREGROUND=black \
    PGPLOT_DEV=/xs \
    PATH=/opt/julia/bin:/software/tempo/bin:/software/tempo2/T2runtime/bin:${PATH} \
    LD_LIBRARY_PATH=/software/tempo2/T2runtime/lib:${LD_LIBRARY_PATH} \
    C_INCLUDE_PATH=/software/tempo2/T2runtime/include:${C_INCLUDE_PATH}

WORKDIR /software

RUN set -euxo pipefail; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        autoconf \
        autogen \
        automake \
        autotools-dev \
        build-essential \
        ca-certificates \
        csh \
        gcc \
        gfortran \
        git \
        gnutls-bin \
        gsl-bin \
        libcfitsio-dev \
        libfftw3-3 \
        libfftw3-bin \
        libfftw3-dev \
        libfftw3-single3 \
        libgsl-dev \
        libltdl-dev \
        libpng-dev \
        libpnglite-dev \
        libtool \
        libx11-dev \
        libxext-dev \
        libxml2 \
        libxml2-dev \
        perl \
        pgplot5 \
        pkg-config \
        unzip \
        wget; \
    rm -rf /var/lib/apt/lists/*

RUN set -euxo pipefail; \
    arch="$(dpkg --print-architecture)"; \
    case "${arch}" in \
        arm64) julia_arch="aarch64" ;; \
        amd64) julia_arch="x86_64" ;; \
        *) echo "Unsupported architecture: ${arch}" >&2; exit 1 ;; \
    esac; \
    julia_minor="$(cut -d. -f1,2 <<<"${JULIA_VERSION}")"; \
    mkdir -p /opt/julia; \
    wget -q "https://julialang-s3.julialang.org/bin/linux/${julia_arch}/${julia_minor}/julia-${JULIA_VERSION}-linux-${julia_arch}.tar.gz" -O /tmp/julia.tar.gz; \
    tar -xzf /tmp/julia.tar.gz -C /opt/julia --strip-components=1; \
    rm -f /tmp/julia.tar.gz; \
    ln -sf /opt/julia/bin/julia /usr/local/bin/julia; \
    julia --version

RUN set -euxo pipefail; \
    printf '%s\n' \
        '#include <float.h>' \
        '#include <stdio.h>' \
        '' \
        'int main(void) {' \
        '    if (sizeof(long double) != 16) {' \
        '        fprintf(stderr, "expected sizeof(long double)=16, got %zu\\n", sizeof(long double));' \
        '        return 1;' \
        '    }' \
        '    if (LDBL_MANT_DIG != 113) {' \
        '        fprintf(stderr, "expected LDBL_MANT_DIG=113, got %d\\n", LDBL_MANT_DIG);' \
        '        return 1;' \
        '    }' \
        '    printf("long double ok: size=%zu mantissa=%d\\n", sizeof(long double), LDBL_MANT_DIG);' \
        '    return 0;' \
        '}' > /tmp/check_long_double.c; \
    gcc /tmp/check_long_double.c -o /tmp/check_long_double; \
    /tmp/check_long_double; \
    rm -f /tmp/check_long_double.c /tmp/check_long_double

RUN set -euxo pipefail; \
    git clone "${TEMPO_REPO}" "${TEMPO}"; \
    git -C "${TEMPO}" checkout "${TEMPO_COMMIT}"; \
    git -C "${TEMPO}" submodule update --init --recursive || true; \
    cp "${TEMPO}"/arm_fix/* "${TEMPO}"/src/; \
    cd "${TEMPO}"; \
    autoreconf --install; \
    ./configure --prefix="${TEMPO}"; \
    make -j"$(nproc)"; \
    make install; \
    find "${TEMPO}" -type f -name '*.zip' -print0 | while IFS= read -r -d '' archive; do \
        unzip -o -d "$(dirname "${archive}")" "${archive}"; \
    done; \
    find "${TEMPO}" -type f -name '*.zip' -delete; \
    test -x "${TEMPO_CMD}"; \
    test -d "${TEMPO}/data_ddstg"

RUN set -euxo pipefail; \
    git clone "${TEMPO2_REPO}" /software/tempo2; \
    git -C /software/tempo2 checkout "${TEMPO2_COMMIT}"; \
    git -C /software/tempo2 submodule update --init --recursive || true; \
    cd /software/tempo2; \
    perl -pi -e 's/chmod \+x/#chmod +x/' bootstrap; \
    ./bootstrap; \
    x_libdir="/usr/lib/$(gcc -print-multiarch)"; \
    ./configure --x-libraries="${x_libdir}" --enable-shared --enable-static --with-pic F77=gfortran; \
    make -j"$(nproc)"; \
    make install; \
    make plugins-install; \
    rm -rf "${TEMPO2}/data_ddstg"; \
    ln -s "${TEMPO}/data_ddstg" "${TEMPO2}/data_ddstg"; \
    test -x "${TEMPO2_CMD}"; \
    test -L "${TEMPO2}/data_ddstg"; \
    [[ "$(readlink -f "${TEMPO2}/data_ddstg")" == "$(readlink -f "${TEMPO}/data_ddstg")" ]]

RUN set -euxo pipefail; \
    julia -e 'println("Julia ready: ", VERSION)'; \
    test -x "${TEMPO_CMD}"; \
    "${TEMPO2_CMD}" -help >/dev/null 2>&1 || true

CMD ["bash"]
