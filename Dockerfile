# Builds gfabase and SQLite in older OS to maximize compatibility
FROM centos:7

ARG target_cpu=ivybridge
ARG zstd_version=1.4.8

ENV RUSTFLAGS="-C target-cpu=${target_cpu} -C link-args=-Wl,-rpath,\$ORIGIN"

# yum
RUN yum install -q -y gcc make unzip git wget perl-Test-Harness epel-release
RUN yum install -q -y aria2 htslib-tools

# rustup
ADD https://sh.rustup.rs /usr/local/bin/rustup-init.sh
RUN chmod +x /usr/local/bin/rustup-init.sh && rustup-init.sh -y
ENV PATH=${PATH}:/root/.cargo/bin

# add source tree
RUN mkdir -p /work/gfabase
ADD . /work/gfabase

# install prebuilt libsqlite3.so.0 from GenomicSQLite releases
WORKDIR /work
RUN wget -nv "https://github.com/mlin/GenomicSQLite/releases/download/v$(cat gfabase/Cargo.toml.in | grep genomicsqlite | sed s/genomicsqlite// | tr -d ' ="')/libsqlite3.so.0.zip"
RUN unzip libsqlite3.so.0.zip && mv libsqlite3.so.0 /usr/local/lib/
RUN ln -s libsqlite3.so.0 /usr/local/lib/libsqlite3.so
ENV LD_LIBRARY_PATH=/usr/local/lib

# cargo build
RUN gfabase/cargo clean && gfabase/cargo build --release

# Zstandard (only needed for tests)
WORKDIR /work
RUN wget -nv -O - https://github.com/facebook/zstd/releases/download/v${zstd_version}/zstd-${zstd_version}.tar.gz | tar zx
WORKDIR /work/zstd-${zstd_version}
RUN make install
        # ^ don't use `make -j $(nproc)` due to buggy old make

# Poise to run a subset of the tests (others aren't convenient to run here as they need modern python)
WORKDIR /work/gfabase
CMD bash -c "sha256sum target/release/gfabase && target/release/gfabase version && prove -v test/{1,2}*.t"
