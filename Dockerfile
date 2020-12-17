# Builds gfabase in older environment to maximize portability
FROM ubuntu:16.04

ARG sqlite_version=3340000
ARG zstd_version=1.4.5
ARG target_cpu=haswell

ENV CFLAGS="-march=${target_cpu} -O3"
ENV RUSTFLAGS="-C target-cpu=${target_cpu} -C link-args=-Wl,-rpath,\$ORIGIN"

# apt
RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y \
        build-essential zip wget aria2 tabix libtest-harness-perl

RUN mkdir -p /work/gfabase

# SQLite
WORKDIR /work
RUN wget -nv https://www.sqlite.org/2020/sqlite-amalgamation-${sqlite_version}.zip \
        && unzip -o sqlite-amalgamation-${sqlite_version}.zip
WORKDIR /work/sqlite-amalgamation-${sqlite_version}
RUN gcc -shared -o libsqlite3.so.0 -fPIC -shared -Wl,-soname,libsqlite3.so.0 -g ${CFLAGS} sqlite3.c
RUN gcc -o sqlite3 -g ${CFLAGS} sqlite3.c shell.c -lpthread -ldl
RUN cp libsqlite3.so.0 /usr/local/lib && cp *.h /usr/local/include && cp sqlite3 /usr/local/bin
RUN ln -s /usr/local/lib/libsqlite3.so.0 /usr/local/lib/libsqlite3.so

# Zstandard (only needed for tests)
WORKDIR /work
RUN wget -nv -O - https://github.com/facebook/zstd/releases/download/v${zstd_version}/zstd-${zstd_version}.tar.gz | tar zx
WORKDIR /work/zstd-${zstd_version}
RUN make install -j $(nproc)

# rustup
ADD https://sh.rustup.rs /usr/local/bin/rustup-init.sh
RUN chmod +x /usr/local/bin/rustup-init.sh && rustup-init.sh -y
ENV PATH=${PATH}:/root/.cargo/bin

RUN ldconfig

# cargo build
ADD . /work/gfabase
WORKDIR /work/gfabase
RUN rm -rf target/ && cargo build --release

# some of the tests aren't convenient to run here, as they involve large downloads or python3.6+
CMD bash -c "prove -v test/{1,2}*.t"
