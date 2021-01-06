# Builds gfabase and SQLite in older OS to maximize compatibility
FROM centos:7

ARG sqlite_version=3340000
ARG zstd_version=1.4.8
ARG target_cpu=core-avx-i

ENV CFLAGS="-march=${target_cpu} -O3"
ENV RUSTFLAGS="-C target-cpu=${target_cpu} -C link-args=-Wl,-rpath,\$ORIGIN"
# https://www.sqlite.org/compile.html
ENV SQLITE_CFLAGS="\
        -DSQLITE_ENABLE_LOAD_EXTENSION \
        -DSQLITE_USE_URI \
        -DSQLITE_LIKE_DOESNT_MATCH_BLOBS \
        -DSQLITE_DEFAULT_MEMSTATUS=0 \
        -DSQLITE_MAX_EXPR_DEPTH=0 \
        -DSQLITE_ENABLE_NULL_TRIM \
        -DSQLITE_USE_ALLOCA \
        -DSQLITE_HAVE_ISNAN \
        -DSQLITE_ENABLE_UPDATE_DELETE_LIMIT \
        -DSQLITE_ENABLE_COLUMN_METADATA \
        -DSQLITE_ENABLE_DBSTAT_VTAB \
        -DSQLITE_ENABLE_FTS5 \
        -DSQLITE_ENABLE_RTREE \
        -DSQLITE_ENABLE_PREUPDATE_HOOK \
        -DSQLITE_ENABLE_SESSION \
"

# yum
RUN yum install -q -y gcc make unzip git wget perl-Test-Harness epel-release
RUN yum install -q -y aria2 htslib-tools

# rustup
ADD https://sh.rustup.rs /usr/local/bin/rustup-init.sh
RUN chmod +x /usr/local/bin/rustup-init.sh && rustup-init.sh -y
ENV PATH=${PATH}:/root/.cargo/bin

RUN mkdir -p /work/gfabase

# SQLite
WORKDIR /work
RUN wget -nv https://www.sqlite.org/2020/sqlite-amalgamation-${sqlite_version}.zip \
        && unzip -o sqlite-amalgamation-${sqlite_version}.zip
WORKDIR /work/sqlite-amalgamation-${sqlite_version}
RUN gcc -shared -o libsqlite3.so.0 -fPIC -shared -Wl,-soname,libsqlite3.so.0 -g ${CFLAGS} ${SQLITE_CFLAGS} sqlite3.c
RUN gcc -o sqlite3 -g ${CFLAGS} ${SQLITE_CFLAGS} sqlite3.c shell.c -lpthread -ldl -lm
RUN cp libsqlite3.so.0 /usr/local/lib && cp *.h /usr/local/include && cp sqlite3 /usr/local/bin
RUN ln -s /usr/local/lib/libsqlite3.so.0 /usr/local/lib/libsqlite3.so
ENV LD_LIBRARY_PATH=/usr/local/lib

# Zstandard (only needed for tests)
WORKDIR /work
RUN wget -nv -O - https://github.com/facebook/zstd/releases/download/v${zstd_version}/zstd-${zstd_version}.tar.gz | tar zx
WORKDIR /work/zstd-${zstd_version}
RUN make install
        # ^ don't use `make -j $(nproc)` due to buggy old make

# cargo build
ADD . /work/gfabase
WORKDIR /work/gfabase
RUN ./cargo clean && ./cargo build --release

# some of the tests aren't convenient to run here, as they involve large downloads or python3.6+
CMD bash -c "sha256sum target/release/gfabase && target/release/gfabase version && prove -v test/{1,2}*.t"
