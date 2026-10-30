# syntax=docker/dockerfile:1

# ==============================================================================
# Stage 1: Build HyPhy 3 from source (C++20 + OpenMP)
# ==============================================================================
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src

COPY CMakeLists.txt /src/
COPY cmake/ /src/cmake/
COPY include/ /src/include/
COPY src/ /src/src/
COPY viz/ /src/viz/

RUN cmake -S /src -B /build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DHYPHY3_GIT_COMMIT=docker \
    && cmake --build /build --target \
        hyphy3 \
        hyphy3-fel \
        hyphy3-contrast-fel \
        hyphy3-meme \
        hyphy3-busted \
        hyphy3-absrel \
        hyphy3-relax \
        hyphy3-fubar \
        hyphy3-prime \
        hyphy3-gard \
        hyphy3-slac \
        hyphy3-dfel \
        hyphy3-mss \
        hyphy3-ohmic \
        hyphy3-tree \
        hyphy3-joint-align \
        hyphy3-characterize \
        hyphy3-expcm \
    && strip /build/hyphy3 /build/hyphy3-*

# ==============================================================================
# Stage 2: Minimal runtime image
# ==============================================================================
FROM ubuntu:24.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/hyphy3 /usr/local/bin/
COPY --from=builder /build/hyphy3-* /usr/local/bin/

WORKDIR /data

ENTRYPOINT ["hyphy3"]
CMD ["--help"]
