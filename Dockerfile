FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential=12.9ubuntu3 \
    cmake=3.22.1-1ubuntu1.22.04.2 \
    python3=3.10.6-1~22.04.1 \
    python3-pip=22.0.2+dfsg-1ubuntu0.7 \
    && rm -rf /var/lib/apt/lists/*

# Pin pip itself
RUN python3 -m pip install --no-cache-dir --upgrade "pip==24.2"

# Pin Python packages
RUN python3 -m pip install --no-cache-dir \
    "numpy==2.2.6" \
    "matplotlib==3.10.8" \
    "scipy==1.15.3"

# Set working directory
WORKDIR /app

# Copy files
COPY data /app/data
COPY openfhe-development-1.1.2 /app/openfhe-development-1.1.2
COPY src /app/src
COPY ARTIFACT-APPENDIX.md /app/ARTIFACT-APPENDIX.md
COPY benchmark-realworld.sh /app/benchmark-realworld.sh
COPY benchmark-synthetic.sh /app/benchmark-synthetic.sh
COPY CMakeLists.txt /app/CMakeLists.txt
COPY LICENSE /app/LICENSE
COPY README.md /app/README.md

# Build OpenFHE
WORKDIR /app/openfhe-development-1.1.2
RUN mkdir build && cd build && \
    cmake .. && \
    make -j && \
    make install

# Set library path
ENV LD_LIBRARY_PATH=/usr/local/lib

# Build the project
WORKDIR /app
RUN mkdir build && cd build && \
    cmake .. && \
    make -j

# Set default command
CMD ["/bin/bash"]