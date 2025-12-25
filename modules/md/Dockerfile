FROM ubuntu:20.04
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ENV DEBIAN_FRONTEND=noninteractive
RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|http://mirrors.aliyun.com/ubuntu/|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com/ubuntu/|http://mirrors.aliyun.com/ubuntu/|g' /etc/apt/sources.list && \
    apt-get update && apt-get install -y --no-install-recommends \
    gromacs \
    build-essential \
    wget \
    bzip2 \
    file \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY Miniconda3-latest-Linux-x86_64.sh /tmp/miniconda.sh
RUN bash /tmp/miniconda.sh -b -p /opt/conda && rm /tmp/miniconda.sh
ENV PATH=/opt/conda/envs/md/bin:/opt/conda/bin:$PATH
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda update -n base -c defaults -y conda && \
    conda create -y -n md python=3.12 && \
    conda install -y -n md -c conda-forge ambertools>=24 parmed numpy=1.26 matplotlib && \
    conda clean -afy
RUN pip install --no-cache-dir fastapi uvicorn
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
EXPOSE 8004
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8004"]
