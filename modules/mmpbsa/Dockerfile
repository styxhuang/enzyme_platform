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
ENV PATH=/opt/conda/envs/mmpbsa/bin:/opt/conda/bin:$PATH
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda update -n base -c defaults -y conda && \
    conda create -y -n mmpbsa python=3.11.8 && \
    conda install -y -n mmpbsa -c conda-forge ambertools>=23 parmed numpy=1.26 matplotlib mpi4py=4.0.1 gromacs && \
    conda clean -afy
RUN pip install --no-cache-dir fastapi uvicorn && \
    /opt/conda/bin/conda run -n mmpbsa pip install --no-cache-dir gmx_MMPBSA
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
ENV MMPBSA_CONDA_ENV=mmpbsa
EXPOSE 8008
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s CMD /opt/conda/bin/conda run -n mmpbsa gmx_MMPBSA -h >/dev/null 2>&1 || exit 1
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8008"]
