FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=$PATH:/usr/local/mgltools/bin
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple

# ===== 1. 使用国内源，加速 apt =====
RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|http://mirrors.aliyun.com/ubuntu/|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com/ubuntu/|http://mirrors.aliyun.com/ubuntu/|g' /etc/apt/sources.list && \
    apt update && apt install -y \
    python3 \
    python3-pip \
    wget \
    tar \
    unzip \
    build-essential \
    git \
    autoconf \
    automake \
    libtool \
    libgl1-mesa-glx \
    libxrender1 \
    libxt6 \
    libsm6 \
    libxext6 \
    software-properties-common \
    openbabel \
    libopenbabel-dev \
    && pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple \
        fastapi \
        uvicorn \
        pydantic \
        httpx \
        meeko \
        rdkit \
    && rm -rf /var/lib/apt/lists/*

# ===== 2. 安装 MGLTools =====
WORKDIR /opt
COPY mgltools_x86_64Linux2_1.5.7p1.tar.gz /opt/MGLTools.tar.gz
RUN tar -zxvf MGLTools.tar.gz && \
    cd mgltools_x86_64Linux* && \
    ./install.sh -d /usr/local/mgltools && \
    rm -rf /opt/MGLTools*

# ===== 3. 安装 AutoDock4 =====
COPY autodocksuite-4.2.6-x86_64Linux2.tar /opt/autodock.tar

RUN tar -xvf /opt/autodock.tar && \
    cp x86_64Linux2/autodock4 x86_64Linux2/autogrid4 /usr/local/bin/ && \
    chmod +x /usr/local/bin/autodock4 /usr/local/bin/autogrid4 && \
    rm -rf /opt/x86_64Linux2 /opt/autodock.tar

# ===== 4. 安装 AutoDock Vina =====
RUN wget -c https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/local/bin/vina && \
    chmod +x /usr/local/bin/vina

ENV ENZYME_DATA_DIR=/data
ENV ADT_PREP_RECEPTOR=/usr/local/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py
ENV PATH=$PATH:/usr/local/mgltools/bin

WORKDIR /app
COPY app.py /app/app.py
EXPOSE 8003
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8003"]
