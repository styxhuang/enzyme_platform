FROM mcr.microsoft.com/devcontainers/python:1-3.11-bullseye
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    file \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir fastapi uvicorn
ARG VINA_TARBALL=autodock_vina_1_1_2_linux_x86.tgz
COPY ${VINA_TARBALL} /tmp/autodock_vina.tar.gz
RUN mkdir -p /opt/vina \
    && tar -xzf /tmp/autodock_vina.tar.gz -C /opt/vina \
    && VINA_BIN=$(find /opt/vina -type f -name vina -perm -111 | head -n 1) \
    && if [ -z "$VINA_BIN" ]; then echo "vina binary not found in tarball" && exit 1; fi \
    && file "$VINA_BIN" \
    && cp "$VINA_BIN" /usr/local/bin/vina \
    && VINA_SPLIT=$(find /opt/vina -type f -name vina_split -perm -111 | head -n 1) \
    && if [ -n "$VINA_SPLIT" ]; then cp "$VINA_SPLIT" /usr/local/bin/vina_split; fi \
    && rm -f /tmp/autodock_vina.tar.gz
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
EXPOSE 8006
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8006"]
