FROM mcr.microsoft.com/devcontainers/python:1-3.11-bullseye
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    file xz-utils \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir fastapi uvicorn
RUN mkdir -p /opt/vina
# Prefer patched Vina 1.2.3 from local bundles if available
COPY patches/mod_vina_docking-v1.3-amend+01/bundles/vina_1.2.3_linux_x86_64.tar.xz /tmp/vina_1.2.3.tar.xz
RUN tar -xJf /tmp/vina_1.2.3.tar.xz -C /opt/vina --no-same-owner \
    && VINA_BIN=$(find /opt/vina -type f -name vina -perm -111 | head -n 1) \
    && if [ -z "$VINA_BIN" ]; then echo "vina binary not found in 1.2.3" && exit 1; fi \
    && file "$VINA_BIN" \
    && cp "$VINA_BIN" /usr/local/bin/vina \
    && VINA_SPLIT=$(find /opt/vina -type f -name vina_split -perm -111 | head -n 1) \
    && if [ -n "$VINA_SPLIT" ]; then cp "$VINA_SPLIT" /usr/local/bin/vina_split; fi \
    && rm -f /tmp/vina_1.2.3.tar.xz || true
# Copy AD4 parameter files for metal support
COPY patches/application-patched-vina-1.2.3.amend2/data/AD4_parameters.dat /opt/vina/AD4_parameters.dat
COPY patches/application-patched-vina-1.2.3.amend2/data/AD4Zn.dat /opt/vina/AD4Zn.dat
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
EXPOSE 8006
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8006"]
