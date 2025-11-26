FROM mcr.microsoft.com/devcontainers/python:1-3.11-bullseye
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    openbabel \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir "numpy<2.0.0" fastapi uvicorn rdkit-pypi pillow meeko httpx
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
EXPOSE 8003
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8003"]
