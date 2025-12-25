FROM mcr.microsoft.com/devcontainers/python:1-3.11-bullseye
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
RUN apt-get update && apt-get install -y --no-install-recommends gromacs && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir fastapi uvicorn numpy matplotlib
COPY app.py /app/app.py
ENV ENZYME_DATA_DIR=/data
EXPOSE 8005
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8005"]
