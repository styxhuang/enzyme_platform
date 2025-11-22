FROM mcr.microsoft.com/devcontainers/python:1-3.11-bullseye
WORKDIR /app
ENV PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
COPY dist/ /app/dist/
EXPOSE 80
CMD ["python", "-m", "http.server", "80", "--directory", "/app/dist"]