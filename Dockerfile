FROM docker.io/akshatvolta/m3gnet:new

RUN pip install minio python-dotenv

COPY .env /config/.env

ENV MINIO_ENDPOINT=192.168.9.13:30002
ENV MINIO_ACCESS_KEY=GbhgmzTziWIanCKPaay0
ENV MINIO_SECRET_KEY=miniostorage