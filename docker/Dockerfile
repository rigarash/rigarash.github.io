FROM node:current-slim
RUN npm install -g zenn-cli
WORKDIR /work
ENTRYPOINT ["npx", "zenn"]
CMD ["preview"]

LABEL org.opencontainers.image.authors="Ryo IGARASHI <rigarash@gmail.com>" \
      org.opencontainers.image.url="https://github.com/rigarash/rigarash.github.io" \
      org.opencontainers.image.source="https://github.com/rigarash/rigarash.github.io/blob/zenn/docker/Dockerfile" \
      org.opencontainers.image.title="Zenn CLI Docker image for everyday use"
