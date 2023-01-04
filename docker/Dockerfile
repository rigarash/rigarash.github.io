FROM node:current
RUN npm install -g zenn-cli
WORKDIR /work
ENTRYPOINT ["npx", "zenn"]
CMD ["preview"]

LABEL org.opencontainers.image.source="https://github.com/rigarash/rigarash.github.io" \
      org.opencontainers.image.authors="Ryo IGARASHI <rigarash@gmail.com>"
