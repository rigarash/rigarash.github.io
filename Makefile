CURRENT_UID := $(shell id -u)
CURRENT_GID := $(shell id -g)

ZENN_IMAGE := ghcr.io/rigarash/zenn:latest

all: preview

build:
	docker image build --tag ${ZENN_IMAGE} docker

%:
	docker container run --init --rm --interactive --tty --volume=${PWD}:/work --user=${CURRENT_UID}:${CURRENT_GID} --publish 8000:8000 ${ZENN_IMAGE} $@

.PHONY: build
