# syntax=docker/dockerfile:1.4

# This Dockerfile is derived from the stjudecloud/rtcg Dockerfile (originally
# written by Michael Macias). Pulled and edited by Clay McLeod on Friday, March
# 3, 2023.

FROM python:3.9.16

ENV POETRY_VERSION=1.4.0 \
    POETRY_HOME=/opt/poetry \
    POETRY_VIRTUALENVS_CREATE=false \
    NGSDERIVE_HOME=/opt/ngsderive

ENV PATH=${POETRY_HOME}/bin:${PATH}

RUN curl -sSL https://install.python-poetry.org | python3 -

COPY poetry.lock pyproject.toml ${NGSDERIVE_HOME}/

WORKDIR $NGSDERIVE_HOME

COPY ngsderive/ ${NGSDERIVE_HOME}/ngsderive/
COPY tests/ ${NGSDERIVE_HOME}/tests/
COPY README.md ${NGSDERIVE_HOME}/README.md

RUN poetry install

RUN <<EOS
  useradd --system --create-home --base-dir /var/lib --shell /bin/bash ngsderive
EOS

USER ngsderive

ENTRYPOINT ["ngsderive"]