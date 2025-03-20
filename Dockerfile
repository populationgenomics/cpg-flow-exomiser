# Any analysis-runner driver image must at least include git.
ARG PY_VER=${PY_VER:-3.10}

FROM python:${PY_VER}-slim-bullseye AS basic

ENV PYTHONDONTWRITEBYTECODE=1

RUN apt update && apt install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        gnupg \
        openjdk-17-jdk-headless \
        zip && \
    pip install --no-cache-dir -U pip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# now do some fun stuff, installing ClinvArbitration
WORKDIR /cpg_exomiser

COPY src src/
COPY pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
