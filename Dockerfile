FROM python:3.11-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    GPX6_HOME=/workspace/GPX6 \
    QDIR=/opt/Q6 \
    QTOOLS_HOME=/opt/qtools \
    PYTHONPATH=/opt/qtools/packages \
    PATH=/opt/qtools/qscripts-cli:/opt/qtools/bin:/opt/Q6/bin:/opt/q6/bin:$PATH

WORKDIR /workspace/GPX6

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    ca-certificates \
    curl \
    gfortran \
    git \
    less \
    libopenmpi-dev \
    make \
    openmpi-bin \
    openssh-client \
    pymol \
    rsync \
    slurm-client \
    texlive-fonts-recommended \
    texlive-latex-extra \
    texlive-science \
    latexmk \
    vim-tiny \
    && rm -rf /var/lib/apt/lists/*

COPY requirements-docker.txt /tmp/requirements-docker.txt
RUN python -m pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir -r /tmp/requirements-docker.txt

# Qtools provides Qpyl plus q_mapper.py and q_analysefeps.py used by the FEP scripts.
RUN git clone --depth 1 https://github.com/mpurg/qtools.git "$QTOOLS_HOME" \
    && if [ -f "$QTOOLS_HOME/requirements.txt" ]; then pip install --no-cache-dir -r "$QTOOLS_HOME/requirements.txt"; fi \
    && pip install --no-cache-dir -e "$QTOOLS_HOME" || true \
    && find "$QTOOLS_HOME" -maxdepth 3 -type f -name "*.py" -exec chmod +x {} \;

RUN git clone --depth 1 https://github.com/qusers/Q6.git "$QDIR" \
    && make -C "$QDIR/src" all COMP=gcc \
    && make -C "$QDIR/src" mpi COMP=gcc \
    && cd "$QDIR/bin" \
    && if [ -e Qprep6 ]; then ln -sf Qprep6 qprep6 && ln -sf Qprep6 qprep5; fi \
    && if [ -e Qdyn6 ]; then ln -sf Qdyn6 qdyn6 && ln -sf Qdyn6 qdyn5_r8; fi \
    && if [ -e Qdyn6p ]; then ln -sf Qdyn6p qdyn6p; fi \
    && if [ -e Qfep6 ]; then ln -sf Qfep6 qfep6 && ln -sf Qfep6 qfep5; fi

COPY docker/entrypoint.sh /usr/local/bin/gpx6-entrypoint
RUN chmod +x /usr/local/bin/gpx6-entrypoint

ENTRYPOINT ["gpx6-entrypoint"]
CMD ["bash"]
