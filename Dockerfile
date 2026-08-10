FROM python:3.12-slim

WORKDIR /app

# Install uv for fast package installation
COPY --from=ghcr.io/astral-sh/uv:0.4.30 /uv /usr/local/bin/uv

# Copy package files
COPY pyproject.toml uv.lock README.md LICENSE.md MANIFEST.in ./
COPY rfinder/ rfinder/

# python-casacore ships a self-contained manylinux_2_28 wheel (confirmed working on this base -
# no system casacore/boost/cfitsio needed). If that ever stops being true for a future version,
# fall back to a kernsuite base instead:
#   FROM kernsuite/base:10
#   RUN docker-apt-install python3-casacore
RUN uv sync --frozen

# Ensure the venv-installed CLI is available
ENV PATH="/app/.venv/bin:${PATH}"

# Make CLI available
CMD ["rfinder", "--help"]
