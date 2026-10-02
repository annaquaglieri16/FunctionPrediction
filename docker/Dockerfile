FROM python:3.11-slim

WORKDIR /app

# Install the plm_benchmark package and its core dependencies from
# pyproject.toml. This provides the methods that the Nextflow processes call. 
# Configs and data are passed in by Nextflow at run time, so they are not copied into the image.
COPY pyproject.toml ./
COPY src/ src/
RUN pip install --no-cache-dir .

# For the embedding stage (large: torch + transformers), install the ML extra:
# RUN pip install --no-cache-dir ".[ml]"
