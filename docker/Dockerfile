FROM python:3.11-slim

WORKDIR /app

COPY requirements-core.txt requirements-ml.txt ./
# Core deps by default (fast to build). Uncomment the ML line once you're
# ready to extract embeddings inside the container (large: torch + transformers).
RUN pip install --no-cache-dir -r requirements-core.txt
# RUN pip install --no-cache-dir -r requirements-ml.txt

COPY src/ src/
COPY configs/ configs/
COPY scripts/ scripts/
COPY tests/ tests/

ENV PYTHONPATH=/app/src

CMD ["pytest", "tests/", "-v"]
