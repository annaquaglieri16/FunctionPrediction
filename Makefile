.PHONY: venv venv-ml kernel test docker-build docker-test

# Core env: digestion, tests, data IO, notebooks. Safe to re-run.
venv:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -r requirements-core.txt ipykernel
	$(MAKE) kernel

# Adds torch + transformers on top of `make venv` (large, slow).
venv-ml:
	.venv/bin/pip install -r requirements-ml.txt

# Registers a kernelspec pointing at .venv/bin/python so VSCode / Jupyter
# can select this env as "Python (plm-benchmark)".
kernel:
	.venv/bin/python -m ipykernel install --user \
		--name plm-benchmark --display-name "Python (plm-benchmark)"

test:
	PYTHONPATH=src pytest tests/ -v

docker-build:
	docker compose -f docker/docker-compose.yml build

docker-test:
	docker compose -f docker/docker-compose.yml run --rm benchmark pytest tests/ -v