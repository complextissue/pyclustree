

.DEFAULT_GOAL := help


.PHONY : help
help:
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)


#
# Virtual environment setup
# You may want to replace python3 with the path to your python3 executable,
# e.g. the output of `pyenv which python3`when using pyenv
#
.PHONY : create-venv
create-venv: ## Create virtualenv
	python3 -m virtualenv '.venv'
	echo "Don't forget to activate with 'source .venv/bin/activate'"

#
# Install requirements
#
.PHONY : install
install: ## Install dependencies for production
	python3 -m pip --require-virtualenv install --upgrade pip
	python3 -m pip --require-virtualenv install flit poetry
	python3 -m flit install --deps production

.PHONY : install-dev
install-dev: install ## Install dependencies for development and production
	python3 -m pip --require-virtualenv install '.[dev]'
	python3 -m certifi
	pre-commit install --hook-type pre-commit --hook-type pre-push
	echo "Please also install pandoc to create the documentation."


#
# Build and upload package
#
.PHONY : build
build: ## Build python package
	python3 -m build
	python3 -m twine check --strict dist/*.whl

.PHONY : upload
upload: ## Publish package with flit
	python3 -m flit publish

#
# Testing
#
.PHONY : pytest
pytest: ## Run pytest with coverage
	python3 -m coverage run -m pytest --maxfail=10

.PHONY : coverage-report
coverage-report: pytest ## Generate coverage report
	python3 -m coverage report -m
	python3 -m coverage html
