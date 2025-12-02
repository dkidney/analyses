#!/bin/bash

pyenv --version
#which pyenv

#pyenv install --list | grep -E "^\s*(\d+\.\d+\.\d+)$" | tail -10
#py_version=$(pyenv install --list | grep -E "^\s*(\d+\.\d+\.\d+)$" | tail -1)
#py_version=$(pyenv install --list | grep -E "^\s*(3\.10\.\d+)$" | tail -1)
py_version=3.10.16
echo ${py_version}
pyenv install ${py_version}
$(pyenv root)/versions/${py_version}/bin/python -m venv .venv
source .venv/bin/activate
echo $(which python)
echo $(python --version)
pip install -U pip ipykernel nbformat numpy pandas matplotlib pyarrow igraph leidenalg
pip freeze > requirements.txt