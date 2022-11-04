#!/bin/zsh

$PYTHON = 

cd Rosetta/main/source/src/python/PyRosetta/
python build.py -j24 --create-package $HOME/my_pyrosetta_package
cd $HOME/my_pyrosetta_package/setup
python setup.py install
