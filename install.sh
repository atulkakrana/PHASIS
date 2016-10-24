#!/bin/bash


echo "In case of any issue at any point in phasTER analyses, contact:"
echo "Atul Kakrana: kakrana@udel.edu"
echo "Reza Hammond: rkweku@udel.edu"

#### Compile PHASER #########
if [ -d "./__pycache__" ]; then
  rm -r __pycache__/ 
fi

python3 -m py_compile phaser.py
cd __pycache__ && chmod u+x phaser.cpython-33.pyc 
mv phaser.cpython-33.pyc phaser && cp phaser ../
cd ..

#### Compile Collapser #########
if [ -d "./__pycache__" ]; then
  rm -r __pycache__/ 
fi

python3 -m py_compile collapser.py
cd __pycache__ && chmod u+x collapser.cpython-33.pyc 
mv collapser.cpython-33.pyc collapser && cp collapser ../
cd ..

#### Compile revFerno #########
if [ -d "./__pycache__" ]; then
  rm -r __pycache__/ 
fi

python3 -m py_compile mirferno.py
cd __pycache__ && chmod u+x mirferno.cpython-33.pyc 
mv mirferno.cpython-33.pyc collapser && cp mirferno ../
cd ..

#### Place core scripts ######
if [ -d "~/.phaster" ]; then
  rm -r ~/.phaster 
fi

mv phaster-core.pl ~/.phaster
# mv phaser ~/.phaster          ### Will need modification of environment variable
# mv collapser ~/.phaster       ### Will need modification of environment variable
# mv revferno ~/.phaster        ### Will need modification of environment variable


##### Clean-up ##############

if [ -d "./__pycache__" ]; then
  rm -r ./__pycache__ 
fi
if [ -f "./phaser.py" ]; then
  rm ~/phaser.py 
fi
if [ -f "./collapser.py" ]; then
  rm ~/.phaster 
fi
if [ -f "./revferno.py" ]; then
  rm ~/.phaster 

#############################
echo "phasTER tool-set is ready be used"