#!/bin/bash

## Updated: version-v1.02 01/04/17 

printf "\nIn case of any issue at any point in phasTER analyses, contact authors at:\n"
# printf "Atul Kakrana: kakrana@udel.edu\n"
# printf "Reza Hammond: rkweku@udel.edu\n\n"
printf "https://github.com/atulkakrana/phasTER/issues\n\n"

#### Compile PHASER #########

if [ -d ./__pycache__ ]; then
  rm -r __pycache__/ 
fi

if [ ! -f ./phaser.py ]
then
  printf "phaser source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n"
  printf "Download complete package and rerun script\n\n"
  exit 1
else
  python3 -m py_compile phaser.py
  cd __pycache__ && mv phaser.*.pyc phaser 
  chmod u+x phaser && cp phaser ../
  printf "Done:      phaser\n"
  cd ..
  rm ./phaser.py
fi


#### Compile Collapser #########

if [ -d ./__pycache__ ]; then
  rm -r __pycache__/ 
  # printf "Found\n"
  # exit 1
fi

if [ ! -f ./collapser.py ]
then
  printf "collapser source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n"
  printf "Download complete package and rerun script\n\n"
  exit 1
else
  python3 -m py_compile collapser.py
  cd __pycache__ && mv collapser.*.pyc collapser 
  chmod u+x collapser && cp collapser ../
  printf "Done:      collapser\n"
  cd ..
  rm ./collapser.py
fi

#### Compile revFerno #########

if [ -d ./__pycache__ ]; then
  rm -r __pycache__/ 
fi

if [ ! -f ./revferno.py ]
then
  printf "revferno source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n"
  printf "Download complete package and rerun script\n\n"
  exit 1
else
  python3 -m py_compile revferno.py
  cd __pycache__ && mv revferno.*.pyc revferno 
  chmod u+x revferno && cp revferno ../
  printf "Done:      revferno\n"
  cd ..
  rm ./revferno.py
fi

#### Place core scripts ######

if [ -d ~/phaster ]
then
  rm -rf ~/phaster
  mkdir ~/phaster 
  # printf "file deleted\n"
else
  mkdir ~/phaster
fi

if [ ! -f ./phasclust.genome.v2.pl ]
then
  printf "phaster-core source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n\n"
  printf "Download complete package and rerun script\n"
  exit 1
else
  cp phasclust.genome.v2.pl ~/phaster
  rm phasclust.genome.v2.pl
fi

if [ ! -f ./phasclust.MUL.v2.pl ]
then
  printf "phaster-core source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n"
  printf "Download complete package and rerun script\n\n"
  exit 1
else
  cp phasclust.MUL.v2.pl ~/phaster
  rm phasclust.MUL.v2.pl
fi

if [ ! -f ./sPARTA.py ]
then
  printf "sPARTA source not found in current directory\n"
  printf "Please check that you downloaded complete archive\n"
  printf "Download complete package and rerun script\n\n"
  exit 1
else
  cp sPARTA.py ~/phaster
  rm sPARTA.py
fi

## Remove old folder
if [ -d ~/.phaster ]
then
  rm -r ~/.phaster 
  # printf "file deleted\n"
fi

mv ~/phaster ~/.phaster
printf "Done:      core scripts\n"

#### Clean up
if [ -d ./__pycache__ ]; then
  rm -r __pycache__/ 
fi
printf "Deleted:   source files\n"
#############################
printf "\nNote:'install.sh' cannot be re-used after a successful installation\n"
printf "Note: If installation fails, then recopy all files and try again\n"
printf "\nphasTER tool-set is ready be used\n"
printf "See readme here: https://github.com/atulkakrana/phasTER\n\n"

exit 1


### Change Log 
## V1.0 -> v1.01
## Edited names of phasTER-core scripts
## v1.01 -> v1.02
## Added sPARTA script