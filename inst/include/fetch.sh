#!/bin/bash

set -e
set -u

# Fetches all of the header files from tatami. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

if [ ! -e source-tatami ]
then 
    git clone https://github.com/LTLA/tatami source-tatami
else 
    cd source-tatami
    git pull
    cd -
fi

cd source-tatami
git checkout ca5835b9d7f9fc544a342eff475a9ba645410ab1
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -

# Same process for raticate.

if [ ! -e source-raticate ]
then 
    git clone https://github.com/LTLA/raticate source-raticate
else 
    cd source-raticate
    git pull
    cd -
fi

cd source-raticate
git checkout 7066ee0278682396c708ea4da7fa92b3a5181feb
rm -rf ../raticate
cp -r include/raticate/ ../raticate
git checkout master
cd -

# Same process for byteme.

if [ ! -e source-byteme ]
then 
    git clone https://github.com/LTLA/byteme source-byteme
else 
    cd source-byteme
    git pull
    cd -
fi

cd source-byteme
git checkout f9ed015693c45424aa2835d79fc8e58128a24c93
rm -rf ../byteme
cp -r include/byteme/ ../byteme
git checkout master
cd -
