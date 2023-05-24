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
git checkout fb428b53b8bfd642db09d71ac34c41351e492b70
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -

# Same process for raticate.

if [ ! -e source-tatami_r ]
then 
    git clone https://github.com/tatami-inc/tatami_r source-tatami_r
else 
    cd source-tatami_r
    git pull
    cd -
fi

cd source-tatami_r
git checkout 57d2d982c4e609c448e4d53738ab00feaf3ad40c
rm -rf ../tatami_r
cp -r include/tatami_r/ ../tatami_r
git checkout master
cd -

# Same process for manticore

if [ ! -e source-manticore ]
then 
    git clone https://github.com/tatami-inc/manticore source-manticore
else 
    cd source-manticore
    git pull
    cd -
fi

cd source-manticore
git checkout 313002dcea65ddb0fd334034f0f618626c736037
rm -rf ../manticore
cp -r include/manticore/ ../manticore
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
