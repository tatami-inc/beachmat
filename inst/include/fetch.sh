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
git checkout 898039be5c0f9a4ff5de6f95ba7bdfd8d3940119
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
git checkout ba2ec1690b01f602ad4d95e800d6aa8bf00e29a5
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
