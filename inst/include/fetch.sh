#!/bin/bash

set -e
set -u

# Fetches all of the header files from tatami. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

if [ ! -e source-tatami ]
then 
    git clone https://github.com/tatami-inc/tatami source-tatami
else 
    cd source-tatami
    git pull
    cd -
fi

cd source-tatami
git checkout d6cb42539e434efce7fd24cd5abcbb1c006f74f3
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -

# Same process for tatami_r.

if [ ! -e source-tatami_r ]
then 
    git clone https://github.com/tatami-inc/tatami_r source-tatami_r
else 
    cd source-tatami_r
    git pull
    cd -
fi

cd source-tatami_r
git checkout 62ec62544ee1f309cb5164517a58fc4d0a3218d8
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
git checkout 96a9fc13fa9a4c6e0d921d4dad9cea132256640d
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
git checkout 1baca1ffddfbf0cd71fa7e285554959e3a116911
rm -rf ../byteme
cp -r include/byteme/ ../byteme
git checkout master
cd -
