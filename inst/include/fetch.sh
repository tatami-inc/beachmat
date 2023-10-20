#!/bin/bash

set -e
set -u

# Fetches all of the header files from tatami. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

harvester() {
    local name=$1
    local url=$2
    local version=$3

    local tmpname=source-${name}
    if [ ! -e $tmpname ]
    then 
        git clone $url $tmpname
    else 
        cd $tmpname
        git checkout master
        git pull
        cd -
    fi

    cd $tmpname
    git checkout $version
    rm -rf ../$name
    cp -r include/$name ../$name
    cd -
}

harvester tatami https://github.com/tatami-inc/tatami v2.1.0
harvester tatami_r https://github.com/tatami-inc/tatami_r v1.0.0
harvester manticore https://github.com/tatami-inc/manticore v1.0.1
harvester byteme https://github.com/LTLA/byteme v1.0.1
harvester tatami_chunked https://github.com/tatami-inc/tatami_chunked v1.0.1
