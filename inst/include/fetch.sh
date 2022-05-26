# Fetches all of the header files from tatami and byteme.

if [ ! -e source-tatami ]
then 
    git clone https://github.com/LTLA/tatami source-tatami
else 
    cd source-tatami
    git pull
    cd -
fi
rm -rf tatami
cp -r source-tatami/include/tatami/ tatami

if [ ! -e source-byteme ]
then 
    git clone https://github.com/LTLA/byteme source-byteme
else 
    cd source-byteme
    git pull
    cd -
fi
rm -rf byteme
cp -r source-byteme/include/byteme byteme
