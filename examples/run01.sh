#!/bin/bash
FASTJET_CONFIG="fastjet-config"
FASTJET_FLAGS=`$FASTJET_CONFIG --cxxflags --libs --plugins=yes`

file=example01

if [ -f $file ]
then
rm $file
fi

g++ -std=c++11 $file.cc -o $file $FASTJET_FLAGS -lDynamicRPlugin
./$file

