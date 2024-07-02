#!/bin/bash
FASTJET_CONFIG="fastjet-config" # if fastjet path is not set, use full path: /path/to/fastjet/bin/fastjet-config
FASTJET_FLAGS=`$FASTJET_CONFIG --cxxflags --libs --plugins=yes`

file=03-jetarea-active

if [ -f $file ]
then
rm $file
fi

g++ -std=c++11 $file.cc -o $file $FASTJET_FLAGS -lDynamicRPlugin
./$file

