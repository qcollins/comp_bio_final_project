#!/bin/bash
SOURCE="$1"
TARGET="$2"

echo "Reading from $SOURCE"
echo "Writing to $TARGET"

for file in $SOURCE/*
do
    echo "Running $file"
    name=${file##*/}
    base=${name%.bpm}
    python ../genecentric/genecentric-go -e data/essentials data/yeast_emap.gi $file $TARGET/$base.gobpm 
done