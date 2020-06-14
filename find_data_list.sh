#!/bin/bash


raw=`find . -name features.tsv.gz |grep filtered`

for name in ${raw[@]};
do 
    echo `dirname $name`
done
