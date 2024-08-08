#!/bin/bash

for i in $1
do
	splitSets --name ${i} --threshold 20 --offset 150
	for j in ${i}.*
	do
		sumSets < ${j} > ${j}.sum
		gaussianKernel 2 500 < ${j}.sum > ${j}.smooth
	done
done
