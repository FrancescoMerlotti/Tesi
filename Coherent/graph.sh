#!/bin/bash

for i in {1..12}
do
	for j in 0 1
	do
		./Coherent.exe $i $j
	done
done