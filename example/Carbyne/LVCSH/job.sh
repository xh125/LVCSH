#!/bin/bash
nnode=10
sed -i "s:nnode:nnode = $nnode !:g" LVCSH.in
for i in $(seq 1 1 $nnode)
	do 
	mkdir node$i
	cp ./lvcsh.bsub ./node$i/	
	sed -i "2s/n0/n$i/g" ./node$i/lvcsh.bsub
	cp ./LVCSH.in ./node$i/
	cd ./node$i
	bsub < lvcsh.bsub
	cd ..
	done
