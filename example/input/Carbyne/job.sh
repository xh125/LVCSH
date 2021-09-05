#!/bin/bash
	for i in {1..10}
	do 
	mkdir node$i
	cp ./lvcsh.bsub ./node$i/	
	sed -i "2s/n0/n$i/g" ./node$i/lvcsh.bsub
	cp ./LVCSH.in ./node$i/
	cd ./node$i
	bsub < lvcsh.bsub
	cd ..
	done
