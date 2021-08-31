#!/bin/bash
	for i in {1..10}
	do 
	mkdir node$i
	cp ./lvcsh.bsub ./node$i/	
	sed -i "2s/node/node$i/s" ./node$i/lvcsh.bsub
#	cp ./lvcsh.pbs ./node$i/
	cp ./LVCSH.in ./node$i/
	cd ./node$i
#	qsub lvcsh.pbs
	bsub < lvcsh.bsub
	cd ..
	done
