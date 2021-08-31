#!/bin/bash
	for i in {1..100}
	do 
	mkdir node$i
	cp ./lvcsh.bsub ./node$i/
#	cp ./lvcsh.pbs ./node$i/
	cp ./LVCSH.in ./node$i/
	cd ./node$i
#	qsub lvcsh.pbs
	bsub < lvcsh.bsub
	cd ..
	done
