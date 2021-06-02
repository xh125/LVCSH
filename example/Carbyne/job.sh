#!/bin/bash
	for i in {1..10}
	do 
	mkdir node$i
	cp ./lvcsh.pbs ./node$i/
	cp ./LVCSH.in ./node$i/
	cd ./node$i
	qsub lvcsh.pbs
	cd ..
	done
