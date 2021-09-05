#!/bin/bash
ncore=28
for i in $(seq 160 40 160)
	do
		mkdir epw$i
		mkdir epw$i/QEfiles
		cp ../epw/epw$i.out epw$i/QEfiles/
		cp lvcsh.bsub epw$i
		sed -i "s/ncore/$ncore/g" epw$i/lvcsh.bsub
		sed -i "2s:lvcsh-epw:lvcsh-epw${i}-n0:g" epw$i/lvcsh.bsub
		cp job.sh epw$i
		cp LVCSH.in epw$i
		sed -i "s:epw40:epw$i:g" epw$i/LVCSH.in
		cp LVCSH.in epw$i/QEfiles
		sed -i "2s/low/high/g" epw$i/QEfiles/LVCSH.in
		sed -i "s:../../QEfiles/epw40.out:epw$i.out:g" epw$i/QEfiles/LVCSH.in
		cp lvcsh-test.bsub epw$i/QEfiles
		cd epw$i/QEfiles
		sed -i "2s/lvcsh-epw40/lvcsh-epw$i/g" lvcsh-test.bsub
		bsub < lvcsh-test.bsub
		cd ../..		
	done
