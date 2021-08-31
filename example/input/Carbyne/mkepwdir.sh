#!/bin/bash
for i in $(seq 80 40 200)
	do
#		mkdir epw$i
#		mkdir epw$i/QEfiles
#		cp QEfiles/epw$i.out epw$i/QEfiles/
		cp lvcsh.bsub epw$i
		sed -i "2s/lvcsh-epw40/lvcsh-epw$i-n0/g" epw$i/lvcsh.bsub
		cp job.sh epw$i
		cp LVCSH.in epw$i
		sed -i "19s/40/$i/g" epw$i/LVCSH.in
		sed -i "20s/40/$i/g" epw$i/LVCSH.in
		sed -i "22s/40/$i/g" epw$i/LVCSH.in
		cp LVCSH.in epw$i/QEfiles
		sed -i "2s/low/high/g" epw$i/QEfiles/LVCSH.in
		sed -i "19s/40/$i/g" epw$i/QEfiles/LVCSH.in
		sed -i "20s/40/$i/g" epw$i/QEfiles/LVCSH.in
		sed -i "22s/epw40.out/epw$i.out/g" epw$i/QEfiles/LVCSH.in
		sed -i "23s/1/10/g" epw$i/QEfiles/LVCSH.in
		sed -i "25s/1000/2/g" epw$i/QEfiles/LVCSH.in
		sed -i "46s/10/1/g" epw$i/QEfiles/LVCSH.in
		cp lvcsh-test.bsub epw$i/QEfiles
		cd epw$i/QEfiles
		sed -i "2s/lvcsh-epw40/lvcsh-epw$i/g" lvcsh-test.bsub
		bsub < lvcsh-test.bsub
		cd ../..		
	done
