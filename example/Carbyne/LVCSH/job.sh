#!/bin/bash
nnode=1
ncore=NCORE_
naver=NAVER_
nsnap=NSNAP_
nefre_sh=40
nhfre_sh=40
sed -i "s:LVCSH:lvcsh:g" LVCSH.in
sed -i "s:LOW:low:g" LVCSH.in
sed -i "s:NNODE:$nnode:g" LVCSH.in
sed -i "s:NCORE:$ncore:g" LVCSH.in
sed -i "s:NSNAP:$nsnap:g" LVCSH.in
sed -i "s:NAVER:$((naver / nnode)):g" LVCSH.in
sed -i "s:!nefre_sh:nefre_sh  = $nefre_sh !:g" LVCSH.in
sed -i "s:!nhfre_sh:nhfre_sh  = $nhfre_sh !:g" LVCSH.in
for i in $(seq 1 1 $nnode) ; do
	mkdir node$i
	cp ./lvcsh.bsub ./node$i/	
	sed -i "2s/n0/n$i/g" ./node$i/lvcsh.bsub
	cp ./LVCSH.in ./node$i/
	cd ./node$i
	bsub < lvcsh.bsub
	cd ..
done
