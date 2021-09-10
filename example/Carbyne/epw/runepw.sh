#!/bin/bash
for i in $(seq 20 20 240)
	do
	cp epw.in epw${i}.in
	sed -i "s:epwwrite:epwwrite=.false. ! :g" epw${i}.in
	sed -i "s:epwread:epwread=.true. !:g" epw${i}.in
	sed -i "s:wannierize:wannierize=.false. !:g" epw${i}.in
	sed -i "s:nkf1:nkf1=$i !:g" epw${i}.in
	sed -i "s:nqf1:nqf1=$i !:g" epw${i}.in

	cp qe-epw.bsub qe-epw${i}.bsub
	sed -i "2s:epw:epw${i}:g" qe-epw${i}.bsub
	sed -i "s:epw.in:epw${i}.in:g" qe-epw${i}.bsub
	sed -i "s:epw.out:epw${i}.out:g" qe-epw${i}.bsub	

	bsub < qe-epw${i}.bsub
  wait 1
	done
