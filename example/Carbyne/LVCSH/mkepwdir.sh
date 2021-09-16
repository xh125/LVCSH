#!/bin/bash
ncore=28
MODULEPATH="/share/home/zw/xiehua/opt/modules-4.7.1/modulefiles"
lvcsh_version="0.6.6"
QUEUE_NAME="privateq-zw"
for i in $(seq 80 40 80)
  do
    mkdir epw$i
    mkdir epw$i/QEfiles
    cp ../epw/epw$i.out epw$i/QEfiles/
    
    cp lvcsh.bsub epw$i
    sed -i "s/ncore/$ncore/g" epw$i/lvcsh.bsub
    sed -i "s:JOB_NAME:lvcsh-epw${i}-n0:g" epw$i/lvcsh.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" epw$i/lvcsh.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH:g" epw$i/lvcsh.bsub
    sed -i "s:version:$lvcsh_version:g" epw$i/lvcsh.bsub
    cp job.sh epw$i
    cp LVCSH.in epw$i
    sed -i "s:./epw.out:../../QEfiles/epw$i.out:g" epw$i/LVCSH.in
    sed -i "s:ncore:ncore         = $ncore !:g" epw$i/LVCSH.in    
    cp lvcsh-test.bsub epw$i/lvcsh-plot.bsub
    sed -i "2s/JOB_NAME/lvcsh-epw$i-plot/g" epw$i/lvcsh-plot.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" epw$i/lvcsh-plot.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH :g" epw$i/lvcsh-plot.bsub
    sed -i "s:version:$lvcsh_version:g" epw$i/lvcsh-plot.bsub
    
    cp LVCSH.in epw$i/QEfiles
    sed -i "s:verbosity:verbosity     = "high" !:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:./epw.out:./epw$i.out:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:naver:naver         = 10 !:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:nsnap:nsnap         = 2  !:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:savedsnap:savedsnap     = 2 !:g" epw$i/QEfiles/LVCSH.in
    
    cp lvcsh-test.bsub epw$i/QEfiles
    cd epw$i/QEfiles
    sed -i "2s/JOB_NAME/lvcsh-epw$i-test/g" lvcsh-test.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" lvcsh-test.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH :g" lvcsh-test.bsub
    sed -i "s:version:$lvcsh_version:g" lvcsh-test.bsub
    bsub < lvcsh-test.bsub
    cd ../..    
  
  done
