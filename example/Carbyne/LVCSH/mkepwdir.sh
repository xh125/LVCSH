#!/bin/bash
# wgs 1
#MODULESH="/share/home/ZhuangW/xh/opt/modules/init/profile.sh"
#MODULEPATH="/share/home/ZhuangW/xh/modulefiles"
#lvcsh_version="0.6.8"
#QUEUE_NAME="privateq-zw"
#ncore=32
#naver=300

#wgs 2
MODULESH="/share/home/zw/xiehua/opt/modules-5.0.0/init/profile.sh"
MODULEPATH="/share/home/zw/xiehua/modulefiles"
lvcsh_version="0.6.8"
QUEUE_NAME="chkpnt_rerun_queue"
#QUEUE_NAME="privateq-zw"
ncore=28
naver=10
nsnap=1000
for i in $(seq 10 10 20) ; do
    rm -rf epw$i
    mkdir epw$i
    mkdir epw$i/QEfiles
    cp ../epw/epw$i.out epw$i/QEfiles/
    
    cp lvcsh.bsub epw$i
    sed -i "s/ncore/$ncore/g" epw$i/lvcsh.bsub
    sed -i "s:JOB_NAME:lvcsh-epw${i}-n0:g" epw$i/lvcsh.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" epw$i/lvcsh.bsub
    sed -i "s:SH_MODULE:$MODULESH:g" epw$i/lvcsh.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH:g" epw$i/lvcsh.bsub
    sed -i "s:version:$lvcsh_version:g" epw$i/lvcsh.bsub
    
    cp job.sh epw$i
    sed -i "s/NCORE_/$ncore/g" epw$i/job.sh
    sed -i "s/NAVER_/$naver/g" epw$i/job.sh
    sed -i "s:NSNAP_:$nsnap:g" epw$i/job.sh

    cp LVCSH.in epw$i
    sed -i "s:./epw.out:../../QEfiles/epw$i.out:g" epw$i/LVCSH.in
    sed -i "s:NCORE:$ncore :g" epw$i/LVCSH.in
    
    cp lvcsh-test.bsub epw$i/lvcsh-plot.bsub
    sed -i "s:JOB_NAME:lvcsh-epw$i-plot:g" epw$i/lvcsh-plot.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" epw$i/lvcsh-plot.bsub
    sed -i "s:SH_MODULE:$MODULESH:g" epw$i/lvcsh-plot.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH :g" epw$i/lvcsh-plot.bsub
    sed -i "s:version:$lvcsh_version:g" epw$i/lvcsh-plot.bsub


    cp LVCSH.in epw$i/QEfiles
    sed -i "s:LVCSH:lvcsh:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:LOW:high:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:./epw.out:./epw$i.out:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:NAVER:10:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:NSNAP:2:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:NNODE:1:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:NCORE:1:g" epw$i/QEfiles/LVCSH.in
    sed -i "s:savedsnap:savedsnap     = 2 !:g" epw$i/QEfiles/LVCSH.in
    
    cp epw$i/lvcsh-plot.bsub epw$i/QEfiles/lvcsh-test.bsub
    cd epw$i/QEfiles
    sed -i "2s/lvcsh-epw$i-plot/lvcsh-epw$i-test/g" lvcsh-test.bsub
    sed -i "s:QUEUE_NAME:$QUEUE_NAME:g" lvcsh-test.bsub
    sed -i "s:DIR_MODULEPATH:$MODULEPATH :g" lvcsh-test.bsub
    sed -i "s:version:$lvcsh_version:g" lvcsh-test.bsub
    bsub < lvcsh-test.bsub
    cd ../..    
  done
