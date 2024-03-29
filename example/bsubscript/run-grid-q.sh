#!/bin/bash

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to calculate the phonon dispersion on a GRID"
$ECHO "for Carbyne. Only q-points are split."


PREFIX="carbyne"
nq=21
TMP_DIR="."

# clean TMP_DIR
for q in $(seq 1 1 $nq) ; do
rm -rf $TMP_DIR/$q
done


for q in $(seq 1 1 $nq) ; do
mkdir $TMP_DIR/$q
cp -r ../scf/$PREFIX.* $TMP_DIR/$q
cat > ph.q$q.in << EOF
phonons calculation
&inputph
  tr2_ph=1.0d-16,
  prefix='$PREFIX',
  outdir='$TMP_DIR/$q'
!  epsil=.true. !use for insulators
  ldisp=.true.
  nq1=40, nq2=1, nq3=1
  start_q=$q
  last_q=$q
!  recover=.true.
!  amass(1)=0.0
  fildyn='$PREFIX.dyn'
  fildvscf='dvscf'
/
EOF

$ECHO "  running the phonon calculation for q= "$q

cp qe-grid.bsub qe-q$q.bsub
sed -i "s:JOBNAME:ph-$PREFIX-q$q:g" qe-q$q.bsub
sed -i "s:INPUTNAME:ph.q$q.in:g" qe-q$q.bsub
sed -i "s:OUTPUTNAME:ph.q$q.out:g" qe-q$q.bsub
sed -i "s:nodelist:nodelist.q$q:g" qe-q$q.bsub
bsub < qe-q$q.bsub
sleep 10

done

