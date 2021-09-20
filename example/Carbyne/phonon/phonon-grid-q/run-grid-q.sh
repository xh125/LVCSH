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
nq=6
TMP_DIR="./"

# clean TMP_DIR
for q in $(seq 1 1 $nq) ; do
rm -rf $TMP_DIR/$q
done


for q in $(seq 1 1 $nq) ; do
mkdir $TMP_DIR/$q
cp -r $TMP_DIR/$PREFIX.* $TMP_DIR/$q
cat > $TMP_DIR/$q/ph.$q.in << EOF
phonons calculation
&inputph
  tr2_ph=1.0d-16,
  prefix='$PREFIX',
  outdir='./'
!  epsil=.true. !use for insulators
  ldisp=.true.
  nq1=10, nq2=1, nq3=1
  start_q=$q
  last_q=$q
  recover=.true.
!  amass(1)=0.0
  fildyn='$PREFIX.dyn'
  fildvscf='$PREFIX.dvscf'
/
EOF

$ECHO "  running the phonon calculation for q= "$q

cp qe-grid.bsub $TMP_DIR/$q/qe-$q.bsub
cd $TMP_DIR/$q
sed -i "s:JOBNAME:phonon-q$q:g" qe-$q.bsub
sed -i "s:INPUTNAME:ph.$q.in:g" qe-$q.bsub
sed -i "s:OUTPUTNAME:ph.$q.out:g" qe-$q.bsub
sed -i "s:nodelist:nodelist.$q:g" qe-$q.bsub
bsub < qe-$q.bsub

cd $EXAMPLE_DIR

done

