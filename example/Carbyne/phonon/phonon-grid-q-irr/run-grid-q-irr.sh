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
$ECHO "for Carbyne. Both q-points and irreps are split."

PREFIX="carbyne"
nq=6
represent=16
TMP_DIR='./'

# clean TMP_DIR
for q in $(seq 1 1 $nq) ; do
for irr in $(seq 1 1 $represent) ; do
rm -rf $TMP_DIR/$q.$irr
done
done

for q in $(seq 1 1 $nq) ; do
for irr in $(seq 1 1 $represent) ; do
mkdir $TMP_DIR/$q.$irr
cat > $TMP_DIR/$q.$irr/ph.$q.$irr.in << EOF
phonons calculation
&inputph
  tr2_ph=1.0d-16,
  prefix='$PREFIX',
  outdir="./"
!  epsil=.true. !use for insulators
  ldisp=.true.
  nq1=10, nq2=1, nq3=1
  start_q=$q
  last_q=$q
  start_irr=$irr
  last_irr=$irr
  recover=.true.
!  amass(1)=0.0
  fildyn='$PREFIX.dyn'
  fildvscf='$PREFIX.dvscf'
/
EOF

cp -r $TMP_DIR/$PREFIX.* $TMP_DIR/$q.$irr
mkdir -p $TMP_DIR/$q.$irr/_ph0/$PREFIX.phsave
cp -r $TMP_DIR/_ph0/$PREFIX.phsave/* $TMP_DIR/$q.$irr/_ph0/$PREFIX.phsave


$ECHO "  running the phonon calculation for q= " $q " irr=" $irr

cp qe-grid.bsub $TMP_DIR/$q.$irr/qe-$q.$irr.bsub
cd $TMP_DIR/$q.$irr
sed -i "s:JOBNAME:phonon-$q.$irr:g" qe-$q.$irr.bsub
sed -i "s:INPUTNAME:ph.$q.$irr.in:g" qe-$q.$irr.bsub
sed -i "s:OUTPUTNAME:ph.$q.$irr.out:g" qe-$q.$irr.bsub
sed -i "s:nodelist:nodelist.$q.$irr:g" qe-$q.$irr.bsub
bsub < qe-$q.$irr.bsub
cd $EXAMPLE_DIR
done
done

