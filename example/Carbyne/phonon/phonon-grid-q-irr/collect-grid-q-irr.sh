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
nq=2
represent=2
TMP_DIR='./'

#
#  collect also the representation 0 (contribution to the dynamical 
#  matrix independent from the induced charge).

for q in `seq 1 $nq ` ; do

for irr in `seq 1 $represent` ; do

\cp -f $TMP_DIR/$q.$irr/_ph0/$PREFIX.phsave/dynmat.$q.$irr.xml $TMP_DIR/_ph0/$PREFIX.phsave 2> /dev/null

done 
#
#  collect also the representation 0 (contribution to the dynamical 
#  matrix independent from the induced charge).
#
\cp -f $TMP_DIR/$q.1/_ph0/$PREFIX.phsave/dynmat.$q.0.xml $TMP_DIR/_ph0/$PREFIX.phsave 2> /dev/null

done 
#
# cp electric field part
#
\cp -f $TMP_DIR/1.1/_ph0/$PREFIX.phsave/tensors.xml $TMP_DIR/_ph0/$PREFIX.phsave 

cat > ph.collect.in << EOF
phonons calculation
&inputph
  tr2_ph=1.0d-16,
  prefix='$PREFIX',
  outdir="$TMP_DIR"
!  epsil=.true. !use for insulators
  ldisp=.true.
  nq1=10, nq2=1, nq3=1
  recover=.true.
!  amass(1)=0.0
  fildyn='$PREFIX.dyn'
  fildvscf='$PREFIX.dvscf'
/
EOF
$ECHO "  running the phonon calculation to collect the results...\c"

cp qe-grid.bsub qe-collect.bsub
sed -i "s:JOBNAME:phonon-collect:g" qe-collect.bsub
sed -i "s:INPUTNAME:ph.collect.in:g" qe-collect.bsub
sed -i "s:OUTPUTNAME:ph.collect.out:g" qe-collect.bsub
sed -i "s:nodelist:nodelist.collect:g" qe-collect.bsub
bsub < qe-collect.bsub


