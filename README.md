# This program is use to calculate the dynamic of exciton in two dimentional hitstructure
 the step is as follows
 (1) in the workspace directory mkdir opt and optmizttion the structure
 (2) in the workspace mkdir phonon-gamma and calculate vibration
      (the POSCAR need to use Caracter POSCAR)
 (3) in the workspace mkdir wannier and calculate wannier Wavefunction ,
      and band project on each local WF use wannier package changed by xh
      could save file wannier90_bandproj.dat
      need set 
      write_xyz=.TRUE.
      bands_num_points= 80
      bands_plot_project=22
 (4) in the workspace mkdir wannierinput and in this directory put the wannier input file 
     (INCAR wannier90.win)and bsub script(vasp.bsub) (bsub name called *wannier eg. 
      mpijob.intelmpi vasp_std > C_wannier and in the line=8)
 (5) in workspace directory exection NsPOSCAR.x with shiftinp file
    shiftinp
&input
dQ        = 0.0020
nstep      = 5
lbsub      = .TRUE.
amass(1)   =  95.942
amass(2)   =  32.065
amass(3)   =  183.841
/
 (6) the exection will mkdir workspace/nomashift/noma_imode/shift_ldQ/
    and in this directory calculate wannier90
 (7) in directory workspace add file SHIN as the exciton SH input file 
    and exect SHexciton.x
 (8) the SHIN as follows
    na1site = 500
    na2site = 500
    na3site = 500
    temp    = 320
    gamma   = 0.02
    dt      = 0.01
    nstep   = 20000
    nsnap   = 100
    naver   = 500
    elecB   = 18
    elecK   = 100
    holeB   = 17
    holeK   = 100
    epsr    = 1.0           !tixi jintai jiedian hangs 
    nshiftstep = 5          !nomal shift in +-5  need same as in shiftinp nstep
    stadQ   = 0.0020        !nomal shift step    need same as in shiftinp dQ
    L_hotphonon= .true.
    hot_mode = 1
    hot_scal = 1.5
 (9)
 