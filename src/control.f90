module control
  !-----------------------------------------------------------------------
  !!
  !! Common variables for the SH program
  !!  
  !-----------------------------------------------------------------------
  use kinds,    only    : dp
  use constants,only    : maxlen
  implicit none
  
  !是否跑SH-MD，用于生成相关参数及文件，用于脚本批处理
  logical :: Lrunsh
  !!sh type -- "elec","hole","elechole","exciton"
  character(len=maxlen) :: shtype
  !!system dimention of real and kpoint parameter 
  integer :: dimention
  !体系的维度
  integer :: na1
  integer :: na2
  integer :: na3
  !在三个方向上的重复单元
  integer :: nk1
  integer :: nk2
  integer :: nk3
  !如果采用k空间的bloch态来描述电子态，k空间网格的大小
  integer :: nq1
  integer :: nq2
  integer :: nq3
  !用于描述声子的q的网格大小
  
  !! representation of the system -- "adiabatic","blochstat","wfstat","atomob","molerob"
  character(len=maxlen):: representation
  !体系电子态描述的表象
  !!when use wfstat representation ,use method of nomal shift to calculate elec-phonon coupling
  !!the parameter must == in normal shift input
  integer :: nshiftstep
  !体系随简正坐标简正模方向有限位移的次数
  real(kind=dp) :: dtadq
  !有限位移的步长，单位为A
  real(kind=dp) :: pj_adjust
  !用于描述TB模型转移积分截断的判据，也可以改用wannier的wannier90_hr.cut.dat的数据
  !!
  real(kind=dp) :: temp
  !体系的温度
  real(kind=dp) :: gamma
  !核动力学的阻尼，单位ps^-1,表征体系的热传导
  real(kind=dp) :: dt
  !动力学的步长，单位fs
  integer :: nstep
  !动力学总的步长
  integer :: nsnap
  !每一步计算的帧数
  integer :: naver
  !系综采样的条数
  
  !! initial core state init_normal_stat--"class","quantum" 
  character(len=maxlen) ::  initnmstat
  !简谐振动动力学的计算方法，class对应能量均分定理，quantum对应谐振子的能量分布
  logical       :: Lrandomsita
  logical       :: L_hotphonon
  integer       :: hot_mode
  real(kind=dp) :: hot_scal
  
  !! initial elec-hole state--"BK","WF","EN","ES"
  !! 分别对应 blochstat 下的kpoint-band,
  !! wfstat -- wannier-function
  !! wfstat -- energy
  !! wfstat -- energy-surface
  integer ::  Num_occupied
  !基态电子占据的能带数
  character(len=maxlen) ::  initehstat
  !初始电子和空穴状态的指定方式
  integer :: initek            !init_elec_K
  integer :: initeb            !init_elec_band
  integer :: inithk            !init_hole_K
  integer :: inithb            !init_hole_band
  !采样能带的K，与band来指定
  integer :: initeWF           !init_elec_WF
  integer :: inithWF           !init_hole_WF
  !采用wannier轨道来指定
  real(kind=dp) :: initeEN     !init_elec_en
  real(kind=dp) :: inithEN     !init_hole_en
  !采用电子和空穴分别与humo和lumo之间的能量差
  integer :: initeES           !init_elec_esurface
  integer :: inithES           !init_hole_esurface
  
  !! use the exciton effects
  logical       :: L_exciton
  !是否考虑电子和空穴之间的激子效应。
  real(kind=dp) :: epsr
  
  !! which SH method do we use
  character(len=maxlen) ::  MSH  !"FSSH","LD-FSSH","SC-FSSH","mSC-FSSH","CC_FSSH"
  !采用哪种surface hopping方法来计算电子在势能面之间的跃迁几率。
  !! do we use decoherence method to account for decoherence correction
  logical::Ldecoherece
  !是否采用退相干修正跃迁几率
  !! which methoe of decoherence do we use
  character(len=maxlen) :: Tdecoherence !!"enbased"
  !采用的退相干修正的方法
  
  !! do we account for feedback of elec to core
  logical:: Lfeedback
  !是否考虑电子的运动对于核运动的影响
  
  
  !! MKL parallel
  integer :: mkl_threads
  !在调用mkl库时采用的线程数
  
  !!!输入文件问价夹路径
  character(len=maxlen) :: SHROOT_dir
  
  
end module control
                       