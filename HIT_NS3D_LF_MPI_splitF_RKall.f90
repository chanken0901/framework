!==================================================================================================================================
!                                                                                                                                 !
!  DNS (Implicit LES) Code based on Finite Difference methods for Compressible Temporal HIT Turbulence                            !
!    Temporally evolving HIT Turbulence is simulated in a periodic box.                                                           !
!                                                                                                                                 !
!    === Temporal advancement ===                                                                                                 !
!      5-stage 4th-order Runge-Kutta for Euler terms                                                                              !
!      1st-order Euler method for viscous and molecular diffusion terms                                                           !
!                                                                                                                                 !
!    === Spatial discretization ===                                                                                               !
!      8th order central difference schemes                                                                                       !
!      The accuracy of the central difference schemes is reduced near the boundary unless the periodic BC is used.                !
!                                                                                                                                 !
!                                                                                                                                 !
!    -------------------------------------------------------------------------------------------------------------------------    !
!    This code is based on the finite difference code in                                                                          !
!       Z.  Wang, Y. Lv, P. He, J. Zhou, K. Cen                                                                                   !
!       Fully explicit implementation of direct numerical simulation for a transientnear-field methane/air diffusion jet flame    !
!       Computers & Fluids, 32 (2009)                                                                                             !
!    -------------------------------------------------------------------------------------------------------------------------    !
!                                                                                                                                 !
!                                                                                                                                 !
!    The temporally files are saved in "temp" directory.                                                                          !
!      -> The files in temp are used for restart the simulations.                                                                 !
!                                                                                                                                 !
!    The output files are saved in "output" directory.                                                                            !
!                                                                                                                                 !
!    The 3D filed files are saved in one file for Ityp_Save = 1                                                                   !
!                                 in (n_slab*n_div) files for Ityp_Save = 1                                                       !
!                                                                                                                                 !
!    [*** NOTE ***]                                                                                                               !
!    For Ityp_Save = 1, the number of MPI processes must be larger than (n_slab*6), and                                           !
!                       Nx-1 must be represented by C*n_slab*n_div with integer C.                                                !
!                                                                                                                                 !
!                                                                                  Programed by T.Watanabe                        !
!                                                                                                                                 !
!==================================================================================================================================
!==== [  Main Module  ] ===========================================================================================================
!==================================================================================================================================

      Module Definition

        Implicit none 

        !--- Set computation time for SX9 ---
        Real*4,Parameter :: Time_limit = 23.5d0 ! 68.0d0              !Real Time (hour)
        Real*4,Parameter :: CPU_limit  = 10000000000.05d0     !CPU_time 
        Integer time_1, time_2, time_rate ; Common /time1/ time_1, time_2, time_rate
        Integer time_max, time_diff       ; Common /time2/ time_max, time_diff
        Integer count_round               ; Common /time3/ count_round
        Integer time_pre                  ; Common /time4/ time_pre

        !----------------------------------------------------------------------------------------
        Integer,Parameter :: Nx = 129                
        Integer,Parameter :: Ny = 129                
        Integer,Parameter :: Nz = 129                

        Integer,Parameter :: Nx_old       = 129   !ê¸å`ï‚äÆÇ∑ÇÈëOÇÃäiéqì_êî
        Integer,Parameter :: Ny_old       = 129   !ê¸å`ï‚äÆÇ∑ÇÈëOÇÃäiéqì_êî
        Integer,Parameter :: Nz_old       = 129   !ê¸å`ï‚äÆÇ∑ÇÈëOÇÃäiéqì_êî
        Integer,Parameter :: Num_rank_old = 48    !ê¸å`ï‚äÆÇ∑ÇÈëOÇÃMPIÉvÉçÉZÉX

        Integer,Parameter :: nbc   = 5-1 !ã´äEÇ…égÇ§ÉZÉãÇÃêî-1

        Integer,Parameter :: Ityp         =   0       ! [énÇﬂÇ©ÇÁ:Ityp = 0, ìríÜÇ©ÇÁ:Ityp = 1, ê¸å`ï‚äÆÇµÇƒêVÇµÇ¢DNSÅFItyp = 2]
        Integer,Parameter :: Ityp_in      =   0       ! [èoå˚ó¨ë¨åvéZÇ∑ÇÈ  :Ityp = 0]
                                                      ! [èoå˚ó¨ë¨Çì«çûÇﬁ  :Ityp = 1]

        Integer,Parameter :: Ityp_inG     =   1       ! [DNSÇ∆ìØÇ∂äiéq        :Ityp = 0]
                                                      ! [èâä˙èåèê›íËópÇÃäiéq :Ityp = 1]

        Integer,Parameter :: Ityp_inS     =   1       ! [generated by solution of Diss equtation : Ityp = 0]
                                                      ! [generated by IFFT                       : Ityp = 1]

        Integer,Parameter :: Ityp_Spec    =   1       ! [generated by Johnsen : Ityp = 0]
                                                      ! [generated by Pope    : Ityp = 1]

        Integer,Parameter :: Ityp_Inter   =   1       ! [generated by Johnsen :Ityp = 0]
    
        Integer,Parameter :: Ityp_dt      =   0       ! [CFLèåè âπë¨:Ityp = 0, ó¨ë¨:Ityp = 1]

        Integer,Parameter :: Ityp_RK      =   2       ! [RK5 for Euler terms :Ityp = 1 ]
                                                      ! [RK5 for   all terms :Ityp = 2 ]

        Integer,Parameter :: Ityp_data    =   2       ! [Small data:Ityp = 1                                                ]
                                                      ! [One file is created for one MPI process:Ityp=2                     ]
                                                      ! [ Ityp_data is for subroutine Itypwrite/Inputdata                   ]

        Integer,Parameter :: Ityp_Save    =   2       ! [Small data                                :Ityp = 1                ]
                                                      ! [Slices of data is saved into each file    :Ityp = 2                ]
                                                      ! [ Ityp_Save is for write_xyz                                        ]

        Integer,Parameter :: Ityp_Force   =   2       ! [Without Forcing                           :Ityp = 0                ]
                                                      ! [With Linear Forcing(Meneveau 2005)        :Ityp = 1                ]
                                                      ! [With Linear Forcing(Livescu  2010)        :Ityp = 2                ]

        Integer,Parameter :: Ityp_Press   =   1       ! [Constant Pressure                             :Ityp = 0            ]
                                                      ! [Initial condition has pressure fluctuations   :Ityp = 1            ]


        Integer,Parameter :: n_div        =   4       ! [Data is saved in n_div slices used when Ityp_Save = 3               ]
        Integer,Parameter :: n_slab       =   2       ! [Field is devided into n_div*n_slab of slabs when Ityp_Save = 3      ]

        Integer,Parameter :: Ityp_data_in =   0       ! [Small data:Ityp = 0, Large data:Ityp = 1 for inflow velocity       ]

        !Integer,Parameter :: ifl2      =    1       ! Statistics on the center line are saved every ifl2 steps
        !Integer,Parameter :: ifl3      =   1       ! inst. field on xy plane is saved every ifl3 steps ifl3 = C*ifl2
        !Integer,Parameter :: ifl4      =   1       ! Lateral profiles of statistics are saved every ifl4 steps
        !Integer,Parameter :: ifl5      =   1       ! Lateral profiles of statistics are saved every ifl4 steps
        !Integer,Parameter :: nlast     = 1       ! ëSÉXÉeÉbÉvêî
        Integer,Parameter :: ifl2      =  100         ! Statistics on the center line are saved every ifl2 steps
        Integer,Parameter :: ifl3      =  10000       ! inst. field on xy plane is saved every ifl3 steps ifl3 = C*ifl2
        Integer,Parameter :: ifl4      =  100         ! Lateral profiles of statistics are saved every ifl4 steps
        Integer,Parameter :: ifl5      =  10000       ! Lateral profiles of statistics are saved every ifl4 steps
        Integer,Parameter :: ifl6      =  10000       ! Lateral profiles of statistics are saved every ifl4 steps
        Integer,Parameter :: nlast     = 999999       ! ëSÉXÉeÉbÉvêî


        !--- èuéûílï€ë∂óp ---
        Real*4,Parameter :: T_inst1   = 10000.0d0    !    6.0d0 
        Real*4,Parameter :: T_inst2   = 10000.0d0    !   12.0d0 
        Real*4,Parameter :: T_inst3   =  1.0d0                  
        Real*4,Parameter :: T_inst4   =  2.5d0                  
        Real*4,Parameter :: T_inst5   =  3.0d0                  
        Real*4,Parameter :: T_inst6   =  5.0d0                  
        Real*4,Parameter :: T_last    = 10000.0d0               

        ! --- Non-dimentional parameters --- 
        ! parameters for initial HIT         
        Real*4,Parameter :: Re_ini     = 30.0d0                  ! Reynolds number based on the initial rms velocity and characteristics length
        !Real*4,Parameter :: Re_ini     =  51.99d0                  ! Reynolds number based on the initial rms velocity and characteristics length
        Real*4,Parameter :: Ma_ini     =   0.1d0                  ! Mach number based on the initial rms velocity 
        ! parameters for IFFT with Johnsen's spectrum
        Real*4,Parameter :: L_lambda   =    2.0d0                  ! target ratio of dilatational to solenoidal kinetic energy
        ! parameters for IFFT with Pope's spectrum
        !Real*4,parameter :: K_C   = 3.2d0
        !Real*4,Parameter :: C_L   = 4.6d0                          ! target ratio of dilatational to solenoidal kinetic energy
        !Real*4,Parameter :: C_eta = 0.40d0                         ! target ratio of dilatational to solenoidal kinetic energy
        !Real*4,Parameter :: p_0   = 2.0d0                          ! target ratio of dilatational to solenoidal kinetic energy
        !Real*4,Parameter :: beta  = 5.2d0                          ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,parameter :: K_C   = 1.5d0
        Real*4,Parameter :: C_L   = 6.78d0                         ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,Parameter :: C_eta = 0.40d0                         ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,Parameter :: p_0   = 2.0d0                          ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,Parameter :: beta  = 5.2d0                          ! target ratio of dilatational to solenoidal kinetic energy
        Real*4 :: f_L,f_eta                                        ! target ratio of dilatational to solenoidal kinetic energy
        Real*4 :: C_EP,Re_lambdaP,Re_LP,L_P,lambda_P,eta_P,E_P,K_P ! target ratio of dilatational to solenoidal kinetic energy
        ! parameters for Linear Forcing
        Real*4,Parameter :: Ed_Es      =    0.0d0                  ! target ratio of dilatational to solenoidal kinetic energy
        Real*4           :: Eta_target                             ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,Parameter :: E_target   =    2.0d0                  ! target ratio of dilatational to solenoidal kinetic energy
        Real*4,Parameter :: Coeff_Eta  =    1.3d0                  ! target ratio of dilatational to solenoidal kinetic energy
        Real*4 :: Ed,Es                                            ! dilatational and solenoidal kinetic energy

        Real*4,Parameter :: Pr     =   0.71d0                      ! Prandtl number
        Real*4,Parameter :: Sc     =    1.0d0                      ! Schmidt number
        Real*4,Parameter :: gh     =    1.4d0                      ! î‰îMî‰gamma (heat capacity ratio)

        !=== Physical Parameters ===
        Real*4,Parameter :: P_inf  = 101300.0d0                    ! Pressure in atmosphere (Back Ground) 1.013d0*(10.0d0**(5.0d0)) 
        Real*4,Parameter :: T_inf  = 300.0d0                       ! Temperature (Back Ground)
        Real*4,Parameter :: R_gas  = 287.0d0                       ! Gas constant [J/(kg*K)]
        Real*4,Parameter :: mu_T0  = 0.00001724d0                  ! Shear viscosity at T0 1.724d0*(10.0d0**(-5.0d0)) 

        Real*4 U_ini   ! Initial rms velocity [m/s]
        Real*4 L_ini   ! Initial turbulence lenght scale [m] 

        !--- Reference value---
        Real*4  U_ref    ! Jet velocity at the nozzle ó¨ì¸ïîó¨ë¨
        Real*4  L_ref    ! Nozzle width
        Real*4 Tm_ref    ! Reference time
        Real*4 Rh_ref    ! Density in ambient flows
        Real*4  P_ref    ! 
        Real*4  R_ref    ! Gas constant
        Real*4  T_ref    ! 
        Real*4 Cp_ref    ! íËà≥î‰îM
        Real*4 Cv_ref    ! íËêœî‰îM
        Real*4 mu_ref    ! îSìx
        Real*4 kT_ref    ! îMägéUåWêî

        Real*4  cRT_UU   ! R_ref*T_ref/U_ref/U_ref    ! à≥óÕïœä∑
        Real*4  cUU_RT   ! U_ref*U_ref/R_ref/T_ref    ! à≥óÕïœä∑


        !---à»â∫ *_refÇ≈ñ≥éüå≥âª
        Real*4 Re
        Real*4 Ma     ! Reference Mach number
        Real*4 Re_ini_test  !--- Reynolds number ---
        Real*4 Ma_ini_test  !--- Mach number---

        !=== Parmeters ===
        Real*4,Parameter :: pi   = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825d0
        Real*8,Parameter :: Lx   = 4.0d0
        Real*8,Parameter :: Ly   = 4.0d0
        Real*8,Parameter :: Lz   = 4.0d0

        Real*8,Parameter  :: dtime0  = 0.5d0
        Real*4,Parameter  :: Cour    = 0.6d0

        !--- For initial velocity fields ---
        Real*8,Parameter :: rms_in = 1.0d0                                   ! Nondimemsional rms velocty 
        Real*8,Parameter :: L_in   = 1.0d0                                   ! Nondimemsional characteristics length scale
        !Real*8,Parameter :: L_in   = 0.01299                                 ! Nondimemsional characteristics length scale
        Real*8,Parameter :: D_in   = L_in*L_in/(2.0d0*pi*dtime0*400.0d0)     ! ïœìÆê∂ê¨Ç≈égÇ§ägéUï˚íˆéÆÇÃägéUåWêî

        Integer,Parameter :: Nx_in = 129
        Integer,Parameter :: Ny_in = 129
        Integer,Parameter :: Nz_in = 129


        !=== ï®óùó íËêî ===
        ! SutherlandÇÃéÆ
        Real*4 T0       ! äÓèÄâ∑ìx
        Real*4 mu0      ! îSê´åWêî at T0 
        Real*4 C_su     ! SutherlandÇÃéÆ

        Real*4 P_bg     ! é¸àÕó¨à≥óÕ
        Real*4 T_bg     ! é¸àÕó¨â∑ìx

        Real*4 Ys_max 
        Real*4 Ys_min 

        !=== Data array ===
        !--- ì∆óßïœêî ---
        Real*4, dimension(:,:,:),allocatable ::    Rh           !ñßìx
        Real*4, dimension(:,:,:),allocatable ::    Ur, Vr, Wr   !ñßìx*ë¨ìx
        Real*4, dimension(:,:,:),allocatable ::    Tr           !ñßìx*â∑ìx

        Real*4, dimension(:,:,:),allocatable ::    Ysr          !ñßìx*Passive Scalar

        !---è]ëÆïœêî---
        Real*4, dimension(:,:,:),allocatable ::    U, V, W      !ë¨ìx
        Real*4, dimension(:,:,:),allocatable ::    T            !â∑ìx
        Real*4, dimension(:,:,:),allocatable ::    Ys           !Passive Scalar

        Real*4, dimension(:,:,:),allocatable ::    Ud,Vd,Wd     !ë¨ìx
        Real*4, dimension(:,:,:),allocatable ::    Us,Vs,Ws     !ë¨ìx

        Real*4, dimension(:,:,:),allocatable ::    P            !ê√à≥
        Real*4, dimension(:,:,:),allocatable ::    Mu           !îSê´åWêî

        !---î˜ï™(x-dir)---
        Real*4, dimension(:,:,:),allocatable ::     dE_xdx
        Real*4, dimension(:,:,:),allocatable ::      dUrdx
        Real*4, dimension(:,:,:),allocatable ::      dRhdx
        Real*4, dimension(:,:,:),allocatable ::       dUdx
        Real*4, dimension(:,:,:),allocatable ::      dUddx
        Real*4, dimension(:,:,:),allocatable ::      dUsdx
        Real*4, dimension(:,:,:),allocatable ::       dVdx
        Real*4, dimension(:,:,:),allocatable ::       dWdx
        Real*4, dimension(:,:,:),allocatable ::       dPdx
        Real*4, dimension(:,:,:),allocatable ::       dTdx
        Real*4, dimension(:,:,:),allocatable ::      dYsdx
        !---î˜ï™(y-dir)---
        Real*4, dimension(:,:,:),allocatable ::     dE_ydy
        Real*4, dimension(:,:,:),allocatable ::      dVrdy
        Real*4, dimension(:,:,:),allocatable ::      dRhdy
        Real*4, dimension(:,:,:),allocatable ::       dUdy
        Real*4, dimension(:,:,:),allocatable ::       dVdy
        Real*4, dimension(:,:,:),allocatable ::      dVddy
        Real*4, dimension(:,:,:),allocatable ::      dVsdy
        Real*4, dimension(:,:,:),allocatable ::       dWdy
        Real*4, dimension(:,:,:),allocatable ::       dPdy
        Real*4, dimension(:,:,:),allocatable ::       dTdy
        Real*4, dimension(:,:,:),allocatable ::      dYsdy
        !---î˜ï™(z-dir)---
        Real*4, dimension(:,:,:),allocatable ::     dE_zdz
        Real*4, dimension(:,:,:),allocatable ::      dWrdz
        Real*4, dimension(:,:,:),allocatable ::      dRhdz
        Real*4, dimension(:,:,:),allocatable ::       dUdz
        Real*4, dimension(:,:,:),allocatable ::       dVdz
        Real*4, dimension(:,:,:),allocatable ::       dWdz
        Real*4, dimension(:,:,:),allocatable ::      dWddz
        Real*4, dimension(:,:,:),allocatable ::      dWsdz
        Real*4, dimension(:,:,:),allocatable ::       dPdz
        Real*4, dimension(:,:,:),allocatable ::       dTdz
        Real*4, dimension(:,:,:),allocatable ::      dYsdz

        Real*4, dimension(:,:,:),allocatable ::    c1_Rh           !1.0d0/ñßìx

        !---Flux---
        Real*4, dimension(:,:,:),allocatable ::    E_x, E_y, E_z
        !---tau_ij---
        Real*4, dimension(:,:,:),allocatable ::    T_xx, T_xy, T_xz
        Real*4, dimension(:,:,:),allocatable ::    T_yx, T_yy, T_yz
        Real*4, dimension(:,:,:),allocatable ::    T_zx, T_zy, T_zz

        !===éûä‘êiçs===
        Real*4, dimension(:,:,:),allocatable ::    dRh , dTr
        Real*4, dimension(:,:,:),allocatable ::    dRh0, dTr0

        Real*4, dimension(:,:,:),allocatable ::    dUr , dVr , dWr
        Real*4, dimension(:,:,:),allocatable ::    dUr0, dVr0, dWr0
        !---éûä‘êiçs(îSê´çÄ)---
        Real*4, dimension(:,:,:),allocatable ::    dRhv, dTrv
        Real*4, dimension(:,:,:),allocatable ::    dUrv, dVrv, dWrv

        !===éûä‘êiçs(ÉXÉJÉâÅ[)===
        Real*4, dimension(:,:,:),allocatable ::    dYsr
        Real*4, dimension(:,:,:),allocatable ::    dYsr0
        !---éûä‘êiçs(ägéUçÄ)---
        Real*4, dimension(:,:,:),allocatable ::    dYsrv

        !===ó¨ì¸ã´äEèåèê›íË===
        Real*8, dimension(:,:,:),allocatable :: Uf_in, Vf_in, Wf_in
        Real*8, dimension(:,:,:),allocatable ::   dUf,   dVf,   dWf

        Real*4 Rand_x, Rand_y, Rand_z
        Integer Nt_in
        Integer N_step
        Integer ifl_i


        Real*4 Ma_max !ç≈ëÂÉ}ÉbÉnêî


        !=== Linear Forcing by Rosales and Meneveau in POF (2005) ===
        Real*4 B_Ceff
        Real*4 B1_Ceff,B2_Ceff

        !=== Split Linear Forcing by Livescu in POF (2010) ===
        Real*4 C_s,C_d
 

        !=== Coordinate ===
        Real*4 x (-4:Nx+4)
        Real*4 y (-4:Ny+4)
        Real*4 z (-4:Nz+4)

        Real*4  dx, c1_dx
        Real*4  dy, c1_dy
        Real*4  dz, c1_dz
        Real*4 dxyz

        !=== Coordinate ===
        Real*4 x_in (0:Nx_in)
        Real*4 y_in (0:Ny_in)
        Real*4 z_in (0:Nz_in)

        Real*4  dx_in, c1_dx_in
        Real*4  dy_in, c1_dy_in
        Real*4  dz_in, c1_dz_in

        !=== éûä‘ ===
        Integer nt, nt0
        Real*4  time
        Real*4  time_0
        Real*4  dtime
        Real*4  dtime_c

        !---Runge Kutta---
        Integer  Num_RK
        Real*4  RK_A(5), RK_B(5)
        Real*4,Parameter ::   RK_A1 = 0.0d0
        Real*4,Parameter ::   RK_A2 = -   73838075589.0d0 /  153778244506.0d0
        Real*4,Parameter ::   RK_A3 = - 6194124222391.0d0 / 4410992767914.0d0
        Real*4,Parameter ::   RK_A4 = -  875753879689.0d0 /  434298951106.0d0
        Real*4,Parameter ::   RK_A5 = -  375951423649.0d0 /  355864889808.0d0

        Real*4,Parameter ::   RK_B1 =     69922298306.0d0 /  679754813293.0d0
        Real*4,Parameter ::   RK_B2 =    151102391305.0d0 /  203957027379.0d0
        Real*4,Parameter ::   RK_B3 =    168159696081.0d0 /  226431017777.0d0
        Real*4,Parameter ::   RK_B4 =    189406283338.0d0 /  403426599621.0d0
        Real*4,Parameter ::   RK_B5 =    184809359603.0d0 /  982122979183.0d0


        !=== Data Write on xy plane ===
        Real*4   Rh_xy (0:Nx,0:Ny)
        Real*4    U_xy (0:Nx,0:Ny)
        Real*4    V_xy (0:Nx,0:Ny)
        Real*4    W_xy (0:Nx,0:Ny)
        Real*4   Te_xy (0:Nx,0:Ny)
        Real*4    P_xy (0:Nx,0:Ny)
        Real*4 divU_xy (0:Nx,0:Ny)
        Real*4   Ys_xy (0:Nx,0:Ny)
        Real*4 Enst_xy (0:Nx,0:Ny)
        Real*4   Ma_xy (0:Nx,0:Ny)

        !=== Data Write on xz plane ===
        Real*4   Rh_xz (0:Nx,0:Nz)
        Real*4    U_xz (0:Nx,0:Nz)
        Real*4    V_xz (0:Nx,0:Nz)
        Real*4    W_xz (0:Nx,0:Nz)
        Real*4   Te_xz (0:Nx,0:Nz)
        Real*4    P_xz (0:Nx,0:Nz)
        Real*4 divU_xz (0:Nx,0:Nz)
        Real*4   Ys_xz (0:Nx,0:Nz)
        Real*4 Enst_xz (0:Nx,0:Nz)
        Real*4   Ma_xz (0:Nx,0:Ny)

        !=== Data Write on yz plane ===
        Real*4   Rh_yz (0:Ny,0:Nz)
        Real*4    U_yz (0:Ny,0:Nz)
        Real*4    V_yz (0:Ny,0:Nz)
        Real*4    W_yz (0:Ny,0:Nz)
        Real*4   Te_yz (0:Ny,0:Nz)
        Real*4    P_yz (0:Ny,0:Nz)
        Real*4 divU_yz (0:Ny,0:Nz)
        Real*4   Ys_yz (0:Ny,0:Nz)
        Real*4 Enst_yz (0:Ny,0:Nz)
        Real*4   Ma_yz (0:Nx,0:Ny)


        !=== Data Array_Center ===
        Real*4, dimension(:,:,:),allocatable ::    Rhc
        Real*4, dimension(:,:,:),allocatable ::    Urc
        Real*4, dimension(:,:,:),allocatable ::    Vrc
        Real*4, dimension(:,:,:),allocatable ::    Wrc
        Real*4, dimension(:,:,:),allocatable ::    Trc
        Real*4, dimension(:,:,:),allocatable ::   Ysrc

        !=== [ Statistics for temporally evolving HIT turbulence ]  ===
        !=== Parameters for HIT Turbulence case ===
        Real*8 ,Parameter :: Lyg = 1.0d0
        Real*8 ,Parameter :: Lzg = 1.0d0
        Integer,Parameter :: Nyg = 61
        Integer,Parameter :: Nzg = 61

        Real*4  yg(-2:Nyg+2), zg(-2:Nzg+2)
        Real*4  dyg, dzg
        Integer j_yg(-2:Ny+2), k_zg(-2:Nz+2)
        Integer Num_Gy 
        Integer Num_Gz 
        Integer Num_Gy2
        Integer Num_Gz2

        Integer ig, jg, kg


        !=== Conventional Statistics ===
        Integer Count_yz (0:Nyg,0:Nzg)

        Real*4 RhsumT    (0:Nyg,0:Nzg)
        Real*4 Rh2sumT   (0:Nyg,0:Nzg)
        Real*4 Rh3sumT   (0:Nyg,0:Nzg)
        Real*4 Rh4sumT   (0:Nyg,0:Nzg)

        Real*4 TsumT     (0:Nyg,0:Nzg)
        Real*4 PsumT     (0:Nyg,0:Nzg)
        Real*4 UsumT     (0:Nyg,0:Nzg)
        Real*4 VsumT     (0:Nyg,0:Nzg)
        Real*4 WsumT     (0:Nyg,0:Nzg)
        Real*4 YssumT    (0:Nyg,0:Nzg)

        Real*4 TTsumT    (0:Nyg,0:Nzg)
        Real*4 PPsumT    (0:Nyg,0:Nzg)
        Real*4 UUsumT    (0:Nyg,0:Nzg)
        Real*4 VVsumT    (0:Nyg,0:Nzg)
        Real*4 WWsumT    (0:Nyg,0:Nzg)
        Real*4 Ys2sumT   (0:Nyg,0:Nzg)

        Real*4 T3sumT    (0:Nyg,0:Nzg)
        Real*4 P3sumT    (0:Nyg,0:Nzg)
        Real*4 U3sumT    (0:Nyg,0:Nzg)
        Real*4 V3sumT    (0:Nyg,0:Nzg)
        Real*4 W3sumT    (0:Nyg,0:Nzg)
        Real*4 Ys3sumT   (0:Nyg,0:Nzg)

        Real*4 T4sumT    (0:Nyg,0:Nzg)
        Real*4 P4sumT    (0:Nyg,0:Nzg)
        Real*4 U4sumT    (0:Nyg,0:Nzg)
        Real*4 V4sumT    (0:Nyg,0:Nzg)
        Real*4 W4sumT    (0:Nyg,0:Nzg)
        Real*4 Ys4sumT   (0:Nyg,0:Nzg)


        Real*4 dUdxsumT  (0:Nyg,0:Nzg)
        Real*4 dVdxsumT  (0:Nyg,0:Nzg)
        Real*4 dWdxsumT  (0:Nyg,0:Nzg)
        Real*4 dUdysumT  (0:Nyg,0:Nzg)
        Real*4 dVdysumT  (0:Nyg,0:Nzg)
        Real*4 dWdysumT  (0:Nyg,0:Nzg)
        Real*4 dUdzsumT  (0:Nyg,0:Nzg)
        Real*4 dVdzsumT  (0:Nyg,0:Nzg)
        Real*4 dWdzsumT  (0:Nyg,0:Nzg)

        Real*4 dUdx2sumT (0:Nyg,0:Nzg)
        Real*4 dVdx2sumT (0:Nyg,0:Nzg)
        Real*4 dWdx2sumT (0:Nyg,0:Nzg)
        Real*4 dUdy2sumT (0:Nyg,0:Nzg)
        Real*4 dVdy2sumT (0:Nyg,0:Nzg)
        Real*4 dWdy2sumT (0:Nyg,0:Nzg)
        Real*4 dUdz2sumT (0:Nyg,0:Nzg)
        Real*4 dVdz2sumT (0:Nyg,0:Nzg)
        Real*4 dWdz2sumT (0:Nyg,0:Nzg)

        Real*4 dUdx3sumT (0:Nyg,0:Nzg)
        Real*4 dVdx3sumT (0:Nyg,0:Nzg)
        Real*4 dWdx3sumT (0:Nyg,0:Nzg)
        Real*4 dUdy3sumT (0:Nyg,0:Nzg)
        Real*4 dVdy3sumT (0:Nyg,0:Nzg)
        Real*4 dWdy3sumT (0:Nyg,0:Nzg)
        Real*4 dUdz3sumT (0:Nyg,0:Nzg)
        Real*4 dVdz3sumT (0:Nyg,0:Nzg)
        Real*4 dWdz3sumT (0:Nyg,0:Nzg)

        Real*4 dUdx4sumT (0:Nyg,0:Nzg)
        Real*4 dVdx4sumT (0:Nyg,0:Nzg)
        Real*4 dWdx4sumT (0:Nyg,0:Nzg)
        Real*4 dUdy4sumT (0:Nyg,0:Nzg)
        Real*4 dVdy4sumT (0:Nyg,0:Nzg)
        Real*4 dWdy4sumT (0:Nyg,0:Nzg)
        Real*4 dUdz4sumT (0:Nyg,0:Nzg)
        Real*4 dVdz4sumT (0:Nyg,0:Nzg)
        Real*4 dWdz4sumT (0:Nyg,0:Nzg)

        Real*4 dPdxsumT  (0:Nyg,0:Nzg)
        Real*4 dPdysumT  (0:Nyg,0:Nzg)
        Real*4 dPdzsumT  (0:Nyg,0:Nzg)
        Real*4 dTdxsumT  (0:Nyg,0:Nzg)
        Real*4 dTdysumT  (0:Nyg,0:Nzg)
        Real*4 dTdzsumT  (0:Nyg,0:Nzg)
        Real*4 dYsdxsumT (0:Nyg,0:Nzg)
        Real*4 dYsdysumT (0:Nyg,0:Nzg)
        Real*4 dYsdzsumT (0:Nyg,0:Nzg)

        Real*4 dPdx2sumT  (0:Nyg,0:Nzg)
        Real*4 dPdy2sumT  (0:Nyg,0:Nzg)
        Real*4 dPdz2sumT  (0:Nyg,0:Nzg)
        Real*4 dTdx2sumT  (0:Nyg,0:Nzg)
        Real*4 dTdy2sumT  (0:Nyg,0:Nzg)
        Real*4 dTdz2sumT  (0:Nyg,0:Nzg)
        Real*4 dYsdx2sumT (0:Nyg,0:Nzg)
        Real*4 dYsdy2sumT (0:Nyg,0:Nzg)
        Real*4 dYsdz2sumT (0:Nyg,0:Nzg)

        Real*4 dPdx3sumT  (0:Nyg,0:Nzg)
        Real*4 dPdy3sumT  (0:Nyg,0:Nzg)
        Real*4 dPdz3sumT  (0:Nyg,0:Nzg)
        Real*4 dTdx3sumT  (0:Nyg,0:Nzg)
        Real*4 dTdy3sumT  (0:Nyg,0:Nzg)
        Real*4 dTdz3sumT  (0:Nyg,0:Nzg)
        Real*4 dYsdx3sumT (0:Nyg,0:Nzg)
        Real*4 dYsdy3sumT (0:Nyg,0:Nzg)
        Real*4 dYsdz3sumT (0:Nyg,0:Nzg)

        Real*4 dPdx4sumT  (0:Nyg,0:Nzg)
        Real*4 dPdy4sumT  (0:Nyg,0:Nzg)
        Real*4 dPdz4sumT  (0:Nyg,0:Nzg)
        Real*4 dTdx4sumT  (0:Nyg,0:Nzg)
        Real*4 dTdy4sumT  (0:Nyg,0:Nzg)
        Real*4 dTdz4sumT  (0:Nyg,0:Nzg)
        Real*4 dYsdx4sumT (0:Nyg,0:Nzg)
        Real*4 dYsdy4sumT (0:Nyg,0:Nzg)
        Real*4 dYsdz4sumT (0:Nyg,0:Nzg)

        Real*4 dRhdxsumT  (0:Nyg,0:Nzg)
        Real*4 dRhdysumT  (0:Nyg,0:Nzg)
        Real*4 dRhdzsumT  (0:Nyg,0:Nzg)

        Real*4 dRhdx2sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdy2sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdz2sumT  (0:Nyg,0:Nzg)

        Real*4 dRhdx3sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdy3sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdz3sumT  (0:Nyg,0:Nzg)

        Real*4 dRhdx4sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdy4sumT  (0:Nyg,0:Nzg)
        Real*4 dRhdz4sumT  (0:Nyg,0:Nzg)


        Real*4 EpsisumT     (0:Nyg,0:Nzg)
        Real*4 Epsi_s_sumT  (0:Nyg,0:Nzg) ! \mu^(\Omega_i)(\Omega_i)
        Real*4 Epsi_c_sumT  (0:Nyg,0:Nzg) ! (4/3)\mu^(S_{kk})^2
        Real*4 Epsi_ih_sumT (0:Nyg,0:Nzg) ! 2\mu[(dU_{i}/dx_{j})(dU_{j}/dx_{i})-S_{kk}S_{kk}]


        Real*4 UVsumT       (0:Nyg,0:Nzg)
        Real*4 VWsumT       (0:Nyg,0:Nzg)
        Real*4 WUsumT       (0:Nyg,0:Nzg)

        Real*4 UYssumT      (0:Nyg,0:Nzg)
        Real*4 VYssumT      (0:Nyg,0:Nzg)
        Real*4 WYssumT      (0:Nyg,0:Nzg)

        Real*4 UPsumT       (0:Nyg,0:Nzg)
        Real*4 VPsumT       (0:Nyg,0:Nzg)
        Real*4 WPsumT       (0:Nyg,0:Nzg)

        Real*4 dUdxPsumT    (0:Nyg,0:Nzg)
        Real*4 dVdyPsumT    (0:Nyg,0:Nzg)
        Real*4 dWdzPsumT    (0:Nyg,0:Nzg)

        Real*4    nusumT    (0:Nyg,0:Nzg)
        Real*4    musumT    (0:Nyg,0:Nzg)
        Real*4   mu2sumT    (0:Nyg,0:Nzg)
        Real*4  divUsumT    (0:Nyg,0:Nzg)
        Real*4 divU2sumT    (0:Nyg,0:Nzg)



        !--- Correlation ---
        Real*4   Cuu_C (0:Nx)
        Real*4   Cvv_C (0:Nx)
        Real*4   Cww_C (0:Nx)
        Real*4   Cpp_C (0:Nx)
        Real*4   Crr_C (0:Nx)
        Real*4   Ctt_C (0:Nx)
        Real*4   Cyy_C (0:Nx)

        Real*4   Um_C
        Real*4   Vm_C
        Real*4   Wm_C
        Real*4   Pm_C
        Real*4   Rm_C
        Real*4   Tm_C
        Real*4   Ym_C

        Real*4   Cuu_C0
        Real*4   Cvv_C0
        Real*4   Cww_C0
        Real*4   Cpp_C0
        Real*4   Crr_C0
        Real*4   Ctt_C0
        Real*4   Cyy_C0

        Real*4   Luu_C

        !===xyzïΩãœ===
        Real*4 RhsumT_xyz    
        Real*4 Rh2sumT_xyz   
        Real*4 Rh3sumT_xyz   
        Real*4 Rh4sumT_xyz   

        Real*4 rkensumT_xyz  

        Real*4 TsumT_xyz     
        Real*4 PsumT_xyz     
        Real*4 UsumT_xyz     
        Real*4 VsumT_xyz     
        Real*4 WsumT_xyz     
        Real*4 YssumT_xyz    
        Real*4 UVWsumT_xyz   

        Real*4 TTsumT_xyz    
        Real*4 PPsumT_xyz    
        Real*4 UUsumT_xyz    
        Real*4 VVsumT_xyz    
        Real*4 WWsumT_xyz    
        Real*4 Ys2sumT_xyz   
        Real*4 UVW2sumT_xyz  

        Real*4 T3sumT_xyz    
        Real*4 P3sumT_xyz    
        Real*4 U3sumT_xyz    
        Real*4 V3sumT_xyz    
        Real*4 W3sumT_xyz    
        Real*4 Ys3sumT_xyz   

        Real*4 T4sumT_xyz    
        Real*4 P4sumT_xyz    
        Real*4 U4sumT_xyz    
        Real*4 V4sumT_xyz    
        Real*4 W4sumT_xyz    
        Real*4 Ys4sumT_xyz   

!=======================================
        Real*4 RUUsumT_xyz     
        Real*4 RVVsumT_xyz     
        Real*4 RWWsumT_xyz     
        Real*4 CTKEsumT_xyz    

        Real*4 PDsumT_xyz      
        Real*4 ForcesumT_xyz   
!=======================================

        Real*4 dUdxsumT_xyz  
        Real*4 dVdxsumT_xyz  
        Real*4 dWdxsumT_xyz  
        Real*4 dUdysumT_xyz  
        Real*4 dVdysumT_xyz  
        Real*4 dWdysumT_xyz  
        Real*4 dUdzsumT_xyz  
        Real*4 dVdzsumT_xyz  
        Real*4 dWdzsumT_xyz  

        Real*4 dUdx2sumT_xyz 
        Real*4 dVdx2sumT_xyz 
        Real*4 dWdx2sumT_xyz 
        Real*4 dUdy2sumT_xyz 
        Real*4 dVdy2sumT_xyz 
        Real*4 dWdy2sumT_xyz 
        Real*4 dUdz2sumT_xyz 
        Real*4 dVdz2sumT_xyz 
        Real*4 dWdz2sumT_xyz 

        Real*4 dUdx3sumT_xyz 
        Real*4 dVdx3sumT_xyz 
        Real*4 dWdx3sumT_xyz 
        Real*4 dUdy3sumT_xyz 
        Real*4 dVdy3sumT_xyz 
        Real*4 dWdy3sumT_xyz 
        Real*4 dUdz3sumT_xyz 
        Real*4 dVdz3sumT_xyz 
        Real*4 dWdz3sumT_xyz 

        Real*4 dUdx4sumT_xyz 
        Real*4 dVdx4sumT_xyz 
        Real*4 dWdx4sumT_xyz 
        Real*4 dUdy4sumT_xyz 
        Real*4 dVdy4sumT_xyz 
        Real*4 dWdy4sumT_xyz 
        Real*4 dUdz4sumT_xyz 
        Real*4 dVdz4sumT_xyz 
        Real*4 dWdz4sumT_xyz 

        Real*4 dPdxsumT_xyz  
        Real*4 dPdysumT_xyz  
        Real*4 dPdzsumT_xyz  
        Real*4 dTdxsumT_xyz  
        Real*4 dTdysumT_xyz  
        Real*4 dTdzsumT_xyz  
        Real*4 dYsdxsumT_xyz 
        Real*4 dYsdysumT_xyz 
        Real*4 dYsdzsumT_xyz 

        Real*4 dPdx2sumT_xyz 
        Real*4 dPdy2sumT_xyz 
        Real*4 dPdz2sumT_xyz 
        Real*4 dTdx2sumT_xyz 
        Real*4 dTdy2sumT_xyz 
        Real*4 dTdz2sumT_xyz 
        Real*4 dYsdx2sumT_xyz
        Real*4 dYsdy2sumT_xyz
        Real*4 dYsdz2sumT_xyz

        Real*4 dPdx3sumT_xyz 
        Real*4 dPdy3sumT_xyz 
        Real*4 dPdz3sumT_xyz 
        Real*4 dTdx3sumT_xyz 
        Real*4 dTdy3sumT_xyz 
        Real*4 dTdz3sumT_xyz 
        Real*4 dYsdx3sumT_xyz
        Real*4 dYsdy3sumT_xyz
        Real*4 dYsdz3sumT_xyz

        Real*4 dPdx4sumT_xyz 
        Real*4 dPdy4sumT_xyz 
        Real*4 dPdz4sumT_xyz 
        Real*4 dTdx4sumT_xyz 
        Real*4 dTdy4sumT_xyz 
        Real*4 dTdz4sumT_xyz 
        Real*4 dYsdx4sumT_xyz
        Real*4 dYsdy4sumT_xyz
        Real*4 dYsdz4sumT_xyz

        Real*4 dRhdxsumT_xyz 
        Real*4 dRhdysumT_xyz 
        Real*4 dRhdzsumT_xyz 

        Real*4 dRhdx2sumT_xyz
        Real*4 dRhdy2sumT_xyz
        Real*4 dRhdz2sumT_xyz

        Real*4 dRhdx3sumT_xyz
        Real*4 dRhdy3sumT_xyz
        Real*4 dRhdz3sumT_xyz

        Real*4 dRhdx4sumT_xyz
        Real*4 dRhdy4sumT_xyz
        Real*4 dRhdz4sumT_xyz

        Real*4 EpsisumT_xyz     
        Real*4 Epsi_s_sumT_xyz   ! \mu^(\Omega_i)(\Omega_i)
        Real*4 Epsi_c_sumT_xyz   ! (4/3)\mu^(S_{kk})^2
        Real*4 Epsi_ih_sumT_xyz  ! 2\mu[(dU_{i}/dx_{j})(dU_{j}/dx_{i})-S_{kk}S_{kk}]

        Real*4 UVsumT_xyz       
        Real*4 VWsumT_xyz       
        Real*4 WUsumT_xyz       

        Real*4 UYssumT_xyz      
        Real*4 VYssumT_xyz      
        Real*4 WYssumT_xyz      

        Real*4 UPsumT_xyz       
        Real*4 VPsumT_xyz       
        Real*4 WPsumT_xyz       

        Real*4 dUdxPsumT_xyz    
        Real*4 dVdyPsumT_xyz    
        Real*4 dWdzPsumT_xyz    

        Real*4 UdPdxsumT_xyz    
        Real*4 VdPdysumT_xyz    
        Real*4 WdPdzsumT_xyz    

        Real*4 forcing_term     

        Real*4    nusumT_xyz    
        Real*4    musumT_xyz    
        Real*4   mu2sumT_xyz    
        Real*4  divUsumT_xyz    
        Real*4 divU2sumT_xyz    


        !--- Correlation ---
        Real*4   Cuu_xyz (0:Nx)
        Real*4   Cvv_xyz (0:Nx)
        Real*4   Cww_xyz (0:Nx)
        Real*4   Cpp_xyz (0:Nx)
        Real*4   Crr_xyz (0:Nx)
        Real*4   Ctt_xyz (0:Nx)
        Real*4   Cyy_xyz (0:Nx)

        Real*4   Cuu_xyz0
        Real*4   Cvv_xyz0
        Real*4   Cww_xyz0
        Real*4   Cpp_xyz0
        Real*4   Crr_xyz0
        Real*4   Ctt_xyz0
        Real*4   Cyy_xyz0

        Real*4   Luu_xyz


        !=== Filename ===
        Character( len = 255 ) :: fname_write

        Character( len = 255 ) :: fname_n_write
        Character( len = 255 ) :: fname_r_write
        Character( len = 255 ) :: fname_u_write
        Character( len = 255 ) :: fname_v_write
        Character( len = 255 ) :: fname_w_write
        Character( len = 255 ) :: fname_t_write
        Character( len = 255 ) :: fname_s_write

        Character( len = 255 ) :: fname_geo
        Character( len = 255 ) :: fname_scl
        Character( len = 255 ) :: fname_case


        !--- CPUéûä‘ÇÃåvë™ ---
        Real*4  Clock_Initial,  Clock_P0,  Clock_P1
        Integer  RClock_T1,  RClock_T2, Rt_rate, Rt_max, Rt_diff

        !=== Parameters ===
        Real*4,Parameter :: c1_pi      = 1.0d0 / pi
        Real*4,Parameter :: c1_Lx      = 1.0d0 / Lx
        Real*4,Parameter :: c1_Ly      = 1.0d0 / Ly
        Real*4,Parameter :: c1_Lz      = 1.0d0 / Lz

        Real*4  c1_Nx
        Real*4  c1_Nz
        Real*4  c1_NxNz
        Real*4  c1_NyNz
        Real*4  c1_NxNyNz

        Real*4 c1_T0
        Real*4 c1_Re
        Real*4 c1_Pr
        Real*4 c1_Sc
        Real*4 c1_gh

        Real*4,Parameter :: c1_2       = 1.0d0 /   2.0d0
        Real*4,Parameter :: c1_3       = 1.0d0 /   3.0d0
        Real*4,Parameter :: c1_4       = 1.0d0 /   4.0d0
        Real*4,Parameter :: c1_8       = 1.0d0 /   8.0d0
        Real*4,Parameter :: c3_8       = 3.0d0 /   8.0d0
        Real*4,Parameter :: c6_8       = 6.0d0 /   8.0d0
        Real*4,Parameter :: c9_8       = 9.0d0 /   8.0d0
        Real*4,Parameter :: c1_100     = 1.0d0 / 100.0d0
        Real*4,Parameter :: c1_256     = 1.0d0 / 256.0d0
        Real*4,Parameter :: c5_8       = 5.0d0 /   8.0d0

        Real*4,Parameter :: c2_3       = 2.0d0 /   3.0d0
        Real*4,Parameter :: c4_3       = 4.0d0 /   3.0d0

        Real*4,Parameter ::    c1_16   =    1.0d0 /  16.0d0
        Real*4,Parameter ::    c9_16   =    9.0d0 /  16.0d0
        Real*4,Parameter ::    c1_24   =    1.0d0 /  24.0d0
        Real*4,Parameter ::   c27_24   =   27.0d0 /  24.0d0
        Real*4,Parameter ::    c1_576  =    1.0d0 / 576.0d0
        Real*4,Parameter ::   c54_576  =   54.0d0 / 576.0d0
        Real*4,Parameter ::  c783_576  =  783.0d0 / 576.0d0
        Real*4,Parameter :: c1460_576  = 1460.0d0 / 576.0d0

        Real*4,Parameter ::    c5_9    =    5.0d0 /   9.0d0
        Real*4,Parameter ::    c8_15   =    8.0d0 /  15.0d0
        Real*4,Parameter ::   c15_16   =   15.0d0 /  16.0d0
        Real*4,Parameter ::  c153_128  =  153.0d0 / 128.0d0

        Real*4,Parameter ::    c3_4    =    3.0d0 /   4.0d0
        Real*4,Parameter ::    c2_15   =    2.0d0 /  15.0d0
        Real*4,Parameter ::    c5_12   =    5.0d0 /  12.0d0
        Real*4,Parameter ::   c17_60   =   17.0d0 /  60.0d0

        Real*4,Parameter ::  c150_256  =  150.0d0 / 256.0d0
        Real*4,Parameter ::   c25_256  =   25.0d0 / 256.0d0
        Real*4,Parameter ::    c3_256  =    3.0d0 / 256.0d0

        Real*4,Parameter :: c29_96     = 29.0d0 /   96.0d0
        Real*4,Parameter :: c37_160    = 37.0d0 /  160.0d0
        Real*4,Parameter ::  c3_40     =  3.0d0 /   40.0d0
        Real*4,Parameter ::  c5_24     =  5.0d0 /   24.0d0

        !--- Filter ---
        Real*4,Parameter ::   c1_1024  =    1.0d0 / 1024.0d0
        Real*4,Parameter ::  c10_1024  =   10.0d0 / 1024.0d0
        Real*4,Parameter ::  c45_1024  =   45.0d0 / 1024.0d0
        Real*4,Parameter :: c120_1024  =  120.0d0 / 1024.0d0
        Real*4,Parameter :: c210_1024  =  210.0d0 / 1024.0d0
        Real*4,Parameter :: c252_1024  =  252.0d0 / 1024.0d0

        Real*4,Parameter ::   c5_1024  =    5.0d0 / 1024.0d0
        Real*4,Parameter ::  c26_1024  =   26.0d0 / 1024.0d0
        Real*4,Parameter ::  c55_1024  =   55.0d0 / 1024.0d0
        Real*4,Parameter ::  c60_1024  =   60.0d0 / 1024.0d0
        Real*4,Parameter ::  c35_1024  =   35.0d0 / 1024.0d0
        Real*4,Parameter :: c126_1024  =  126.0d0 / 1024.0d0
        Real*4,Parameter :: c155_1024  =  155.0d0 / 1024.0d0
        Real*4,Parameter :: c110_1024  =  110.0d0 / 1024.0d0
        Real*4,Parameter :: c226_1024  =  226.0d0 / 1024.0d0
        Real*4,Parameter :: c205_1024  =  205.0d0 / 1024.0d0
        Real*4,Parameter :: c251_1024  =  251.0d0 / 1024.0d0


        !--- Filter ---
        Real*4,Parameter ::   c3565_10368  =    3565.0d0 / 10368.0d0
        Real*4,Parameter ::   c3091_12960  =    3091.0d0 / 12960.0d0
        Real*4,Parameter ::   c1997_25920  =    1997.0d0 / 25920.0d0
        Real*4,Parameter ::    c149_12960  =     149.0d0 / 12960.0d0
        Real*4,Parameter ::   c107_103680  =     107.0d0 / 103680.0d0

        !--- ç∑ï™ ---
        Real*4,Parameter ::    c4_5  =  4.0d0 /   5.0d0
        Real*4,Parameter ::    c1_5  =  1.0d0 /   5.0d0
        Real*4,Parameter ::  c4_105  =  4.0d0 / 105.0d0
        Real*4,Parameter ::  c1_280  =  1.0d0 / 280.0d0

        Real*4,Parameter ::   c11_6  = 11.0d0 /   6.0d0
        Real*4,Parameter ::   c18_6  = 18.0d0 /   6.0d0
        Real*4,Parameter ::    c9_6  =  9.0d0 /   6.0d0
        Real*4,Parameter ::    c2_6  =  2.0d0 /   6.0d0

        Real*4,Parameter ::    c3_6  =  3.0d0 /   6.0d0
        Real*4,Parameter ::    c6_6  =  6.0d0 /   6.0d0
        Real*4,Parameter ::    c1_6  =  1.0d0 /   6.0d0

        Real*4,Parameter ::    c1_12 =  1.0d0 /  12.0d0
        Real*4,Parameter ::    c8_12 =  8.0d0 /  12.0d0

        Real*4,Parameter ::    c1_60 =  1.0d0 /  60.0d0
        Real*4,Parameter ::    c9_60 =  9.0d0 /  60.0d0
        Real*4,Parameter ::   c45_60 = 45.0d0 /  60.0d0



        !=== Parmeter for Communication ===
        Real*4, dimension(:,:,:), allocatable :: dumcomxy_r4,dumcomyz_r4,dumcomzx_r4
        Real*4, dimension(:,:,:), allocatable :: dumcomxy_s4,dumcomyz_s4,dumcomzx_s4
        Real*8, dimension(:,:,:), allocatable :: dumcomxy_r8,dumcomyz_r8,dumcomzx_r8
        Real*8, dimension(:,:,:), allocatable :: dumcomxy_s8,dumcomyz_s8,dumcomzx_s8
        Integer (kind=4)                      :: icom,jcom

        !--- The process (rank) number for filenames ---
        Integer fn0_rank, fn1_rank, fn2_rank, fn3_rank

        !--- Data_Save ---
        Integer nsb
        Integer is_sta, is_end
        Integer ib_sta, ib_end
        Integer data_num, rank_write



        Integer Iter_end

      End Module Definition

!****************************************************************************************************************************
!***************    IFFTópÇÃÉÇÉWÉÖÅ[Éã                ***********************************************************************
!****************************************************************************************************************************


   Module Definition_IFFT
     Implicit none

      integer :: n1,n2,n3
 
      Real*8 :: b,c,d,e,f,g,h
      Real*8 :: akre,akk,ak12,ak112,ak212,ak3k,ak12k,nk
      Real*8 :: eek1,eek
      Real*8 :: are,aim,bre,bim
      Real*8 :: rel,le,uvari,viscor,dounen,urms
!   ---Pope's spectrum parameter
      Real*8 :: fl1,fl2,fl,fet1,fet
      Real*8 :: cam,be,clf,ce,p0

!   ---IFFTópÇÃïœêîÇíËã`(Johnsen)
      Real*8 :: umt0f,relam0f,akre0f,urms0f,lam0f
      Real*8 :: ccf,tcf,viscof,rcf,pcf,ucf,vcf,wcf,tauf
      Real*8 :: dounen0f,Asp,diss0f,rel0f,el0f,eta0f
      Real*8 :: length
      Real*8, parameter :: coasp=1.0



!    ---óêêî---
      Real*8, dimension(:,:,:), allocatable :: psi1,psi2,psi3
      Real*8, dimension(:,:,:), allocatable :: theta1, theta2, theta3
      Real*8, dimension(:,:,:), allocatable :: Fr, Fi     
      Real*8, dimension(:), allocatable :: a

   End module Definition_IFFT



!==== [  End Module  ] ============================================================================================================

!==================================================================================================================================
!****** [  Program main  ] ********************************************************************************************************
!==================================================================================================================================

      Program main

        Use module_mpi
        Use Definition
        Implicit none 
        Integer i, j, k
        !--- MPI ---
        integer js,je,ks,ke,ierror,ierr,ierf
        !-----------

        !--- MPI ---
        Call mpi_init(ierror)
        Call mpi_comm_size(mpi_comm_world,nprocs,ierr)
        Call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
        Call mp_setup_division( Nx-1, Ny-1, Nz-1)
        !-----------

        Call Write_SubDomain

        If(my_rank==root) then
          Call system_clock(time_1) 
        End If
        Call mp_bcast_i( time_1 )
        time_pre = time_1
        count_round = 0


       ! Call Cpu_time( Clock_Initial ) ! CPUéûä‘ÇÃåvë™
       ! Clock_P0 = Clock_Initial

       !TEST Call System_clock( RClock_T1 )   ! äJénéûÇãLò^

        !--- Allocate Data Array ---
        !Call Sub_Allocate_MPI_Global
        !Call Sub_Allocate_MPI_Local

        If( Ityp_inG==0) Call Sub_Allocate_MPI_Local_IN
        If( Ityp_inG==1) Call Sub_Allocate_MPI_Local_IN_CG

        !If(my_rank==root) then
        !  Open(60,File = 'output/CPJET_Log.txt'        , Status = 'unknown' )
        !  Open(70,File = 'output/PJET_nt_Time.txt'    , Access = 'append' )
        !  write(70,*) 'nt Time dtime, B_Ceff'
        !End If

        If(my_rank==root) then
          Open(60,File = 'output/HIT_Log.txt'        , Status = 'unknown' )
          Open(70,File = 'output/HIT_nt_Time.txt'    , Access = 'append'  )
          Open(90,File = 'output/HIT_Log_write_time.txt'        , Status = 'unknown' )
          write(70,*) 'nt Time dtime'
          write(90,*) 'nt write_time Mt Re'

          !--- For Statistics ---
          If( Ityp==0) then 
            Open(81,File = 'output/HIT_TimeSt0_Mean.txt'       , Status = 'unknown' )
            write(81,*) 'nt Time U(3) V(4) W(5) Rh(6) T(7) P(8) divU(9) mu(10) nu(11)'

            Open(83,File = 'output/HIT_TimeSt0_var.txt'        , Status = 'unknown' )
            write(83,*) 'nt Time U2(3) V2(4) W2(5) Mat1(6) Mat2(7) Rh2(8) T2(9) P2(10) divU2(11) mu2(12) Uskew(13) Uflat(14)'

            Open(85,File = 'output/HIT_TimeSt0_Scales.txt'     , Status = 'unknown' )
            write(85,*) 'nt Time Luu(3) lx(4) ly(5) lz(6) Rel(7) eta1(8) eta2(9) dudx_skew(10) dudx_flat(11)'

            Open(87,File = 'output/HIT_TimeSt0_Diss.txt'       , Status = 'unknown' )
            write(87,*) 'nt Time Epsi_tot(3) C_Epsi(4) Epsi_s(5) Epsi_c(6) Epsi_ih(7) pdivU(8)'

            Open(89,File = 'output/HIT_TimeSt0_Energy.txt'      , Status = 'unknown' )
            write(89,*) 'nt Time RU2(3) RV2(4) RW2(5) CTKE(6)'

          Else
            Open(81,File = 'output/HIT_TimeSt0_Mean.txt'       , Access = 'append'  )
            Open(83,File = 'output/HIT_TimeSt0_var.txt'        , Access = 'append'  )
            Open(85,File = 'output/HIT_TimeSt0_Scales.txt'     , Access = 'append'  )
            Open(87,File = 'output/HIT_TimeSt0_Diss.txt'       , Access = 'append'  )
            Open(89,File = 'output/HIT_TimeSt0_Energy.txt'     , Access = 'append'  )

          End If

        End If

        !--- Set Coordinates and Parameters ---
        Call Set_Physical_Parameter
        Call Set_Parameter
        Call Set_Coordinate_DNS
        Call Set_Coordinate_DNS_Inflow

        If(my_rank==root) then
          Call Write_Parameter
          Write(60,*) ' End: Set_Coordinate_DNS'
        End If

        !=== Velocity Fluctuations in the jet nozzle ===
        !--- Calc. Velocity Fluctuations ---
        If( Ityp == 0 ) then
          If( Ityp_in .EQ. 0 ) then
            If( Ityp_inG==0) Call Inflow_Data_Calc
            If( Ityp_inG==1) Call Inflow_Data_Calc_CG
          End If

          !--- Read Velocity Fluctuations ---
          If( Ityp_in == 1 ) then
            If( Ityp_inG==0) Call Read_UVWf_in
            If( Ityp_inG==1) Call Read_UVWf_in_CG
          End If
        End If

        !=== [ Set Initial Field ] ===
        If( Ityp == 0 ) Call Inputdata_Ityp0_T      !èâä˙èåèçÏê¨
        If( Ityp == 1 ) Call Inputdata_Ityp1        !ì«Ç›çûÇÒÇ≈ë±Ç´Ç©ÇÁ
        !If( Ityp == 2 ) Call Inputdata_Ityp2        !ê¸å`ï‚äÆÇµÇƒêVÇµÇ¢DNSÇ÷
        If( Ityp == 2 ) Call Inputdata_Ityp3        !ê¸å`ï‚äÆÇµÇƒêVÇµÇ¢DNSÇ÷

        Call write_xy
        Call write_xz
        Call write_yz

        !--------------------------------------------------------
      !  Call Cpu_time( Clock_P1 ) ! CPUéûä‘ÇÃåvë™

        Call System_clock(RClock_T2, Rt_rate, Rt_max) ! é¿éûä‘ÇÃåvë™
        If ( RClock_T2 < RClock_T1 ) then
          Rt_diff = (Rt_max - RClock_T1) + RClock_T2 + 1
        else
          Rt_diff = RClock_T2 - RClock_T1
        End If

       !TEST If(my_rank==root) then
       !TEST   Write(60,*) 'Time [s] for setting initial fields CPU Real'
       !TEST   Write(60,*) Rt_diff/Dble(Rt_rate)
       !TEST   Clock_P0 = Clock_P1
       !TEST   Write(60,*) 'Time [s] for computing 1step Real'
       !TEST End If
        !--------------------------------------------------------

        !----------------------------------------------------------------------------
        !=== [ Temporal Advancement Runge-Kutta ] ===================================
        !----------------------------------------------------------------------------
        Do nt = nt0, nlast ,1

         ! !--- Computation Time ---
         ! Call System_clock( RClock_T1 ) 

          !=== dtime ===
          Call Calc_dtime

          time = time + dtime

          !If(my_rank==root) then
          !  write(70,*) nt, Time, dtime, B_Ceff,B1_Ceff,B2_Ceff
          !  write( 70,*) nt, Time, dtime, EpsisumT_xyz,  UdPdxsumT_xyz+VdPdysumT_xyz+WdPdzsumT_xyz, rkensumT_xyz
          !  write(170,*) nt, Time, dtime, rkensumT_xyz, EpsisumT_xyz ,UdPdxsumT_xyz, forcing_term,       &
          !                                              EpsisumT_xyz+UdPdxsumT_xyz+forcing_term
          !End If

          !=== éûä‘êœï™ ===
          !Call RK_5EU_1VIS
          !=== éûä‘êœï™ ===
          If(Ityp_RK == 1) Call RK_5EU_1VIS
          If(Ityp_RK == 2) Call RK_EU_VIS

          !--- Computation Time ---
         ! Call Cpu_time( Clock_P1 )
         !TEST Call System_clock(RClock_T2, Rt_rate, Rt_max) ! é¿éûä‘ÇÃåvë™
         !TEST If ( RClock_T2 < RClock_T1 ) then
         !TEST   Rt_diff = (Rt_max - RClock_T1) + RClock_T2 + 1
         !TEST else
         !TEST   Rt_diff = RClock_T2 - RClock_T1
         !TEST End If

         !TEST If(my_rank==root) then
         !TEST   Write(60,777) Dble(Rt_diff)/Dble(Rt_rate)
         !TEST End If

           If( Rh(1,j_sta,k_sta) /= Rh(1,j_sta,k_sta) ) then
             Write(*,*) 'Error', nt
             Stop
           End If

          !=== Statistics ===
          If( mod(Nt,ifl2) == 0 .or. Nt == 1) then
            !--- For output data ---
            Call Calc_Var
            Call BC_Periodic
            Call Calc_Derivative
            Call Calc_Tau
            !-----------------------

            If( mod(Nt,ifl3) == 0 .or. Nt == 1 ) then
              Call write_xy
              Call write_xz
              Call write_yz
              !Call write_xy_force

              !If(my_rank==root) then
              !  Open(75,File = 'output/HIT_VIS_Time.txt'    , Access = 'append' )
              !  Open(76,File = 'output/HIT_VIS_nt.txt'      , Access = 'append' )
              !  write(75,*) Time
              !  write(76,*) nt
              !  Close(75)
              !  Close(76)
              !End If
            End If

           ! !--- Statistics ---
           Call Sum_xyz_T
           Call Cor_xyz

           If( my_rank==root ) then
             Call Sum_Time_development_xyz

             !--- Data write ---
             If( mod(Nt,ifl4) == 0 .or. Nt == 1) Call Sumwrite_T_xyz
             If( mod(Nt,ifl4) == 0 .or. Nt == 1) Call Corwrite_T

           End If
          End If

          If( mod(Nt,ifl6) == 0 ) then
            Call Itypwrite_STI    !  3D data for shock turbulence interaction 
            Call write_xyz        !  3D data for shock turbulence interaction 
            Call Sum_xyz_T
            write(90,*) nt,                                                &
                        Time,                                              &
                        Sqrt(   ( UUsumT_xyz - UsumT_xyz*UsumT_xyz )       &  !Mt
                              + ( VVsumT_xyz - VsumT_xyz*VsumT_xyz )       &  !Mt
                              + ( WWsumT_xyz - WsumT_xyz*WsumT_xyz )  )    &  !Mt
                        /Sqrt(gh*TsumT_xyz*cRT_UU),                        &  !Mt
                          Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )            &  !Rel
                        /(Sqrt( dUdx2sumT_xyz - dUdxsumT_xyz*dUdxsumT_xyz ) ) &  !Rel
                        *Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )             &  !Rel
                        *RhsumT_xyz/musumT_xyz * Re                              !Rel
          End If
          !If( mod(Nt,ifl6) == 0 ) write(90,*) nt,Time   !  3D data for shock turbulence interaction 

          !Call Filter_F1(  Rh , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
          !Call Filter_F1(  Ur , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
          !Call Filter_F1(  Vr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
          !Call Filter_F1(  Wr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
          !Call Filter_F1(  Tr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
          !Call Filter_F1( Ysr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )

          !=== Save 3D Field ===
          !If(  time - dtime < T_inst1 .and. T_inst1 <= time ) Call write_xyz
          !If(  time - dtime < T_inst2 .and. T_inst2 <= time ) Call write_xyz
          !If(  time - dtime < T_inst3 .and. T_inst3 <= time ) Call write_xyz
          !If(  time - dtime < T_inst4 .and. T_inst4 <= time ) Call write_xyz
          !If(  time - dtime < T_inst5 .and. T_inst5 <= time ) Call write_xyz
          !If(  time - dtime < T_inst6 .and. T_inst6 <= time ) Call write_xyz

          !=== Computation RealTime ===
          If(my_rank==root) then
            Call system_clock(time_2, time_rate, time_max)

            If( time_pre > time_2 ) then 
              count_round = count_round + 1
            End If

            time_pre = time_2
            time_2 = time_max*count_round+time_2
          End If

          Call mp_bcast_i( time_2    )
          Call mp_bcast_i( time_rate )
          Call mp_bcast_i( time_max  )

          time_diff = time_2 - time_1


          !--- End Computation ---
          If( (dble(time_diff)/dble(time_rate)/3600.0d0) >= Time_limit ) then
            Call Itypwrite

            !--- For output data ---
            Call mp_send_recv_pre(  U , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            Call mp_send_recv_pre(  V , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            Call mp_send_recv_pre(  W , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            !--- x dir ---
            Call Sub_Derivative_MPI( 1  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dUdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !--- y dir ---
            Call Sub_Derivative_MPI( 2  ,   Vr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dVrdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !--- z dir ---
            Call Sub_Derivative_MPI( 3  ,   Wr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dWrdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !------------------------

            Call write_xy
            Call write_xz
            Call write_yz

            If(my_rank==root) then
              Close( 60)
              Close( 70)
              Close( 77)
              Close(170)
            End If

            Call mpi_finalize(ierror)

            Stop
          End If

          !--- End Computation ---
          If( time > T_last ) then

            !--- For output data ---
            Call mp_send_recv_pre(  U , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            Call mp_send_recv_pre(  V , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            Call mp_send_recv_pre(  W , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
            !--- x dir ---
            Call Sub_Derivative_MPI( 1  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dUdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !--- y dir ---
            Call Sub_Derivative_MPI( 2  ,   Vr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dVrdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !--- z dir ---
            Call Sub_Derivative_MPI( 3  ,   Wr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                           dWrdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            !------------------------

            Call write_xy
            Call Itypwrite

            If(my_rank==root) then
              Close(60)
              Close(70)
              Close(77)
            End If

            Call mpi_finalize(ierror)

            Stop
          End If


        End Do
        !----------------------------------------------------------------------------
        !=== [ End Temporal Advancement] ============================================
        !----------------------------------------------------------------------------


        close(71)
        close(73)

        If(my_rank==root) then
          Write(60,*) '  '
          Write(60,*) ' time : ', time
          Write(60,*) ' OK| Calculation was completed.'

          Close(60)
          Close(70)

          Close(77)
        End If


  777 Format(E13.6, 200(1x, E13.6)) 

        Call mpi_finalize(ierror)

        Stop
      End Program

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!==================================================================================================================================
!========= [ Subroutine Write_SubDomain ] =========================================================================================
!==================================================================================================================================
      Subroutine Write_SubDomain

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        !=== fname ===
        fn3_rank =   my_rank       / 1000
        fn2_rank = ( my_rank - fn3_rank * 1000     ) / 100
        fn1_rank = ( my_rank - fn3_rank * 1000 - fn2_rank * 100     ) / 10
        fn0_rank = ( my_rank - fn3_rank * 1000 - fn2_rank * 100 - fn1_rank * 10     ) / 1

        fname_write = 'output/SubDomain_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48)//'.txt'
        Open(80,File = fname_write   , Status = 'unknown' ) 
        Write(80,*) my_rank
        Write(80,*) 'my_rank='//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48)
        Write(80,*) j_sta
        Write(80,*) j_end
        Write(80,*) k_sta
        Write(80,*) k_end
        Close(80)

        Return
      End Subroutine

!==================================================================================================================================
!========= [ Subroutine Sub_Allocate_MPI_Global ] =================================================================================
!==================================================================================================================================
      Subroutine Sub_Allocate_MPI_Global

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k

        Return
      End Subroutine

!==================================================================================================================================
!========= [ Subroutine Sub_Allocate_MPI_Local ] ==================================================================================
!==================================================================================================================================
      Subroutine Sub_Allocate_MPI_Local

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        allocate(    Rh   (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(    Ur   (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(    Vr   (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(    Wr   (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(    Tr   (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(    Ysr  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )

        !allocate(    U  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        !allocate(    V  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        !allocate(    W  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Ys (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    P  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Mu (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    Ud (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Vd (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Wd (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Us (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Vs (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    Ws (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        !---î˜ï™(x-dir)---
        allocate(    dE_xdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dUrdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dRhdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dUdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dUddx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dUsdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dVdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dWdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dPdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dTdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dYsdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        !---î˜ï™(y-dir)---
        allocate(    dE_ydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dVrdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dRhdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dUdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dVdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dVddy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dVsdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dWdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dPdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dTdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dYsdy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        !---î˜ï™(z-dir)---
        allocate(    dE_zdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dWrdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dRhdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dUdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dVdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dWdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dWddz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dWsdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dPdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(      dTdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(     dYsdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(     c1_Rh (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )


        !---Flux---
        allocate(    E_x (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    E_y (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    E_z (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        !---tau_ij---
        allocate(    T_xx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_xy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_xz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    T_yx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_yy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_yz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    T_zx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_zy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    T_zz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        !===éûä‘êiçs===
        allocate(    dRh (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dTr (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dRh0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dTr0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    dUr (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dVr (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dWr (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        allocate(    dUr0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dVr0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dWr0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        !---éûä‘êiçs(îSê´çÄ)---
        allocate(    dRhv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dTrv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dUrv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dVrv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dWrv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        !---éûä‘êiçs(ÉXÉJÉâÅ[)---
        allocate(    dYsr (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dYsr0(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        !---éûä‘êiçs(ägéUçÄ)---
        allocate(    dYsrv(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )


        Rh  = 0.0d0
        Ur  = 0.0d0
        Vr  = 0.0d0
        Wr  = 0.0d0
        Tr  = 0.0d0
        Ysr = 0.0d0

        !U  = 0.0d0
        !V  = 0.0d0
        !W  = 0.0d0
        T  = 0.0d0
        Ys = 0.0d0

        P  = 0.0d0
        Mu = 0.0d0

        dE_xdx = 0.0d0
         dUrdx = 0.0d0
         dRhdx = 0.0d0
          dUdx = 0.0d0
          dVdx = 0.0d0
          dWdx = 0.0d0
          dPdx = 0.0d0
          dTdx = 0.0d0
         dYsdx = 0.0d0

        dE_ydy = 0.0d0
         dVrdy = 0.0d0
         dRhdy = 0.0d0
          dUdy = 0.0d0
          dVdy = 0.0d0
          dWdy = 0.0d0
          dPdy = 0.0d0
          dTdy = 0.0d0
         dYsdy = 0.0d0

        dE_zdz = 0.0d0
         dWrdz = 0.0d0
         dRhdz = 0.0d0
          dUdz = 0.0d0
          dVdz = 0.0d0
          dWdz = 0.0d0
          dPdz = 0.0d0
          dTdz = 0.0d0
         dYsdz = 0.0d0

         c1_Rh = 0.0d0

        E_x = 0.0d0
        E_y = 0.0d0
        E_z = 0.0d0


        T_xx = 0.0d0
        T_xy = 0.0d0
        T_xz = 0.0d0

        T_yx = 0.0d0
        T_yy = 0.0d0
        T_yz = 0.0d0

        T_zx = 0.0d0
        T_zy = 0.0d0
        T_zz = 0.0d0


        dRh  = 0.0d0
        dTr  = 0.0d0
        dRh0 = 0.0d0
        dTr0 = 0.0d0

        dUr  = 0.0d0
        dVr  = 0.0d0
        dWr  = 0.0d0

        dUr0 = 0.0d0
        dVr0 = 0.0d0
        dWr0 = 0.0d0


        dRhv = 0.0d0
        dTrv = 0.0d0
        dUrv = 0.0d0
        dVrv = 0.0d0
        dWrv = 0.0d0


        dYsr  = 0.0d0
        dYsr0 = 0.0d0

        dYsrv = 0.0d0


        Return
      End Subroutine
!==================================================================================================================================
!========= [ Subroutine Sub_Allocate_MPI_Local ] ==================================================================================
!==================================================================================================================================
      Subroutine Sub_Allocate_MPI_Local_UVW

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        allocate(    U  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    V  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    W  (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

        U  = 0.0d0
        V  = 0.0d0
        W  = 0.0d0

        Return
      End Subroutine
!==================================================================================================================================
!========= [ Subroutine Sub_Allocate_MPI_Local_IN ] ===============================================================================
!==================================================================================================================================
      Subroutine Sub_Allocate_MPI_Local_IN

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        !---èâä˙èåè---
        allocate(    Uf_in(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)     )
        allocate(    Vf_in(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)     )
        allocate(    Wf_in(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)     )

        allocate(    dUf (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dVf (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
        allocate(    dWf (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )



        Return
      End Subroutine
!==================================================================================================================================
!========= [ Subroutine Sub_Allocate_MPI_Local_IN ] ===============================================================================
!==================================================================================================================================
      Subroutine Sub_Allocate_MPI_Local_IN_CG

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        !---èâä˙èåè---
        allocate(    Uf_in(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )
        allocate(    Vf_in(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )
        allocate(    Wf_in(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )

        allocate(    dUf  (-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )
        allocate(    dVf  (-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )
        allocate(    dWf  (-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)     )



        Return
      End Subroutine
!==================================================================================================================================
!========= [ Subroutine Sub_Deallocate_MPI_Local_IN_CG ] ==========================================================================
!==================================================================================================================================
      Subroutine Sub_Deallocate_MPI_Local_IN_CG

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k


        !---èâä˙èåè---
        Deallocate(    Uf_in  )
        Deallocate(    Vf_in  )
        Deallocate(    Wf_in  )

        Deallocate(    dUf    )
        Deallocate(    dVf    )
        Deallocate(    dWf    )


        Return
      End Subroutine
!==================================================================================================================================
!****************Å@ç¿ïWåvéZÅ@******************************************************************************************************
!==================================================================================================================================
!========= [ Subroutine Set_Coordinate_DNS ] ======================================================================================
      Subroutine Set_Coordinate_DNS

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k, n ! ïœêî
        Integer Num_G, check

        !=== x coordinate ===
           dx =  Lx / (Dble(Nx-1))
        c1_dx =  1.0d0 / dx

        Do i = -4, Nx+4, 1
          x(i) = Dble(i-1) * dx
        End Do


        !=== y coordinate ===
           dy =  Ly / (Dble(Ny-1))
        c1_dy =  1.0d0 / dy

        Do j = -4, Ny+4, 1
          y(j) = Dble(j-1) * dy
        End Do


        !=== z coordinate ===
           dz =  Lz / (Dble(Nz-1))
        c1_dz =  1.0d0 / dz

        Do k = -4, Nz+4, 1
          z(k) = Dble(k-1) * dz
        End Do

        !--- local HIT spacing ---
        Do j = 1, Ny-1, 1
          dxyz = (dx*dy*dz)**c1_3
        End Do


    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!=============================================================================================================================
!========= [ Subroutine Set_Coordinate_DNS_Inflow ] ==========================================================================
      Subroutine Set_Coordinate_DNS_Inflow

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k, n ! ïœêî

        !=== x coordinate ===
           dx_in =  Lx / (Dble(Nx_in-1))
        c1_dx_in =  1.0d0 / dx_in

        Do i = 0, Nx_in, 1
          x_in(i) = Dble(i-1) * dx_in
        End Do


        !=== y coordinate ===
           dy_in =  Ly / (Dble(Ny_in-1))
        c1_dy_in =  1.0d0 / dy_in

        Do j = 0, Ny_in, 1
          y_in(j) = Dble(j-1) * dy_in
        End Do


        !=== z coordinate ===
           dz_in =  Lz / (Dble(Nz_in-1))
        c1_dz_in =  1.0d0 / dz_in

        Do k = 0, Nz_in, 1
          z_in(k) = Dble(k-1) * dz_in
        End Do


        Return
      End Subroutine
!=============================================================================================================================
!========= [ Subroutine Set_Physical_Parameter ] =============================================================================
      Subroutine Set_Physical_Parameter

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k ! ïœêî
        Real*8  ReL     ! ïœêî

        ReL = (3.0d0/20.0d0)*Re_ini*Re_ini

        !--- The speed of the towed HIT ---
        U_ini = Ma_ini*sqrt(gh*R_gas*T_inf)/sqrt(3.0d0)       ! Initial rms velocity [m/s]
        L_ini = ReL   *(  mu_T0 * ( (T_inf/273.0d0 )**(1.5d0)       )       &
                                * (273.0d0+110.4d0)/(T_inf+110.4d0) )       &
                         / ( sqrt(1.5d0)*U_ini*(P_inf/R_gas/T_inf) )                 ! Initial turbulence lenght scale [m] 

        !--- Reference value---
         U_ref   = U_ini                      ! Initial rms velocity [m/s]
         L_ref   = L_ini                      ! Initial turbulence lenght scale [m] 
        Tm_ref   = L_ref/U_ref                ! Reference time
        Rh_ref   = P_inf/R_gas/T_inf          ! Density in initial flows
         P_ref   = Rh_ref*U_ref*U_ref         ! 
         R_ref   = R_gas                      ! Gas constant
         T_ref   = P_ref/R_ref/Rh_ref         ! 
        Cp_ref   = gh*R_ref / (gh-1.0d0)      ! íËà≥î‰îM
        Cv_ref   =    R_ref / (gh-1.0d0)      ! íËêœî‰îM
        mu_ref   = mu_T0                      ! îSìx
        kT_ref   = mu_ref*Cp_ref / Pr         ! îMägéUåWêî

         cRT_UU  = 1.0d0   !R_ref*T_ref/U_ref/U_ref    ! à≥óÕïœä∑
         cUU_RT  = 1.0d0   !U_ref*U_ref/R_ref/T_ref    ! à≥óÕïœä∑

        !=== ï®óùó íËêî ===
        ! SutherlandÇÃéÆ
        T0    = 273.0d0 / T_ref             ! äÓèÄâ∑ìx
        mu0   = 0.00001724d0 / mu_ref       ! îSê´åWêî at T0 
        C_su  = 110.4d0 / T_ref             ! SutherlandÇÃéÆ

        P_bg  = P_inf /P_ref                ! é¸àÕó¨à≥óÕ
        T_bg  = T_inf /T_ref                ! é¸àÕó¨â∑ìx

        Ys_max = 1.0d0
        Ys_min = 0.0d0

        !--- Reference nondimentional parameters ---
        Re = Rh_ref*U_ref*L_ref/mu_ref  !Reynolds numbe
        !Ma = U_ref/sqrt(R_ref*T_ref)     !Mach number
        Ma = U_ref/sqrt(gh*R_ref*T_ref)     !Mach number

        !--- Mesh Reynolds number ---
        Re_ini_test   =   U_ini                                                  &
                         *L_ini                                                  &
                         *(P_inf/R_gas/T_inf)                                   &
                          /(  mu_T0 * ( (T_inf/273.0d0 )**(1.5d0)       )       &
                                 * (273.0d0+110.4d0)/(T_inf+110.4d0) )  

        !--- Mach number ---
        Ma_ini_test   = U_ini/sqrt(gh*R_gas*T_inf)

        !--- Target Kolmogorov Scale ---
        Eta_target    = Coeff_eta*dx/(Nx-1)

        c1_T0      = 1.0d0 / T0
        c1_Re      = 1.0d0 / Re
        c1_Pr      = 1.0d0 / Pr
        c1_Sc      = 1.0d0 / Sc
        c1_gh      = 1.0d0 / gh


    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!=============================================================================================================================
!========= [ Subroutine Set_Parameter ] ======================================================================================
      Subroutine Set_Parameter

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî

        !---Runge-Kutta---
        RK_A(1) = RK_A1
        RK_A(2) = RK_A2
        RK_A(3) = RK_A3
        RK_A(4) = RK_A4
        RK_A(5) = RK_A5

        RK_B(1) = RK_B1
        RK_B(2) = RK_B2
        RK_B(3) = RK_B3
        RK_B(4) = RK_B4
        RK_B(5) = RK_B5


        !---ï™êîÇÃê›íË---
        c1_Nz     = 1.0d0 / (  Dble(Nx-1)                       )
        c1_Nz     = 1.0d0 / (                        Dble(Nz-1) )
        c1_NxNz   = 1.0d0 / (  Dble(Nx-1)           *Dble(Nz-1) )
        c1_NxNyNz = 1.0d0 / (  Dble(Nx-1)*Dble(Ny-1)*Dble(Nz-1) )
        c1_NyNz   = 1.0d0 / (             Dble(Ny-1)*Dble(Nz-1) )


    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!===============================================================================================================================
!========= [ Subroutine Write_Parameter ] ======================================================================================
      Subroutine Write_Parameter

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî

        Open(80,File = 'output/HIT_Para.txt'        , Status = 'unknown' ) 

        Write(80,*) 'Nx     = ', Nx
        Write(80,*) 'Ny     = ', Ny
        Write(80,*) 'Nz     = ', Nz

        Write(80,*) 'dx     = ', dx
        Write(80,*) 'dy     = ', dy
        Write(80,*) 'dz     = ', dz

        Write(80,*) 'Re_ini       = ', Re_ini
        Write(80,*) 'Ma_ini       = ', Ma_ini
        Write(80,*) 'Pr           = ', Pr
        Write(80,*) 'Sc           = ', Sc
        Write(80,*) 'gh           = ', gh
        Write(80,*) '               '
        Write(80,*) 'Re_ini_test  = ', Re_ini_test
        Write(80,*) 'Ma_ini_test  = ', Ma_ini_test
        Write(80,*) '         '
        Write(80,*) '         '

        Write(80,*) 'Parameters of flow fields with dimensions m, s, K'
        Write(80,*) 'P_inf  = ', P_inf
        Write(80,*) 'T_inf  = ', T_inf
        Write(80,*) 'T_inf = ', T_inf

        Write(80,*) 'R_gas  = ', R_gas
        Write(80,*) 'mu_T0  = ', mu_T0

        Write(80,*) 'U_ini = ', U_ini
        Write(80,*) 'L_ini = ', L_ini


        Write(80,*) '         '
        Write(80,*) '         '

        Write(80,*) 'Reference values used for normalizing the governing equations'
        Write(80,*) 'Re     = ', Re
        Write(80,*) 'Ma     = ', Ma
        Write(80,*) '         '
        Write(80,*) 'U_ref  = ', U_ref
        Write(80,*) 'L_ref  = ', L_ref
        Write(80,*) 'Tm_ref = ', Tm_ref
        Write(80,*) 'P_ref  = ', P_ref
        Write(80,*) 'Rh_ref = ', Rh_ref
        Write(80,*) 'R_ref  = ', R_ref
        Write(80,*) 'Cp_ref = ', Cp_ref
        Write(80,*) 'Cv_ref = ', Cv_ref
        Write(80,*) 'T_ref  = ', T_ref
        Write(80,*) 'mu_ref = ', mu_ref
        Write(80,*) 'kT_ref = ', kT_ref
        Write(80,*) '         '
        Write(80,*) 'cRT_UU = ', cRT_UU
        Write(80,*) 'cUU_RT = ', cUU_RT

        Write(80,*) '         '
        Write(80,*) '         '
        Write(80,*) 'Sutherland mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)'
        Write(80,*) 'c1_T0  = ', c1_T0
        Write(80,*) 'T0     = ', T0
        Write(80,*) 'C_su   = ', C_su

        Write(80,*) '         '
        Write(80,*) '         '
        Write(80,*) 'Computational Parameter'

        Write(80,*) 'Ityp         = ', Ityp   
        Write(80,*) 'Ityp_in      = ', Ityp_in

        Write(80,*) 'Ityp_inG     = ', Ityp_inG
        Write(80,*) 'Ityp_dt      = ', Ityp_dt
        Write(80,*) 'Ityp_data    = ', Ityp_data
        Write(80,*) 'Ityp_save    = ', Ityp_save
        Write(80,*) 'n_div        = ', n_div
        Write(80,*) 'Ityp_data_in = ', Ityp_data_in

        Write(80,*) 'ifl2         = ', ifl2 
        Write(80,*) 'ifl3         = ', ifl3 
        Write(80,*) 'ifl4         = ', ifl4 
        Write(80,*) 'nlast        = ', nlast

        Write(80,*) ' '
        Write(80,*) 'Initial Condition '
        Write(80,*) 'Nx_in =  ', Nx_in
        Write(80,*) 'Ny_in =  ', Ny_in
        Write(80,*) 'Nz_in =  ', Nz_in

        Write(80,*) ' '
        Write(80,*) 'LF target '
        Write(80,*) 'E_target =  ', E_target
        Write(80,*) 'Ed_Es    =  ', Ed_Es

        Close(80)



        !--- Write cod. ---
        Open(80, File = 'output/xyz.cod' )
        Do i = 1, Nx-1, 1
          Write(80,777) x(i)
        End Do
        Do j = 1, Ny-1, 1
          Write(80,777) y(j)
        End Do
        Do k = 1, Nz-1, 1
          Write(80,777) z(k)
        End Do
        Close(80)



    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!==================================================================================================================================
!*** [ Initial Field for Temporally Evolving Jet ] ********************************************************************************
!==================================================================================================================================
!========= [ Subroutine Ityp_0 ] ==================================================================================================
      Subroutine Inputdata_Ityp0_T

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Integer n
        Real*4 sum_UV, sum_UV2
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        nt0    = 1
        nt     = 0
        time   = 0.0d0
        time_0 = 0.0d0

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- Allocate Data Array ---
        Call Sub_Allocate_MPI_Global
        Call Sub_Allocate_MPI_Local_UVW

        If( Ityp_inG==0 ) then
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1, Nx-1, 1
            U(i,j,k) = Uf_in(i,j,k)
            V(i,j,k) = Vf_in(i,j,k)
            W(i,j,k) = Wf_in(i,j,k)
          End Do; End Do; End Do
        End If

        If( Ityp_inG==1 ) then
          Call Interpolate_DNS_Inflow
          Call Sub_Deallocate_MPI_Local_IN_CG
          Call Sub_Allocate_MPI_Local
        End If

        If( Ityp_inG==1 ) then
           If( Ityp_inS==1 ) then
            Call HD_init
            Call Make_Dispersion_one
           End If
        End If

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          P (i,j,k) = P_bg
          Rh(i,j,k) = P_bg / T_bg
          Tr(i,j,k) = T_bg * Rh(i,j,k)

          Ur(i,j,k) = (  U(i,j,k)*rms_in ) * Rh(i,j,k)
          Vr(i,j,k) = (  V(i,j,k)*rms_in ) * Rh(i,j,k)
          Wr(i,j,k) = (  W(i,j,k)*rms_in ) * Rh(i,j,k)

         Ysr(i,j,k) = 0.0d0* Rh(i,j,k)
        End Do; End Do; End Do

        If( Ityp_Press==1 ) then
          Call Calc_initial_pressure_fluctuations
        End If

        Call Calc_Var
        Call BC_Periodic
        Call Calc_Derivative
        Call Calc_Tau

        ! !--- Statistics ---
        Call Sum_xyz_T
        Call Cor_xyz

        If( my_rank==root ) then
          Call Sum_Time_development_xyz
        End If

        Return
      End Subroutine
!==============================================================================================================================
!========= [ Subroutine Ityp_1 ] ==============================================================================================
      Subroutine Inputdata_Ityp1

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        Integer ns, is, ie

        !--- MPI ---
        Integer js,je,ks,ke
        Integer j_sta_dammy,j_end_dammy,k_sta_dammy,k_end_dammy
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call Sub_Allocate_MPI_Global
        Call Sub_Allocate_MPI_Local_UVW
        Call Sub_Deallocate_MPI_Local_IN_CG
        Call Sub_Allocate_MPI_Local



        !==========================================================================================
        If(Ityp_data .eq. 1) then 
          Open(40, File = 'temp/Inst_Data_Nt', Form = 'unFormatted', Status = 'unknown')
          Read(40) Nt, time
          Close(40)

          !---------------------
          Open(40, File = 'temp/Inst_Data_Rh', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Rh(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

          !---------------------
          Open(40, File = 'temp/Inst_Data_Ur', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Ur(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

          !---------------------
          Open(40, File = 'temp/Inst_Data_Vr', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Vr(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

          !---------------------
          Open(40, File = 'temp/Inst_Data_Wr', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Wr(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

          !---------------------
          Open(40, File = 'temp/Inst_Data_Tr', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Tr(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

          !---------------------
          Open(40, File = 'temp/Inst_Data_Ysr', Form = 'unFormatted', Status = 'unknown')

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do i = -4, Nx+4, 1
            Read(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
          End Do
          Close(40)

          Do k = ks-5,ke+5,1
          Do j = js-5,je+5,1
          Do i = -4,Nx+4,1
            Ysr(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate (Data_r4 )

        End If

        !==========================================================================================
        If(Ityp_data .eq. 2) then 
          !Open(40, File = 'temp/Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48),&
          !     Form = 'unFormatted', Status = 'unknown')           !ñºëÂÉXÉpÉRÉì
          Open(40, File = 'Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48),&
               Form = 'unFormatted', Status = 'unknown')            !   ínãÖÉVÉ~ÉÖÉåÅ[É^

write(*,*) 'Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48)



          Read(40) Nt, time
          time_0 = 0.0d0
          Read(40) j_sta_dammy,j_end_dammy,k_sta_dammy,k_end_dammy
          Read(40) Rh
          Read(40) Ur
          Read(40) Vr
          Read(40) Wr
          Read(40) Tr
          Read(40) Ysr
          Close(40)
        End If



        !--- BC ---
        !Call Calc_Var
        !Call BC_Periodic
        Call Calc_Var
        Call BC_Periodic
        Call Calc_Derivative

        nt0 = 1 + nt

  999   Continue

        Return
      End Subroutine
!==============================================================================================================================
!========= [ Subroutine Ityp_1 ] ==============================================================================================
      Subroutine Inputdata_Ityp2

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî 
        Real*4, dimension(:,:,:),allocatable :: Data_r4

        Integer ns, is, ie
        Integer xi_old, yj_old, zk_old
        Integer js_old,je_old,ks_old,ke_old

        Real*8  dx_old, dy_old, dz_old
        Real*8  ax, by, cz    !ï‚ä‘éÆåWêî
        Real*8  x_old (-4:Nx_old+4)
        Real*8  y_old (-4:Ny_old+4)
        Real*8  z_old (-4:Nz_old+4)

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call Sub_Allocate_MPI_Global
        Call Sub_Allocate_MPI_Local_UVW
        Call Sub_Deallocate_MPI_Local_IN_CG
        Call Sub_Allocate_MPI_Local

        !=== x coordinate ===
        dx_old =  Lx / (Dble(Nx_old-1))
        Do i = -4, Nx_old+4, 1
          x_old(i) = Dble(i-1) * dx_old
        End Do

        !=== y coordinate ===
        dy_old =  Ly / (Dble(Ny_old-1))
        Do j = -4, Ny_old+4, 1
          y_old(j) = Dble(j-1) * dy_old
        End Do

        !=== z coordinate ===
        dz_old =  Lz / (Dble(Nz_old-1))
        Do k = -4, Nz_old+4, 1
          z_old(k) = Dble(k-1) * dz_old
        End Do

        !==========================================================================================
          Open(40, File = 'temp/Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48),&
               Form = 'unFormatted', Status = 'unknown')
          Read(40) Nt, time
          time_0 = time
          Read(40) js_old,je_old,ks_old,ke_old
          Nt   = 0
          time = 0.0d0

        !=====  Rh interpolation  ============================================================================
          Allocate (Data_r4(-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4


          !!====== ã´äEèåèÇÃê›íË ======
          !! --- x direction ---
          !Do k = ks_old, ke_old, 1
          !Do j = js_old, je_old, 1
          !  Data_r4(Nx_old+0,j,k) =  Data_r4(       1+0,j,k)
          !  Data_r4(Nx_old+1,j,k) =  Data_r4(       1+1,j,k)
          !  Data_r4(Nx_old+2,j,k) =  Data_r4(       1+2,j,k)
          !  Data_r4(Nx_old+3,j,k) =  Data_r4(       1+3,j,k)
          !  Data_r4(Nx_old+4,j,k) =  Data_r4(       1+4,j,k)
          !  Data_r4(     0-0,j,k) =  Data_r4(Nx_old-1-0,j,k)
          !  Data_r4(     0-1,j,k) =  Data_r4(Nx_old-1-1,j,k)
          !  Data_r4(     0-2,j,k) =  Data_r4(Nx_old-1-2,j,k)
          !  Data_r4(     0-3,j,k) =  Data_r4(Nx_old-1-3,j,k)
          !  Data_r4(     0-4,j,k) =  Data_r4(Nx_old-1-4,j,k)
          !End Do; End Do
          !! --- y,z direction ---
          !Call BC_Periodic_y_dir_MPI_R4_init( -4,Nx_old+4,js_old-5,je_old+5,ks_old-5,ke_old+5, Data_r4 )
          !Call BC_Periodic_z_dir_MPI_R4_init( -4,Nx_old+4,js_old-5,je_old+5,ks_old-5,ke_old+5, Data_r4 )


          !============================

          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Rh(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

!If(xi_old == Nx_old .or. yj_old == Ny_old .or. zk_old == Nz_old ) then
!If(Data_r4(xi_old, yj_old  ,zk_old) == 0.0d0) then
!write(*,*) i,j,k,xi_old, yj_old  ,zk_old,Data_r4(xi_old, yj_old  ,zk_old)
!End If
!End If
          End Do; End Do; End Do
          Deallocate(  Data_r4  )

          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1, Nx-1
          !  x_int = int( x(i) / dx_in + 1.0d0 ) + 1
          !  y_int = int( y(j) / dy_in + 1.0d0 ) + 1
          !  z_int = int( z(k) / dz_in + 1.0d0 ) + 1
          !  ax = ( x_in(x_int)  - (x(i)) ) / dx_in
          !  by = ( y_in(y_int)  - (y(j)) ) / dy_in
          !  cz = ( z_in(z_int)  - (z(k)) ) / dz_in
          !  U(i,j,k) =             cz * (      by  * (            ax  * Uf_in(x_int-1, y_int-1,z_int-1)       &
          !                                             + (1.0d0 - ax) * Uf_in(x_int  , y_int-1,z_int-1)  )    &
          !                            + (1.0d0-by) * (            ax  * Uf_in(x_int-1, y_int  ,z_int-1)       &
          !                                             + (1.0d0 - ax) * Uf_in(x_int  , y_int  ,z_int-1)  )  ) &
          !           +    (1.0d0 - cz)* (      by  * (            ax  * Uf_in(x_int-1, y_int-1,z_int  )       &
          !                                             + (1.0d0 - ax) * Uf_in(x_int  , y_int-1,z_int  )  )    &
          !                            + (1.0d0-by) * (            ax  * Uf_in(x_int-1, y_int  ,z_int  )       &
          !                                             + (1.0d0 - ax) * Uf_in(x_int  , y_int  ,z_int  )  )  )

          !End Do; End Do; End Do


            !xi_old = int( x(i) / dx_old)! + 1
            !yj_old = int( y(j) / dy_old)! + 1
            !zk_old = int( z(k) / dz_old)! + 1

            !ax = ( x(i) -  x_old(xi_old)) / dx_old
            !by = ( y(j) -  y_old(yj_old)) / dy_old
            !cz = ( z(k) -  z_old(zk_old)) / dz_old

            !Rh(i,j,k) = + (1.0d0 - cz) * (   (1.0d0-by) * (   (1.0d0 - ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
            !                                                + (        ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
            !                               + (      by) * (   (1.0d0 - ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
            !                                                + (        ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
            !            + (        cz) * (   (1.0d0-by) * (   (1.0d0 - ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
            !                                                + (        ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
            !                               + (      by) * (   (1.0d0 - ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
            !                                                + (        ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  








        !=====  Ur interpolation  ============================================================================
          Allocate (Data_r4 (-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4


!aaa
        If(my_rank==root) then
          Open (50, File = 'test_Nxoldm2.txt', Status = 'unknown' ) 
          i = Nx_old-2
          Do j = js_old-5,je_old+5, 1
          Do k = ks_old-5,ke_old+5, 1
            Write(50,*) y_old(j), z_old(k), Data_r4(i,j,k)
          End Do; End Do
          Close(50) 

          Open (50, File = 'test_Nxoldm1.txt', Status = 'unknown' ) 
          i = Nx_old-1
          Do j = js_old-5,je_old+5, 1
          Do k = ks_old-5,ke_old+5, 1
            Write(50,*) y_old(j), z_old(k), Data_r4(i,j,k)
          End Do; End Do
          Close(50) 

          Open (50, File = 'test_Nxold.txt', Status = 'unknown' ) 
          i = Nx_old
          Do j = js_old-5,je_old+5, 1
          Do k = ks_old-5,ke_old+5, 1
            Write(50,*) y_old(j), z_old(k), Data_r4(i,j,k)
          End Do; End Do
          Close(50) 

          Open (50, File = 'test_Nxoldp1.txt', Status = 'unknown' ) 
          i = Nx_old+1
          Do j = js_old-5,je_old+5, 1
          Do k = ks_old-5,ke_old+5, 1
            Write(50,*) y_old(j), z_old(k), Data_r4(i,j,k)
          End Do; End Do
          Close(50) 

          Open (50, File = 'test_Nxoldp2.txt', Status = 'unknown' ) 
          i = Nx_old+2
          Do j = js_old-5,je_old+5, 1
          Do k = ks_old-5,ke_old+5, 1
            Write(50,*) y_old(j), z_old(k), Data_r4(i,j,k)
          End Do; End Do
          Close(50) 
        End If
!aaa
stop



          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1

            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old

            Ur(i,j,k) =   (        cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

          End Do; End Do; End Do
          Deallocate(  Data_r4  )

        !=====  Vr interpolation  ============================================================================
          Allocate (Data_r4 (-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1

            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old

            Vr(i,j,k) =   (        cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

          End Do; End Do; End Do
          Deallocate(  Data_r4  )

        !=====  Wr interpolation  ============================================================================
          Allocate (Data_r4 (-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1

            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old

            Wr(i,j,k) =   (        cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

          End Do; End Do; End Do
          Deallocate(  Data_r4  )

        !=====  Tr interpolation  ============================================================================
          Allocate (Data_r4 (-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4

          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1

            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old

            Tr(i,j,k) =   (        cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        !=====  Ysr interpolation  ============================================================================
          Allocate (Data_r4 (-4:Nx_old+4,js_old-5:je_old+5,ks_old-5:ke_old+5))
          Data_r4 = 0.0d0
          Read(40) Data_r4

          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1

            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old

            Ysr(i,j,k) =  (        cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old-1)  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old-1)       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz)* (   (      by) * (   (        ax) * Data_r4(xi_old-1, yj_old-1,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old-1,zk_old  )  )    &
                                          + (1.0d0-by) * (   (        ax) * Data_r4(xi_old-1, yj_old  ,zk_old  )       &
                                                           + (1.0d0 - ax) * Data_r4(xi_old  , yj_old  ,zk_old  )  )  )  

          End Do; End Do; End Do
          Deallocate(  Data_r4  )

          Close(40)

        ! --- BC x direction --- 
        Do k = ks, ke, 1
        Do j = js, je, 1

          !--- Rh ---
          Rh(Nx+0,j,k) =  Rh(0+1,j,k)
          Rh(Nx+1,j,k) =  Rh(1+1,j,k)
          Rh(Nx+2,j,k) =  Rh(2+1,j,k)
          Rh(Nx+3,j,k) =  Rh(3+1,j,k)
          Rh(Nx+4,j,k) =  Rh(4+1,j,k)

          Rh( 0,j,k)   =  Rh(Nx-1-0,j,k)
          Rh(-1,j,k)   =  Rh(Nx-1-1,j,k)
          Rh(-2,j,k)   =  Rh(Nx-1-2,j,k)
          Rh(-3,j,k)   =  Rh(Nx-1-3,j,k)
          Rh(-4,j,k)   =  Rh(Nx-1-4,j,k)
          !--- Ur ---
          Ur(Nx+0,j,k) =  Ur(0+1,j,k)
          Ur(Nx+1,j,k) =  Ur(1+1,j,k)
          Ur(Nx+2,j,k) =  Ur(2+1,j,k)
          Ur(Nx+3,j,k) =  Ur(3+1,j,k)
          Ur(Nx+4,j,k) =  Ur(4+1,j,k)
                          
          Ur( 0,j,k)   =  Ur(Nx-1-0,j,k)
          Ur(-1,j,k)   =  Ur(Nx-1-1,j,k)
          Ur(-2,j,k)   =  Ur(Nx-1-2,j,k)
          Ur(-3,j,k)   =  Ur(Nx-1-3,j,k)
          Ur(-4,j,k)   =  Ur(Nx-1-4,j,k)
          !--- Vr ---
          Vr(Nx+0,j,k) =  Vr(0+1,j,k)
          Vr(Nx+1,j,k) =  Vr(1+1,j,k)
          Vr(Nx+2,j,k) =  Vr(2+1,j,k)
          Vr(Nx+3,j,k) =  Vr(3+1,j,k)
          Vr(Nx+4,j,k) =  Vr(4+1,j,k)
                          
          Vr( 0,j,k)   =  Vr(Nx-1-0,j,k)
          Vr(-1,j,k)   =  Vr(Nx-1-1,j,k)
          Vr(-2,j,k)   =  Vr(Nx-1-2,j,k)
          Vr(-3,j,k)   =  Vr(Nx-1-3,j,k)
          Vr(-4,j,k)   =  Vr(Nx-1-4,j,k)
          !--- Wr ---
          Wr(Nx+0,j,k) =  Wr(0+1,j,k)
          Wr(Nx+1,j,k) =  Wr(1+1,j,k)
          Wr(Nx+2,j,k) =  Wr(2+1,j,k)
          Wr(Nx+3,j,k) =  Wr(3+1,j,k)
          Wr(Nx+4,j,k) =  Wr(4+1,j,k)
                          
          Wr( 0,j,k)   =  Wr(Nx-1-0,j,k)
          Wr(-1,j,k)   =  Wr(Nx-1-1,j,k)
          Wr(-2,j,k)   =  Wr(Nx-1-2,j,k)
          Wr(-3,j,k)   =  Wr(Nx-1-3,j,k)
          Wr(-4,j,k)   =  Wr(Nx-1-4,j,k)
          !--- Tr ---
          Tr(Nx+0,j,k) =  Tr(0+1,j,k)
          Tr(Nx+1,j,k) =  Tr(1+1,j,k)
          Tr(Nx+2,j,k) =  Tr(2+1,j,k)
          Tr(Nx+3,j,k) =  Tr(3+1,j,k)
          Tr(Nx+4,j,k) =  Tr(4+1,j,k)
                          
          Tr( 0,j,k)   =  Tr(Nx-1-0,j,k)
          Tr(-1,j,k)   =  Tr(Nx-1-1,j,k)
          Tr(-2,j,k)   =  Tr(Nx-1-2,j,k)
          Tr(-3,j,k)   =  Tr(Nx-1-3,j,k)
          Tr(-4,j,k)   =  Tr(Nx-1-4,j,k)
          !--- Ysr ---
          Ysr(Nx+0,j,k) =  Ysr(0+1,j,k)
          Ysr(Nx+1,j,k) =  Ysr(1+1,j,k)
          Ysr(Nx+2,j,k) =  Ysr(2+1,j,k)
          Ysr(Nx+3,j,k) =  Ysr(3+1,j,k)
          Ysr(Nx+4,j,k) =  Ysr(4+1,j,k)
                           
          Ysr( 0,j,k)   =  Ysr(Nx-1-0,j,k)
          Ysr(-1,j,k)   =  Ysr(Nx-1-1,j,k)
          Ysr(-2,j,k)   =  Ysr(Nx-1-2,j,k)
          Ysr(-3,j,k)   =  Ysr(Nx-1-3,j,k)
          Ysr(-4,j,k)   =  Ysr(Nx-1-4,j,k)
        End Do ; End Do

        ! --- BC y direction --- 
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )

        ! --- BC z direction --- 
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )

        ! --- ë≥óÃàÊ --- 
        Call mp_send_recv_pre( Rh , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Ur , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Vr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Wr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Tr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre(Ysr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

        !--- BC ---
        !Call Calc_Var
        !Call BC_Periodic
        Call Calc_Var
        Call BC_Periodic
        Call Calc_Derivative

        nt0 = 1 + nt

  999   Continue

        Return
      End Subroutine
!==============================================================================================================================
!========= [ Subroutine Ityp_1 ] ==============================================================================================
      Subroutine Inputdata_Ityp3

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî 
        Real*4, dimension(:,:,:),allocatable :: Data_r4
        Real*4, dimension(:,:,:),allocatable :: Data_r_all

        Integer n_rank

        Integer fn3_rank_r
        Integer fn2_rank_r
        Integer fn1_rank_r
        Integer fn0_rank_r

        Integer ns, is, ie
        Integer xi_old, yj_old, zk_old
        Integer js_old(0:Num_rank_old-1)
        Integer je_old(0:Num_rank_old-1)
        Integer ks_old(0:Num_rank_old-1)
        Integer ke_old(0:Num_rank_old-1)

        Real*8  dx_old, dy_old, dz_old
        Real*8  ax, by, cz    !ï‚ä‘éÆåWêî
        Real*8  x_old (-4:Nx_old+4)
        Real*8  y_old (-4:Ny_old+4)
        Real*8  z_old (-4:Nz_old+4)

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call Sub_Allocate_MPI_Global
        Call Sub_Allocate_MPI_Local_UVW
        Call Sub_Deallocate_MPI_Local_IN_CG
        Call Sub_Allocate_MPI_Local

        !=== x coordinate ===
        dx_old =  Lx / (Dble(Nx_old-1))
        Do i = -4, Nx_old+4, 1
          x_old(i) = Dble(i-1) * dx_old
        End Do

        !=== y coordinate ===
        dy_old =  Ly / (Dble(Ny_old-1))
        Do j = -4, Ny_old+4, 1
          y_old(j) = Dble(j-1) * dy_old
        End Do

        !=== z coordinate ===
        dz_old =  Lz / (Dble(Nz_old-1))
        Do k = -4, Nz_old+4, 1
          z_old(k) = Dble(k-1) * dz_old
        End Do

        Allocate (Data_r_all(-4:Nx_old+4,-4:Ny_old+4,-4:Nz_old+4))
        Data_r_all = 0.0d0

        Do n_rank = 0, Num_rank_old-1, 1
          !=== fname ===
          fn3_rank_r =   n_rank              / 1000
          fn2_rank_r = ( n_rank - fn3_rank_r * 1000            ) / 100
          fn1_rank_r = ( n_rank - fn3_rank_r * 1000 - fn2_rank_r * 100            ) / 10
          fn0_rank_r = ( n_rank - fn3_rank_r * 1000 - fn2_rank_r * 100 - fn1_rank_r * 10     ) / 1

          Open(40+n_rank, File = 'temp/Inst_Data_myrank'                         &
                                 //char(fn3_rank_r+48)//char(fn2_rank_r+48)      &
                                 //char(fn1_rank_r+48)//char(fn0_rank_r+48),     &
               Form = 'unFormatted', Status = 'unknown')
          Read(40+n_rank) Nt, time
          time_0 = time
          Nt     = 0
          time   = 0.0d0
          Read(40+n_rank) js_old(n_rank),  &
                          je_old(n_rank),  &
                          ks_old(n_rank),  &
                          ke_old(n_rank)
        End Do

        !===  Rh ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Rh(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do





        !===  Ur  ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Ur(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do








        !===  Vr  ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Vr(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do









        !===  Wr  ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Wr(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do








        !===  Tr  ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Tr(i,j,k) =   (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do

        !===  Ysr  ===
        Do n_rank = 0, Num_rank_old-1, 1
          Allocate (Data_r4(-4:Nx_old+4,js_old(n_rank)-5:je_old(n_rank)+5,ks_old(n_rank)-5:ke_old(n_rank)+5))
          Data_r4 = 0.0d0
          Read(40+n_rank) Data_r4

          Do k = ks_old(n_rank),   ke_old(n_rank), 1
          Do j = js_old(n_rank),   je_old(n_rank), 1
          Do i =  1, Nx_old-1, 1
            Data_r_all(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do
          Deallocate(  Data_r4  )
        End Do

        !====== ã´äEèåèÇÃê›íË ======
        Do k = 1, Nz_old-1, 1
        Do j = 1, Ny_old-1, 1
          Data_r_all(Nx_old+0,j,k) =  Data_r_all(       1+0,j,k)
          Data_r_all(Nx_old+1,j,k) =  Data_r_all(       1+1,j,k)
          Data_r_all(Nx_old+2,j,k) =  Data_r_all(       1+2,j,k)
          Data_r_all(Nx_old+3,j,k) =  Data_r_all(       1+3,j,k)
          Data_r_all(Nx_old+4,j,k) =  Data_r_all(       1+4,j,k)
          Data_r_all(     0-0,j,k) =  Data_r_all(Nx_old-1-0,j,k)
          Data_r_all(     0-1,j,k) =  Data_r_all(Nx_old-1-1,j,k)
          Data_r_all(     0-2,j,k) =  Data_r_all(Nx_old-1-2,j,k)
          Data_r_all(     0-3,j,k) =  Data_r_all(Nx_old-1-3,j,k)
          Data_r_all(     0-4,j,k) =  Data_r_all(Nx_old-1-4,j,k)
        End Do; End Do

        Do k =  1, Nz_old-1, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,Ny_old+0,k) =  Data_r_all(i,       1+0,k)
          Data_r_all(i,Ny_old+1,k) =  Data_r_all(i,       1+1,k)
          Data_r_all(i,Ny_old+2,k) =  Data_r_all(i,       1+2,k)
          Data_r_all(i,Ny_old+3,k) =  Data_r_all(i,       1+3,k)
          Data_r_all(i,Ny_old+4,k) =  Data_r_all(i,       1+4,k)
          Data_r_all(i,     0-0,k) =  Data_r_all(i,Ny_old-1-0,k)
          Data_r_all(i,     0-1,k) =  Data_r_all(i,Ny_old-1-1,k)
          Data_r_all(i,     0-2,k) =  Data_r_all(i,Ny_old-1-2,k)
          Data_r_all(i,     0-3,k) =  Data_r_all(i,Ny_old-1-3,k)
          Data_r_all(i,     0-4,k) =  Data_r_all(i,Ny_old-1-4,k)
        End Do; End Do

        Do j = -4, Ny_old+4, 1
        Do i = -4, Nx_old+4, 1
          Data_r_all(i,j,Nz_old+0) =  Data_r_all(i,j,       1+0)
          Data_r_all(i,j,Nz_old+1) =  Data_r_all(i,j,       1+1)
          Data_r_all(i,j,Nz_old+2) =  Data_r_all(i,j,       1+2)
          Data_r_all(i,j,Nz_old+3) =  Data_r_all(i,j,       1+3)
          Data_r_all(i,j,Nz_old+4) =  Data_r_all(i,j,       1+4)
          Data_r_all(i,j,     0-0) =  Data_r_all(i,j,Nz_old-1-0)
          Data_r_all(i,j,     0-1) =  Data_r_all(i,j,Nz_old-1-1)
          Data_r_all(i,j,     0-2) =  Data_r_all(i,j,Nz_old-1-2)
          Data_r_all(i,j,     0-3) =  Data_r_all(i,j,Nz_old-1-3)
          Data_r_all(i,j,     0-4) =  Data_r_all(i,j,Nz_old-1-4)
        End Do; End Do


          !============================
          Do k = ks,   ke, 1
          Do j = js,   je, 1
          Do i =  1, Nx-1, 1
            xi_old = int( x(i) / dx_old + 1.0d0 ) + 1
            yj_old = int( y(j) / dy_old + 1.0d0 ) + 1
            zk_old = int( z(k) / dz_old + 1.0d0 ) + 1
            ax = ( x_old(xi_old)  - (x(i)) ) / dx_old
            by = ( y_old(yj_old)  - (y(j)) ) / dy_old
            cz = ( z_old(zk_old)  - (z(k)) ) / dz_old
            Ysr(i,j,k) =  (        cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old-1)  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old-1)       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old-1)  )  ) &
                        + (1.0d0 - cz) * (   (      by) * (   (        ax) * Data_r_all(xi_old-1, yj_old-1,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old-1,zk_old  )  )    &
                                           + (1.0d0-by) * (   (        ax) * Data_r_all(xi_old-1, yj_old  ,zk_old  )       &
                                                            + (1.0d0 - ax) * Data_r_all(xi_old  , yj_old  ,zk_old  )  )  )  
          End Do; End Do; End Do




        Do n_rank = 0, Num_rank_old-1, 1
          Close(40+n_rank)
        End Do

        Deallocate(  Data_r_all  )




        ! --- BC x direction --- 
        Do k = ks, ke, 1
        Do j = js, je, 1

          !--- Rh ---
          Rh(Nx+0,j,k) =  Rh(0+1,j,k)
          Rh(Nx+1,j,k) =  Rh(1+1,j,k)
          Rh(Nx+2,j,k) =  Rh(2+1,j,k)
          Rh(Nx+3,j,k) =  Rh(3+1,j,k)
          Rh(Nx+4,j,k) =  Rh(4+1,j,k)

          Rh( 0,j,k)   =  Rh(Nx-1-0,j,k)
          Rh(-1,j,k)   =  Rh(Nx-1-1,j,k)
          Rh(-2,j,k)   =  Rh(Nx-1-2,j,k)
          Rh(-3,j,k)   =  Rh(Nx-1-3,j,k)
          Rh(-4,j,k)   =  Rh(Nx-1-4,j,k)
          !--- Ur ---
          Ur(Nx+0,j,k) =  Ur(0+1,j,k)
          Ur(Nx+1,j,k) =  Ur(1+1,j,k)
          Ur(Nx+2,j,k) =  Ur(2+1,j,k)
          Ur(Nx+3,j,k) =  Ur(3+1,j,k)
          Ur(Nx+4,j,k) =  Ur(4+1,j,k)
                          
          Ur( 0,j,k)   =  Ur(Nx-1-0,j,k)
          Ur(-1,j,k)   =  Ur(Nx-1-1,j,k)
          Ur(-2,j,k)   =  Ur(Nx-1-2,j,k)
          Ur(-3,j,k)   =  Ur(Nx-1-3,j,k)
          Ur(-4,j,k)   =  Ur(Nx-1-4,j,k)
          !--- Vr ---
          Vr(Nx+0,j,k) =  Vr(0+1,j,k)
          Vr(Nx+1,j,k) =  Vr(1+1,j,k)
          Vr(Nx+2,j,k) =  Vr(2+1,j,k)
          Vr(Nx+3,j,k) =  Vr(3+1,j,k)
          Vr(Nx+4,j,k) =  Vr(4+1,j,k)
                          
          Vr( 0,j,k)   =  Vr(Nx-1-0,j,k)
          Vr(-1,j,k)   =  Vr(Nx-1-1,j,k)
          Vr(-2,j,k)   =  Vr(Nx-1-2,j,k)
          Vr(-3,j,k)   =  Vr(Nx-1-3,j,k)
          Vr(-4,j,k)   =  Vr(Nx-1-4,j,k)
          !--- Wr ---
          Wr(Nx+0,j,k) =  Wr(0+1,j,k)
          Wr(Nx+1,j,k) =  Wr(1+1,j,k)
          Wr(Nx+2,j,k) =  Wr(2+1,j,k)
          Wr(Nx+3,j,k) =  Wr(3+1,j,k)
          Wr(Nx+4,j,k) =  Wr(4+1,j,k)
                          
          Wr( 0,j,k)   =  Wr(Nx-1-0,j,k)
          Wr(-1,j,k)   =  Wr(Nx-1-1,j,k)
          Wr(-2,j,k)   =  Wr(Nx-1-2,j,k)
          Wr(-3,j,k)   =  Wr(Nx-1-3,j,k)
          Wr(-4,j,k)   =  Wr(Nx-1-4,j,k)
          !--- Tr ---
          Tr(Nx+0,j,k) =  Tr(0+1,j,k)
          Tr(Nx+1,j,k) =  Tr(1+1,j,k)
          Tr(Nx+2,j,k) =  Tr(2+1,j,k)
          Tr(Nx+3,j,k) =  Tr(3+1,j,k)
          Tr(Nx+4,j,k) =  Tr(4+1,j,k)
                          
          Tr( 0,j,k)   =  Tr(Nx-1-0,j,k)
          Tr(-1,j,k)   =  Tr(Nx-1-1,j,k)
          Tr(-2,j,k)   =  Tr(Nx-1-2,j,k)
          Tr(-3,j,k)   =  Tr(Nx-1-3,j,k)
          Tr(-4,j,k)   =  Tr(Nx-1-4,j,k)
          !--- Ysr ---
          Ysr(Nx+0,j,k) =  Ysr(0+1,j,k)
          Ysr(Nx+1,j,k) =  Ysr(1+1,j,k)
          Ysr(Nx+2,j,k) =  Ysr(2+1,j,k)
          Ysr(Nx+3,j,k) =  Ysr(3+1,j,k)
          Ysr(Nx+4,j,k) =  Ysr(4+1,j,k)
                           
          Ysr( 0,j,k)   =  Ysr(Nx-1-0,j,k)
          Ysr(-1,j,k)   =  Ysr(Nx-1-1,j,k)
          Ysr(-2,j,k)   =  Ysr(Nx-1-2,j,k)
          Ysr(-3,j,k)   =  Ysr(Nx-1-3,j,k)
          Ysr(-4,j,k)   =  Ysr(Nx-1-4,j,k)
        End Do ; End Do

        ! --- BC y direction --- 
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr  )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )

        ! --- BC z direction --- 
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr  )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )

        ! --- ë≥óÃàÊ --- 
        Call mp_send_recv_pre( Rh , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Ur , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Vr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Wr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Tr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre(Ysr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

          Nt   = 0
          time = 0.0d0

        !--- BC ---
        !Call Calc_Var
        !Call BC_Periodic
        Call Calc_Var
        Call BC_Periodic
        Call Calc_Derivative

        Do j = js, je
          Do k = ks, ke
          Do i = 1, Nx-1

             UsumT_xyz =    UsumT_xyz + (    U(i,j,k)              ) * c1_NxNyNz
             VsumT_xyz =    VsumT_xyz + (    V(i,j,k)              ) * c1_NxNyNz
             WsumT_xyz =    WsumT_xyz + (    W(i,j,k)              ) * c1_NxNyNz
            UUsumT_xyz =   UUsumT_xyz + (    U(i,j,k) *  U(i,j,k)  ) * c1_NxNyNz
            VVsumT_xyz =   VVsumT_xyz + (    V(i,j,k) *  V(i,j,k)  ) * c1_NxNyNz
            WWsumT_xyz =   WWsumT_xyz + (    W(i,j,k) *  W(i,j,k)  ) * c1_NxNyNz

          End Do; End Do
        End Do

        Call mp_allsumr(    UsumT_xyz  )
        Call mp_allsumr(    VsumT_xyz  )
        Call mp_allsumr(    WsumT_xyz  )
        Call mp_allsumr(    UUsumT_xyz )
        Call mp_allsumr(    VVsumT_xyz )
        Call mp_allsumr(    WWsumT_xyz )

        Do j = js, je
          Do k = ks, ke
          Do i = 1, Nx-1

             U(i,j,k) = U(i,j,k)/(sqrt(UUsumT_xyz - UsumT_xyz*UsumT_xyz))
             V(i,j,k) = V(i,j,k)/(sqrt(VVsumT_xyz - VsumT_xyz*VsumT_xyz))
             W(i,j,k) = W(i,j,k)/(sqrt(WWsumT_xyz - WsumT_xyz*WsumT_xyz))

          End Do; End Do
        End Do

        nt0 = 1 + nt

  999   Continue

        Return
      End Subroutine
!===============================================================================================================================
!****************Å@éûä‘êœï™Å@***************************************************************************************************
!===============================================================================================================================
!========= [ Subroutine RK_5EU_1VIS ] ==========================================================================================
      Subroutine RK_5EU_1VIS

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        !---îÒîSê´çÄRunge Kutta 5-stage 4th-order---
        Do Num_RK = 1, 5, 1
          Call Calc_Var
          Call BC_Periodic
          Call Calc_Derivative

          Call NS3D_SkewSym      !===îÒîSê´çÄ ( ë¨ìx )===
          Call ScalarEq_SkewSym  !===îÒîSê´çÄ (Scalar)===

          !---ñßìx, ë¨ìx, â∑ìx îÒîSê´çÄÇÃéûä‘êiçs---
          !---Scalar îÒîSê´çÄÇÃéûä‘êiçs---
          Call RK5_All

        End Do

        !---îSê´çÄ ägéUçÄ 1st-order---
        Call Calc_Var
        Call BC_Periodic
        Call Calc_Derivative

        Call Calc_Tau
        Call Linear_Forcing

        Call NS3D_RHO_T_Vis  !===îSê´çÄ===
        Call ScalarEq_Vis    !===ägéUçÄ===

        Call Eu1_All

        Call BC_Periodic

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
!========= [ Subroutine RK_EU_VIS ] ==========================================================================================
      Subroutine RK_EU_VIS

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        !---îÒîSê´çÄRunge Kutta 5-stage 4th-order---
        Do Num_RK = 1, 5, 1
          Call Calc_Var
          Call BC_Periodic
          Call Calc_Derivative
          Call Calc_Tau

          If (Num_RK == 1) Call Linear_Forcing

          Call NS3D_SkewSym      !===îÒîSê´çÄ ( ë¨ìx )===
          Call ScalarEq_SkewSym  !===îÒîSê´çÄ (Scalar)===
          Call NS3D_RHO_T_Vis    !===îSê´çÄ===
          Call ScalarEq_Vis      !===ägéUçÄ===

          !---ñßìx, ë¨ìx, â∑ìxÇÃéûä‘êiçs---
          Call RK5_All_NS3D

        End Do

        !=== Sponge Zone ===
        Call Filter_F1(  Rh , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
        Call Filter_F1(  Ur , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
        Call Filter_F1(  Vr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
        Call Filter_F1(  Wr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
        Call Filter_F1(  Tr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )
        Call Filter_F1( Ysr , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5 )

        Call BC_Periodic


        Return
      End Subroutine

!=============================================================================================================================
!=====[ Runge Kutta 5 ]=======================================================================================================
!=============================================================================================================================
      Subroutine RK5_All_NS3D

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 :: f1,f2,f3
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- Viscous terms ---
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1

           If( Ityp_Force == 0 ) then
             f1 = 0.0d0
             f2 = 0.0d0
             f3 = 0.0d0
           End If

           If( Ityp_Force == 1 ) then
             f1 = B_Ceff*Rh(i,j,k)*U(i,j,k)
             f2 = B_Ceff*Rh(i,j,k)*V(i,j,k)
             f3 = B_Ceff*Rh(i,j,k)*W(i,j,k)
           End If

           If( Ityp_Force == 2 ) then
             f1 = + c_d*sqrt( Rh(i,j,k) )*Ud(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Us(i,j,k)   
             f2 = + c_d*sqrt( Rh(i,j,k) )*Vd(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Vs(i,j,k)   
             f3 = + c_d*sqrt( Rh(i,j,k) )*Wd(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Ws(i,j,k)   
           End If

           dUr(i,j,k) =  dUr(i,j,k) +  dUrv(i,j,k) + f1
           dVr(i,j,k) =  dVr(i,j,k) +  dVrv(i,j,k) + f2
           dWr(i,j,k) =  dWr(i,j,k) +  dWrv(i,j,k) + f3
           dTr(i,j,k) =  dTr(i,j,k) +  dTrv(i,j,k) - (gh-1.0d0)*(f1*U(i,j,k)+f2*V(i,j,k)+f3*W(i,j,k))

          dYsr(i,j,k) = dYsr(i,j,k) + dYsrv(i,j,k)
        End Do; End Do; End Do


        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          Rh(i,j,k) =  Rh(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dRh0(i,j,k) + dtime * dRh(i,j,k) )
          Ur(i,j,k) =  Ur(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dUr0(i,j,k) + dtime * dUr(i,j,k) )
          Vr(i,j,k) =  Vr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dVr0(i,j,k) + dtime * dVr(i,j,k) )
          Wr(i,j,k) =  Wr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dWr0(i,j,k) + dtime * dWr(i,j,k) )
          Tr(i,j,k) =  Tr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dTr0(i,j,k) + dtime * dTr(i,j,k) )

          !--- Scalar ---
          Ysr(i,j,k) = Ysr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dYsr0(i,j,k) + dtime * dYsr(i,j,k) )
        End Do; End Do; End Do



        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          dRh0(i,j,k) =  ( RK_A(Num_RK) * dRh0(i,j,k) + dtime * dRh(i,j,k) )
          dUr0(i,j,k) =  ( RK_A(Num_RK) * dUr0(i,j,k) + dtime * dUr(i,j,k) )
          dVr0(i,j,k) =  ( RK_A(Num_RK) * dVr0(i,j,k) + dtime * dVr(i,j,k) )
          dWr0(i,j,k) =  ( RK_A(Num_RK) * dWr0(i,j,k) + dtime * dWr(i,j,k) )
          dTr0(i,j,k) =  ( RK_A(Num_RK) * dTr0(i,j,k) + dtime * dTr(i,j,k) )

          !--- Scalar ---
          dYsr0(i,j,k) =  ( RK_A(Num_RK) * dYsr0(i,j,k) + dtime * dYsr(i,j,k) )

        End Do; End Do; End Do




        Return
      End Subroutine

!=============================================================================================================================
!=====[ Runge Kutta 5 ]=======================================================================================================
!=============================================================================================================================
      Subroutine RK5_All

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          Rh(i,j,k) =  Rh(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dRh0(i,j,k) + dtime * dRh(i,j,k) )
          Ur(i,j,k) =  Ur(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dUr0(i,j,k) + dtime * dUr(i,j,k) )
          Vr(i,j,k) =  Vr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dVr0(i,j,k) + dtime * dVr(i,j,k) )
          Wr(i,j,k) =  Wr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dWr0(i,j,k) + dtime * dWr(i,j,k) )
          Tr(i,j,k) =  Tr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dTr0(i,j,k) + dtime * dTr(i,j,k) )

          !--- Scalar ---
          Ysr(i,j,k) = Ysr(i,j,k) + RK_B(Num_RK) * ( RK_A(Num_RK) * dYsr0(i,j,k) + dtime * dYsr(i,j,k) )

        End Do; End Do; End Do

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          dRh0(i,j,k) =  ( RK_A(Num_RK) * dRh0(i,j,k) + dtime * dRh(i,j,k) )
          dUr0(i,j,k) =  ( RK_A(Num_RK) * dUr0(i,j,k) + dtime * dUr(i,j,k) )
          dVr0(i,j,k) =  ( RK_A(Num_RK) * dVr0(i,j,k) + dtime * dVr(i,j,k) )
          dWr0(i,j,k) =  ( RK_A(Num_RK) * dWr0(i,j,k) + dtime * dWr(i,j,k) )
          dTr0(i,j,k) =  ( RK_A(Num_RK) * dTr0(i,j,k) + dtime * dTr(i,j,k) )

          !--- Scalar ---
          dYsr0(i,j,k) =  ( RK_A(Num_RK) * dYsr0(i,j,k) + dtime * dYsr(i,j,k) )

        End Do; End Do; End Do



        Return
      End Subroutine
!=============================================================================================================================
!=====[ Eu1_All ]=============================================================================================================
!=============================================================================================================================
      Subroutine Eu1_All

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Real*4 :: f1,f2,f3

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1

           If( Ityp_Force == 0 ) then
             f1 = 0.0d0
             f2 = 0.0d0
             f3 = 0.0d0
           End If

           If( Ityp_Force == 1 ) then
             f1 = B_Ceff*Rh(i,j,k)*U(i,j,k)
             f2 = B_Ceff*Rh(i,j,k)*V(i,j,k)
             f3 = B_Ceff*Rh(i,j,k)*W(i,j,k)
           End If

           If( Ityp_Force == 2 ) then
             f1 = + c_d*sqrt( Rh(i,j,k) )*Ud(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Us(i,j,k)   
             f2 = + c_d*sqrt( Rh(i,j,k) )*Vd(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Vs(i,j,k)   
             f3 = + c_d*sqrt( Rh(i,j,k) )*Wd(i,j,k)  &
                  + c_s*sqrt( Rh(i,j,k) )*Ws(i,j,k)   
           End If

           dUr(i,j,k) =  dUr(i,j,k) +  dUrv(i,j,k) + f1
           dVr(i,j,k) =  dVr(i,j,k) +  dVrv(i,j,k) + f2
           dWr(i,j,k) =  dWr(i,j,k) +  dWrv(i,j,k) + f3
           dTr(i,j,k) =  dTr(i,j,k) +  dTrv(i,j,k) - (gh-1.0d0)*(f1*U(i,j,k)+f2*V(i,j,k)+f3*W(i,j,k))

          dYsr(i,j,k) = dYsr(i,j,k) + dYsrv(i,j,k)
        End Do; End Do; End Do

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          Ur(i,j,k) =  Ur(i,j,k) + dtime * dUrv(i,j,k)
          Vr(i,j,k) =  Vr(i,j,k) + dtime * dVrv(i,j,k)
          Wr(i,j,k) =  Wr(i,j,k) + dtime * dWrv(i,j,k)
          Tr(i,j,k) =  Tr(i,j,k) + dtime * dTrv(i,j,k)

         !--- Scalar ---
         Ysr(i,j,k) = Ysr(i,j,k) + dtime * dYsrv(i,j,k)

        End Do; End Do; End Do


        Return
      End Subroutine
!=============================================================================================================================
!=====[ dRh, dU, dT Calc. ]===================================================================================================
!=============================================================================================================================
      Subroutine NS3D_SkewSym

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 dF1dx, dF1dy, dcF1_F2dx, dcF1_F2dy
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

       !---dRh---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dRh(i,j,k) = - c1_2 * ( dUrdx(i,j,k) + dVrdy(i,j,k) + dWrdz(i,j,k) )    &
         !              - c1_2 * (  (  U (i,j,k) *dRhdx(i,j,k)                     &
         !                           + V (i,j,k) *dRhdy(i,j,k)                     &
         !                           + W (i,j,k) *dRhdz(i,j,k)   )                 &
         !                        + (  Rh(i,j,k) *dUdx (i,j,k)                     &
         !                           + Rh(i,j,k) *dVdy (i,j,k)                     &
         !                           + Rh(i,j,k) *dWdz (i,j,k)   )   )

          dRh(i,j,k) = - c1_2 * ( ( dUrdx(i,j,k) + dVrdy(i,j,k) + dWrdz(i,j,k) )    &
                                 +(  (  U (i,j,k) *dRhdx(i,j,k)                     &
                                      + V (i,j,k) *dRhdy(i,j,k)                     &
                                      + W (i,j,k) *dRhdz(i,j,k)   )                 &
                                   + Rh(i,j,k) *(  dUdx (i,j,k)                     &
                                                 + dVdy (i,j,k)                     &
                                                 + dWdz (i,j,k)   )   )  )
        End Do; End Do; End Do


       !---dUr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = -4+ks, ke+4, 1
        Do j = -4+js, je+4, 1
        Do i = -4+1 , Nx-1+4
          E_x(i,j,k) =  Ur(i,j,k)*U (i,j,k)
          E_y(i,j,k) =  Ur(i,j,k)*V (i,j,k)
          E_z(i,j,k) =  Ur(i,j,k)*W (i,j,k)
        End Do; End Do; End Do

       ! Call BC_Periodic_E_xyz

       ! Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )


       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dUr(i,j,k) = - cRT_UU* dPdx(i,j,k)                                       &
         !              - c1_2 * ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
         !              - c1_2 * (  (  U (i,j,k) *U (i,j,k) *dRhdx(i,j,k)           &
         !                           + U (i,j,k) *V (i,j,k) *dRhdy(i,j,k)           &
         !                           + U (i,j,k) *W (i,j,k) *dRhdz(i,j,k) )         &

         !                        + (  Rh(i,j,k) *U (i,j,k) *dUdx (i,j,k)           &
         !                           + Rh(i,j,k) *V (i,j,k) *dUdy (i,j,k)           &
         !                           + Rh(i,j,k) *W (i,j,k) *dUdz (i,j,k) )         &

         !                        + (  Rh(i,j,k) *U (i,j,k) *dUdx (i,j,k)           &
         !                           + Rh(i,j,k) *U (i,j,k) *dVdy (i,j,k)           &
         !                           + Rh(i,j,k) *U (i,j,k) *dWdz (i,j,k) )   )

          dUr(i,j,k) = - cRT_UU* dPdx(i,j,k)                                       &
                       - c1_2*( ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
                               +(  U (i,j,k) *(  U (i,j,k) *dRhdx(i,j,k)           &
                                                +V (i,j,k) *dRhdy(i,j,k)           &
                                                +W (i,j,k) *dRhdz(i,j,k) )         &

                                 + Rh(i,j,k) *(  U (i,j,k) *( dUdx (i,j,k)         &
                                                             +dUdx (i,j,k)         &
                                                             +dVdy (i,j,k)         &
                                                             +dWdz (i,j,k) )       &
                                                +V (i,j,k) *dUdy (i,j,k)           &
                                                +W (i,j,k) *dUdz (i,j,k)    )  ) ) 
                       !!--- Linear Forcing ---
                       !+ c_d*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*          Ud(i,j,k) )       &
                       !+ c_s*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*(U(i,j,k)-Ud(i,j,k)))        
                       !+ B_Ceff *Ur(i,j,k)
                       !+ B1_Ceff *Ur(i,j,k)

        End Do; End Do; End Do


       !---dVr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = -4+ks, ke+4, 1
        Do j = -4+js, je+4, 1
        Do i = -4+1 , Nx-1+4
          E_x(i,j,k) = Vr(i,j,k)*U (i,j,k)
          E_y(i,j,k) = Vr(i,j,k)*V (i,j,k)
          E_z(i,j,k) = Vr(i,j,k)*W (i,j,k)
        End Do; End Do; End Do

       ! Call BC_Periodic_E_xyz

       ! Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )


       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dVr(i,j,k) = - cRT_UU* dPdy(i,j,k)                                       &
         !              - c1_2 * ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
         !              - c1_2 * (  (  V (i,j,k) *U (i,j,k) *dRhdx(i,j,k)           &
         !                           + V (i,j,k) *V (i,j,k) *dRhdy(i,j,k)           &
         !                           + V (i,j,k) *W (i,j,k) *dRhdz(i,j,k) )         &

         !                        + (  Rh(i,j,k) *U (i,j,k) *dVdx (i,j,k)           &
         !                           + Rh(i,j,k) *V (i,j,k) *dVdy (i,j,k)           &
         !                           + Rh(i,j,k) *W (i,j,k) *dVdz (i,j,k) )         &

         !                        + (  Rh(i,j,k) *V (i,j,k) *dUdx (i,j,k)           &
         !                           + Rh(i,j,k) *V (i,j,k) *dVdy (i,j,k)           &
         !                           + Rh(i,j,k) *V (i,j,k) *dWdz (i,j,k) )   )

          dVr(i,j,k) = - cRT_UU* dPdy(i,j,k)                                       &
                       - c1_2*( ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
                               +(  V (i,j,k) *(  U (i,j,k) *dRhdx(i,j,k)           &
                                                +V (i,j,k) *dRhdy(i,j,k)           &
                                                +W (i,j,k) *dRhdz(i,j,k) )         &

                                 + Rh(i,j,k) *(  U (i,j,k) *dVdx (i,j,k)           &
                                                +V (i,j,k) *( dVdy (i,j,k)         &
                                                             +dUdx (i,j,k)         &
                                                             +dVdy (i,j,k)         &
                                                             +dWdz (i,j,k) )       &
                                                +W (i,j,k) *dVdz (i,j,k)    )  ) ) 
                       !!--- Linear Forcing ---
                       !+ c_d*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*          Vd(i,j,k) )       &
                       !+ c_s*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*(V(i,j,k)-Vd(i,j,k)))        

                       !+ B_Ceff *Vr(i,j,k)
                       !+ B1_Ceff *Vr(i,j,k)

        End Do; End Do; End Do


       !---dWr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = -4+ks, ke+4, 1
        Do j = -4+js, je+4, 1
        Do i = -4+1 , Nx-1+4
          E_x(i,j,k) = Wr(i,j,k)*U (i,j,k)
          E_y(i,j,k) = Wr(i,j,k)*V (i,j,k)
          E_z(i,j,k) = Wr(i,j,k)*W (i,j,k)
        End Do; End Do; End Do

       ! Call BC_Periodic_E_xyz

       ! Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dWr(i,j,k) = - cRT_UU* dPdz(i,j,k)                                       &
         !              - c1_2 * ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
         !              - c1_2 * (  (  W (i,j,k) *U (i,j,k) *dRhdx(i,j,k)           &
         !                           + W (i,j,k) *V (i,j,k) *dRhdy(i,j,k)           &
         !                           + W (i,j,k) *W (i,j,k) *dRhdz(i,j,k) )         &

         !                        + (  Rh(i,j,k) *U (i,j,k) *dWdx (i,j,k)           &
         !                           + Rh(i,j,k) *V (i,j,k) *dWdy (i,j,k)           &
         !                           + Rh(i,j,k) *W (i,j,k) *dWdz (i,j,k) )         &

         !                        + (  Rh(i,j,k) *W (i,j,k) *dUdx (i,j,k)           &
         !                           + Rh(i,j,k) *W (i,j,k) *dVdy (i,j,k)           &
         !                           + Rh(i,j,k) *W (i,j,k) *dWdz (i,j,k) )   )

          dWr(i,j,k) = - cRT_UU* dPdz(i,j,k)                                       &
                       - c1_2*( ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )  &
                               +(  W (i,j,k) *(  U (i,j,k) *dRhdx(i,j,k)           &
                                                +V (i,j,k) *dRhdy(i,j,k)           &
                                                +W (i,j,k) *dRhdz(i,j,k) )         &

                                 + Rh(i,j,k) *(  U (i,j,k) *dWdx (i,j,k)           &
                                                +V (i,j,k) *dWdy (i,j,k)           &
                                                +W (i,j,k) *( dWdz (i,j,k)         &
                                                             +dUdx (i,j,k)         &
                                                             +dVdy (i,j,k)         &
                                                             +dWdz (i,j,k) ))  ) ) 
                       !!--- Linear Forcing ---
                       !+ c_d*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*          Wd(i,j,k) )       &
                       !+ c_s*sqrt( Rh(i,j,k) )*(sqrt( Rh(i,j,k) )*(W(i,j,k)-Wd(i,j,k)))        

                       !+ B_Ceff *Wr(i,j,k)
                       !+ B1_Ceff *Wr(i,j,k)

        End Do; End Do; End Do


       !---dTr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = -4+ks, ke+4, 1
        Do j = -4+js, je+4, 1
        Do i = -4+1 , Nx-1+4
          E_x(i,j,k) = Tr(i,j,k)*U (i,j,k)
          E_y(i,j,k) = Tr(i,j,k)*V (i,j,k)
          E_z(i,j,k) = Tr(i,j,k)*W (i,j,k)
        End Do; End Do; End Do

       ! Call BC_Periodic_E_xyz

       ! Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dTr(i,j,k) = -  (gh-1.0d0)*P(i,j,k)*( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )   &
         !              - c1_2 * ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )           &
         !              - c1_2 * (  (  T (i,j,k) *U (i,j,k) *dRhdx(i,j,k)                    &
         !                           + T (i,j,k) *V (i,j,k) *dRhdy(i,j,k)                    &
         !                           + T (i,j,k) *W (i,j,k) *dRhdz(i,j,k)   )                &

         !                        + (  Rh(i,j,k) *U (i,j,k) *dTdx (i,j,k)                    &
         !                           + Rh(i,j,k) *V (i,j,k) *dTdy (i,j,k)                    &
         !                           + Rh(i,j,k) *W (i,j,k) *dTdz (i,j,k)   )                &

         !                        + (  Rh(i,j,k) *T (i,j,k) *dUdx (i,j,k)                    &
         !                           + Rh(i,j,k) *T (i,j,k) *dVdy (i,j,k)                    &
         !                           + Rh(i,j,k) *T (i,j,k) *dWdz (i,j,k)   )   )

          dTr(i,j,k) = -  (gh-1.0d0)*P(i,j,k)*( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )   &
                       - c1_2*( ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )           &
                               +(  T (i,j,k) *(  U (i,j,k) *dRhdx(i,j,k)                    &
                                                +V (i,j,k) *dRhdy(i,j,k)                    &
                                                +W (i,j,k) *dRhdz(i,j,k)   )                &

                                 + Rh(i,j,k) *(  U (i,j,k) *dTdx (i,j,k)                    &
                                                +V (i,j,k) *dTdy (i,j,k)                    &
                                                +W (i,j,k) *dTdz (i,j,k)                    &

                                                +T (i,j,k)*(dUdx (i,j,k)                    &
                                                           +dVdy (i,j,k)                    &
                                                           +dWdz (i,j,k)   ))    )  )
        End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!=====[ dYsr Calc. ]==========================================================================================================
!=============================================================================================================================
      Subroutine ScalarEq_SkewSym

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 dF1dx, dF1dy, dcF1_F2dx, dcF1_F2dy
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

       !---dYsr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = -4+ks, ke+4, 1
        Do j = -4+js, je+4, 1
        Do i = -4+1 , Nx-1+4
          E_x(i,j,k) =  Ysr(i,j,k)*U (i,j,k)
          E_y(i,j,k) =  Ysr(i,j,k)*V (i,j,k)
          E_z(i,j,k) =  Ysr(i,j,k)*W (i,j,k)
        End Do; End Do; End Do

       ! Call BC_Periodic_E_xyz

       ! Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
       ! Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dYsr(i,j,k) = - c1_2 * ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )      &
         !               - c1_2 * (  (  Ys (i,j,k) *U  (i,j,k) *dRhdx(i,j,k)             &
         !                            + Ys (i,j,k) *V  (i,j,k) *dRhdy(i,j,k)             &
         !                            + Ys (i,j,k) *W  (i,j,k) *dRhdz(i,j,k)    )        &
         !                         + (  Rh (i,j,k) *U  (i,j,k) *dYsdx(i,j,k)             &
         !                            + Rh (i,j,k) *V  (i,j,k) *dYsdy(i,j,k)             &
         !                            + Rh (i,j,k) *W  (i,j,k) *dYsdz(i,j,k)    )        &
         !                         + (  Rh (i,j,k) *Ys (i,j,k) *dUdx (i,j,k)             &
         !                            + Rh (i,j,k) *Ys (i,j,k) *dVdy (i,j,k)             &
         !                            + Rh (i,j,k) *Ys (i,j,k) *dWdz (i,j,k)    )   )

          dYsr(i,j,k) = - c1_2*( ( dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) )      &
                                +(  Ys (i,j,k) *(  U  (i,j,k) *dRhdx(i,j,k)             &
                                                  +V  (i,j,k) *dRhdy(i,j,k)             &
                                                  +W  (i,j,k) *dRhdz(i,j,k)     )       &
                                  + Rh (i,j,k) *(  U  (i,j,k) *dYsdx(i,j,k)             &
                                                  +V  (i,j,k) *dYsdy(i,j,k)             &
                                                  +W  (i,j,k) *dYsdz(i,j,k)             &
                                                  +Ys (i,j,k) *( dUdx (i,j,k)           &
                                                                +dVdy (i,j,k)           &
                                                                +dWdz (i,j,k) ) )  )   )

        End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!=====[ dRhv, dUv, dTv Calc. ]================================================================================================
!=============================================================================================================================
      Subroutine NS3D_RHO_T_Vis

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 dF1dx, dF1dy, dcF1_F2dx, dcF1_F2dy
        Real*4 f1,f2,f3
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


         !---dRh---
        !  Do k = 1, Nz-1
        !  Do j = 1, Ny-1
        !  Do i = 1, Nx-1
        !    dRhv(i,j,k) = 0.0d0
        !  End Do; End Do; End Do

       !---dUr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          E_x(i,j,k) =  T_xx(i,j,k)*c1_Re
          E_y(i,j,k) =  T_xy(i,j,k)*c1_Re
          E_z(i,j,k) =  T_xz(i,j,k)*c1_Re
        End Do; End Do; End Do

        Call BC_Periodic_E_xyz

        Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          !f1 = + c_d*sqrt( Rh(i,j,k) )*Ud(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Us(i,j,k)        
          !f1 = B_Ceff*Ur(i,j,k)

          dUrv(i,j,k) =   dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k)
         ! dUrv(i,j,k) =   dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) + f1

                        !+ B_Ceff *Ur(i,j,k)      !--- Linear Forcing ---

        End Do; End Do; End Do

       !---dVr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          E_x(i,j,k) =  T_yx(i,j,k)*c1_Re
          E_y(i,j,k) =  T_yy(i,j,k)*c1_Re
          E_z(i,j,k) =  T_yz(i,j,k)*c1_Re
        End Do; End Do; End Do

        Call BC_Periodic_E_xyz

        Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          !f2 = + c_d*sqrt( Rh(i,j,k) )*Vd(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Vs(i,j,k)        
          !f2 = B_Ceff*Vr(i,j,k)

          dVrv(i,j,k) =   dE_xdx(i,j,k)+ dE_ydy(i,j,k) + dE_zdz(i,j,k)
          !dVrv(i,j,k) =   dE_xdx(i,j,k)+ dE_ydy(i,j,k) + dE_zdz(i,j,k) + f2

                        !+ B_Ceff *Vr(i,j,k)      !--- Linear Forcing ---
        End Do; End Do; End Do

       !---dWr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          E_x(i,j,k) =  T_zx(i,j,k)*c1_Re
          E_y(i,j,k) =  T_zy(i,j,k)*c1_Re
          E_z(i,j,k) =  T_zz(i,j,k)*c1_Re
        End Do; End Do; End Do

        Call BC_Periodic_E_xyz

        Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          !f3 = + c_d*sqrt( Rh(i,j,k) )*Wd(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Ws(i,j,k)        
          !f3 = B_Ceff*Wr(i,j,k)

          !dWrv(i,j,k) =   dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k) + f3
          dWrv(i,j,k) =   dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k)

                        !+ B_Ceff *Wr(i,j,k)      !--- Linear Forcing ---
        End Do; End Do; End Do

        !---dTr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          E_x(i,j,k) =  Mu(i,j,k)*dTdx(i,j,k)
          E_y(i,j,k) =  Mu(i,j,k)*dTdy(i,j,k)
          E_z(i,j,k) =  Mu(i,j,k)*dTdz(i,j,k)

         ! E_x(i,j,k) =  kT(i,j,k)*dTdx(i,j,k)
         ! E_y(i,j,k) =  kT(i,j,k)*dTdy(i,j,k)
         ! E_z(i,j,k) =  kT(i,j,k)*dTdz(i,j,k)
        End Do; End Do; End Do

        Call BC_Periodic_E_xyz

        Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
         ! dTrv(i,j,k) =  + cUU_RT*(gh-1.0d0)*c1_Re*( + T_xx(i,j,k)*dUdx(i,j,k)         &
         !                                            + T_yx(i,j,k)*dVdx(i,j,k)         &
         !                                            + T_zx(i,j,k)*dWdx(i,j,k)   )     &

         !                + cUU_RT*(gh-1.0d0)*c1_Re*( + T_xy(i,j,k)*dUdy(i,j,k)         &
         !                                            + T_yy(i,j,k)*dVdy(i,j,k)         &
         !                                            + T_zy(i,j,k)*dWdy(i,j,k)   )     &

         !                + cUU_RT*(gh-1.0d0)*c1_Re*( + T_xz(i,j,k)*dUdz(i,j,k)         &
         !                                            + T_yz(i,j,k)*dVdz(i,j,k)         &
         !                                            + T_zz(i,j,k)*dWdz(i,j,k)   )     &

         !                           +  gh*c1_Re*c1_Pr*( dE_xdx(i,j,k) )                &
         !                           +  gh*c1_Re*c1_Pr*( dE_zdz(i,j,k) )

          !f1 = + c_d*sqrt( Rh(i,j,k) )*Ud(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Us(i,j,k)        
          !f2 = + c_d*sqrt( Rh(i,j,k) )*Vd(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Vs(i,j,k)        
          !f3 = + c_d*sqrt( Rh(i,j,k) )*Wd(i,j,k)       &
          !     + c_s*sqrt( Rh(i,j,k) )*Ws(i,j,k)        

          !f1 = B_ceff*Ur(i,j,k)
          !f2 = B_ceff*Vr(i,j,k)
          !f3 = B_ceff*Wr(i,j,k)

          dTrv(i,j,k) =  + cUU_RT*(gh-1.0d0)*c1_Re*(  ( + T_xx(i,j,k)*dUdx(i,j,k)         &
                                                        + T_yx(i,j,k)*dVdx(i,j,k)         &
                                                        + T_zx(i,j,k)*dWdx(i,j,k)   )     &

                                                    + ( + T_xy(i,j,k)*dUdy(i,j,k)         &
                                                        + T_yy(i,j,k)*dVdy(i,j,k)         &
                                                        + T_zy(i,j,k)*dWdy(i,j,k)   )     &

                                                    + ( + T_xz(i,j,k)*dUdz(i,j,k)         &
                                                        + T_yz(i,j,k)*dVdz(i,j,k)         &
                                                        + T_zz(i,j,k)*dWdz(i,j,k)   ) )   &

                                    +  gh*c1_Re*c1_Pr*(  dE_xdx(i,j,k)                    &
                                                        +dE_ydy(i,j,k)                    &
                                                        +dE_zdz(i,j,k) )                   
                         !- (gh-1.0d0)*( f1*U(i,j,k) + f2*V(i,j,k) + f3*W(i,j,k) )


        End Do; End Do; End Do



        Return
      End Subroutine
!=============================================================================================================================
!=====[ dYsrv Calc. ]=========================================================================================================
!=============================================================================================================================
      Subroutine ScalarEq_Vis

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 dF1dx, dF1dy, dcF1_F2dx, dcF1_F2dy
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        !---dYsr---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          E_x(i,j,k) =  Mu(i,j,k)*c1_Sc*c1_Re*dYsdx(i,j,k)
          E_y(i,j,k) =  Mu(i,j,k)*c1_Sc*c1_Re*dYsdy(i,j,k)
          E_z(i,j,k) =  Mu(i,j,k)*c1_Sc*c1_Re*dYsdz(i,j,k)

         ! E_x(i,j,k) =  Dys(i,j,k)*c1_Sc*c1_Re*dYsdx(i,j,k)
         ! E_y(i,j,k) =  Dys(i,j,k)*c1_Sc*c1_Re*dYsdy(i,j,k)
         ! E_z(i,j,k) =  Dys(i,j,k)*c1_Sc*c1_Re*dYsdz(i,j,k)
        End Do; End Do; End Do

        Call BC_Periodic_E_xyz

        Call mp_send_recv_pre(E_x,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_y,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre(E_z,4,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

       ! Call Sub_dF1dx_Local_MPI ( E_x, dE_xdx )
       ! Call Sub_dF1dy_Local_MPI ( E_y, dE_ydy )
       ! Call Sub_dF1dz_Local_MPI ( E_z, dE_zdz )

        Call Sub_Derivative_MPI( 1  ,   E_x   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_xdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 2  ,   E_y   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_ydy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
        Call Sub_Derivative_MPI( 3  ,   E_z   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                       dE_zdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          dYsrv(i,j,k) =  dE_xdx(i,j,k) + dE_ydy(i,j,k) + dE_zdz(i,j,k)
        End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!*** [Periodic BC] ***********************************************************************************************************
!=============================================================================================================================
      Subroutine BC_Periodic_Tau

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call BC_Periodic_Tau_x_para
        Call BC_Periodic_Tau_y_para
        Call BC_Periodic_Tau_z_para

        Call mp_send_recv_pre( T_xx , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_xy , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_xz , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

        Call mp_send_recv_pre( T_yx , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_yy , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_yz , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

        Call mp_send_recv_pre( T_zx , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_zy , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( T_zz , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

        Return
      End Subroutine
!=============================================================================================================================
!*** [Periodic BC] ***********************************************************************************************************
!=============================================================================================================================
      Subroutine BC_Periodic

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call BC_Periodic_x_para
        Call BC_Periodic_y_para
        Call BC_Periodic_z_para

        Call mp_send_recv_pre( Rh , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Ur , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Vr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Wr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre( Tr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )
        Call mp_send_recv_pre(Ysr , 4, -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5 )

        !--- Calculate U, V, W, P, T, Ys in the boundary region between mpi processes ---
        Call Calc_Var_BC_MPI_r4


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluid variables x direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_x_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        ! Rh 
        ! Ur 
        ! Vr 
        ! Wr 
        ! Tr 

        ! U  
        ! V  
        ! W  
        ! T  

        ! P  
        ! Mu 
        ! kT 

        ! Ysr
        ! Ys 


        !====== ã´äEèåèÇÃê›íË ======
        Do k = ks, ke, 1
        Do j = js, je, 1
          !Do i = 0, nbc, 1
          !  F1(Nx+i,j,k) =  F1(i+1,j,k)
          !End Do
          !Do i = -nbc, 0, 1
          !  F1(i,j,k) =  F1(Nx-1+i,j,k)
          !End Do

          !--- Rh ---
          Rh(Nx+0,j,k) =  Rh(0+1,j,k)
          Rh(Nx+1,j,k) =  Rh(1+1,j,k)
          Rh(Nx+2,j,k) =  Rh(2+1,j,k)
          Rh(Nx+3,j,k) =  Rh(3+1,j,k)
          Rh(Nx+4,j,k) =  Rh(4+1,j,k)

          Rh( 0,j,k)   =  Rh(Nx-1-0,j,k)
          Rh(-1,j,k)   =  Rh(Nx-1-1,j,k)
          Rh(-2,j,k)   =  Rh(Nx-1-2,j,k)
          Rh(-3,j,k)   =  Rh(Nx-1-3,j,k)
          Rh(-4,j,k)   =  Rh(Nx-1-4,j,k)

          !--- Ur ---
          Ur(Nx+0,j,k) =  Ur(0+1,j,k)
          Ur(Nx+1,j,k) =  Ur(1+1,j,k)
          Ur(Nx+2,j,k) =  Ur(2+1,j,k)
          Ur(Nx+3,j,k) =  Ur(3+1,j,k)
          Ur(Nx+4,j,k) =  Ur(4+1,j,k)
                          
          Ur( 0,j,k)   =  Ur(Nx-1-0,j,k)
          Ur(-1,j,k)   =  Ur(Nx-1-1,j,k)
          Ur(-2,j,k)   =  Ur(Nx-1-2,j,k)
          Ur(-3,j,k)   =  Ur(Nx-1-3,j,k)
          Ur(-4,j,k)   =  Ur(Nx-1-4,j,k)

          !--- Vr ---
          Vr(Nx+0,j,k) =  Vr(0+1,j,k)
          Vr(Nx+1,j,k) =  Vr(1+1,j,k)
          Vr(Nx+2,j,k) =  Vr(2+1,j,k)
          Vr(Nx+3,j,k) =  Vr(3+1,j,k)
          Vr(Nx+4,j,k) =  Vr(4+1,j,k)

          Vr( 0,j,k)   =  Vr(Nx-1-0,j,k)
          Vr(-1,j,k)   =  Vr(Nx-1-1,j,k)
          Vr(-2,j,k)   =  Vr(Nx-1-2,j,k)
          Vr(-3,j,k)   =  Vr(Nx-1-3,j,k)
          Vr(-4,j,k)   =  Vr(Nx-1-4,j,k)

          !--- Wr ---
          Wr(Nx+0,j,k) =  Wr(0+1,j,k)
          Wr(Nx+1,j,k) =  Wr(1+1,j,k)
          Wr(Nx+2,j,k) =  Wr(2+1,j,k)
          Wr(Nx+3,j,k) =  Wr(3+1,j,k)
          Wr(Nx+4,j,k) =  Wr(4+1,j,k)

          Wr( 0,j,k)   =  Wr(Nx-1-0,j,k)
          Wr(-1,j,k)   =  Wr(Nx-1-1,j,k)
          Wr(-2,j,k)   =  Wr(Nx-1-2,j,k)
          Wr(-3,j,k)   =  Wr(Nx-1-3,j,k)
          Wr(-4,j,k)   =  Wr(Nx-1-4,j,k)

          !--- Tr ---
          Tr(Nx+0,j,k) =  Tr(0+1,j,k)
          Tr(Nx+1,j,k) =  Tr(1+1,j,k)
          Tr(Nx+2,j,k) =  Tr(2+1,j,k)
          Tr(Nx+3,j,k) =  Tr(3+1,j,k)
          Tr(Nx+4,j,k) =  Tr(4+1,j,k)

          Tr( 0,j,k)   =  Tr(Nx-1-0,j,k)
          Tr(-1,j,k)   =  Tr(Nx-1-1,j,k)
          Tr(-2,j,k)   =  Tr(Nx-1-2,j,k)
          Tr(-3,j,k)   =  Tr(Nx-1-3,j,k)
          Tr(-4,j,k)   =  Tr(Nx-1-4,j,k)

          !--- U ---
          U (Nx+0,j,k) =  U (0+1,j,k)
          U (Nx+1,j,k) =  U (1+1,j,k)
          U (Nx+2,j,k) =  U (2+1,j,k)
          U (Nx+3,j,k) =  U (3+1,j,k)
          U (Nx+4,j,k) =  U (4+1,j,k)

          U ( 0,j,k)   =  U (Nx-1-0,j,k)
          U (-1,j,k)   =  U (Nx-1-1,j,k)
          U (-2,j,k)   =  U (Nx-1-2,j,k)
          U (-3,j,k)   =  U (Nx-1-3,j,k)
          U (-4,j,k)   =  U (Nx-1-4,j,k)

          !--- V ---
          V (Nx+0,j,k) =  V (0+1,j,k)
          V (Nx+1,j,k) =  V (1+1,j,k)
          V (Nx+2,j,k) =  V (2+1,j,k)
          V (Nx+3,j,k) =  V (3+1,j,k)
          V (Nx+4,j,k) =  V (4+1,j,k)

          V ( 0,j,k)   =  V (Nx-1-0,j,k)
          V (-1,j,k)   =  V (Nx-1-1,j,k)
          V (-2,j,k)   =  V (Nx-1-2,j,k)
          V (-3,j,k)   =  V (Nx-1-3,j,k)
          V (-4,j,k)   =  V (Nx-1-4,j,k)

          !--- W ---
          W (Nx+0,j,k) =  W (0+1,j,k)
          W (Nx+1,j,k) =  W (1+1,j,k)
          W (Nx+2,j,k) =  W (2+1,j,k)
          W (Nx+3,j,k) =  W (3+1,j,k)
          W (Nx+4,j,k) =  W (4+1,j,k)

          W ( 0,j,k)   =  W (Nx-1-0,j,k)
          W (-1,j,k)   =  W (Nx-1-1,j,k)
          W (-2,j,k)   =  W (Nx-1-2,j,k)
          W (-3,j,k)   =  W (Nx-1-3,j,k)
          W (-4,j,k)   =  W (Nx-1-4,j,k)

          !--- T ---
          T (Nx+0,j,k) =  T (0+1,j,k)
          T (Nx+1,j,k) =  T (1+1,j,k)
          T (Nx+2,j,k) =  T (2+1,j,k)
          T (Nx+3,j,k) =  T (3+1,j,k)
          T (Nx+4,j,k) =  T (4+1,j,k)

          T ( 0,j,k)   =  T (Nx-1-0,j,k)
          T (-1,j,k)   =  T (Nx-1-1,j,k)
          T (-2,j,k)   =  T (Nx-1-2,j,k)
          T (-3,j,k)   =  T (Nx-1-3,j,k)
          T (-4,j,k)   =  T (Nx-1-4,j,k)

          !--- P ---
          P (Nx+0,j,k) =  P (0+1,j,k)
          P (Nx+1,j,k) =  P (1+1,j,k)
          P (Nx+2,j,k) =  P (2+1,j,k)
          P (Nx+3,j,k) =  P (3+1,j,k)
          P (Nx+4,j,k) =  P (4+1,j,k)

          P ( 0,j,k)   =  P (Nx-1-0,j,k)
          P (-1,j,k)   =  P (Nx-1-1,j,k)
          P (-2,j,k)   =  P (Nx-1-2,j,k)
          P (-3,j,k)   =  P (Nx-1-3,j,k)
          P (-4,j,k)   =  P (Nx-1-4,j,k)

          !--- Mu ---
          Mu(Nx+0,j,k) =  Mu(0+1,j,k)
          Mu(Nx+1,j,k) =  Mu(1+1,j,k)
          Mu(Nx+2,j,k) =  Mu(2+1,j,k)
          Mu(Nx+3,j,k) =  Mu(3+1,j,k)
          Mu(Nx+4,j,k) =  Mu(4+1,j,k)
                          
          Mu( 0,j,k)   =  Mu(Nx-1-0,j,k)
          Mu(-1,j,k)   =  Mu(Nx-1-1,j,k)
          Mu(-2,j,k)   =  Mu(Nx-1-2,j,k)
          Mu(-3,j,k)   =  Mu(Nx-1-3,j,k)
          Mu(-4,j,k)   =  Mu(Nx-1-4,j,k)

          !--- kT ---
         ! kT(Nx+0,j,k) =  kT(0+1,j,k)
         ! kT(Nx+1,j,k) =  kT(1+1,j,k)
         ! kT(Nx+2,j,k) =  kT(2+1,j,k)
         ! kT(Nx+3,j,k) =  kT(3+1,j,k)
         ! kT(Nx+4,j,k) =  kT(4+1,j,k)
         !                 
         ! kT( 0,j,k)   =  kT(Nx-1-0,j,k)
         ! kT(-1,j,k)   =  kT(Nx-1-1,j,k)
         ! kT(-2,j,k)   =  kT(Nx-1-2,j,k)
         ! kT(-3,j,k)   =  kT(Nx-1-3,j,k)
         ! kT(-4,j,k)   =  kT(Nx-1-4,j,k)

          !--- Ysr ---
          Ysr(Nx+0,j,k) =  Ysr(0+1,j,k)
          Ysr(Nx+1,j,k) =  Ysr(1+1,j,k)
          Ysr(Nx+2,j,k) =  Ysr(2+1,j,k)
          Ysr(Nx+3,j,k) =  Ysr(3+1,j,k)
          Ysr(Nx+4,j,k) =  Ysr(4+1,j,k)

          Ysr( 0,j,k)   =  Ysr(Nx-1-0,j,k)
          Ysr(-1,j,k)   =  Ysr(Nx-1-1,j,k)
          Ysr(-2,j,k)   =  Ysr(Nx-1-2,j,k)
          Ysr(-3,j,k)   =  Ysr(Nx-1-3,j,k)
          Ysr(-4,j,k)   =  Ysr(Nx-1-4,j,k)

          !--- Ys ---
          Ys(Nx+0,j,k) =  Ys(0+1,j,k)
          Ys(Nx+1,j,k) =  Ys(1+1,j,k)
          Ys(Nx+2,j,k) =  Ys(2+1,j,k)
          Ys(Nx+3,j,k) =  Ys(3+1,j,k)
          Ys(Nx+4,j,k) =  Ys(4+1,j,k)

          Ys( 0,j,k)   =  Ys(Nx-1-0,j,k)
          Ys(-1,j,k)   =  Ys(Nx-1-1,j,k)
          Ys(-2,j,k)   =  Ys(Nx-1-2,j,k)
          Ys(-3,j,k)   =  Ys(Nx-1-3,j,k)
          Ys(-4,j,k)   =  Ys(Nx-1-4,j,k)

        End Do; End Do

        Return
      End Subroutine


!=============================================================================================================================
!=====[ Fluid variables z direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_z_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr )

        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluid variables y direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_y_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Rh )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ur )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vr )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wr )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Tr )

        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ysr )


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluid variables x direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_Tau_x_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !====== ã´äEèåèÇÃê›íË ======
        Do k = ks, ke, 1
        Do j = js, je, 1
          !--- T_xx ---
          T_xx(Nx+0,j,k) =  T_xx(0+1,j,k)
          T_xx(Nx+1,j,k) =  T_xx(1+1,j,k)
          T_xx(Nx+2,j,k) =  T_xx(2+1,j,k)
          T_xx(Nx+3,j,k) =  T_xx(3+1,j,k)
          T_xx(Nx+4,j,k) =  T_xx(4+1,j,k)

          T_xx( 0,j,k)   =  T_xx(Nx-1-0,j,k)
          T_xx(-1,j,k)   =  T_xx(Nx-1-1,j,k)
          T_xx(-2,j,k)   =  T_xx(Nx-1-2,j,k)
          T_xx(-3,j,k)   =  T_xx(Nx-1-3,j,k)
          T_xx(-4,j,k)   =  T_xx(Nx-1-4,j,k)

          !--- T_xy ---
          T_xy(Nx+0,j,k) =  T_xy(0+1,j,k)
          T_xy(Nx+1,j,k) =  T_xy(1+1,j,k)
          T_xy(Nx+2,j,k) =  T_xy(2+1,j,k)
          T_xy(Nx+3,j,k) =  T_xy(3+1,j,k)
          T_xy(Nx+4,j,k) =  T_xy(4+1,j,k)

          T_xy( 0,j,k)   =  T_xy(Nx-1-0,j,k)
          T_xy(-1,j,k)   =  T_xy(Nx-1-1,j,k)
          T_xy(-2,j,k)   =  T_xy(Nx-1-2,j,k)
          T_xy(-3,j,k)   =  T_xy(Nx-1-3,j,k)
          T_xy(-4,j,k)   =  T_xy(Nx-1-4,j,k)

          !--- T_xz ---
          T_xz(Nx+0,j,k) =  T_xz(0+1,j,k)
          T_xz(Nx+1,j,k) =  T_xz(1+1,j,k)
          T_xz(Nx+2,j,k) =  T_xz(2+1,j,k)
          T_xz(Nx+3,j,k) =  T_xz(3+1,j,k)
          T_xz(Nx+4,j,k) =  T_xz(4+1,j,k)

          T_xz( 0,j,k)   =  T_xz(Nx-1-0,j,k)
          T_xz(-1,j,k)   =  T_xz(Nx-1-1,j,k)
          T_xz(-2,j,k)   =  T_xz(Nx-1-2,j,k)
          T_xz(-3,j,k)   =  T_xz(Nx-1-3,j,k)
          T_xz(-4,j,k)   =  T_xz(Nx-1-4,j,k)

          !--- T_yx ---
          T_yx(Nx+0,j,k) =  T_yx(0+1,j,k)
          T_yx(Nx+1,j,k) =  T_yx(1+1,j,k)
          T_yx(Nx+2,j,k) =  T_yx(2+1,j,k)
          T_yx(Nx+3,j,k) =  T_yx(3+1,j,k)
          T_yx(Nx+4,j,k) =  T_yx(4+1,j,k)

          T_yx( 0,j,k)   =  T_yx(Nx-1-0,j,k)
          T_yx(-1,j,k)   =  T_yx(Nx-1-1,j,k)
          T_yx(-2,j,k)   =  T_yx(Nx-1-2,j,k)
          T_yx(-3,j,k)   =  T_yx(Nx-1-3,j,k)
          T_yx(-4,j,k)   =  T_yx(Nx-1-4,j,k)

          !--- T_yy ---
          T_yy(Nx+0,j,k) =  T_yy(0+1,j,k)
          T_yy(Nx+1,j,k) =  T_yy(1+1,j,k)
          T_yy(Nx+2,j,k) =  T_yy(2+1,j,k)
          T_yy(Nx+3,j,k) =  T_yy(3+1,j,k)
          T_yy(Nx+4,j,k) =  T_yy(4+1,j,k)

          T_yy( 0,j,k)   =  T_yy(Nx-1-0,j,k)
          T_yy(-1,j,k)   =  T_yy(Nx-1-1,j,k)
          T_yy(-2,j,k)   =  T_yy(Nx-1-2,j,k)
          T_yy(-3,j,k)   =  T_yy(Nx-1-3,j,k)
          T_yy(-4,j,k)   =  T_yy(Nx-1-4,j,k)

          !--- T_yz ---
          T_yz(Nx+0,j,k) =  T_yz(0+1,j,k)
          T_yz(Nx+1,j,k) =  T_yz(1+1,j,k)
          T_yz(Nx+2,j,k) =  T_yz(2+1,j,k)
          T_yz(Nx+3,j,k) =  T_yz(3+1,j,k)
          T_yz(Nx+4,j,k) =  T_yz(4+1,j,k)

          T_yz( 0,j,k)   =  T_yz(Nx-1-0,j,k)
          T_yz(-1,j,k)   =  T_yz(Nx-1-1,j,k)
          T_yz(-2,j,k)   =  T_yz(Nx-1-2,j,k)
          T_yz(-3,j,k)   =  T_yz(Nx-1-3,j,k)
          T_yz(-4,j,k)   =  T_yz(Nx-1-4,j,k)

          !--- T_zx ---
          T_zx(Nx+0,j,k) =  T_zx(0+1,j,k)
          T_zx(Nx+1,j,k) =  T_zx(1+1,j,k)
          T_zx(Nx+2,j,k) =  T_zx(2+1,j,k)
          T_zx(Nx+3,j,k) =  T_zx(3+1,j,k)
          T_zx(Nx+4,j,k) =  T_zx(4+1,j,k)

          T_zx( 0,j,k)   =  T_zx(Nx-1-0,j,k)
          T_zx(-1,j,k)   =  T_zx(Nx-1-1,j,k)
          T_zx(-2,j,k)   =  T_zx(Nx-1-2,j,k)
          T_zx(-3,j,k)   =  T_zx(Nx-1-3,j,k)
          T_zx(-4,j,k)   =  T_zx(Nx-1-4,j,k)

          !--- T_zy ---
          T_zy(Nx+0,j,k) =  T_zy(0+1,j,k)
          T_zy(Nx+1,j,k) =  T_zy(1+1,j,k)
          T_zy(Nx+2,j,k) =  T_zy(2+1,j,k)
          T_zy(Nx+3,j,k) =  T_zy(3+1,j,k)
          T_zy(Nx+4,j,k) =  T_zy(4+1,j,k)

          T_zy( 0,j,k)   =  T_zy(Nx-1-0,j,k)
          T_zy(-1,j,k)   =  T_zy(Nx-1-1,j,k)
          T_zy(-2,j,k)   =  T_zy(Nx-1-2,j,k)
          T_zy(-3,j,k)   =  T_zy(Nx-1-3,j,k)
          T_zy(-4,j,k)   =  T_zy(Nx-1-4,j,k)

          !--- T_zz ---
          T_zz(Nx+0,j,k) =  T_zz(0+1,j,k)
          T_zz(Nx+1,j,k) =  T_zz(1+1,j,k)
          T_zz(Nx+2,j,k) =  T_zz(2+1,j,k)
          T_zz(Nx+3,j,k) =  T_zz(3+1,j,k)
          T_zz(Nx+4,j,k) =  T_zz(4+1,j,k)

          T_zz( 0,j,k)   =  T_zz(Nx-1-0,j,k)
          T_zz(-1,j,k)   =  T_zz(Nx-1-1,j,k)
          T_zz(-2,j,k)   =  T_zz(Nx-1-2,j,k)
          T_zz(-3,j,k)   =  T_zz(Nx-1-3,j,k)
          T_zz(-4,j,k)   =  T_zz(Nx-1-4,j,k)

        End Do; End Do

        Return
      End Subroutine


!=============================================================================================================================
!=====[ Fluid variables z direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_Tau_z_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_xx )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_xy )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_xz )

        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_yx )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_yy )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_yz )

        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_zx )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_zy )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, T_zz )


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluid variables y direction ]=========================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_Tau_y_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_xx )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_xy )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_xz )

        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_yx )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_yy )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_yz )

        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_zx )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_zy )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5,  T_zz )


        Return
      End Subroutine
!=============================================================================================================================
!=====[ BC_Periodic_y_dir_MPI_R8 ]============================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_y_dir_MPI_R8( n11 , n12 , n21 , n22 , n31 , n32 , F1 )

        Use module_mpi
        Use Definition
        Implicit none
        Integer,intent(in) :: n11 , n21 , n31
        Integer,intent(in) :: n12 , n22 , n32
        Real*8,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
	integer ierror,ierr
        integer dum_len
        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        If(j_sta==1.or.j_end==Ny-1) then
          Allocate(dumcomzx_r8(1:Nx-1,5,ks:ke))
          Allocate(dumcomzx_s8(1:Nx-1,5,ks:ke))
          Do icom=0,Ndiv_Nz-1
            If(k_sta==kksta(icom)) then

              If(j_end==Ny-1) then
                Do k = ks,ke
                Do i = 1, Nx-1, 1
                  dumcomzx_s8(i,1,k) = F1(i,Ny-5,k)
                  dumcomzx_s8(i,2,k) = F1(i,Ny-4,k)
                  dumcomzx_s8(i,3,k) = F1(i,Ny-3,k)
                  dumcomzx_s8(i,4,k) = F1(i,Ny-2,k)
                  dumcomzx_s8(i,5,k) = F1(i,Ny-1,k)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
                Call mpi_isend(dumcomzx_s8(1,1,ks),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(0,icom),1,MPI_COMM_WORLD,isend(1),ierr)
                Call mpi_irecv(dumcomzx_r8(1,1,ks),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(0,icom),1,MPI_COMM_WORLD,isend(2),ierr)
              End If

              If(j_sta==1) then
                Do k = ks, ke
                Do i = 1, Nx-1, 1
                  dumcomzx_s8(i,1,k) = F1(i,1,k)
                  dumcomzx_s8(i,2,k) = F1(i,2,k)
                  dumcomzx_s8(i,3,k) = F1(i,3,k)
                  dumcomzx_s8(i,4,k) = F1(i,4,k)
                  dumcomzx_s8(i,5,k) = F1(i,5,k)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
                Call mpi_isend(dumcomzx_s8(1,1,ks),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(2),ierr)
                Call mpi_irecv(dumcomzx_r8(1,1,ks),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(1),ierr)
              End if

              Call mpi_wait(isend(1),istatus,ierr)
              Call mpi_wait(isend(2),istatus,ierr)

            End if
          End Do


          If(j_end==Ny-1) then
            Do k = ks, ke
            Do i = 1, Nx-1, 1
              F1(i,Ny  ,k)=dumcomzx_r8(i,1,k)
              F1(i,Ny+1,k)=dumcomzx_r8(i,2,k)
              F1(i,Ny+2,k)=dumcomzx_r8(i,3,k)
              F1(i,Ny+3,k)=dumcomzx_r8(i,4,k)
              F1(i,Ny+4,k)=dumcomzx_r8(i,5,k)
            End Do; End Do 
          End If

          If(j_sta==1) then
            Do k = ks, ke
            Do i = 1, Nx-1, 1
              F1(i, 0,k)=dumcomzx_r8(i,5,k)
              F1(i,-1,k)=dumcomzx_r8(i,4,k)
              F1(i,-2,k)=dumcomzx_r8(i,3,k)
              F1(i,-3,k)=dumcomzx_r8(i,2,k)
              F1(i,-4,k)=dumcomzx_r8(i,1,k)
            End Do; End Do
          End if

          Deallocate(dumcomzx_r8)
          Deallocate(dumcomzx_s8)
        End If


        Return
      End Subroutine

!=============================================================================================================================
!=====[ BC_Periodic_z_dir_MPI_R8 ]============================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_z_dir_MPI_R8( n11 , n12 , n21 , n22 , n31 , n32 , F1 )

        Use module_mpi
        Use Definition
        Implicit none
        Integer,intent(in) :: n11 , n21 , n31
        Integer,intent(in) :: n12 , n22 , n32
        Real*8,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
	integer ierror,ierr
        integer dum_len
        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        If(k_sta==1.or.k_end==Nz-1) then
          Allocate(dumcomxy_r8(1:Nx-1,js:je,5))
          Allocate(dumcomxy_s8(1:Nx-1,js:je,5))
          Do icom=0,Ndiv_Ny-1
            If(j_sta==jjsta(icom)) then

              If(k_end==Nz-1) then
                Do j = js,je
                Do i = 1, Nx-1, 1
                  dumcomxy_s8(i,j,1) = F1(i,j,Nz-5)
                  dumcomxy_s8(i,j,2) = F1(i,j,Nz-4)
                  dumcomxy_s8(i,j,3) = F1(i,j,Nz-3)
                  dumcomxy_s8(i,j,4) = F1(i,j,Nz-2)
                  dumcomxy_s8(i,j,5) = F1(i,j,Nz-1)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(je-js+1)*5
                Call mpi_isend(dumcomxy_s8(1,js,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(icom,0),1,MPI_COMM_WORLD,isend(1),ierr)
                Call mpi_irecv(dumcomxy_r8(1,js,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(icom,0),1,MPI_COMM_WORLD,isend(2),ierr)
              End If

              If(k_sta==1) then
                Do j = js, je
                Do i = 1, Nx-1, 1
                  dumcomxy_s8(i,j,1) = F1(i,j,1)
                  dumcomxy_s8(i,j,2) = F1(i,j,2)
                  dumcomxy_s8(i,j,3) = F1(i,j,3)
                  dumcomxy_s8(i,j,4) = F1(i,j,4)
                  dumcomxy_s8(i,j,5) = F1(i,j,5)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(je-js+1)*5
                Call mpi_isend(dumcomxy_s8(1,js,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(2),ierr)
                Call mpi_irecv(dumcomxy_r8(1,js,1),dum_len,MPI_DOUBLE_PRECISION,                        &
                                                  itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(1),ierr)
              End if

              Call mpi_wait(isend(1),istatus,ierr)
              Call mpi_wait(isend(2),istatus,ierr)

            End if
          End Do


          If(k_end==Nz-1) then
            Do j = js, je
            Do i = 1, Nx-1, 1
              F1(i,j,Nz  )=dumcomxy_r8(i,j,1)
              F1(i,j,Nz+1)=dumcomxy_r8(i,j,2)
              F1(i,j,Nz+2)=dumcomxy_r8(i,j,3)
              F1(i,j,Nz+3)=dumcomxy_r8(i,j,4)
              F1(i,j,Nz+4)=dumcomxy_r8(i,j,5)
            End Do; End Do
          End If

          If(k_sta==1) then
            Do j = js, je
            Do i = 1, Nx-1, 1
              F1(i,j, 0)=dumcomxy_r8(i,j,5)
              F1(i,j,-1)=dumcomxy_r8(i,j,4)
              F1(i,j,-2)=dumcomxy_r8(i,j,3)
              F1(i,j,-3)=dumcomxy_r8(i,j,2)
              F1(i,j,-4)=dumcomxy_r8(i,j,1)
            End Do; End Do
          End if

          Deallocate(dumcomxy_r8)
          Deallocate(dumcomxy_s8)
        End If

        Return
      End Subroutine
!=============================================================================================================================
!=====[ BC_Periodic_z_dir_MPI_R4 ]============================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_z_dir_MPI_R4( n11 , n12 , n21 , n22 , n31 , n32 , F1 )

        Use module_mpi
        Use Definition
        Implicit none
        Integer,intent(in) :: n11 , n21 , n31
        Integer,intent(in) :: n12 , n22 , n32
        Real*4,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
	integer ierror,ierr
        integer dum_len
        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        If(k_sta==1.or.k_end==Nz-1) then
          Allocate(dumcomxy_r4(1:Nx-1,js:je,5))
          Allocate(dumcomxy_s4(1:Nx-1,js:je,5))
          Do icom=0,Ndiv_Ny-1
            If(j_sta==jjsta(icom)) then

              If(k_end==Nz-1) then
                Do j = js,je
                Do i = 1, Nx-1, 1
                  dumcomxy_s4(i,j,1) = F1(i,j,Nz-5)
                  dumcomxy_s4(i,j,2) = F1(i,j,Nz-4)
                  dumcomxy_s4(i,j,3) = F1(i,j,Nz-3)
                  dumcomxy_s4(i,j,4) = F1(i,j,Nz-2)
                  dumcomxy_s4(i,j,5) = F1(i,j,Nz-1)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(je-js+1)*5
                Call mpi_isend(dumcomxy_s4(1,js,1),dum_len,MPI_REAL,itable(icom,0),1,MPI_COMM_WORLD,isend(1),ierr)
                Call mpi_irecv(dumcomxy_r4(1,js,1),dum_len,MPI_REAL,itable(icom,0),1,MPI_COMM_WORLD,isend(2),ierr)
              End If

              If(k_sta==1) then
                Do j = js, je
                Do i = 1, Nx-1, 1
                  dumcomxy_s4(i,j,1) = F1(i,j,1)
                  dumcomxy_s4(i,j,2) = F1(i,j,2)
                  dumcomxy_s4(i,j,3) = F1(i,j,3)
                  dumcomxy_s4(i,j,4) = F1(i,j,4)
                  dumcomxy_s4(i,j,5) = F1(i,j,5)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(je-js+1)*5
                Call mpi_isend(dumcomxy_s4(1,js,1),dum_len,MPI_REAL,itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(2),ierr)
                Call mpi_irecv(dumcomxy_r4(1,js,1),dum_len,MPI_REAL,itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(1),ierr)
              End if

              Call mpi_wait(isend(1),istatus,ierr)
              Call mpi_wait(isend(2),istatus,ierr)

            End if
          End Do


          If(k_end==Nz-1) then
            Do j = js, je
            Do i = 1, Nx-1, 1
              F1(i,j,Nz  )=dumcomxy_r4(i,j,1)
              F1(i,j,Nz+1)=dumcomxy_r4(i,j,2)
              F1(i,j,Nz+2)=dumcomxy_r4(i,j,3)
              F1(i,j,Nz+3)=dumcomxy_r4(i,j,4)
              F1(i,j,Nz+4)=dumcomxy_r4(i,j,5)
            End Do; End Do
          End If

          If(k_sta==1) then
            Do j = js, je
            Do i = 1, Nx-1, 1
              F1(i,j, 0)=dumcomxy_r4(i,j,5)
              F1(i,j,-1)=dumcomxy_r4(i,j,4)
              F1(i,j,-2)=dumcomxy_r4(i,j,3)
              F1(i,j,-3)=dumcomxy_r4(i,j,2)
              F1(i,j,-4)=dumcomxy_r4(i,j,1)
            End Do; End Do
          End if

          Deallocate(dumcomxy_r4)
          Deallocate(dumcomxy_s4)
        End If

        Return
      End Subroutine
!=============================================================================================================================
!=====[ BC_Periodic_y_dir_MPI_R4 ]============================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_y_dir_MPI_R4( n11 , n12 , n21 , n22 , n31 , n32 , F1 )

        Use module_mpi
        Use Definition
        Implicit none
        Integer,intent(in) :: n11 , n21 , n31
        Integer,intent(in) :: n12 , n22 , n32
        Real*4,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
	integer ierror,ierr
        integer dum_len
        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        If(j_sta==1.or.j_end==Ny-1) then
          Allocate(dumcomxy_r4(1:Nx-1,5,ks:ke))
          Allocate(dumcomxy_s4(1:Nx-1,5,ks:ke))
          Do icom=0,Ndiv_Nz-1
            If(k_sta==kksta(icom)) then

              If(j_end==Ny-1) then
                Do k = ks,ke
                Do i = 1, Nx-1, 1
                  dumcomxy_s4(i,1,k) = F1(i,Ny-5,k)
                  dumcomxy_s4(i,2,k) = F1(i,Ny-4,k)
                  dumcomxy_s4(i,3,k) = F1(i,Ny-3,k)
                  dumcomxy_s4(i,4,k) = F1(i,Ny-2,k)
                  dumcomxy_s4(i,5,k) = F1(i,Ny-1,k)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
                Call mpi_isend(dumcomxy_s4(1,1,ks),dum_len,MPI_REAL,itable(0,icom),1,MPI_COMM_WORLD,isend(1),ierr)
                Call mpi_irecv(dumcomxy_r4(1,1,ks),dum_len,MPI_REAL,itable(0,icom),1,MPI_COMM_WORLD,isend(2),ierr)
              End If

              If(j_sta==1) then
                Do k = ks, ke
                Do i = 1, Nx-1, 1
                  dumcomxy_s4(i,1,k) = F1(i,1,k)
                  dumcomxy_s4(i,2,k) = F1(i,2,k)
                  dumcomxy_s4(i,3,k) = F1(i,3,k)
                  dumcomxy_s4(i,4,k) = F1(i,4,k)
                  dumcomxy_s4(i,5,k) = F1(i,5,k)
                End Do; End Do
                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
                Call mpi_isend(dumcomxy_s4(1,1,ks),dum_len,MPI_REAL,itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(2),ierr)
                Call mpi_irecv(dumcomxy_r4(1,1,ks),dum_len,MPI_REAL,itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(1),ierr)
              End if

              Call mpi_wait(isend(1),istatus,ierr)
              Call mpi_wait(isend(2),istatus,ierr)

            End if
          End Do


          If(j_end==Ny-1) then
            Do k = ks, ke
            Do i = 1, Nx-1, 1
              F1(i,Ny  ,k)=dumcomxy_r4(i,1,k)
              F1(i,Ny+1,k)=dumcomxy_r4(i,2,k)
              F1(i,Ny+2,k)=dumcomxy_r4(i,3,k)
              F1(i,Ny+3,k)=dumcomxy_r4(i,4,k)
              F1(i,Ny+4,k)=dumcomxy_r4(i,5,k)
            End Do; End Do
          End If

          If(j_sta==1) then
            Do k = ks, ke
            Do i = 1, Nx-1, 1
              F1(i, 0,k)=dumcomxy_r4(i,5,k)
              F1(i,-1,k)=dumcomxy_r4(i,4,k)
              F1(i,-2,k)=dumcomxy_r4(i,3,k)
              F1(i,-3,k)=dumcomxy_r4(i,2,k)
              F1(i,-4,k)=dumcomxy_r4(i,1,k)
            End Do; End Do
          End if

          Deallocate(dumcomxy_r4)
          Deallocate(dumcomxy_s4)
        End If

        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluxes x direction ]==================================================================================================
!=============================================================================================================================
      Subroutine BC_Periodic_E_xyz

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî

        Call BC_E_xyz_x_para
        Call BC_E_xyz_y_para
        Call BC_E_xyz_z_para


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluxes x direction ]==================================================================================================
!=============================================================================================================================
      Subroutine BC_E_xyz_x_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !====== ã´äEèåèÇÃê›íË ======
        Do k = ks, ke, 1
        Do j = js, je, 1
          !--- E_x ---
          E_x(Nx+0,j,k) =  E_x(0+1,j,k)
          E_x(Nx+1,j,k) =  E_x(1+1,j,k)
          E_x(Nx+2,j,k) =  E_x(2+1,j,k)
          E_x(Nx+3,j,k) =  E_x(3+1,j,k)
          E_x(Nx+4,j,k) =  E_x(4+1,j,k)

          E_x( 0,j,k)   =  E_x(Nx-1-0,j,k)
          E_x(-1,j,k)   =  E_x(Nx-1-1,j,k)
          E_x(-2,j,k)   =  E_x(Nx-1-2,j,k)
          E_x(-3,j,k)   =  E_x(Nx-1-3,j,k)
          E_x(-4,j,k)   =  E_x(Nx-1-4,j,k)

          !--- E_y ---
          E_y(Nx+0,j,k) =  E_y(0+1,j,k)
          E_y(Nx+1,j,k) =  E_y(1+1,j,k)
          E_y(Nx+2,j,k) =  E_y(2+1,j,k)
          E_y(Nx+3,j,k) =  E_y(3+1,j,k)
          E_y(Nx+4,j,k) =  E_y(4+1,j,k)

          E_y( 0,j,k)   =  E_y(Nx-1-0,j,k)
          E_y(-1,j,k)   =  E_y(Nx-1-1,j,k)
          E_y(-2,j,k)   =  E_y(Nx-1-2,j,k)
          E_y(-3,j,k)   =  E_y(Nx-1-3,j,k)
          E_y(-4,j,k)   =  E_y(Nx-1-4,j,k)

          !--- E_z ---
          E_z(Nx+0,j,k) =  E_z(0+1,j,k)
          E_z(Nx+1,j,k) =  E_z(1+1,j,k)
          E_z(Nx+2,j,k) =  E_z(2+1,j,k)
          E_z(Nx+3,j,k) =  E_z(3+1,j,k)
          E_z(Nx+4,j,k) =  E_z(4+1,j,k)

          E_z( 0,j,k)   =  E_z(Nx-1-0,j,k)
          E_z(-1,j,k)   =  E_z(Nx-1-1,j,k)
          E_z(-2,j,k)   =  E_z(Nx-1-2,j,k)
          E_z(-3,j,k)   =  E_z(Nx-1-3,j,k)
          E_z(-4,j,k)   =  E_z(Nx-1-4,j,k)
        End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluxes y direction ]==================================================================================================
!=============================================================================================================================
      Subroutine BC_E_xyz_y_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî


        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_x )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_y )
        Call BC_Periodic_y_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_z )


        Return
      End Subroutine
!=============================================================================================================================
!=====[ Fluxes z direction ]==================================================================================================
!=============================================================================================================================
      Subroutine BC_E_xyz_z_para

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî


        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_x )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_y )
        Call BC_Periodic_z_dir_MPI_R4( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, E_z )


        Return
      End Subroutine
!=============================================================================================================================
!*************  ïœêîåvéZ  ****************************************************************************************************
!=============================================================================================================================
      Subroutine Calc_Var

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          c1_Rh(i,j,k) = 1.0d0 / Rh(i,j,k)
        End Do; End Do; End Do

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          U (i,j,k) = Ur (i,j,k) * c1_Rh(i,j,k)
          V (i,j,k) = Vr (i,j,k) * c1_Rh(i,j,k)
          W (i,j,k) = Wr (i,j,k) * c1_Rh(i,j,k)
          T (i,j,k) = Tr (i,j,k) * c1_Rh(i,j,k)

          Ys(i,j,k) = Ysr(i,j,k) * c1_Rh(i,j,k)

          P  (i,j,k) = Tr(i,j,k)
          Mu (i,j,k) = mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)
         ! kT (i,j,k) = Mu(i,j,k)
         ! Dys(i,j,k) = Mu(i,j,k)
        End Do; End Do; End Do


        Return
      End Subroutine
!=============================================================================================================================
!*************  ïœêîåvéZ  ****************************************************************************************************
!=============================================================================================================================
      Subroutine Calc_Var_BC_MPI_r4

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !---------------------------------------------
        Do k = ks-4, ks-1, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          c1_Rh(i,j,k) = 1.0d0 / Rh(i,j,k)
        End Do; End Do; End Do

        Do k = ks-4, ks-1, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          U (i,j,k) = Ur (i,j,k) * c1_Rh(i,j,k)
          V (i,j,k) = Vr (i,j,k) * c1_Rh(i,j,k)
          W (i,j,k) = Wr (i,j,k) * c1_Rh(i,j,k)
          T (i,j,k) = Tr (i,j,k) * c1_Rh(i,j,k)

          Ys(i,j,k) = Ysr(i,j,k) * c1_Rh(i,j,k)

          P  (i,j,k) = Tr(i,j,k)
          Mu (i,j,k) = mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)
         ! kT (i,j,k) = Mu(i,j,k)
         ! Dys(i,j,k) = Mu(i,j,k)
        End Do; End Do; End Do

        !---------------------------------------------
        Do k = ke+1, ke+4, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          c1_Rh(i,j,k) = 1.0d0 / Rh(i,j,k)
        End Do; End Do; End Do

        Do k = ke+1, ke+4, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          U (i,j,k) = Ur (i,j,k) * c1_Rh(i,j,k)
          V (i,j,k) = Vr (i,j,k) * c1_Rh(i,j,k)
          W (i,j,k) = Wr (i,j,k) * c1_Rh(i,j,k)
          T (i,j,k) = Tr (i,j,k) * c1_Rh(i,j,k)

          Ys(i,j,k) = Ysr(i,j,k) * c1_Rh(i,j,k)

          P  (i,j,k) = Tr(i,j,k)
          Mu (i,j,k) = mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)
         ! kT (i,j,k) = Mu(i,j,k)
         ! Dys(i,j,k) = Mu(i,j,k)
        End Do; End Do; End Do

        !---------------------------------------------
        Do k = ks, ke, 1
        Do j = js-4, js-1, 1
        Do i = 1, Nx-1, 1
          c1_Rh(i,j,k) = 1.0d0 / Rh(i,j,k)
        End Do; End Do; End Do

        Do k = ks, ke, 1
        Do j = js-4, js-1, 1
        Do i = 1, Nx-1, 1
          U (i,j,k) = Ur (i,j,k) * c1_Rh(i,j,k)
          V (i,j,k) = Vr (i,j,k) * c1_Rh(i,j,k)
          W (i,j,k) = Wr (i,j,k) * c1_Rh(i,j,k)
          T (i,j,k) = Tr (i,j,k) * c1_Rh(i,j,k)

          Ys(i,j,k) = Ysr(i,j,k) * c1_Rh(i,j,k)

          P  (i,j,k) = Tr(i,j,k)
          Mu (i,j,k) = mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)
         ! kT (i,j,k) = Mu(i,j,k)
         ! Dys(i,j,k) = Mu(i,j,k)
        End Do; End Do; End Do

        !---------------------------------------------
        Do k = ks, ke, 1
        Do j = je+1, je+4, 1
        Do i = 1, Nx-1, 1
          c1_Rh(i,j,k) = 1.0d0 / Rh(i,j,k)
        End Do; End Do; End Do

        Do k = ks, ke, 1
        Do j = je+1, je+4, 1
        Do i = 1, Nx-1, 1
          U (i,j,k) = Ur (i,j,k) * c1_Rh(i,j,k)
          V (i,j,k) = Vr (i,j,k) * c1_Rh(i,j,k)
          W (i,j,k) = Wr (i,j,k) * c1_Rh(i,j,k)
          T (i,j,k) = Tr (i,j,k) * c1_Rh(i,j,k)

          Ys(i,j,k) = Ysr(i,j,k) * c1_Rh(i,j,k)

          P  (i,j,k) = Tr(i,j,k)
          Mu (i,j,k) = mu0 * ( ( T(i,j,k)*c1_T0 )**(1.5d0) ) * (T0+C_su)/(T(i,j,k)+C_su)
         ! kT (i,j,k) = Mu(i,j,k)
         ! Dys(i,j,k) = Mu(i,j,k)
        End Do; End Do; End Do
        !---------------------------------------------


        Return
      End Subroutine
!=============================================================================================================================
!*************  î˜ï™åvéZ  ****************************************************************************************************
!=============================================================================================================================
      Subroutine Calc_Derivative

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî




       !--- x dir ---
       Call Sub_Derivative_MPI( 1  ,   Ur  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dUrdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   Rh  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dRhdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       Call Sub_Derivative_MPI( 1  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dUdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   V   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dVdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   W   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dWdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   P   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dPdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   T   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dTdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 1  ,   Ys  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dYsdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )


       !--- y dir ---
       Call Sub_Derivative_MPI( 2  ,   Vr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dVrdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   Rh  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dRhdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       Call Sub_Derivative_MPI( 2  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dUdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   V   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dVdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   W   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dWdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   P   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dPdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   T   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dTdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 2  ,   Ys  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dYsdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )


       !--- z dir ---
       Call Sub_Derivative_MPI( 3  ,   Wr  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dWrdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   Rh  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dRhdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

       Call Sub_Derivative_MPI( 3  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dUdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   V   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dVdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   W   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dWdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   P   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dPdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   T   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dTdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
       Call Sub_Derivative_MPI( 3  ,   Ys  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,      &
                                      dYsdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )


       ! Call Sub_dF1dx_Global_MPI(  Ur,  dUrdx )
       ! Call Sub_dF1dx_Global_MPI(  Rh,  dRhdx )
       ! Call Sub_dF1dx_Local_MPI (  U ,  dUdx  )
       ! Call Sub_dF1dx_Local_MPI (  V ,  dVdx  )
       ! Call Sub_dF1dx_Local_MPI (  W ,  dWdx  )
       ! Call Sub_dF1dx_Local_MPI (  P ,  dPdx  )
       ! Call Sub_dF1dx_Local_MPI (  T ,  dTdx  )
       ! Call Sub_dF1dx_Local_MPI ( Ys ,  dYsdx )

       ! Call Sub_dF1dy_Global_MPI(  Vr,  dVrdy )
       ! Call Sub_dF1dy_Global_MPI(  Rh,  dRhdy )
       ! Call Sub_dF1dy_Local_MPI (  U ,  dUdy  )
       ! Call Sub_dF1dy_Local_MPI (  V ,  dVdy  )
       ! Call Sub_dF1dy_Local_MPI (  W ,  dWdy  )
       ! Call Sub_dF1dy_Local_MPI (  P ,  dPdy  )
       ! Call Sub_dF1dy_Local_MPI (  T ,  dTdy  )
       ! Call Sub_dF1dy_Local_MPI ( Ys ,  dYsdy )

       ! Call Sub_dF1dz_Global_MPI(  Wr,  dWrdz )
       ! Call Sub_dF1dz_Global_MPI(  Rh,  dRhdz )
       ! Call Sub_dF1dz_Local_MPI (  U ,  dUdz  )
       ! Call Sub_dF1dz_Local_MPI (  V ,  dVdz  )
       ! Call Sub_dF1dz_Local_MPI (  W ,  dWdz  )
       ! Call Sub_dF1dz_Local_MPI (  P ,  dPdz  )
       ! Call Sub_dF1dz_Local_MPI (  T ,  dTdz  )
       ! Call Sub_dF1dz_Local_MPI ( Ys ,  dYsdz )

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Filter_F1( F1 , n11 , n12 , n21 , n22 , n31 , n32  )
        Use module_mpi
        Use Definition
        Implicit none   
        Integer,intent(in)   :: n11 , n21 , n31
        Integer,intent(in)   :: n12 , n22 , n32
        Real*4,intent(inout) ::  F1(n11 :n12 ,n21 :n22 ,n31 :n32 )
        Real*4                   F2(n11 :n12 ,n21 :n22 ,n31 :n32 )
        Real*4                   F3(n11 :n12 ,n21 :n22 ,n31 :n32 )

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- BC in x dir ---
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          F1(Nx+0,j,k) =  F1(0+1,j,k)
          F1(Nx+1,j,k) =  F1(1+1,j,k)
          F1(Nx+2,j,k) =  F1(2+1,j,k)
          F1(Nx+3,j,k) =  F1(3+1,j,k)
          F1(Nx+4,j,k) =  F1(4+1,j,k)

          F1( 0,j,k)   =  F1(Nx-1-0,j,k)
          F1(-1,j,k)   =  F1(Nx-1-1,j,k)
          F1(-2,j,k)   =  F1(Nx-1-2,j,k)
          F1(-3,j,k)   =  F1(Nx-1-3,j,k)
          F1(-4,j,k)   =  F1(Nx-1-4,j,k)
        End Do; End Do; End Do

        !--- x dir ---
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          F2(i,j,k) = F1(i,j,k)                                    &
              - (  + c252_1024* ( F1(i  ,j,k)               )      &
                   - c210_1024* ( F1(i-1,j,k) + F1(i+1,j,k) )      &
                   + c120_1024* ( F1(i-2,j,k) + F1(i+2,j,k) )      &
                   -  c45_1024* ( F1(i-3,j,k) + F1(i+3,j,k) )      &
                   +  c10_1024* ( F1(i-4,j,k) + F1(i+4,j,k) )      &
                   -   c1_1024* ( F1(i-5,j,k) + F1(i+5,j,k) )  )
        End Do; End Do; End Do



        !--- z dir ---
        Call BC_Periodic_z_dir_MPI_R4(n11,n12,n21,n22,n31,n32, F2)
        Call mp_send_recv_pre_jk(F2,5,3, n11,n12,n21,n22,n31,n32 )

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          F3(i,j,k) = F2(i,j,k)                                      &
                - (  + c252_1024* ( F2(i,j,k  )               )      &
                     - c210_1024* ( F2(i,j,k-1) + F2(i,j,k+1) )      &
                     + c120_1024* ( F2(i,j,k-2) + F2(i,j,k+2) )      &
                     -  c45_1024* ( F2(i,j,k-3) + F2(i,j,k+3) )      &
                     +  c10_1024* ( F2(i,j,k-4) + F2(i,j,k+4) )      &
                     -   c1_1024* ( F2(i,j,k-5) + F2(i,j,k+5) )  )
        End Do; End Do; End Do





        !--- y dir ---
        Call BC_Periodic_y_dir_MPI_R4(n11,n12,n21,n22,n31,n32, F3)
        Call mp_send_recv_pre_jk(F3,5,2, n11,n12,n21,n22,n31,n32 )

        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          F1(i,j,k) = F3(i,j,k)                                    &
              - (  + c252_1024* ( F3(i,j  ,k)               )      &
                   - c210_1024* ( F3(i,j-1,k) + F3(i,j+1,k) )      &
                   + c120_1024* ( F3(i,j-2,k) + F3(i,j+2,k) )      &
                   -  c45_1024* ( F3(i,j-3,k) + F3(i,j+3,k) )      &
                   +  c10_1024* ( F3(i,j-4,k) + F3(i,j+4,k) )      &
                   -   c1_1024* ( F3(i,j-5,k) + F3(i,j+5,k) )  )
        End Do; End Do; End Do


        !--- F2 -> F1 ---
 !      Do k = 1, Nz-1, 1
 !      Do j = 1, Ny-1, 1
 !      Do i = 1, Nx-1, 1
 !        F1(i,j,k) = F2(i,j,k)
 !      End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Sub_Derivative_MPI(idir,   F1, n11 , n12 , n21 , n22 , n31 , n32 , &
                                           dF1, n11d, n12d, n21d, n22d, n31d, n32d )

        Use module_mpi
        Use Definition
        Implicit none   
        Integer,intent(in) :: idir
        Integer,intent(in) :: n11 , n21 , n31
        Integer,intent(in) :: n12 , n22 , n32
        Integer,intent(in) :: n11d, n21d, n31d
        Integer,intent(in) :: n12d, n22d, n32d
        Real*4,intent( in) ::  F1(n11 :n12 ,n21 :n22 ,n31 :n32 )
        Real*4,intent(out) :: dF1(n11d:n12d,n21d:n22d,n31d:n32d)

        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !=== x dir ==========================================================================
        If(idir==1) then
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1, Nx-1, 1
             dF1(i,j,k)     = (     c4_5*( F1(i+1,j,k) - F1(i-1,j,k) )           &
                                -   c1_5*( F1(i+2,j,k) - F1(i-2,j,k) )           &
                                + c4_105*( F1(i+3,j,k) - F1(i-3,j,k) )           &
                                - c1_280*( F1(i+4,j,k) - F1(i-4,j,k) )  )*c1_dx
          End Do; End Do; End Do

        End If



        !=== y dir ==========================================================================
        If(idir==2) then
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1, Nx-1, 1
             dF1  (i,j,k)   = (     c4_5*( F1(i,j+1,k) - F1(i,j-1,k) )           &
                                -   c1_5*( F1(i,j+2,k) - F1(i,j-2,k) )           &
                                + c4_105*( F1(i,j+3,k) - F1(i,j-3,k) )           &
                                - c1_280*( F1(i,j+4,k) - F1(i,j-4,k) )  )*c1_dy
          End Do; End Do; End Do

        End If


        !=== z dir ==========================================================================
        If(idir==3) then
          Do j = js, je, 1
          Do k = ks, ke, 1
          Do i = 1, Nx-1, 1
             dF1  (i,j,k)   = (     c4_5*( F1(i,j,k+1) - F1(i,j,k-1) )           &
                                -   c1_5*( F1(i,j,k+2) - F1(i,j,k-2) )           &
                                + c4_105*( F1(i,j,k+3) - F1(i,j,k-3) )           &
                                - c1_280*( F1(i,j,k+4) - F1(i,j,k-4) )  )*c1_dz
          End Do; End Do; End Do

        End If



        Return
      End Subroutine

!=============================================================================================================================
!=============================================================================================================================
      Subroutine Calc_Tau

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Real*4 dF1dx, dF1dy, dcF1_F2dx, dcF1_F2dy

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

       ! Do k = 1, Nz-1, 1
       ! Do j = 1, Ny-1, 1
       ! Do i = 1, Nx-1, 1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1

          T_xx(i,j,k) = c2_3*mu(i,j,k)*( 2.0d0*dUdx(i,j,k) - dVdy(i,j,k) - dWdz(i,j,k) )
          T_xy(i,j,k) =      mu(i,j,k)*(       dUdy(i,j,k) + dVdx(i,j,k) )
          T_xz(i,j,k) =      mu(i,j,k)*(       dUdz(i,j,k) + dWdx(i,j,k) )

          T_yx(i,j,k) =      mu(i,j,k)*(       dVdx(i,j,k) + dUdy(i,j,k) )
          T_yy(i,j,k) = c2_3*mu(i,j,k)*( 2.0d0*dVdy(i,j,k) - dWdz(i,j,k) - dUdx(i,j,k) )
          T_yz(i,j,k) =      mu(i,j,k)*(       dVdz(i,j,k) + dWdy(i,j,k) )

          T_zx(i,j,k) =      mu(i,j,k)*(       dWdx(i,j,k) + dUdz(i,j,k) )
          T_zy(i,j,k) =      mu(i,j,k)*(       dWdy(i,j,k) + dVdz(i,j,k) )
          T_zz(i,j,k) = c2_3*mu(i,j,k)*( 2.0d0*dWdz(i,j,k) - dUdx(i,j,k) - dVdy(i,j,k) )

        End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Calc_dtime

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        dtime_c = dtime0

        !--- dt is determined by the velocity and the sound velocity ---
        If( Ityp_dt .eq. 0 ) then
          Do k = ks, ke, 1
            If( k == (Nz-1)/2 ) then
              Do j = js, je, 1
                If( j == (Ny-1)/2 ) then
                  Do i = 1, Nx-1, 1


          !  dtime_c = min(   dtime_c,                                                              &
          !                 1.0d0/(    abs(U(i,j,k)) / dx                                           &
          !                          + abs(V(i,j,k)) / dy                                           &
          !                          + abs(W(i,j,k)) / dz                                           &
          !                        + Sqrt( gh*R_gas*T_HIT/(U_ref*U_ref)                             &

            dtime_c = min(   dtime_c,                                                                  &
                           1.0d0/(    abs(U(i,j,k)) / dx                                               &
                                    + abs(V(i,j,k)) / dy                                               &
                                    + abs(W(i,j,k)) / dz                                               &
                                  + Sqrt( gh*Tr(i,j,k)/Rh(i,j,k)*cRT_UU                                 &
                                         *( c1_dx*c1_dx + c1_dy    *c1_dy     + c1_dz*c1_dz) ) )  )

                  End Do
                End If
              End Do
            End If
          End Do
        End If

        If( Ityp_dt .eq. 1 ) then
          Do k = ks, ke, 1
            If( k == (Nz-1)/2 ) then
              Do j = js, je, 1
                If( j == (Ny-1)/2 ) then
                  Do i = 1, Nx-1, 1
            dtime_c = min(   dtime_c,                                                              &
                           1.0d0/(    abs(U(i,j,k)) / dx                                           &
                                    + abs(V(i,j,k)) / dy                                           &
                                    + abs(W(i,j,k)) / dz   )  )

                  End Do
                End If
              End Do
            End If
          End Do
        End If


        Call mp_allminr( dtime_c )

        dtime = Cour*dtime_c

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Linear_Forcing

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Real*4 Sxx_i, Sxy_i, Sxz_i
        Real*4 Syx_i, Syy_i, Syz_i
        Real*4 Szx_i, Szy_i, Szz_i
        Real*4 Epsi_i

        Real*4 Epsi_d,Epsi_s
        Real*4 Ken_d,Ken_s
        Real*4 PD_term
        Real*4 f1,f2,f3

        Real*4, dimension(:,:,:),allocatable :: dT_xxdx,dT_xydy,dT_xzdz
        Real*4, dimension(:,:,:),allocatable :: dT_yxdx,dT_yydy,dT_yzdz
        Real*4, dimension(:,:,:),allocatable :: dT_zxdx,dT_zydy,dT_zzdz

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        If( Ityp_Force == 0 ) then
          B_Ceff = 0.0d0
          C_d    = 0.0d0
          C_s    = 0.0d0
        End If

        If( Ityp_Force == 1 ) then
          !--- For output data ---
          !Call Calc_Var
          !Call BC_Periodic
          !Call Calc_Derivative
          !Call Calc_Tau
          !Call BC_Periodic_Tau
          !-----------------------
          !allocate(     dT_xxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_xydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_xzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

          !Call Sub_Derivative_MPI( 1  ,   T_xx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_xy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_xz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 1  ,   T_yx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_yy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_yz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 1  ,   T_zx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_zy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_zz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

          !===ïΩãœ===
          UdPdxsumT_xyz   = 0.0d0
          VdPdysumT_xyz   = 0.0d0
          WdPdzsumT_xyz   = 0.0d0
          EpsisumT_xyz    = 0.0d0
          rkensumT_xyz    = 0.0d0

          Do j = js, je
            Do k = ks, ke
            Do i = 1, Nx-1

               !Epsi_i = - ( + U(i,j,k)*( + dT_xxdx(i,j,k)     &
               !                          + dT_xydy(i,j,k)     &
               !                          + dT_xzdz(i,j,k)  )  &
               !             + V(i,j,k)*( + dT_yxdx(i,j,k)     &
               !                          + dT_yydy(i,j,k)     &
               !                          + dT_yzdz(i,j,k)  )  &
               !             + W(i,j,k)*( + dT_zxdx(i,j,k)     &
               !                          + dT_zydy(i,j,k)     &
               !                          + dT_zzdz(i,j,k)  )  ) / Re 
               !EpsisumT_xyz = EpsisumT_xyz + (   Epsi_i   ) * c1_NxNyNz

              !===ïΩãœ===
            UdPdxsumT_xyz =   UdPdxsumT_xyz + U(i,j,k)*( dPdx(i,j,k) ) * c1_NxNyNz
            VdPdysumT_xyz =   VdPdysumT_xyz + V(i,j,k)*( dPdy(i,j,k) ) * c1_NxNyNz
            WdPdzsumT_xyz =   WdPdzsumT_xyz + W(i,j,k)*( dPdz(i,j,k) ) * c1_NxNyNz

              rkensumT_xyz = rkensumT_xyz + (     Rh(i,j,k)*U(i,j,k)*U(i,j,k)                      &
                                               +  Rh(i,j,k)*V(i,j,k)*V(i,j,k)                      &
                                               +  Rh(i,j,k)*W(i,j,k)*W(i,j,k)   ) * c1_NxNyNz

            End Do; End Do
          End Do


          !---/* MPI */---
          !---------------------------------------------------
          !Call mp_allsumr(   EpsisumT_xyz )
          Call mp_allsumr(  UdPdxsumT_xyz )
          Call mp_allsumr(  VdPdysumT_xyz )
          Call mp_allsumr(  WdPdzsumT_xyz )
          Call mp_allsumr(   rkensumT_xyz )

          B_Ceff = (E_target + UdPdxsumT_xyz+VdPdysumT_xyz+WdPdzsumT_xyz)/rkensumT_xyz

        End If

        If( Ityp_Force == 2 ) then

         Call Helmholtz_decomposition

          !--- For output data ---
          !Call Calc_Var
          !Call BC_Periodic
          !Call Calc_Derivative
          !Call Calc_Tau
          !Call BC_Periodic_Tau
          !-----------------------

          !===dissipation calculation===
          !allocate(     dT_xxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_xydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_xzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_yzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          !allocate(     dT_zzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

          !Call Sub_Derivative_MPI( 1  ,   T_xx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_xy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_xz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_xzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 1  ,   T_yx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_yy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_yz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_yzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 1  ,   T_zx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 2  ,   T_zy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          !Call Sub_Derivative_MPI( 3  ,   T_zz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
          !                               dT_zzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

          !EpsisumT_xyz    = 0.0d0

          !Do j = js, je
          !  Do k = ks, ke
          !  Do i = 1, Nx-1

          !     Epsi_i = - ( + U(i,j,k)*( + dT_xxdx(i,j,k) + dT_xydy(i,j,k) + dT_xzdz(i,j,k)  )  &
          !                  + V(i,j,k)*( + dT_yxdx(i,j,k) + dT_yydy(i,j,k) + dT_yzdz(i,j,k)  )  &
          !                  + W(i,j,k)*( + dT_zxdx(i,j,k) + dT_zydy(i,j,k) + dT_zzdz(i,j,k)  )  ) / Re

          !    EpsisumT_xyz = EpsisumT_xyz + (   Epsi_i   ) * c1_NxNyNz

          !  End Do; End Do
          !End Do

          !Call mp_allsumr(   EpsisumT_xyz )

          Epsi_d = E_target*Ed_Es/(1.0d0+Ed_Es)
          Epsi_s = E_target*1.0d0/(1.0d0+Ed_Es)
          !---------------------------------------------------

          !===Kinetic energy===
          rkensumT_xyz    = 0.0d0
          Ken_d           = 0.0d0
          Ken_s           = 0.0d0

          Do j = js, je
            Do k = ks, ke
            Do i = 1, Nx-1

              Ken_d = Ken_d + (     Ud(i,j,k)*Ud(i,j,k)                  &
                                  + Vd(i,j,k)*Vd(i,j,k)                  &
                                  + Wd(i,j,k)*Wd(i,j,k)  ) * c1_NxNyNz    

              Ken_s = Ken_s + (     Us(i,j,k)*Us(i,j,k)                  &
                                  + Vs(i,j,k)*Vs(i,j,k)                  &
                                  + Ws(i,j,k)*Ws(i,j,k)  ) * c1_NxNyNz    

            End Do; End Do
          End Do

          Call mp_allsumr( Ken_d )
          Call mp_allsumr( Ken_s )

          Ken_d = 0.5d0*Ken_d
          Ken_s = 0.5d0*Ken_s

          !===plessure dilatation correlation===
          UdPdxsumT_xyz   = 0.0d0
          VdPdysumT_xyz   = 0.0d0
          WdPdzsumT_xyz   = 0.0d0

          Do j = js, je
            Do k = ks, ke
            Do i = 1, Nx-1

            UdPdxsumT_xyz = UdPdxsumT_xyz + U(i,j,k)*( dPdx(i,j,k) ) * c1_NxNyNz
            VdPdysumT_xyz = VdPdysumT_xyz + V(i,j,k)*( dPdy(i,j,k) ) * c1_NxNyNz
            WdPdzsumT_xyz = WdPdzsumT_xyz + W(i,j,k)*( dPdz(i,j,k) ) * c1_NxNyNz

            End Do; End Do
          End Do

          Call mp_allsumr(  UdPdxsumT_xyz )
          Call mp_allsumr(  VdPdysumT_xyz )
          Call mp_allsumr(  WdPdzsumT_xyz )

          PD_term = -(UdPdxsumT_xyz + VdPdysumT_xyz + WdPdzsumT_xyz)

          !---------------------------------------------------

          C_d  = (Epsi_d - PD_term)/(2.0d0*Ken_d)
          C_s  =  Epsi_s           /(2.0d0*Ken_s)

          !---------------------------------------------------

          !    PDsumT_xyz = 0.0d0
          ! ForcesumT_xyz = 0.0d0
          ! UdPdxsumT_xyz = 0.0d0
          ! VdPdysumT_xyz = 0.0d0
          ! WdPdzsumT_xyz = 0.0d0
          !  rkensumT_xyz = 0.0d0
          !  EpsisumT_xyz = 0.0d0

          ! Do k = ks, ke
          ! Do j = js, je
          !   Do i = 1, Nx-1

          !      Epsi_i = - ( + U(i,j,k)*( + dT_xxdx(i,j,k) + dT_xydy(i,j,k) + dT_xzdz(i,j,k)  )  &
          !                   + V(i,j,k)*( + dT_yxdx(i,j,k) + dT_yydy(i,j,k) + dT_yzdz(i,j,k)  )  &
          !                   + W(i,j,k)*( + dT_zxdx(i,j,k) + dT_zydy(i,j,k) + dT_zzdz(i,j,k)  )  ) / Re

          !      f1 = + c_d*sqrt(Rh(i,j,k))*Ud(i,j,k)   &
          !           + c_s*sqrt(Rh(i,j,k))*Us(i,j,k)    
          !      f2 = + c_d*sqrt(Rh(i,j,k))*Vd(i,j,k)   &
          !           + c_s*sqrt(Rh(i,j,k))*Vs(i,j,k)    
          !      f3 = + c_d*sqrt(Rh(i,j,k))*Wd(i,j,k)   &
          !           + c_s*sqrt(Rh(i,j,k))*Ws(i,j,k)    

          !      dUddx(i,j,k) = f1*U(i,j,k)
          !      dVddy(i,j,k) = f2*V(i,j,k)
          !      dWddz(i,j,k) = f3*W(i,j,k)
          !      dUsdx(i,j,k) = - (gh-1.0d0)*(f1*U(i,j,k)+f2*V(i,j,k)+f3*W(i,j,k))

          !   !===ïΩãœ===
          !      PDsumT_xyz =    PDsumT_xyz + ( + P(i,j,k)*( + dUdx(i,j,k)                  &
          !                                               + dVdy(i,j,k)                  &
          !                                               + dWdz(i,j,k) )  )* c1_NxNyNz   

          !   ForcesumT_xyz = ForcesumT_xyz + ( f1*U(i,j,k) + f2*V(i,j,k) + f3*W(i,j,k) ) * c1_NxNyNz   

          !   UdPdxsumT_xyz = UdPdxsumT_xyz + U(i,j,k)*( dPdx(i,j,k) ) * c1_NxNyNz
          !   VdPdysumT_xyz = VdPdysumT_xyz + V(i,j,k)*( dPdy(i,j,k) ) * c1_NxNyNz
          !   WdPdzsumT_xyz = WdPdzsumT_xyz + W(i,j,k)*( dPdz(i,j,k) ) * c1_NxNyNz
          !    EpsisumT_xyz =  EpsisumT_xyz + (   Epsi_i   ) * c1_NxNyNz


          !    rkensumT_xyz =  rkensumT_xyz + (  +  0.5d0*Rh(i,j,k)*U(i,j,k)*U(i,j,k)                      &
          !                                      +  0.5d0*Rh(i,j,k)*V(i,j,k)*V(i,j,k)                      &
          !                                      +  0.5d0*Rh(i,j,k)*W(i,j,k)*W(i,j,k)   ) * c1_NxNyNz
          !   End Do; End Do
          ! End Do

          ! Call mp_allsumr(     PDsumT_xyz )
          ! Call mp_allsumr(  ForcesumT_xyz )
          ! Call mp_allsumr(   EpsisumT_xyz )
          ! Call mp_allsumr(  UdPdxsumT_xyz )
          ! Call mp_allsumr(  VdPdysumT_xyz )
          ! Call mp_allsumr(  WdPdzsumT_xyz )
          ! Call mp_allsumr(   rkensumT_xyz )

          ! If( my_rank==root ) then

          !   Write(89,777)                                          &  
          !   Dble(nt),                                              &!1
          !   time,                                                  &!2
          !   !---Epsi---                                               
          !        PDsumT_xyz                                     ,  &!3
          !    -( UdPdxsumT_xyz + VdPdysumT_xyz + WdPdzsumT_xyz ) ,  &!4
          !        EpsisumT_xyz                                   ,  &!5
          !        ForcesumT_xyz                                  ,  &!6
          !        PDsumT_xyz-EpsisumT_xyz+ForcesumT_xyz          ,  &!7
          !        rkensumT_xyz                                   ,  &!8
          !        Ken_d                                          ,  &!9
          !        Ken_s                                          ,  &!9
          !        (PDsumT_xyz-EpsisumT_xyz+ForcesumT_xyz)*dtime      !10
          ! End If

        End If


        Return
      End Subroutine
!=============================================================================================================================
!****************Å@Helmholtz decomposition  **********************************************************************************
!=============================================================================================================================
      Subroutine Helmholtz_decomposition

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k         ! ïœêî
        Integer n_HD            ! ïœêî
        Real*4 :: Rh_sqrt                                            !ñßìx
        Real*4 :: Usd,Vsd,Wsd                                        !ñßìx
        Real*4, dimension(:,:,:),allocatable :: Potential,Source     !ñßìx

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        allocate(    Potential(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(       Source(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- ä÷êîÇÃèÄîı ---
        !---------------------------------------------

        !---É”ÇÃèÄîı---
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Potential (i,j,k) = 0.0d0
         End Do; End Do; End Do

        !---source term(âEï”)ÇÃèÄîı---

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           U(i,j,k) = sqrt(Rh(i,j,k))*U(i,j,k)      !U(i,j,k)=Å„Éœ*U
           V(i,j,k) = sqrt(Rh(i,j,k))*V(i,j,k)      !V(i,j,k)=Å„Éœ*V
           W(i,j,k) = sqrt(Rh(i,j,k))*W(i,j,k)      !W(i,j,k)=Å„Éœ*W
         End Do; End Do; End Do

         UsumT_xyz     = 0.0d0
         VsumT_xyz     = 0.0d0
         WsumT_xyz     = 0.0d0

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           UsumT_xyz =   UsumT_xyz + ( U(i,j,k) ) * c1_NxNyNz
           VsumT_xyz =   VsumT_xyz + ( V(i,j,k) ) * c1_NxNyNz
           WsumT_xyz =   WsumT_xyz + ( W(i,j,k) ) * c1_NxNyNz
         End Do; End Do; End Do

         Call mp_allsumr( UsumT_xyz )
         Call mp_allsumr( VsumT_xyz )
         Call mp_allsumr( WsumT_xyz )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           U(i,j,k) =  U(i,j,k) - UsumT_xyz
           V(i,j,k) =  V(i,j,k) - VsumT_xyz
           W(i,j,k) =  W(i,j,k) - WsumT_xyz
         End Do; End Do; End Do

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, U )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, V )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, W )

         Call Sub_Derivative_MPI( 1  ,   U  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dUdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Call Sub_Derivative_MPI( 2  ,   V  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dVdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Call Sub_Derivative_MPI( 3  ,   W  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dWdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k)
         End Do; End Do; End Do

        !Call Poisson_SOR(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call Poisson_EQ_CG(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

        !--- ë¨ìxÇÃî≠éUê¨ï™åvéZ ---
        !---------------------------------------------

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Potential )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Ud(i,j,k) = - ( Potential(i+1,j  ,k  ) - Potential(i-1,j  ,k  ) )/(2.0d0*dx)   !Ud(i,j,k)=(Å„Éœ*U)d
           Vd(i,j,k) = - ( Potential(i  ,j+1,k  ) - Potential(i  ,j-1,k  ) )/(2.0d0*dy)   !Vd(i,j,k)=(Å„Éœ*V)d
           Wd(i,j,k) = - ( Potential(i  ,j  ,k+1) - Potential(i  ,j  ,k-1) )/(2.0d0*dz)   !Wd(i,j,k)=(Å„Éœ*W)d
           Us(i,j,k) = U(i,j,k) - Ud(i,j,k)                                               !Us(i,j,k)=(Å„Éœ*U)s
           Vs(i,j,k) = V(i,j,k) - Vd(i,j,k)                                               !Vs(i,j,k)=(Å„Éœ*V)s
           Ws(i,j,k) = W(i,j,k) - Wd(i,j,k)                                               !Ws(i,j,k)=(Å„Éœ*W)s
         End Do; End Do; End Do

        !--- ã´äEÇÃílåvéZ ---
        !---------------------------------------------

          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ud )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Us )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vs )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ws )

        !---------------------------------------------

        !---É”sÇÃèÄîı---
         !Do k = ks, ke, 1
         !Do j = js, je, 1
         !Do i = 1, Nx-1, 1
         !  Potential (i,j,k) = 0.0d0
         !End Do; End Do; End Do

         !Do n_HD = 1, 10

         !  !---source term(âEï”)ÇÃèÄîı---

         !   Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Us )
         !   Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vs )
         !   Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ws )

         !   Call Sub_Derivative_MPI( 1  ,   Us  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                  dUsdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !   Call Sub_Derivative_MPI( 2  ,   Vs  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                  dVsdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !   Call Sub_Derivative_MPI( 3  ,   Ws  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                  dWsdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

         !   Do k = ks, ke, 1
         !   Do j = js, je, 1
         !   Do i = 1, Nx-1, 1
         !     Source    (i,j,k) = dUsdx(i,j,k) + dVsdy(i,j,k) + dWsdz(i,j,k)
         !   End Do; End Do; End Do

         !  !Call Poisson_SOR(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
         !  Call Poisson_EQ_CG(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

         !  !--- ë¨ìxÇÃî≠éUê¨ï™åvéZ ---
         !  !---------------------------------------------

         !   Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Potential )

         !   Do k = ks, ke, 1
         !   Do j = js, je, 1
         !   Do i = 1, Nx-1, 1
         !     Usd =  - ( Potential(i+1,j  ,k  ) - Potential(i-1,j  ,k  ) )/(2.0d0*dx)
         !     Vsd =  - ( Potential(i  ,j+1,k  ) - Potential(i  ,j-1,k  ) )/(2.0d0*dy)
         !     Wsd =  - ( Potential(i  ,j  ,k+1) - Potential(i  ,j  ,k-1) )/(2.0d0*dz)
         !     Ud(i,j,k) = Ud(i,j,k) + Usd
         !     Vd(i,j,k) = Vd(i,j,k) + Vsd
         !     Wd(i,j,k) = Wd(i,j,k) + Wsd
         !     Us(i,j,k) = Us(i,j,k) - Usd
         !     Vs(i,j,k) = Vs(i,j,k) - Vsd
         !     Ws(i,j,k) = Ws(i,j,k) - Wsd
         !   End Do; End Do; End Do

         !End Do
         !---------------------------------------------

         !--- ã´äEÇÃílåvéZ ---
         !---------------------------------------------

          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ud )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Us )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vs )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ws )

         !---------------------------------------------
         !Call Sub_Derivative_MPI( 1  ,  Ud  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dUddx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !Call Sub_Derivative_MPI( 2  ,  Vd  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dVddy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !Call Sub_Derivative_MPI( 3  ,  Wd  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dWddz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !Call Sub_Derivative_MPI( 1  ,  Us  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dUsdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !Call Sub_Derivative_MPI( 2  ,  Vs  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dVsdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !Call Sub_Derivative_MPI( 3  ,  Ws  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,    &
         !                              dWsdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5     )
         !---------------------------------------------

         !If( mod(Nt,ifl5) == 0 ) then

         !  Call Sub_Derivative_MPI( 1  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                 dUdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 2  ,   V   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                 dVdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 3  ,   W   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                 dWdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 1  ,  Ud   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dUddx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 2  ,  Vd   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dVddy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 3  ,  Wd   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dWddz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 1  ,  Us   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dUsdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 2  ,  Vs   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dVsdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         !  Call Sub_Derivative_MPI( 3  ,  Ws   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
         !                                dWsdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

         !  Call write_xy_P
         !  Call write_yz_P

         !End If

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Rh_sqrt = 1.0d0/sqrt(Rh(i,j,k))
           U(i,j,k) = U(i,j,k)*Rh_sqrt             !U(i,j,k)=U
           V(i,j,k) = V(i,j,k)*Rh_sqrt             !V(i,j,k)=V
           W(i,j,k) = W(i,j,k)*Rh_sqrt             !W(i,j,k)=W
         End Do; End Do; End Do

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, U )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, V )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, W )

         Call Sub_Derivative_MPI( 1  ,   U   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                        dUdx , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         Call Sub_Derivative_MPI( 2  ,   V   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                        dVdy , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
         Call Sub_Derivative_MPI( 3  ,   W   , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                        dWdz , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

         deallocate( Potential )
         deallocate( Source    )

        Return
      End Subroutine
!==================================================================================================================================
!==================================================================================================================================
!=========================================================================================================================
!========= [ Poisson Solver based on the BiCGStab method  ] ==============================================================
!=========================================================================================================================
!  BiCGSTAB method
! 
! [Algorithm] Solve Ax = b with BiCGSTAB 
! 1) Compute r0 = b - Ax (Choose r0* arbitary, r0* = r0, (r0, r0*) != 0)
!            p = r = r0
!    (c1 = (r0*, r0) = (r0, r0))
! 2) Loop...
!    c2    = (r0*, Ap)
!    alpha = c1 / c2
!    e     = r - alpha*Ap
!    c3    = (e, Ae) / (Ae, Ae)
! 
!    x_new = x + alpha*p + c3*e
!    r_new = e - c3*Ae
!    if(EPS > ||r||/||b||) break
! 
!    c1    = (r_new, r0*) = (r_new, r0)
!    beta  = c1 / (c2*c3)
!    p_new = r_new + beta*(p - c3*Ap)
! 
      Subroutine Poisson_EQ_CG(Phi, source, n11,n12,n21,n22,n31,n32)

        Use module_mpi
        Use Definition
        Implicit none
        Real*4 ,intent(inout) :: Phi   (n11:n12,n21:n22,n31:n32)
        Real*4 ,intent(inout) :: source(n11:n12,n21:n22,n31:n32)
        Integer,intent(in   ) :: n11,n12,n21,n22,n31,n32
        Integer i, j, k, iter 

        !--- BiCGStab ---
        Integer,Parameter :: IterLimit =  50
        Real*4 :: r_n,s_n,t_n
        Real*4 :: r_d,s_d,t_d
        Real*4 :: f_sum
        !Real*4 :: a_CG, a_CG_n, a_CG_d,a_CG_u
        !Real*4 :: b_CG, b_CG_n, b_CG_d,b_CG_u
        !Real*4 :: w_CG, w_CG_n, w_CG_d,w_CG_u
        Real*8 :: a_CG, a_CG_n, a_CG_d,a_CG_u
        Real*8 :: b_CG, b_CG_n, b_CG_d,b_CG_u
        Real*8 :: w_CG, w_CG_n, w_CG_d,w_CG_u
!NECs
        Real*8 :: w_dUr(-2:Nx+2,j_sta-3:j_end+3,k_sta-3:k_end+3)
        Real*8 :: w_dTr(-2:Nx+2,j_sta-3:j_end+3,k_sta-3:k_end+3)
!NECe

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !======== Poissonï˚íˆéÆ  Åﬁ^2É” = -g  ÇÃâñ@ =============================
        !========Å@Å¶É”:ãÅÇﬂÇΩÇ¢ä÷êî(target)  g:source term ===================

        !--- dW = b = Div U / dt --- 
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          dWr(i,j,k) = -source(i,j,k)
          !dW(i,j,k) = (   ( c1_24*U(i-2,j,k) - c9_8*U(i-1,j,k) + c9_8*U(i,j,k) - c1_24*U(i+1,j,k) )*c1_dx                 &
          !              + (                  -      V(i,j-1,k) +      V(i,j,k)                    )*c1_dy(j)              &
          !              + ( c1_24*W(i,j,k-2) - c9_8*W(i,j,k-1) + c9_8*W(i,j,k) - c1_24*W(i,j,k+1) )*c1_dz    )*c1_dtime
        End Do; End Do; End Do

        !Call BC_Preiodic_xz1_Grad0_y_R4( -2,Nx+2,j_sta-3,j_end+3,k_sta-3,k_end+3, P )
        Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Phi )

        !--- dV = r0 = b - Ax --- 
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
            dVr(i,j,k) = dWr(i,j,k)                                                                       &
                         - ( + ( Phi(i+1,j  ,k  ) - 2.0d0*Phi(i,j,k) + Phi(i-1,j  ,k  )) * (c1_dx*c1_dx)  &
                             + ( Phi(i  ,j+1,k  ) - 2.0d0*Phi(i,j,k) + Phi(i  ,j-1,k  )) * (c1_dy*c1_dy)  &
                             + ( Phi(i  ,j  ,k+1) - 2.0d0*Phi(i,j,k) + Phi(i  ,j  ,k-1)) * (c1_dz*c1_dz)  )
        End Do; End Do; End Do

        !--- dYsr = r0*   = r0 --- 
        !--- dUr  = p     = r0 --- 
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1,Nx-1, 1
          dUr (i,j,k) = dVr(i,j,k)
          dYsr(i,j,k) = dVr(i,j,k)
        End Do; End Do; End Do

        a_CG = 0.0d0 ; w_CG = 0.0d0 ; b_CG = 0.0d0
        Iter_end = IterLimit

        !=== Computing \fai ===
        Do iter = 1, IterLimit, 1

          a_CG_n = 0.0d0 ; a_CG_d = 0.0d0
          w_CG_n = 0.0d0 ; w_CG_d = 0.0d0
          b_CG_n = 0.0d0 ; b_CG_d = 0.0d0

          r_n = 0.0d0 ; s_n = 0.0d0
          r_d = 0.0d0 ; s_d = 0.0d0

          !=== Computing a_CG ===
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, dUr )

          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1,Nx-1, 1
            a_CG_n = a_CG_n + ( dYsr(i,j,k) * dVr(i,j,k) )
!NECs
            w_dUr(i,j,k)=    ( dUr(i+1,j  ,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
                          +  ( dUr(i  ,j+1,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
                          +  ( dUr(i  ,j  ,k+1) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j  ,k-1)) * (c1_dz*c1_dz)     
            a_CG_d = a_CG_d + dYsr(i,j,k) * w_dUr(i,j,k)
!NECe
          End Do; End Do; End Do

!NECs
          Call mp_allsumr8_2( a_CG_n,a_CG_d )
!NECe

          !--- a_CG_n = alpha = c1 / c2 ---
          a_CG = a_CG_n / a_CG_d

          !=== Computing t_k+1 ===
          !--- dTr = dV - a_CG* A(dUr)  e = r - alpha*Ap ---
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1,Nx-1, 1
!NECs
            dTr(i,j,k) = dVr(i,j,k) - a_CG * w_dUr(i,j,k)
!NECe
          End Do; End Do; End Do

          !=== Computing w_CG ===
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, dTr )

          !--- w_CG_n = ( e, Ae) ---
          !--- w_CG_d = (Ae, Ae) ---
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1,Nx-1, 1
            w_dTr(i,j,k) =     ( dTr(i+1,j  ,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
                            +  ( dTr(i  ,j+1,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
                            +  ( dTr(i  ,j  ,k+1) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j  ,k-1)) * (c1_dz*c1_dz)     
            w_CG_n = w_CG_n + dTr(i,j,k) * w_dTr(i,j,k)
            w_CG_d = w_CG_d + w_dTr(i,j,k)**2
          End Do; End Do; End Do
!NECs
          Call mp_allsumr8_2( w_CG_n,w_CG_d )
!NECe

          !--- w_CG = c3 = ( e, Ae)/(Ae, Ae) ---
          w_CG = w_CG_n / w_CG_d

          !=== Computing X_k+1 ===
          !--- P = x_new = x + alpha*p + c3*e ---
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1,Nx-1, 1
            Phi(i,j,k) = Phi(i,j,k) + a_CG * dUr(i,j,k) + w_CG * dTr(i,j,k)
          !=== Computing Çí_k+1 ===
          !---  dV = r_new  = e - c3*Ae ---
            dVr(i,j,k) = dTr(i,j,k) - w_CG*w_dTr(i,j,k)
          !=== Computing b_CG ===
          !---  b_CG_n =  c1 = (r0*,r_new)  ---
            b_CG_n = b_CG_n + (dYsr(i,j,k) * dVr(i,j,k))
          End Do; End Do; End Do

          Call mp_allsumr8( b_CG_n )

          !---  b_CG = c1 / (c2*c3)--- 
          b_CG = ( b_CG_n ) / ( a_CG_d * w_CG )

          !=== Computing P_k+1 ===
          !--- (dTr = p_new) p_new = r_new + beta*(p - c3*Ap) ---
          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1,Nx-1, 1
!NECs
            dUr(i,j,k) =   dVr(i,j,k) + b_CG * ( dUr(i,j,k) - w_CG*w_dUr(i,j,k) )
!NECe
          End Do; End Do; End Do

          !--- a_CG_n = c1 = (r0*, r0)  inner product 
          !--- a_CG_d = c2 = (r0*, Ap)  inner product 
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  f_sum  = dYsr(i,j,k)*dVr(i,j,k)
          !  t_n    = a_CG_n + ( f_sum + r_n )
          !  r_n    = ( f_sum + r_n ) - ( t_n - a_CG_n )
          !  a_CG_n = t_n
          !  f_sum  = dYsr(i,j,k) * ( +  ( dUr(i+1,j  ,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i-1,j  ,k  )) * (c1_dx*c1_dx)  &
          !                           +  ( dUr(i  ,j+1,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j-1,k  )) * (c1_dy*c1_dy)  &
          !                           +  ( dUr(i  ,j  ,k+1) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j  ,k-1)) * (c1_dz*c1_dz)  )
          !  t_d    = a_CG_d + ( f_sum + r_d )
          !  r_d    = ( f_sum + r_d ) - ( t_d - a_CG_d )
          !  a_CG_d = t_d
          !End Do; End Do; End Do
          !Call mp_allsumr( a_CG_n )
          !Call mp_allsumr( a_CG_d )

          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  a_CG_n = a_CG_n + ( dYsr(i,j,k) * dVr(i,j,k) )
          !  a_CG_d = a_CG_d + dYsr(i,j,k)    &
          !         * ( +  ( dUr(i+1,j  ,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
          !             +  ( dUr(i  ,j+1,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
          !             +  ( dUr(i  ,j  ,k+1) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j  ,k-1)) * (c1_dz*c1_dz)  ) 
          !End Do; End Do; End Do
          !Call mp_allsumr8( a_CG_n )
          !Call mp_allsumr8( a_CG_d )

          !--- w_CG_n = ( e, Ae) ---
          !--- w_CG_d = (Ae, Ae) ---
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  f_sum = dTr(i,j,k)                                                                         &
          !          * ( +  ( dTr(i+1,j  ,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
          !              +  ( dTr(i  ,j+1,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
          !              +  ( dTr(i  ,j  ,k+1) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j  ,k-1)) * (c1_dz*c1_dz)    )
          !  t_n    = w_CG_n + ( f_sum + r_n )
          !  r_n    = ( f_sum + r_n ) - ( t_n - w_CG_n )
          !  w_CG_n = t_n
          !  f_sum =  ( + ( dTr(i+1,j  ,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
          !             + ( dTr(i  ,j+1,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
          !             + ( dTr(i  ,j  ,k+1) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j  ,k-1)) * (c1_dz*c1_dz) )  &
          !          *( + ( dTr(i+1,j  ,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i-1,j  ,k  )) * (c1_dx*c1_dx)    &
          !             + ( dTr(i  ,j+1,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j-1,k  )) * (c1_dy*c1_dy)    &
          !             + ( dTr(i  ,j  ,k+1) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j  ,k-1)) * (c1_dz*c1_dz) )   
          !  t_d    = w_CG_d + ( f_sum + r_d )
          !  r_d    = ( f_sum + r_d ) - ( t_d - w_CG_d )
          !  w_CG_d = t_d
          !End Do; End Do; End Do
          !Call mp_allsumr( w_CG_n )
          !Call mp_allsumr( w_CG_d )
          !Call mp_allsumr8( w_CG_n )
          !Call mp_allsumr8( w_CG_d )

          !Call mp_allsumr( w_CG_n )
          !Call mp_allsumr( w_CG_d )
          !Call mp_allsumr8( w_CG_n )
          !Call mp_allsumr8( w_CG_d )

          !=== Computing Çí_k+1 ===
          !---  dV = r_new  = e - c3*Ae ---
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  dVr(i,j,k) =   dTr(i,j,k)    &
          !               - w_CG*( + ( dTr(i+1,j  ,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i-1,j  ,k  )) * (c1_dx*c1_dx) &
          !                        + ( dTr(i  ,j+1,k  ) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j-1,k  )) * (c1_dy*c1_dy) &
          !                        + ( dTr(i  ,j  ,k+1) - 2.0d0*dTr(i,j,k)+ dTr(i  ,j  ,k-1)) * (c1_dz*c1_dz) )
          !End Do; End Do; End Do

          !=== Computing b_CG ===
          !---  b_CG_n =  c1 = (r0*,r_new)  ---
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  f_sum  = (dYsr(i,j,k) * dVr(i,j,k))
          !  t_n    = b_CG_n + ( f_sum + r_n )
          !  r_n    = ( f_sum + r_n ) - ( t_n - b_CG_n )
          !  b_CG_n = t_n
          !End Do; End Do; End Do

          !Call mp_allsumr( b_CG_n )
          !Call mp_allsumr8( b_CG_n )



          !=== Computing P_k+1 ===
          !--- (dTr = p_new) p_new = r_new + beta*(p - c3*Ap) ---
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  dTr(i,j,k) =   dVr(i,j,k) + b_CG * ( dUr(i,j,k) - w_CG   &
          !              * ( +  ( dUr(i+1,j  ,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i-1,j  ,k  )) * (c1_dx*c1_dx)   &
          !                  +  ( dUr(i  ,j+1,k  ) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j-1,k  )) * (c1_dy*c1_dy)   &
          !                  +  ( dUr(i  ,j  ,k+1) - 2.0d0*dUr(i,j,k)+ dUr(i  ,j  ,k-1)) * (c1_dz*c1_dz) ) )
          !End Do; End Do; End Do
          !Do k = ks, ke, 1
          !Do j = js, je, 1
          !Do i = 1,Nx-1, 1
          !  dUr(i,j,k) = dTr(i,j,k)
          !End Do; End Do; End Do

        End Do

        !=== Correct Velocity ===
        Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Phi )


  777 Format(50e25.16e3, 200(1x, 50e25.16e3)) 

      Return
      End Subroutine

!============================================================================================================================
!=============================================================================================================================
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!=============================================================================================================================
!=============================================================================================================================

!=============================================================================================================================
!*** [Periodic BC] ***********************************************************************************************************
!=============================================================================================================================
      Subroutine BC_Periodic_xyz(n11,n12,n21,n22,n31,n32,F_BC)

        Use module_mpi
        Use Definition
        Implicit none   
        Integer,intent(in   ) :: n11,n12,n21,n22,n31,n32
        Real*4 ,intent(inout) :: F_BC(n11:n12,n21:n22,n31:n32)
        Integer i, j, k ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !====== ã´äEèåèÇÃê›íË ======
        Do k = ks, ke, 1
        Do j = js, je, 1
          !Do i = 0, nbc, 1
          !  F1(Nx+i,j,k) =  F1(i+1,j,k)
          !End Do
          !Do i = -nbc, 0, 1
          !  F1(i,j,k) =  F1(Nx-1+i,j,k)
          !End Do

          !--- Rh ---
          F_BC(Nx+0,j,k) =  F_BC(0+1,j,k)
          F_BC(Nx+1,j,k) =  F_BC(1+1,j,k)
          F_BC(Nx+2,j,k) =  F_BC(2+1,j,k)
          F_BC(Nx+3,j,k) =  F_BC(3+1,j,k)
          F_BC(Nx+4,j,k) =  F_BC(4+1,j,k)

          F_BC( 0,j,k)   =  F_BC(Nx-1-0,j,k)
          F_BC(-1,j,k)   =  F_BC(Nx-1-1,j,k)
          F_BC(-2,j,k)   =  F_BC(Nx-1-2,j,k)
          F_BC(-3,j,k)   =  F_BC(Nx-1-3,j,k)
          F_BC(-4,j,k)   =  F_BC(Nx-1-4,j,k)

        End do; End do

        Call BC_Periodic_y_dir_MPI_R4( n11,n12,n21,n22,n31,n32, F_BC )
        Call BC_Periodic_z_dir_MPI_R4( n11,n12,n21,n22,n31,n32, F_BC )

        Call mp_send_recv_pre( F_BC , 4, n11,n12,n21,n22,n31,n32 )

        Return
      End Subroutine

!=============================================================================================================================
!========= [ Subroutine write_xy ] ===========================================================================================
      Subroutine write_xy_P

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
          Rh_xy (i,j) = 0.0d0
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
           W_xy (i,j) = 0.0d0
          Te_xy (i,j) = 0.0d0
           P_xy (i,j) = 0.0d0
        divU_xy (i,j) = 0.0d0
          Ys_xy (i,j) = 0.0d0
        Enst_xy (i,j) = 0.0d0

        End Do; End Do

        Do k = ks, ke, 1
          If(k==(Nz-1)/2) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
              Rh_xy (i,j) =   U (i,j,k)
               U_xy (i,j) =   V (i,j,k)
               V_xy (i,j) =   W (i,j,k)
               W_xy (i,j) =   Ud(i,j,k)
              Te_xy (i,j) =   Vd(i,j,k)
               P_xy (i,j) =   Wd(i,j,k)
            !divU_xy (i,j) =   Us(i,j,k)  !div U 
            !  Ys_xy (i,j) =   Vs(i,j,k)  !div Ud
            !Enst_xy (i,j) =   Ws(i,j,k)  !div Us
            !divU_xy (i,j) = dUdx(i,j,k)  !div U 
            !  Ys_xy (i,j) = dVdy(i,j,k)  !div Ud
            !Enst_xy (i,j) = dWdz(i,j,k)  !div Us
            divU_xy (i,j) = dUdx (i,j,k)+dVdy (i,j,k)+dWdz (i,j,k)  !div U 
              Ys_xy (i,j) = dUddx(i,j,k)+dVddy(i,j,k)+dWddz(i,j,k)  !div Ud
            Enst_xy (i,j) = dUsdx(i,j,k)+dVsdy(i,j,k)+dWsdz(i,j,k)  !div Us
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Rh_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    W_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Te_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    P_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, divU_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ys_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, Enst_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then

          fname_geo = 'output/HIT_Poisson_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_Poisson_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'


          i_xyz_min =  1   ;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  Nz-1;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)


          fname_geo = 'HIT_Poisson_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'


            !--- Inst Data ---
            desc = 'UP'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'VP'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'WP'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ud'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Vd'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Wd'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'divUt'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'divUd'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'divUs'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) (((  Real(Enst_xy(i,j))  , i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


          Close(50)

        End If

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!=============================================================================================================================
!========= [ Subroutine write_xy ] ===========================================================================================
      Subroutine write_xy_force

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz
        Real*4 Sxx_i,Sxy_i,Sxz_i
        Real*4 Syx_i,Syy_i,Syz_i
        Real*4 Szx_i,Szy_i,Szz_i
        Real*4 Epsi_i

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
          Rh_xy (i,j) = 0.0d0
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
           W_xy (i,j) = 0.0d0
          Te_xy (i,j) = 0.0d0
           P_xy (i,j) = 0.0d0
        divU_xy (i,j) = 0.0d0
          Ys_xy (i,j) = 0.0d0
        Enst_xy (i,j) = 0.0d0

        End Do; End Do

        Do k = ks, ke, 1
          If(k==(Nz-1)/2) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
             Sxx_i = 0.5d0*( dUdx(i,j,k) + dUdx(i,j,k) )
             Sxy_i = 0.5d0*( dUdy(i,j,k) + dVdx(i,j,k) )
             Sxz_i = 0.5d0*( dUdz(i,j,k) + dWdx(i,j,k) )
             Syx_i = 0.5d0*( dVdx(i,j,k) + dUdy(i,j,k) )
             Syy_i = 0.5d0*( dVdy(i,j,k) + dVdy(i,j,k) )
             Syz_i = 0.5d0*( dVdz(i,j,k) + dWdy(i,j,k) )
             Szx_i = 0.5d0*( dWdx(i,j,k) + dUdz(i,j,k) )
             Szy_i = 0.5d0*( dWdy(i,j,k) + dVdz(i,j,k) )
             Szz_i = 0.5d0*( dWdz(i,j,k) + dWdz(i,j,k) )
             Epsi_i = (   T_xx(i,j,k)*Sxx_i + T_xy(i,j,k)*Sxy_i + T_xz(i,j,k)*Sxz_i        &
                        + T_yx(i,j,k)*Syx_i + T_yy(i,j,k)*Syy_i + T_yz(i,j,k)*Syz_i        &
                        + T_zx(i,j,k)*Szx_i + T_zy(i,j,k)*Szy_i + T_zz(i,j,k)*Szz_i )/ Re

              Rh_xy (i,j) =   dUddx(i,j,k)
               U_xy (i,j) =   dVddy(i,j,k)
               V_xy (i,j) =   dWddz(i,j,k)
               W_xy (i,j) =   dUsdx(i,j,k)
              Te_xy (i,j) =   U(i,j,k)*dPdx(i,j,k) + V(i,j,k)*dPdy(i,j,k) + W(i,j,k)*dPdz(i,j,k)
               P_xy (i,j) =   P(i,j,k)*(dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k))
            divU_xy (i,j) =   Epsi_i  !div U 
            !  Ys_xy (i,j) =   Vs(i,j,k)  !div Ud
            !Enst_xy (i,j) =   Ws(i,j,k)  !div Us
            !divU_xy (i,j) = dUdx(i,j,k)  !div U 
            !  Ys_xy (i,j) = dVdy(i,j,k)  !div Ud
            !Enst_xy (i,j) = dWdz(i,j,k)  !div Us
            !divU_xy (i,j) = dUdx (i,j,k)+dVdy (i,j,k)+dWdz (i,j,k)  !div U 
            !  Ys_xy (i,j) = dUddx(i,j,k)+dVddy(i,j,k)+dWddz(i,j,k)  !div Ud
            !Enst_xy (i,j) = dUsdx(i,j,k)+dVsdy(i,j,k)+dWsdz(i,j,k)  !div Us
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Rh_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    W_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Te_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    P_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, divU_xy(1:Nx-1,1:Ny-1) )
        !Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ys_xy(1:Nx-1,1:Ny-1) )
        !Call mp_sum_r22( 1, Nx-1, 1, Ny-1, Enst_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then

          fname_geo = 'output/HIT_force_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_force_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'


          i_xyz_min =  1   ;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  Nz-1;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)

          fname_geo = 'HIT_force_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'


            !--- Inst Data ---
            desc = 'f1'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'f2'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'f3'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'fe'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'PG_V'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'PD'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'dissipation'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

          Close(50)

        End If

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!============================================================================================================================
!========= [ Subroutine write_yz ] ==========================================================================================
!============================================================================================================================
      Subroutine write_yz_P

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1


        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Do k = 1, Nz-1, 1
        Do j = 1, Ny-1, 1
            Rh_yz (j,k) = 0.0d0
             U_yz (j,k) = 0.0d0
             V_yz (j,k) = 0.0d0
             W_yz (j,k) = 0.0d0
            Te_yz (j,k) = 0.0d0
             P_yz (j,k) = 0.0d0
          divU_yz (j,k) = 0.0d0
            Ys_yz (j,k) = 0.0d0
          Enst_yz (j,k) = 0.0d0
        End Do; End Do


        i=Nx-1
        Do k = ks, ke, 1
        Do j = js, je, 1
              Rh_yz (i,j) =   U (i,j,k)
               U_yz (i,j) =   V (i,j,k)
               V_yz (i,j) =   W (i,j,k)
               W_yz (i,j) =   Ud(i,j,k)
              Te_yz (i,j) =   Vd(i,j,k)
               P_yz (i,j) =   Wd(i,j,k)
            !divU_yz (i,j) =   Us(i,j,k)  !div U 
            !  Ys_yz (i,j) =   Vs(i,j,k)  !div Ud
            !Enst_yz (i,j) =   Ws(i,j,k)  !div Us
            !divU_xy (i,j) = dUdx(i,j,k)  !div U 
            !  Ys_xy (i,j) = dVdy(i,j,k)  !div Ud
            !Enst_xy (i,j) = dWdz(i,j,k)  !div Us
            divU_xy (i,j) = dUdx (i,j,k)+dVdy (i,j,k)+dWdz (i,j,k)  !div U 
              Ys_xy (i,j) = dUddx(i,j,k)+dVddy(i,j,k)+dWddz(i,j,k)  !div Ud
            Enst_xy (i,j) = dUsdx(i,j,k)+dVsdy(i,j,k)+dWsdz(i,j,k)  !div Us

        End Do; End Do


        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Rh_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    U_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    V_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    W_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Te_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    P_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1, divU_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Ys_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1, Enst_yz (1:Ny-1,1:Nz-1) )


        If(my_rank==root) then

          fname_geo = 'output/HIT_Poisson_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_Poisson_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'


          i_xyz_min =  Nx-1;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  1   ;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)




          fname_geo = 'HIT_Poisson_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'



            !--- Inst Data ---
            desc = 'UP'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'VP'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'WP'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)



            !--- Inst Data ---
            desc = 'Ud'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Vd'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Wd'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)




            !--- Inst Data ---
            desc = 'divUt'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'divUd'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)



            !--- Inst Data ---
            desc = 'divUs'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) (((  Real(Enst_yz(j,k))  , i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

          Close(50)

        End If


   ! 777 Format(E13.6, 200(1x, E13.6)) 
    777 Format(50e25.16e3, 200(1x, 50e25.16e3)) 

        Return
      End Subroutine


!=============================================================================================================================
!=============================================================================================================================
      Subroutine Calc_Ken

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Real*4 Sxx_i, Sxy_i, Sxz_i
        Real*4 Syx_i, Syy_i, Syz_i
        Real*4 Szx_i, Szy_i, Szz_i
        Real*4 Epsi_i
        Real*4 Rhosum,Urms_Favre
        Real*4, dimension(:,:,:),allocatable :: dT_xxdx,dT_xydy,dT_xzdz
        Real*4, dimension(:,:,:),allocatable :: dT_yxdx,dT_yydy,dT_yzdz
        Real*4, dimension(:,:,:),allocatable :: dT_zxdx,dT_zydy,dT_zzdz

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        If( Ityp_Force == 0 ) then
          B2_Ceff = 0.0d0
        End If

        If( Ityp_Force == 1 ) then
          !--- For output data ---
          Call Calc_Var
          Call BC_Periodic
          Call Calc_Derivative
          Call Calc_Tau
          !-----------------------
          allocate(     dT_xxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_xydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_xzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_yxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_yydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_yzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_zxdx (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_zydy (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )
          allocate(     dT_zzdz (-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)      )

          Call Sub_Derivative_MPI( 1  ,   T_xx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_xxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 2  ,   T_xy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_xydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 3  ,   T_xz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_xzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 1  ,   T_yx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_yxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 2  ,   T_yy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_yydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 3  ,   T_yz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_yzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 1  ,   T_zx  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_zxdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 2  ,   T_zy  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_zydy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
          Call Sub_Derivative_MPI( 3  ,   T_zz  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                         dT_zzdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

          !===ïΩãœ===
           rkensumT_xyz    = 0.0d0
           EpsisumT_xyz    = 0.0d0
          UdPdxsumT_xyz    = 0.0d0

          Do j = js, je
            Do k = ks, ke
            Do i = 1, Nx-1

               Epsi_i =   ( + U(i,j,k)*( + dT_xxdx(i,j,k) + dT_xydy(i,j,k) + dT_xzdz(i,j,k)  )  &
                            + V(i,j,k)*( + dT_yxdx(i,j,k) + dT_yydy(i,j,k) + dT_yzdz(i,j,k)  )  &
                            + W(i,j,k)*( + dT_zxdx(i,j,k) + dT_zydy(i,j,k) + dT_zzdz(i,j,k)  )  ) / Re 

              rkensumT_xyz = rkensumT_xyz + (  +  0.5d0*Rh(i,j,k)*U(i,j,k)*U(i,j,k)                      &
                                               +  0.5d0*Rh(i,j,k)*V(i,j,k)*V(i,j,k)                      &
                                               +  0.5d0*Rh(i,j,k)*W(i,j,k)*W(i,j,k)   ) * c1_NxNyNz

              EpsisumT_xyz =  EpsisumT_xyz + (   Epsi_i   ) * c1_NxNyNz
             UdPdxsumT_xyz = UdPdxsumT_xyz + ( -U(i,j,k)*dPdx(i,j,k) -V(i,j,k)*dPdy(i,j,k) -W(i,j,k)*dPdz(i,j,k)) * c1_NxNyNz

            End Do; End Do
          End Do


          !---/* MPI */---
          !---------------------------------------------------
          Call mp_allsumr(  EpsisumT_xyz  )
          Call mp_allsumr(  rkensumT_xyz  )
          Call mp_allsumr(  UdPdxsumT_xyz )

          Urms_Favre = rkensumT_xyz/Rhosum

          forcing_term = (-EpsisumT_xyz - UdPdxsumT_xyz)/(2.0d0*rkensumT_xyz)*2.0d0*rkensumT_xyz

          !If(my_rank==root) then
          !  writqe(70,*) 'B2',B2_Ceff
          !End If

        End If

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
!========= [ Subroutine Itypwrite ] =========================================================================================
      Subroutine Itypwrite

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k

        Integer ns, is, ie

        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        If(Ityp_data .eq. 1) then 
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Nt', Form = 'unFormatted', Status = 'unknown')
            Write(40) Nt, time
            Close(40)
          End If

          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Rh', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Rh(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Ur', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Ur(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Vr', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Vr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )


          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Wr', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Wr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Tr', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Tr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = 'temp/Inst_Data_Ysr', Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Ysr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

        End If





        !============================================================================
        If(Ityp_data .eq. 2) then 
          Open(40, File = 'temp/Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48),&
               Form = 'unFormatted', Status = 'unknown')
          Write(40) Nt, time
          Write(40) j_sta,j_end,k_sta,k_end
          Write(40) Rh
          Write(40) Ur
          Write(40) Vr
          Write(40) Wr
          Write(40) Tr
          Write(40) Ysr
          Close(40)
        End If




  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
!========= [ Subroutine Itypwrite ] =========================================================================================
      Subroutine Itypwrite_STI

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        Integer ns, is, ie
        Integer fn0, fn1, fn2, fn3, fn4, fn5, fn6

        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        fn6 =   nt       / 1000000
        fn5 = ( nt - fn6 * 1000000     ) / 100000
        fn4 = ( nt - fn6 * 1000000 - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn6 * 1000000 - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn6 * 1000000 - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn6 * 1000000 - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn6 * 1000000 - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1

        !============================================================================
        !If(Ityp_data .eq. 2) then 
          !Open(40, File = 'temp/Inst_Data_myrank'//char(fn3_rank+48)//char(fn2_rank+48)//char(fn1_rank+48)//char(fn0_rank+48),&
          !     Form = 'unFormatted', Status = 'unknown')
          Open(40, File = 'output/Inst_Data_myrank'//char(fn3_rank+48)              &
                                                   //char(fn2_rank+48)              &
                                                   //char(fn1_rank+48)              &
                                                   //char(fn0_rank+48)              &
                                              //'_'//char(fn6     +48)              &
                                                   //char(fn5     +48)              &
                                                   //char(fn4     +48)              &
                                                   //char(fn3     +48)              &
                                                   //char(fn2     +48)              &
                                                   //char(fn1     +48)              &
                                                   //char(fn0     +48)//'step.txt' ,&
               Status = 'unknown')
          Write(40,*) Nt, time
          Write(40,*) j_sta,j_end,k_sta,k_end
          Do k = ks,  ke
          Do j = js,  je
          Do i =  1,Nx-1
            Write(40,777) Rh(i,j,k)
          End Do ; End Do ; End Do

          Do k = ks,  ke
          Do j = js,  je
          Do i =  1,Nx-1
            Write(40,777) Ur(i,j,k)
          End Do ; End Do ; End Do

          Do k = ks,  ke
          Do j = js,  je
          Do i =  1,Nx-1
            Write(40,777) Vr(i,j,k)
          End Do ; End Do ; End Do

          Do k = ks,  ke
          Do j = js,  je
          Do i =  1,Nx-1
            Write(40,777) Wr(i,j,k)
          End Do ; End Do ; End Do

          Do k = ks,  ke
          Do j = js,  je
          Do i =  1,Nx-1
            Write(40,777) Tr(i,j,k)
          End Do ; End Do ; End Do


          !Write(40,777) Rh
          !Write(40,777) Ur
          !Write(40,777) Vr
          !Write(40,777) Wr
          !Write(40,777) Tr
          !Write(40,777) Ysr
          Close(40)
        !End If




  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!=============================================================================================================================
!========= [ Subroutine write_xyz ] ==========================================================================================
!=============================================================================================================================
      Subroutine write_xyz

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k

        Integer ns, is, ie
        Integer nsb_l

        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end




        !=======================================================================
        If(Ityp_Save .eq. 1) then 

          If( my_rank==root) then
            Call subfname_write_xyz_1
          End If

          If( my_rank==root) then
            Open(40, File = fname_n_write, Form = 'unFormatted', Status = 'unknown')
            Write(40) Nt, time
            Close(40)
          End If

          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_r_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Rh(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_u_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Ur(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_v_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Vr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )


          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_w_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Wr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_t_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Tr(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

          !---------------------
          If( my_rank==root) then
            Open(40, File = fname_s_write, Form = 'unFormatted', Status = 'unknown')
          End If

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Ys(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Do i = -4, Nx+4, 1
              Write(40) ((Data_r4 (i,j,k),j=-4,Ny+4),k=-4,Nz+4)
            End Do
            Close(40)
          End If

          Deallocate (Data_r4 )

        End If




        !===========================================================================
        If(Ityp_Save .eq. 2) then
          If( my_rank==root) then
            Call subfname_write_xyz_2
          End If

          If( my_rank==root) then
            Open(40, File = fname_n_write, Form = 'unFormatted', Status = 'unknown')
            Write(40) Nt, time
            Close(40)
          End If


          !---------------------
          !--- Rh ---
          data_num = n_slab*0
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Rh(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Rh No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_r_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do


          !---------------------
          !--- Ur ---
          data_num = n_slab*1
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Ur(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Ur No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_u_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do


          !---------------------
          !--- Vr ---
          data_num = n_slab*2
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Vr(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Vr No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_v_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do


          !---------------------
          !--- Wr ---
          data_num = n_slab*3
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Wr(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Wr No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_w_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do


          !---------------------
          !--- Tr ---
          data_num = n_slab*4
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Tr(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Tr No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_t_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do


          !---------------------
          !--- Ysr ---
          data_num = n_slab*5
          rank_write = my_rank - data_num

          Do ns = 1, n_div, 1
            is = 1 + (ns-1) * (Nx-1)/n_div
            ie =     (ns  ) * (Nx-1)/n_div
            Allocate ( Data_r4(is:ie,1:Ny-1,1:Nz-1) )
            Data_r4 = 0.0d0

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = is, ie, 1
              Data_r4(i,j,k) = Ysr(i,j,k)
            End Do; End Do; End Do

            Do i = is, ie, 1
              Call mp_allsum_r32( i, i,  1, Ny-1,  1, Nz-1 , Data_r4(i:i,1:Ny-1,1:Nz-1) )
            End Do

            If( 1<=rank_write .and. rank_write<=n_slab ) then
              nsb_l = rank_write
              nsb = (ns-1)*(n_slab) + (nsb_l)
              ib_sta = 1 + (ns-1)*(Nx-1)/(n_div) + (nsb_l-1)*(Nx-1)/(n_div*n_slab)
              ib_end =     (ns-1)*(Nx-1)/(n_div) + (nsb_l  )*(Nx-1)/(n_div*n_slab)
              Write(*,*) 'Ysr No.', nsb, ib_sta, ib_end, -ib_sta+ib_end+1, (Nx-1)/(n_div*n_slab)

              Call subfname_write_xyz_2
              Open(40, File = fname_s_write, Form = 'unFormatted', Status = 'unknown')
              Do i = ib_sta, ib_end, 1
                Write(40) ((Data_r4 (i,j,k),j=1,Ny-1),k=1,Nz-1)
              End Do
              Close(40)
            End If

            Deallocate (Data_r4 )
          End Do

        End If





        !===========================================================================================
        !--- Data on centerline ---
        Allocate(   Rhc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )
        Allocate(   Urc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )
        Allocate(   Vrc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )
        Allocate(   Wrc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )
        Allocate(   Trc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )
        Allocate(  Ysrc  (-4:Nx+4, (Ny-1)/2-4:(Ny-1)/2+4, -4:Nz+4)    )

         Rhc = 0.0d0
         Urc = 0.0d0
         Vrc = 0.0d0
         Wrc = 0.0d0
         Trc = 0.0d0
        Ysrc = 0.0d0


        Do j = js, je,1
          If( (Ny-1)/2-4<=j .and. j<=(Ny-1)/2+4 ) then
            Do k = ks, ke,1
            Do i = -4, Nx+4
               Rhc(i,j,k) =  Rh(i,j,k)
               Urc(i,j,k) =  Ur(i,j,k)
               Vrc(i,j,k) =  Vr(i,j,k)
               Wrc(i,j,k) =  Wr(i,j,k)
               Trc(i,j,k) =  Tr(i,j,k)
              Ysrc(i,j,k) = Ysr(i,j,k)
            End Do; End Do
          End If
        End Do


        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,   Rhc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )
        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,   Urc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )
        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,   Vrc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )
        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,   Wrc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )
        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,   Trc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )
        Call mp_allsum_r32( -4, Nx+4,  (Ny-1)/2-4, (Ny-1)/2+4,  1, Nz-1 ,  Ysrc(-4:Nx+4,  (Ny-1)/2-4:(Ny-1)/2+4,  1:Nz-1) )


        Do j = (Ny-1)/2-4,(Ny-1)/2+4
        Do i =-4,Nx+4
           Rhc (i,j,   0) =  Rhc(i,j,Nz-1)
           Rhc (i,j,  -1) =  Rhc(i,j,Nz-2)
           Rhc (i,j,  -2) =  Rhc(i,j,Nz-3)
           Rhc (i,j,  -3) =  Rhc(i,j,Nz-4)
           Rhc (i,j,  -4) =  Rhc(i,j,Nz-5)
           Rhc (i,j,Nz+0) =  Rhc(i,j,   1)
           Rhc (i,j,Nz+1) =  Rhc(i,j,   2)
           Rhc (i,j,Nz+2) =  Rhc(i,j,   3)
           Rhc (i,j,Nz+3) =  Rhc(i,j,   4)
           Rhc (i,j,Nz+4) =  Rhc(i,j,   5)

           Urc (i,j,   0) =  Urc(i,j,Nz-1)
           Urc (i,j,  -1) =  Urc(i,j,Nz-2)
           Urc (i,j,  -2) =  Urc(i,j,Nz-3)
           Urc (i,j,  -3) =  Urc(i,j,Nz-4)
           Urc (i,j,  -4) =  Urc(i,j,Nz-5)
           Urc (i,j,Nz+0) =  Urc(i,j,   1)
           Urc (i,j,Nz+1) =  Urc(i,j,   2)
           Urc (i,j,Nz+2) =  Urc(i,j,   3)
           Urc (i,j,Nz+3) =  Urc(i,j,   4)
           Urc (i,j,Nz+4) =  Urc(i,j,   5)

           Vrc (i,j,   0) =  Vrc(i,j,Nz-1)
           Vrc (i,j,  -1) =  Vrc(i,j,Nz-2)
           Vrc (i,j,  -2) =  Vrc(i,j,Nz-3)
           Vrc (i,j,  -3) =  Vrc(i,j,Nz-4)
           Vrc (i,j,  -4) =  Vrc(i,j,Nz-5)
           Vrc (i,j,Nz+0) =  Vrc(i,j,   1)
           Vrc (i,j,Nz+1) =  Vrc(i,j,   2)
           Vrc (i,j,Nz+2) =  Vrc(i,j,   3)
           Vrc (i,j,Nz+3) =  Vrc(i,j,   4)
           Vrc (i,j,Nz+4) =  Vrc(i,j,   5)

           Wrc (i,j,   0) =  Wrc(i,j,Nz-1)
           Wrc (i,j,  -1) =  Wrc(i,j,Nz-2)
           Wrc (i,j,  -2) =  Wrc(i,j,Nz-3)
           Wrc (i,j,  -3) =  Wrc(i,j,Nz-4)
           Wrc (i,j,  -4) =  Wrc(i,j,Nz-5)
           Wrc (i,j,Nz+0) =  Wrc(i,j,   1)
           Wrc (i,j,Nz+1) =  Wrc(i,j,   2)
           Wrc (i,j,Nz+2) =  Wrc(i,j,   3)
           Wrc (i,j,Nz+3) =  Wrc(i,j,   4)
           Wrc (i,j,Nz+4) =  Wrc(i,j,   5)

           Trc (i,j,   0) =  Trc(i,j,Nz-1)
           Trc (i,j,  -1) =  Trc(i,j,Nz-2)
           Trc (i,j,  -2) =  Trc(i,j,Nz-3)
           Trc (i,j,  -3) =  Trc(i,j,Nz-4)
           Trc (i,j,  -4) =  Trc(i,j,Nz-5)
           Trc (i,j,Nz+0) =  Trc(i,j,   1)
           Trc (i,j,Nz+1) =  Trc(i,j,   2)
           Trc (i,j,Nz+2) =  Trc(i,j,   3)
           Trc (i,j,Nz+3) =  Trc(i,j,   4)
           Trc (i,j,Nz+4) =  Trc(i,j,   5)

          Ysrc (i,j,   0) =  Ysrc(i,j,Nz-1)
          Ysrc (i,j,  -1) =  Ysrc(i,j,Nz-2)
          Ysrc (i,j,  -2) =  Ysrc(i,j,Nz-3)
          Ysrc (i,j,  -3) =  Ysrc(i,j,Nz-4)
          Ysrc (i,j,  -4) =  Ysrc(i,j,Nz-5)
          Ysrc (i,j,Nz+0) =  Ysrc(i,j,   1)
          Ysrc (i,j,Nz+1) =  Ysrc(i,j,   2)
          Ysrc (i,j,Nz+2) =  Ysrc(i,j,   3)
          Ysrc (i,j,Nz+3) =  Ysrc(i,j,   4)
          Ysrc (i,j,Nz+4) =  Ysrc(i,j,   5)

        End Do; End Do




        If( my_rank == root ) then
          Call subfname_write_xyz_center
          Open(File = fname_write, Form='unFormatted',Status='unknown', Unit=40)
          Write(40)  Rhc
          Write(40)  Urc
          Write(40)  Vrc
          Write(40)  Wrc
          Write(40)  Trc
          Write(40) Ysrc
          Close(40)
        End If

        Deallocate(   Rhc    )
        Deallocate(   Urc    )
        Deallocate(   Vrc    )
        Deallocate(   Wrc    )
        Deallocate(   Trc    )
        Deallocate(  Ysrc    )

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine subfname_write_xyz_center

        Use Definition
        Implicit none
        Integer fn0 , fn1 , fn2 , fn3 , fn4 , fn5
        Integer fsb0, fsb1, fsb2, fsb3

       !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1


       !=== fname ===
        fsb3 = ( nsb      ) / 1000
        fsb2 = ( nsb - fsb3 * 1000     ) / 100
        fsb1 = ( nsb - fsb3 * 1000 - fsb2 * 100     ) / 10
        fsb0 = ( nsb - fsb3 * 1000 - fsb2 * 100 - fsb1 * 10 ) / 1

     fname_write = &
     'output/Inst_Center_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

        Return
      End Subroutine

!=============================================================================================================================
!=============================================================================================================================
      Subroutine subfname_write_xyz_1

        Use module_mpi
        Use Definition
        Implicit none
        Integer fn0, fn1, fn2, fn3, fn4, fn5

       !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1

     fname_n_write = &
     'output/Inst_Nt_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_r_write = &
     'output/Inst_Rh_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_u_write = &
     'output/Inst_Ur_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_v_write = &
     'output/Inst_Vr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_w_write = &
     'output/Inst_Wr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_t_write = &
     'output/Inst_Tr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_s_write = &
     'output/Inst_S_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'


        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine subfname_write_xyz_2

        Use Definition
        Implicit none
        Integer fn0 , fn1 , fn2 , fn3 , fn4 , fn5
        Integer fsb0, fsb1, fsb2, fsb3

       !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1


       !=== fname ===
        fsb3 = ( nsb      ) / 1000
        fsb2 = ( nsb - fsb3 * 1000     ) / 100
        fsb1 = ( nsb - fsb3 * 1000 - fsb2 * 100     ) / 10
        fsb0 = ( nsb - fsb3 * 1000 - fsb2 * 100 - fsb1 * 10 ) / 1

     fname_n_write = &
     'output/Inst_Nt_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'

     fname_r_write = &
     'output/Inst_Rh_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'

     fname_u_write = &
     'output/Inst_Ur_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'

     fname_v_write = &
     'output/Inst_Vr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'

     fname_w_write = &
     'output/Inst_Wr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'

     fname_t_write = &
     'output/Inst_Tr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'

     fname_s_write = &
     'output/Inst_Ysr_xyz_'//char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step_'&
                         //char(fsb3+48)//char(fsb2+48)//char(fsb1+48)//char(fsb0+48)//'.bin'


        Return
      End Subroutine

!=============================================================================================================================
!========= [ Subroutine write_xy ] ===========================================================================================
      Subroutine write_xy

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
          Rh_xy (i,j) = 0.0d0
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
           W_xy (i,j) = 0.0d0
          Te_xy (i,j) = 0.0d0
           P_xy (i,j) = 0.0d0
        divU_xy (i,j) = 0.0d0
          Ys_xy (i,j) = 0.0d0
        Enst_xy (i,j) = 0.0d0
          Ma_xy (i,j) = 0.0d0

        End Do; End Do

        Do k = ks, ke, 1
          If(k==Nz-1) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
              Rh_xy (i,j) =   Rh(i,j,k)
               U_xy (i,j) =    U(i,j,k)
               V_xy (i,j) =    V(i,j,k)
               W_xy (i,j) =    W(i,j,k)
              Te_xy (i,j) =    T(i,j,k)
               P_xy (i,j) =    P(i,j,k)
            divU_xy (i,j) = dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k)
              Ys_xy (i,j) =   Ys(i,j,k)
            Enst_xy (i,j) =   c1_2*( - dVdz(i,j,k) + dWdy(i,j,k) )*( - dVdz(i,j,k) + dWdy(i,j,k) )  &
                            + c1_2*( - dWdx(i,j,k) + dUdz(i,j,k) )*( - dWdx(i,j,k) + dUdz(i,j,k) )  &
                            + c1_2*( - dUdy(i,j,k) + dVdx(i,j,k) )*( - dUdy(i,j,k) + dVdx(i,j,k) )
              Ma_xy (i,j) =  sqrt(U(i,j,k)**2.0d0 + V(i,j,k)**2.0d0 + W(i,j,k)**2.0d0)/Sqrt(gh*T(i,j,k)*cRT_UU)

            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Rh_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    W_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Te_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    P_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, divU_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ys_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, Enst_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ma_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then

          fname_geo = 'output/HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'

          i_xyz_min =  1   ;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  Nz-1;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)




          fname_geo = 'HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'


            !--- Inst Data ---
            desc = 'Rh'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'U'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'V'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'W'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'T'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'P'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'divU'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ys'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'logEnst'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) (((  Real(  dlog10( dmax1( 0.0000001d0, Enst_xy(i,j))) )  &
                          , i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Ma'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) ((( Real(Ma_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


          Close(50)




        End If

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!=============================================================================================================================
!========= [ Subroutine write_xz ] ===========================================================================================
      Subroutine write_xz

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1


        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do k = 1, Nz-1, 1
        Do i = 1, Nx-1, 1
          Rh_xz (i,k) = 0.0d0
           U_xz (i,k) = 0.0d0
           V_xz (i,k) = 0.0d0
           W_xz (i,k) = 0.0d0
          Te_xz (i,k) = 0.0d0
           P_xz (i,k) = 0.0d0
        divU_xz (i,k) = 0.0d0
          Ys_xz (i,k) = 0.0d0
        Enst_xz (i,k) = 0.0d0
          Ma_xz (i,k) = 0.0d0

        End Do; End Do

        Do j = js, je, 1
          If(j==1) then
            Do k = ks, ke, 1
            Do i = 1, Nx-1, 1
              Rh_xz(i,k) =   Rh(i,j,k)
               U_xz(i,k) =    U(i,j,k)
               V_xz(i,k) =    V(i,j,k)
               W_xz(i,k) =    W(i,j,k)
              Te_xz(i,k) =    T(i,j,k)
               P_xz(i,k) =    P(i,j,k)
            divU_xz(i,k) = dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k)
              Ys_xz(i,k) =   Ys(i,j,k)
            Enst_xz(i,k) =   c1_2*( - dVdz(i,j,k) + dWdy(i,j,k) )*( - dVdz(i,j,k) + dWdy(i,j,k) )  &
                           + c1_2*( - dWdx(i,j,k) + dUdz(i,j,k) )*( - dWdx(i,j,k) + dUdz(i,j,k) )  &
                           + c1_2*( - dUdy(i,j,k) + dVdx(i,j,k) )*( - dUdy(i,j,k) + dVdx(i,j,k) )
              Ma_xz(i,k) =  sqrt(U(i,j,k)**2.0d0 + V(i,j,k)**2.0d0 + W(i,j,k)**2.0d0)/Sqrt(gh*T(i,j,k)*cRT_UU)

            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,   Rh_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    U_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    V_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    W_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,   Te_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    P_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1, divU_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,   Ys_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1, Enst_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,   Ma_xz(1:Nx-1,1:Nz-1) )

        If(my_rank==root) then

          fname_geo = 'output/HIT_xz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_xz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'


          i_xyz_min =  1   ;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  1
          k_xyz_min =  1   ;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)




          fname_geo = 'HIT_xz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'


            !--- Inst Data ---
            desc = 'Rh'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'U'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'V'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'W'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'T'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'P'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'divU'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ys'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'logEnst'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) (((  Real(  dlog10( dmax1( 0.0000001d0, Enst_xz(i,k))) )  &
                          , i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ma'
            fname_scl = 'output/HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ma_xz(i,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

          Close(50)




        End If

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!============================================================================================================================
!========= [ Subroutine write_yz ] ==========================================================================================
!============================================================================================================================
      Subroutine write_yz

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1


        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Do k = 1, Nz-1, 1
        Do j = 1, Ny-1, 1
            Rh_yz (j,k) = 0.0d0
             U_yz (j,k) = 0.0d0
             V_yz (j,k) = 0.0d0
             W_yz (j,k) = 0.0d0
            Te_yz (j,k) = 0.0d0
             P_yz (j,k) = 0.0d0
          divU_yz (j,k) = 0.0d0
            Ys_yz (j,k) = 0.0d0
          Enst_yz (j,k) = 0.0d0
        End Do; End Do


        i=Nx
        Do k = ks, ke, 1
        Do j = js, je, 1
            Rh_yz (j,k) =    Rh(i,j,k)
             U_yz (j,k) =     U(i,j,k)
             V_yz (j,k) =     V(i,j,k)
             W_yz (j,k) =     W(i,j,k)
            Te_yz (j,k) =     T(i,j,k)
             P_yz (j,k) =     P(i,j,k)
          divU_yz (j,k) =  dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k)
            Ys_yz (j,k) =    Ys(i,j,k)
          Enst_yz (j,k) =   c1_2*( - dVdz(i,j,k) + dWdy(i,j,k) )*( - dVdz(i,j,k) + dWdy(i,j,k) )  &
                          + c1_2*( - dWdx(i,j,k) + dUdz(i,j,k) )*( - dWdx(i,j,k) + dUdz(i,j,k) )  &
                          + c1_2*( - dUdy(i,j,k) + dVdx(i,j,k) )*( - dUdy(i,j,k) + dVdx(i,j,k) )
            Ma_yz (j,k) =  sqrt(U(i,j,k)**2.0d0 + V(i,j,k)**2.0d0 + W(i,j,k)**2.0d0)/Sqrt(gh*T(i,j,k)*cRT_UU)

        End Do; End Do


        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Rh_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    U_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    V_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    W_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Te_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,    P_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1, divU_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1,   Ys_yz (1:Ny-1,1:Nz-1) )
        Call mp_sum_r22( 1, Ny-1, 1, Nz-1, Enst_yz (1:Ny-1,1:Nz-1) )


        If(my_rank==root) then

          fname_geo = 'output/HIT_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'


          i_xyz_min =  Nx-1;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  1   ;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)




          fname_geo = 'HIT_yz_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'



            !--- Inst Data ---
            desc = 'Rh'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'U'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'V'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)



            !--- Inst Data ---
            desc = 'W'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'T'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'P'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)




            !--- Inst Data ---
            desc = 'divU'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Ys'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)



            !--- Inst Data ---
            desc = 'logEnst'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) (((  Real(  dlog10( dmax1( 0.0000001d0, Enst_yz(j,k))) )  &
                          , i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ma'
            fname_scl = 'output/HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ma_yz(j,k)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_yz_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

          Close(50)

        End If


   ! 777 Format(E13.6, 200(1x, E13.6)) 
    777 Format(50e25.16e3, 200(1x, 50e25.16e3)) 

        Return
      End Subroutine
!================================================================================================================================
!*** [ Statistics subroutine      ] *********************************************************************************************
!================================================================================================================================
!      Subroutine Sum_yz_T
!
!        Use module_mpi
!        Use Definition
!        Implicit none
!        Integer i, j, k
!        Real*4 Sxx_i, Sxy_i, Sxz_i
!        Real*4 Syx_i, Syy_i, Syz_i
!        Real*4 Szx_i, Szy_i, Szz_i
!        Real*4 Omega_x_i, Omega_y_i, Omega_z_i
!        Real*4 Epsi_i
!        Real*4 Epsi_s_i
!        Real*4 Epsi_c_i
!        Real*4 Epsi_ih_i
!
!        !--- MPI ---
!        Integer js,je,ks,ke
!        !-----------
!
!        js=j_sta
!        je=j_end
!        ks=k_sta
!        ke=k_end
!
!
!        Do k = 0, Nzg, 1
!        Do j = 0, Nyg, 1
!             Count_yz(j,k) = 0
!
!             RhsumT    (j,k) = 0.0d0
!             Rh2sumT   (j,k) = 0.0d0
!             Rh3sumT   (j,k) = 0.0d0
!             Rh4sumT   (j,k) = 0.0d0
!
!             TsumT     (j,k) = 0.0d0
!             PsumT     (j,k) = 0.0d0
!             UsumT     (j,k) = 0.0d0
!             VsumT     (j,k) = 0.0d0
!             WsumT     (j,k) = 0.0d0
!             YssumT    (j,k) = 0.0d0
!
!             TTsumT    (j,k) = 0.0d0
!             PPsumT    (j,k) = 0.0d0
!             UUsumT    (j,k) = 0.0d0
!             VVsumT    (j,k) = 0.0d0
!             WWsumT    (j,k) = 0.0d0
!             Ys2sumT   (j,k) = 0.0d0
!
!             T3sumT    (j,k) = 0.0d0
!             P3sumT    (j,k) = 0.0d0
!             U3sumT    (j,k) = 0.0d0
!             V3sumT    (j,k) = 0.0d0
!             W3sumT    (j,k) = 0.0d0
!             Ys3sumT   (j,k) = 0.0d0
!
!             T4sumT    (j,k) = 0.0d0
!             P4sumT    (j,k) = 0.0d0
!             U4sumT    (j,k) = 0.0d0
!             V4sumT    (j,k) = 0.0d0
!             W4sumT    (j,k) = 0.0d0
!             Ys4sumT   (j,k) = 0.0d0
!
!             dUdxsumT  (j,k) = 0.0d0
!             dVdxsumT  (j,k) = 0.0d0
!             dWdxsumT  (j,k) = 0.0d0
!             dUdysumT  (j,k) = 0.0d0
!             dVdysumT  (j,k) = 0.0d0
!             dWdysumT  (j,k) = 0.0d0
!             dUdzsumT  (j,k) = 0.0d0
!             dVdzsumT  (j,k) = 0.0d0
!             dWdzsumT  (j,k) = 0.0d0
!
!             dUdx2sumT (j,k) = 0.0d0
!             dVdx2sumT (j,k) = 0.0d0
!             dWdx2sumT (j,k) = 0.0d0
!             dUdy2sumT (j,k) = 0.0d0
!             dVdy2sumT (j,k) = 0.0d0
!             dWdy2sumT (j,k) = 0.0d0
!             dUdz2sumT (j,k) = 0.0d0
!             dVdz2sumT (j,k) = 0.0d0
!             dWdz2sumT (j,k) = 0.0d0
!
!             dUdx3sumT (j,k) = 0.0d0
!             dVdx3sumT (j,k) = 0.0d0
!             dWdx3sumT (j,k) = 0.0d0
!             dUdy3sumT (j,k) = 0.0d0
!             dVdy3sumT (j,k) = 0.0d0
!             dWdy3sumT (j,k) = 0.0d0
!             dUdz3sumT (j,k) = 0.0d0
!             dVdz3sumT (j,k) = 0.0d0
!             dWdz3sumT (j,k) = 0.0d0
!
!             dUdx4sumT (j,k) = 0.0d0
!             dVdx4sumT (j,k) = 0.0d0
!             dWdx4sumT (j,k) = 0.0d0
!             dUdy4sumT (j,k) = 0.0d0
!             dVdy4sumT (j,k) = 0.0d0
!             dWdy4sumT (j,k) = 0.0d0
!             dUdz4sumT (j,k) = 0.0d0
!             dVdz4sumT (j,k) = 0.0d0
!             dWdz4sumT (j,k) = 0.0d0
!
!             dPdxsumT  (j,k) = 0.0d0
!             dPdysumT  (j,k) = 0.0d0
!             dPdzsumT  (j,k) = 0.0d0
!             dTdxsumT  (j,k) = 0.0d0
!             dTdysumT  (j,k) = 0.0d0
!             dTdzsumT  (j,k) = 0.0d0
!             dYsdxsumT (j,k) = 0.0d0
!             dYsdysumT (j,k) = 0.0d0
!             dYsdzsumT (j,k) = 0.0d0
!
!             dPdx2sumT  (j,k) = 0.0d0
!             dPdy2sumT  (j,k) = 0.0d0
!             dPdz2sumT  (j,k) = 0.0d0
!             dTdx2sumT  (j,k) = 0.0d0
!             dTdy2sumT  (j,k) = 0.0d0
!             dTdz2sumT  (j,k) = 0.0d0
!             dYsdx2sumT (j,k) = 0.0d0
!             dYsdy2sumT (j,k) = 0.0d0
!             dYsdz2sumT (j,k) = 0.0d0
!
!             dPdx3sumT  (j,k) = 0.0d0
!             dPdy3sumT  (j,k) = 0.0d0
!             dPdz3sumT  (j,k) = 0.0d0
!             dTdx3sumT  (j,k) = 0.0d0
!             dTdy3sumT  (j,k) = 0.0d0
!             dTdz3sumT  (j,k) = 0.0d0
!             dYsdx3sumT (j,k) = 0.0d0
!             dYsdy3sumT (j,k) = 0.0d0
!             dYsdz3sumT (j,k) = 0.0d0
!
!             dPdx4sumT  (j,k) = 0.0d0
!             dPdy4sumT  (j,k) = 0.0d0
!             dPdz4sumT  (j,k) = 0.0d0
!             dTdx4sumT  (j,k) = 0.0d0
!             dTdy4sumT  (j,k) = 0.0d0
!             dTdz4sumT  (j,k) = 0.0d0
!             dYsdx4sumT (j,k) = 0.0d0
!             dYsdy4sumT (j,k) = 0.0d0
!             dYsdz4sumT (j,k) = 0.0d0
!
!                 EpsisumT (j,k) = 0.0d0
!              Epsi_s_sumT (j,k) = 0.0d0 ! \mu^(\Omega_i)(\Omega_i)
!              Epsi_c_sumT (j,k) = 0.0d0 ! (4/3)\mu^(S_{kk})^2
!             Epsi_ih_sumT (j,k) = 0.0d0 ! 2\mu[(dU_{i}/dx_{j})(dU_{j}/dx_{i})-S_{kk}S_{kk}]
!
!             UVsumT       (j,k) = 0.0d0
!             VWsumT       (j,k) = 0.0d0
!             WUsumT       (j,k) = 0.0d0
!
!             UYssumT      (j,k) = 0.0d0
!             VYssumT      (j,k) = 0.0d0
!             WYssumT      (j,k) = 0.0d0
!
!             UPsumT       (j,k) = 0.0d0
!             VPsumT       (j,k) = 0.0d0
!             WPsumT       (j,k) = 0.0d0
!
!             dUdxPsumT    (j,k) = 0.0d0
!             dVdyPsumT    (j,k) = 0.0d0
!             dWdzPsumT    (j,k) = 0.0d0
!
!
!                nusumT     (j,k) = 0.0d0
!                musumT     (j,k) = 0.0d0
!               mu2sumT     (j,k) = 0.0d0
!              divUsumT     (j,k) = 0.0d0
!             divU2sumT     (j,k) = 0.0d0
!
!        End Do; End Do
!
!
!        Do j = js, je
!          Do k = ks, ke
!          Do i = 1, Nx-1
!             jg = j_yg(j)
!             kg = k_zg(k)
!
!             Sxx_i = 0.5d0*( dUdx(i,j,k) + dUdx(i,j,k) )
!             Sxy_i = 0.5d0*( dUdy(i,j,k) + dVdx(i,j,k) )
!             Sxz_i = 0.5d0*( dUdz(i,j,k) + dWdx(i,j,k) )
!
!             Syx_i = 0.5d0*( dVdx(i,j,k) + dUdy(i,j,k) )
!             Syy_i = 0.5d0*( dVdy(i,j,k) + dVdy(i,j,k) )
!             Syz_i = 0.5d0*( dVdz(i,j,k) + dWdy(i,j,k) )
!
!             Szx_i = 0.5d0*( dWdx(i,j,k) + dUdz(i,j,k) )
!             Szy_i = 0.5d0*( dWdy(i,j,k) + dVdz(i,j,k) )
!             Szz_i = 0.5d0*( dWdz(i,j,k) + dWdz(i,j,k) )
!
!             Omega_x_i =  - dVdz(i,j,k) + dWdy(i,j,k)
!             Omega_y_i =  - dWdx(i,j,k) + dUdz(i,j,k)
!             Omega_z_i =  - dUdy(i,j,k) + dVdx(i,j,k)
!
!             Epsi_i = (   T_xx(i,j,k)*Sxx_i + T_xy(i,j,k)*Sxy_i + T_xz(i,j,k)*Sxz_i        &
!                        + T_yx(i,j,k)*Syx_i + T_yy(i,j,k)*Syy_i + T_yz(i,j,k)*Syz_i        &
!                        + T_zx(i,j,k)*Szx_i + T_zy(i,j,k)*Szy_i + T_zz(i,j,k)*Szz_i )*c1_Re
!
!             Epsi_s_i =       mu(i,j,k)*( Omega_x_i*Omega_x_i + Omega_y_i*Omega_y_i + Omega_z_i*Omega_z_i )*c1_Re
!             Epsi_c_i =  c4_3*mu(i,j,k)*( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )          &
!                                       *( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )*c1_Re
!
!             Epsi_ih_i = 2.0d0*mu(i,j,k)*(   dUdx(i,j,k)*dUdx(i,j,k) + dVdx(i,j,k)*dUdy(i,j,k) + dWdx(i,j,k)*dUdz(i,j,k)      &
!                                           + dUdy(i,j,k)*dVdx(i,j,k) + dVdy(i,j,k)*dVdy(i,j,k) + dWdy(i,j,k)*dVdz(i,j,k)      &
!                                           + dUdz(i,j,k)*dWdx(i,j,k) + dVdz(i,j,k)*dWdy(i,j,k) + dWdz(i,j,k)*dWdz(i,j,k)      &
!
!                                           - ( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )                                      &
!                                            *( dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k) )  )*c1_Re
!
!
!
!            !===ïΩãœ===
!            Count_yz (jg,kg) = Count_yz (jg,kg) + 1
!
!            RhsumT (jg,kg) = RhsumT (jg,kg) + (   Rh(i,j,k)     ) 
!            Rh2sumT(jg,kg) = Rh2sumT(jg,kg) + (   Rh(i,j,k)*Rh(i,j,k)     )
!            Rh3sumT(jg,kg) = Rh3sumT(jg,kg) + (   Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)     )
!            Rh4sumT(jg,kg) = Rh4sumT(jg,kg) + (   Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)     )
!
!             TsumT(jg,kg) =  TsumT(jg,kg) + (    T(i,j,k)     )
!             PsumT(jg,kg) =  PsumT(jg,kg) + (    P(i,j,k)     )
!             UsumT(jg,kg) =  UsumT(jg,kg) + (    U(i,j,k)     )
!             VsumT(jg,kg) =  VsumT(jg,kg) + (    V(i,j,k)     )
!             WsumT(jg,kg) =  WsumT(jg,kg) + (    W(i,j,k)     )
!            YssumT(jg,kg) = YssumT(jg,kg) + (   Ys(i,j,k)     )
!
!            TTsumT(jg,kg) =  TTsumT(jg,kg) + (    T(i,j,k) *  T(i,j,k)  )
!            PPsumT(jg,kg) =  PPsumT(jg,kg) + (    P(i,j,k) *  P(i,j,k)  )
!            UUsumT(jg,kg) =  UUsumT(jg,kg) + (    U(i,j,k) *  U(i,j,k)  )
!            VVsumT(jg,kg) =  VVsumT(jg,kg) + (    V(i,j,k) *  V(i,j,k)  )
!            WWsumT(jg,kg) =  WWsumT(jg,kg) + (    W(i,j,k) *  W(i,j,k)  )
!           Ys2sumT(jg,kg) = Ys2sumT(jg,kg) + (   Ys(i,j,k) * Ys(i,j,k)  )
!
!            T3sumT(jg,kg) =  T3sumT(jg,kg) + (    T(i,j,k) *  T(i,j,k) *  T(i,j,k) )
!            P3sumT(jg,kg) =  P3sumT(jg,kg) + (    P(i,j,k) *  P(i,j,k) *  P(i,j,k) )
!            U3sumT(jg,kg) =  U3sumT(jg,kg) + (    U(i,j,k) *  U(i,j,k) *  U(i,j,k) )
!            V3sumT(jg,kg) =  V3sumT(jg,kg) + (    V(i,j,k) *  V(i,j,k) *  V(i,j,k) )
!            W3sumT(jg,kg) =  W3sumT(jg,kg) + (    W(i,j,k) *  W(i,j,k) *  W(i,j,k) )
!           Ys3sumT(jg,kg) = Ys3sumT(jg,kg) + (   Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) )
!
!            T4sumT(jg,kg) =  T4sumT(jg,kg) + (    T(i,j,k) *  T(i,j,k) *  T(i,j,k) *  T(i,j,k) )
!            P4sumT(jg,kg) =  P4sumT(jg,kg) + (    P(i,j,k) *  P(i,j,k) *  P(i,j,k) *  P(i,j,k) )
!            U4sumT(jg,kg) =  U4sumT(jg,kg) + (    U(i,j,k) *  U(i,j,k) *  U(i,j,k) *  U(i,j,k) )
!            V4sumT(jg,kg) =  V4sumT(jg,kg) + (    V(i,j,k) *  V(i,j,k) *  V(i,j,k) *  V(i,j,k) )
!            W4sumT(jg,kg) =  W4sumT(jg,kg) + (    W(i,j,k) *  W(i,j,k) *  W(i,j,k) *  W(i,j,k) )
!           Ys4sumT(jg,kg) = Ys4sumT(jg,kg) + (   Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) )
!
!
!          dUdxsumT(jg,kg) =   dUdxsumT(jg,kg) + ( dUdx(i,j,k) )
!          dVdxsumT(jg,kg) =   dVdxsumT(jg,kg) + ( dVdx(i,j,k) )
!          dWdxsumT(jg,kg) =   dWdxsumT(jg,kg) + ( dWdx(i,j,k) )
!          dUdysumT(jg,kg) =   dUdysumT(jg,kg) + ( dUdy(i,j,k) )
!          dVdysumT(jg,kg) =   dVdysumT(jg,kg) + ( dVdy(i,j,k) )
!          dWdysumT(jg,kg) =   dWdysumT(jg,kg) + ( dWdy(i,j,k) )
!          dUdzsumT(jg,kg) =   dUdzsumT(jg,kg) + ( dUdz(i,j,k) )
!          dVdzsumT(jg,kg) =   dVdzsumT(jg,kg) + ( dVdz(i,j,k) )
!          dWdzsumT(jg,kg) =   dWdzsumT(jg,kg) + ( dWdz(i,j,k) )
!
!         dUdx2sumT(jg,kg) =  dUdx2sumT(jg,kg) + ( dUdx(i,j,k)*dUdx(i,j,k) )
!         dVdx2sumT(jg,kg) =  dVdx2sumT(jg,kg) + ( dVdx(i,j,k)*dVdx(i,j,k) )
!         dWdx2sumT(jg,kg) =  dWdx2sumT(jg,kg) + ( dWdx(i,j,k)*dWdx(i,j,k) )
!         dUdy2sumT(jg,kg) =  dUdy2sumT(jg,kg) + ( dUdy(i,j,k)*dUdy(i,j,k) )
!         dVdy2sumT(jg,kg) =  dVdy2sumT(jg,kg) + ( dVdy(i,j,k)*dVdy(i,j,k) )
!         dWdy2sumT(jg,kg) =  dWdy2sumT(jg,kg) + ( dWdy(i,j,k)*dWdy(i,j,k) )
!         dUdz2sumT(jg,kg) =  dUdz2sumT(jg,kg) + ( dUdz(i,j,k)*dUdz(i,j,k) )
!         dVdz2sumT(jg,kg) =  dVdz2sumT(jg,kg) + ( dVdz(i,j,k)*dVdz(i,j,k) )
!         dWdz2sumT(jg,kg) =  dWdz2sumT(jg,kg) + ( dWdz(i,j,k)*dWdz(i,j,k) )
!
!         dUdx3sumT(jg,kg) =  dUdx3sumT(jg,kg) + ( dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k) )
!         dVdx3sumT(jg,kg) =  dVdx3sumT(jg,kg) + ( dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k) )
!         dWdx3sumT(jg,kg) =  dWdx3sumT(jg,kg) + ( dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k) )
!         dUdy3sumT(jg,kg) =  dUdy3sumT(jg,kg) + ( dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k) )
!         dVdy3sumT(jg,kg) =  dVdy3sumT(jg,kg) + ( dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k) )
!         dWdy3sumT(jg,kg) =  dWdy3sumT(jg,kg) + ( dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k) )
!         dUdz3sumT(jg,kg) =  dUdz3sumT(jg,kg) + ( dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k) )
!         dVdz3sumT(jg,kg) =  dVdz3sumT(jg,kg) + ( dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k) )
!         dWdz3sumT(jg,kg) =  dWdz3sumT(jg,kg) + ( dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k) )
!
!         dUdx4sumT(jg,kg) =  dUdx4sumT(jg,kg) + ( dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k) )
!         dVdx4sumT(jg,kg) =  dVdx4sumT(jg,kg) + ( dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k) )
!         dWdx4sumT(jg,kg) =  dWdx4sumT(jg,kg) + ( dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k) )
!         dUdy4sumT(jg,kg) =  dUdy4sumT(jg,kg) + ( dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k) )
!         dVdy4sumT(jg,kg) =  dVdy4sumT(jg,kg) + ( dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k) )
!         dWdy4sumT(jg,kg) =  dWdy4sumT(jg,kg) + ( dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k) )
!         dUdz4sumT(jg,kg) =  dUdz4sumT(jg,kg) + ( dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k) )
!         dVdz4sumT(jg,kg) =  dVdz4sumT(jg,kg) + ( dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k) )
!         dWdz4sumT(jg,kg) =  dWdz4sumT(jg,kg) + ( dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k) )
!
!          dPdxsumT(jg,kg) =   dPdxsumT(jg,kg) + ( dPdx(i,j,k) )
!          dPdysumT(jg,kg) =   dPdysumT(jg,kg) + ( dPdy(i,j,k) )
!          dPdzsumT(jg,kg) =   dPdzsumT(jg,kg) + ( dPdz(i,j,k) )
!          dTdxsumT(jg,kg) =   dTdxsumT(jg,kg) + ( dTdx(i,j,k) )
!          dTdysumT(jg,kg) =   dTdysumT(jg,kg) + ( dTdy(i,j,k) )
!          dTdzsumT(jg,kg) =   dTdzsumT(jg,kg) + ( dTdz(i,j,k) )
!
!         dPdx2sumT(jg,kg) =   dPdx2sumT(jg,kg) + ( dPdx(i,j,k)*dPdx(i,j,k) )
!         dPdy2sumT(jg,kg) =   dPdy2sumT(jg,kg) + ( dPdy(i,j,k)*dPdy(i,j,k) )
!         dPdz2sumT(jg,kg) =   dPdz2sumT(jg,kg) + ( dPdz(i,j,k)*dPdz(i,j,k) )
!         dTdx2sumT(jg,kg) =   dTdx2sumT(jg,kg) + ( dTdx(i,j,k)*dTdx(i,j,k) )
!         dTdy2sumT(jg,kg) =   dTdy2sumT(jg,kg) + ( dTdy(i,j,k)*dTdy(i,j,k) )
!         dTdz2sumT(jg,kg) =   dTdz2sumT(jg,kg) + ( dTdz(i,j,k)*dTdz(i,j,k) )
!
!         dPdx3sumT(jg,kg) =   dPdx3sumT(jg,kg) + ( dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k) )
!         dPdy3sumT(jg,kg) =   dPdy3sumT(jg,kg) + ( dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k) )
!         dPdz3sumT(jg,kg) =   dPdz3sumT(jg,kg) + ( dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k) )
!         dTdx3sumT(jg,kg) =   dTdx3sumT(jg,kg) + ( dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k) )
!         dTdy3sumT(jg,kg) =   dTdy3sumT(jg,kg) + ( dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k) )
!         dTdz3sumT(jg,kg) =   dTdz3sumT(jg,kg) + ( dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k) )
!
!         dPdx4sumT(jg,kg) =   dPdx4sumT(jg,kg) + ( dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k) )
!         dPdy4sumT(jg,kg) =   dPdy4sumT(jg,kg) + ( dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k) )
!         dPdz4sumT(jg,kg) =   dPdz4sumT(jg,kg) + ( dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k) )
!         dTdx4sumT(jg,kg) =   dTdx4sumT(jg,kg) + ( dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k) )
!         dTdy4sumT(jg,kg) =   dTdy4sumT(jg,kg) + ( dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k) )
!         dTdz4sumT(jg,kg) =   dTdz4sumT(jg,kg) + ( dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k) )
!
!          dRhdxsumT(jg,kg) =   dRhdxsumT(jg,kg) + ( dRhdx(i,j,k) )
!          dRhdysumT(jg,kg) =   dRhdysumT(jg,kg) + ( dRhdy(i,j,k) )
!          dRhdzsumT(jg,kg) =   dRhdzsumT(jg,kg) + ( dRhdz(i,j,k) )
!
!         dRhdx2sumT(jg,kg) =   dRhdx2sumT(jg,kg) + ( dRhdx(i,j,k)*dRhdx(i,j,k) )
!         dRhdy2sumT(jg,kg) =   dRhdy2sumT(jg,kg) + ( dRhdy(i,j,k)*dRhdy(i,j,k) )
!         dRhdz2sumT(jg,kg) =   dRhdz2sumT(jg,kg) + ( dRhdz(i,j,k)*dRhdz(i,j,k) )
!
!         dRhdx3sumT(jg,kg) =   dRhdx3sumT(jg,kg) + ( dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k) )
!         dRhdy3sumT(jg,kg) =   dRhdy3sumT(jg,kg) + ( dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k) )
!         dRhdz3sumT(jg,kg) =   dRhdz3sumT(jg,kg) + ( dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k) )
!
!         dRhdx4sumT(jg,kg) =   dRhdx4sumT(jg,kg) + ( dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k) )
!         dRhdy4sumT(jg,kg) =   dRhdy4sumT(jg,kg) + ( dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k) )
!         dRhdz4sumT(jg,kg) =   dRhdz4sumT(jg,kg) + ( dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k) )
!
!         dYsdxsumT(jg,kg) =   dYsdxsumT(jg,kg) + ( dYsdx(i,j,k) )
!         dYsdysumT(jg,kg) =   dYsdysumT(jg,kg) + ( dYsdy(i,j,k) )
!         dYsdzsumT(jg,kg) =   dYsdzsumT(jg,kg) + ( dYsdz(i,j,k) )
!
!        dYsdx2sumT(jg,kg) =   dYsdx2sumT(jg,kg) + ( dYsdx(i,j,k)*dYsdx(i,j,k) )
!        dYsdy2sumT(jg,kg) =   dYsdy2sumT(jg,kg) + ( dYsdy(i,j,k)*dYsdy(i,j,k) )
!        dYsdz2sumT(jg,kg) =   dYsdz2sumT(jg,kg) + ( dYsdz(i,j,k)*dYsdz(i,j,k) )
!
!        dYsdx3sumT(jg,kg) =   dYsdx3sumT(jg,kg) + ( dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k) )
!        dYsdy3sumT(jg,kg) =   dYsdy3sumT(jg,kg) + ( dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k) )
!        dYsdz3sumT(jg,kg) =   dYsdz3sumT(jg,kg) + ( dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k) )
!
!        dYsdx4sumT(jg,kg) =   dYsdx4sumT(jg,kg) + ( dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k) )
!        dYsdy4sumT(jg,kg) =   dYsdy4sumT(jg,kg) + ( dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k) )
!        dYsdz4sumT(jg,kg) =   dYsdz4sumT(jg,kg) + ( dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k) )
!
!
!
!             EpsisumT(jg,kg) =     EpsisumT(jg,kg) + (    Epsi_i  )
!          Epsi_s_sumT(jg,kg) =  Epsi_s_sumT(jg,kg) + (  Epsi_s_i  )
!          Epsi_c_sumT(jg,kg) =  Epsi_c_sumT(jg,kg) + (  Epsi_c_i  )
!         Epsi_ih_sumT(jg,kg) = Epsi_ih_sumT(jg,kg) + ( Epsi_ih_i  )
!
!
!            UVsumT(jg,kg) =    UVsumT(jg,kg) + (     U(i,j,k)* V(i,j,k)     )
!            VWsumT(jg,kg) =    VWsumT(jg,kg) + (     V(i,j,k)* W(i,j,k)     )
!            WUsumT(jg,kg) =    WUsumT(jg,kg) + (     W(i,j,k)* U(i,j,k)     )
!
!           UYssumT(jg,kg) =   UYssumT(jg,kg) + (     U(i,j,k)*Ys(i,j,k)     )
!           VYssumT(jg,kg) =   VYssumT(jg,kg) + (     V(i,j,k)*Ys(i,j,k)     )
!           WYssumT(jg,kg) =   WYssumT(jg,kg) + (     W(i,j,k)*Ys(i,j,k)     )
!
!            UPsumT(jg,kg) =    UPsumT(jg,kg) + (     U(i,j,k)* P(i,j,k)     )
!            VPsumT(jg,kg) =    VPsumT(jg,kg) + (     V(i,j,k)* P(i,j,k)     )
!            WPsumT(jg,kg) =    WPsumT(jg,kg) + (     W(i,j,k)* P(i,j,k)     )
!
!         dUdxPsumT(jg,kg) = dUdxPsumT(jg,kg) + (  dUdx(i,j,k)* P(i,j,k)     )
!         dVdyPsumT(jg,kg) = dVdyPsumT(jg,kg) + (  dVdy(i,j,k)* P(i,j,k)     )
!         dWdzPsumT(jg,kg) = dWdzPsumT(jg,kg) + (  dWdz(i,j,k)* P(i,j,k)     )
!
!            nusumT   (jg,kg) =     nusumT  (jg,kg) +   ( mu(i,j,k)/Rh(i,j,k)                 )
!            musumT   (jg,kg) =     musumT  (jg,kg) +   ( mu(i,j,k)                           )
!           mu2sumT   (jg,kg) =    mu2sumT  (jg,kg) +   ( mu(i,j,k)*mu(i,j,k)                 )
!          divUsumT   (jg,kg) =   divUsumT  (jg,kg) +   ( dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k) )
!         divU2sumT   (jg,kg) =  divU2sumT  (jg,kg) +   ( dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k) )    &
!                                                     * ( dUdx(i,j,k)+dVdy(i,j,k)+dWdz(i,j,k) )
!
!
!          End Do; End Do
!        End Do
!
!
!
!        !---/* MPI */---
!        !---------------------------------------------------
!        Call mp_allsum_i22( 0,Nyg, 0,Nzg, Count_yz(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   RhsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Rh2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Rh3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Rh4sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    TsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    PsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    UsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    VsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    WsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   YssumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   TTsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   PPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   UUsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   VVsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   WWsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Ys2sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   T3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   P3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   U3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   V3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   W3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Ys3sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   T4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   P4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   U4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   V4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   W4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  Ys4sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdzsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdzsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdzsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdz2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdz2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdz2sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdz3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdz3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdz3sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dUdz4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dVdz4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dWdz4sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dPdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dPdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dPdzsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dTdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dTdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dTdzsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdz2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdz2sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdz3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdz3sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dPdz4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dTdz4sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dRhdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dRhdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dRhdzsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdz2sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdz3sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dRhdz4sumT(0:Nyg,0:Nzg) )
!
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dYsdxsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dYsdysumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   dYsdzsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdx2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdy2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdz2sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdx3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdy3sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdz3sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdx4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdy4sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,  dYsdz4sumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,       EpsisumT (0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    Epsi_s_sumT (0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,    Epsi_c_sumT (0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,   Epsi_ih_sumT (0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         UVsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         VWsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         WUsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,        UYssumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,        VYssumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,        WYssumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         UPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         VPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         WPsumT(0:Nyg,0:Nzg) )
!
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,      dUdxPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,      dVdyPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,      dWdzPsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         nusumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,         musumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,        mu2sumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,       divUsumT(0:Nyg,0:Nzg) )
!        Call mp_allsum_r22( 0,Nyg, 0,Nzg,      divU2sumT(0:Nyg,0:Nzg) )
!
!
!        !--- Average ---
!        Do k = 1, Nzg-1
!        Do j = 1, Nyg-1
!
!               RhsumT(j,k) =             RhsumT(j,k) /Dble( Count_yz(j,k) )
!              Rh2sumT(j,k) =            Rh2sumT(j,k) /Dble( Count_yz(j,k) )
!              Rh3sumT(j,k) =            Rh3sumT(j,k) /Dble( Count_yz(j,k) )
!              Rh4sumT(j,k) =            Rh4sumT(j,k) /Dble( Count_yz(j,k) )
!                TsumT(j,k) =              TsumT(j,k) /Dble( Count_yz(j,k) )
!                PsumT(j,k) =              PsumT(j,k) /Dble( Count_yz(j,k) )
!                UsumT(j,k) =              UsumT(j,k) /Dble( Count_yz(j,k) )
!                VsumT(j,k) =              VsumT(j,k) /Dble( Count_yz(j,k) )
!                WsumT(j,k) =              WsumT(j,k) /Dble( Count_yz(j,k) )
!               YssumT(j,k) =             YssumT(j,k) /Dble( Count_yz(j,k) )
!                             
!               TTsumT(j,k) =             TTsumT(j,k) /Dble( Count_yz(j,k) )
!               PPsumT(j,k) =             PPsumT(j,k) /Dble( Count_yz(j,k) )
!               UUsumT(j,k) =             UUsumT(j,k) /Dble( Count_yz(j,k) )
!               VVsumT(j,k) =             VVsumT(j,k) /Dble( Count_yz(j,k) )
!               WWsumT(j,k) =             WWsumT(j,k) /Dble( Count_yz(j,k) )
!              Ys2sumT(j,k) =            Ys2sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!               T3sumT(j,k) =             T3sumT(j,k) /Dble( Count_yz(j,k) )
!               P3sumT(j,k) =             P3sumT(j,k) /Dble( Count_yz(j,k) )
!               U3sumT(j,k) =             U3sumT(j,k) /Dble( Count_yz(j,k) )
!               V3sumT(j,k) =             V3sumT(j,k) /Dble( Count_yz(j,k) )
!               W3sumT(j,k) =             W3sumT(j,k) /Dble( Count_yz(j,k) )
!              Ys3sumT(j,k) =            Ys3sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!               T4sumT(j,k) =             T4sumT(j,k) /Dble( Count_yz(j,k) )
!               P4sumT(j,k) =             P4sumT(j,k) /Dble( Count_yz(j,k) )
!               U4sumT(j,k) =             U4sumT(j,k) /Dble( Count_yz(j,k) )
!               V4sumT(j,k) =             V4sumT(j,k) /Dble( Count_yz(j,k) )
!               W4sumT(j,k) =             W4sumT(j,k) /Dble( Count_yz(j,k) )
!              Ys4sumT(j,k) =            Ys4sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!             dUdxsumT(j,k) =           dUdxsumT(j,k) /Dble( Count_yz(j,k) )
!             dVdxsumT(j,k) =           dVdxsumT(j,k) /Dble( Count_yz(j,k) )
!             dWdxsumT(j,k) =           dWdxsumT(j,k) /Dble( Count_yz(j,k) )
!             dUdysumT(j,k) =           dUdysumT(j,k) /Dble( Count_yz(j,k) )
!             dVdysumT(j,k) =           dVdysumT(j,k) /Dble( Count_yz(j,k) )
!             dWdysumT(j,k) =           dWdysumT(j,k) /Dble( Count_yz(j,k) )
!             dUdzsumT(j,k) =           dUdzsumT(j,k) /Dble( Count_yz(j,k) )
!             dVdzsumT(j,k) =           dVdzsumT(j,k) /Dble( Count_yz(j,k) )
!             dWdzsumT(j,k) =           dWdzsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dUdx2sumT(j,k) =          dUdx2sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdx2sumT(j,k) =          dVdx2sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdx2sumT(j,k) =          dWdx2sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdy2sumT(j,k) =          dUdy2sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdy2sumT(j,k) =          dVdy2sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdy2sumT(j,k) =          dWdy2sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdz2sumT(j,k) =          dUdz2sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdz2sumT(j,k) =          dVdz2sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdz2sumT(j,k) =          dWdz2sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dUdx3sumT(j,k) =          dUdx3sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdx3sumT(j,k) =          dVdx3sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdx3sumT(j,k) =          dWdx3sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdy3sumT(j,k) =          dUdy3sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdy3sumT(j,k) =          dVdy3sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdy3sumT(j,k) =          dWdy3sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdz3sumT(j,k) =          dUdz3sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdz3sumT(j,k) =          dVdz3sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdz3sumT(j,k) =          dWdz3sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dUdx4sumT(j,k) =          dUdx4sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdx4sumT(j,k) =          dVdx4sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdx4sumT(j,k) =          dWdx4sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdy4sumT(j,k) =          dUdy4sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdy4sumT(j,k) =          dVdy4sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdy4sumT(j,k) =          dWdy4sumT(j,k) /Dble( Count_yz(j,k) )
!            dUdz4sumT(j,k) =          dUdz4sumT(j,k) /Dble( Count_yz(j,k) )
!            dVdz4sumT(j,k) =          dVdz4sumT(j,k) /Dble( Count_yz(j,k) )
!            dWdz4sumT(j,k) =          dWdz4sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!             dPdxsumT(j,k) =           dPdxsumT(j,k) /Dble( Count_yz(j,k) )
!             dPdysumT(j,k) =           dPdysumT(j,k) /Dble( Count_yz(j,k) )
!             dPdzsumT(j,k) =           dPdzsumT(j,k) /Dble( Count_yz(j,k) )
!             dTdxsumT(j,k) =           dTdxsumT(j,k) /Dble( Count_yz(j,k) )
!             dTdysumT(j,k) =           dTdysumT(j,k) /Dble( Count_yz(j,k) )
!             dTdzsumT(j,k) =           dTdzsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dPdx2sumT(j,k) =          dPdx2sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdy2sumT(j,k) =          dPdy2sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdz2sumT(j,k) =          dPdz2sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdx2sumT(j,k) =          dTdx2sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdy2sumT(j,k) =          dTdy2sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdz2sumT(j,k) =          dTdz2sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dPdx3sumT(j,k) =          dPdx3sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdy3sumT(j,k) =          dPdy3sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdz3sumT(j,k) =          dPdz3sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdx3sumT(j,k) =          dTdx3sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdy3sumT(j,k) =          dTdy3sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdz3sumT(j,k) =          dTdz3sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dPdx4sumT(j,k) =          dPdx4sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdy4sumT(j,k) =          dPdy4sumT(j,k) /Dble( Count_yz(j,k) )
!            dPdz4sumT(j,k) =          dPdz4sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdx4sumT(j,k) =          dTdx4sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdy4sumT(j,k) =          dTdy4sumT(j,k) /Dble( Count_yz(j,k) )
!            dTdz4sumT(j,k) =          dTdz4sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dRhdxsumT(j,k) =          dRhdxsumT(j,k) /Dble( Count_yz(j,k) )
!            dRhdysumT(j,k) =          dRhdysumT(j,k) /Dble( Count_yz(j,k) )
!            dRhdzsumT(j,k) =          dRhdzsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dRhdx2sumT(j,k) =         dRhdx2sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdy2sumT(j,k) =         dRhdy2sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdz2sumT(j,k) =         dRhdz2sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dRhdx3sumT(j,k) =         dRhdx3sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdy3sumT(j,k) =         dRhdy3sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdz3sumT(j,k) =         dRhdz3sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dRhdx4sumT(j,k) =         dRhdx4sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdy4sumT(j,k) =         dRhdy4sumT(j,k) /Dble( Count_yz(j,k) )
!           dRhdz4sumT(j,k) =         dRhdz4sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!                             
!            dYsdxsumT(j,k) =          dYsdxsumT(j,k) /Dble( Count_yz(j,k) )
!            dYsdysumT(j,k) =          dYsdysumT(j,k) /Dble( Count_yz(j,k) )
!            dYsdzsumT(j,k) =          dYsdzsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dYsdx2sumT(j,k) =         dYsdx2sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdy2sumT(j,k) =         dYsdy2sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdz2sumT(j,k) =         dYsdz2sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dYsdx3sumT(j,k) =         dYsdx3sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdy3sumT(j,k) =         dYsdy3sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdz3sumT(j,k) =         dYsdz3sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!           dYsdx4sumT(j,k) =         dYsdx4sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdy4sumT(j,k) =         dYsdy4sumT(j,k) /Dble( Count_yz(j,k) )
!           dYsdz4sumT(j,k) =         dYsdz4sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!             EpsisumT(j,k) =           EpsisumT(j,k) /Dble( Count_yz(j,k) )
!          Epsi_s_sumT(j,k) =        Epsi_s_sumT(j,k) /Dble( Count_yz(j,k) )
!          Epsi_c_sumT(j,k) =        Epsi_c_sumT(j,k) /Dble( Count_yz(j,k) )
!         Epsi_ih_sumT(j,k) =       Epsi_ih_sumT(j,k) /Dble( Count_yz(j,k) )
!                             
!               UVsumT(j,k) =             UVsumT(j,k) /Dble( Count_yz(j,k) )
!               VWsumT(j,k) =             VWsumT(j,k) /Dble( Count_yz(j,k) )
!               WUsumT(j,k) =             WUsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!              UYssumT(j,k) =            UYssumT(j,k) /Dble( Count_yz(j,k) )
!              VYssumT(j,k) =            VYssumT(j,k) /Dble( Count_yz(j,k) )
!              WYssumT(j,k) =            WYssumT(j,k) /Dble( Count_yz(j,k) )
!                             
!               UPsumT(j,k) =             UPsumT(j,k) /Dble( Count_yz(j,k) )
!               VPsumT(j,k) =             VPsumT(j,k) /Dble( Count_yz(j,k) )
!               WPsumT(j,k) =             WPsumT(j,k) /Dble( Count_yz(j,k) )
!                             
!            dUdxPsumT(j,k) =          dUdxPsumT(j,k) /Dble( Count_yz(j,k) )
!            dVdyPsumT(j,k) =          dVdyPsumT(j,k) /Dble( Count_yz(j,k) )
!            dWdzPsumT(j,k) =          dWdzPsumT(j,k) /Dble( Count_yz(j,k) )
!               nusumT(j,k) =             nusumT(j,k) /Dble( Count_yz(j,k) )
!               musumT(j,k) =             musumT(j,k) /Dble( Count_yz(j,k) )
!              mu2sumT(j,k) =            mu2sumT(j,k) /Dble( Count_yz(j,k) )
!             divUsumT(j,k) =           divUsumT(j,k) /Dble( Count_yz(j,k) )
!            divU2sumT(j,k) =          divU2sumT(j,k) /Dble( Count_yz(j,k) )
!
!        End Do; End Do
!
!        Return
!      End Subroutine
!================================================================================================================================
!*** [ Statistics subroutine      ] *********************************************************************************************
!================================================================================================================================
      Subroutine Cor_Center


        Use module_mpi
        Use Definition
        Implicit none 
        Integer i, j, k
        Integer n, i_np
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        Integer fn0, fn1, fn2, fn3, fn4, fn5

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end



        Do n = 0, Nx, 1
          Cuu_C(n) = 0.0d0
          Cvv_C(n) = 0.0d0
          Cww_C(n) = 0.0d0
          Cpp_C(n) = 0.0d0
          Crr_C(n) = 0.0d0
          Ctt_C(n) = 0.0d0
          Cyy_C(n) = 0.0d0
        End Do


        Do k = ks, ke, 1
        Do j = js, je, 1
          If( j_yg(j) == (Nyg-1)/2 .and. k_zg(k) == (Nzg-1)/2 ) then
            Um_C = 0.0d0 
            Vm_C = 0.0d0 
            Wm_C = 0.0d0 
            Pm_C = 0.0d0 
            Rm_C = 0.0d0 
            Tm_C = 0.0d0 
            Ym_C = 0.0d0 


            Do i = 1,Nx-1, 1
              Um_C = Um_C + U (i,j,k) * c1_Nx
              Vm_C = Vm_C + V (i,j,k) * c1_Nx
              Wm_C = Wm_C + W (i,j,k) * c1_Nx
              Pm_C = Pm_C + P (i,j,k) * c1_Nx
              Rm_C = Rm_C + Rh(i,j,k) * c1_Nx
              Tm_C = Tm_C + T (i,j,k) * c1_Nx
              Ym_C = Ym_C + Ys(i,j,k) * c1_Nx
            End Do


            Do i = 1,Nx-1, 1
              Do n = 0, (Nx-1)/2, 1
                i_np = i + n
                i_np = max(1   , i_np)
                i_np = min(Nx-1, i_np)

                Cuu_C(n) = Cuu_C(n)  +  (U (i,j,k)-Um_C)*(U (i_np,j,k)-Um_C)
                Cvv_C(n) = Cvv_C(n)  +  (V (i,j,k)-Vm_C)*(V (i_np,j,k)-Vm_C)
                Cww_C(n) = Cww_C(n)  +  (W (i,j,k)-Wm_C)*(W (i_np,j,k)-Wm_C)
                Cpp_C(n) = Cpp_C(n)  +  (P (i,j,k)-Pm_C)*(P (i_np,j,k)-Pm_C)
                Crr_C(n) = Crr_C(n)  +  (Rh(i,j,k)-Rm_C)*(Rh(i_np,j,k)-Rm_C)
                Ctt_C(n) = Ctt_C(n)  +  (T (i,j,k)-Tm_C)*(T (i_np,j,k)-Tm_C)
                Cyy_C(n) = Cyy_C(n)  +  (Ys(i,j,k)-Ym_C)*(Ys(i_np,j,k)-Ym_C)
              End Do
            End Do

          End If
        End Do; End Do


        Call mp_allsum_r12( 0,  Nx, Cuu_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cvv_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cww_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cpp_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Crr_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Ctt_C(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cyy_C(0:Nx) )


        Cuu_C0 =  Cuu_C(0)
        Cvv_C0 =  Cvv_C(0)
        Cww_C0 =  Cww_C(0)
        Cpp_C0 =  Cpp_C(0)
        Crr_C0 =  Crr_C(0)
        Ctt_C0 =  Ctt_C(0)
        Cyy_C0 =  Cyy_C(0)


        Do n = 0, Nx, 1
          Cuu_C(n) = Cuu_C(n) / Cuu_C0
          Cvv_C(n) = Cvv_C(n) / Cvv_C0
          Cww_C(n) = Cww_C(n) / Cww_C0
          Cpp_C(n) = Cpp_C(n) / Cpp_C0
          Crr_C(n) = Crr_C(n) / Crr_C0
          Ctt_C(n) = Ctt_C(n) / Ctt_C0
          Cyy_C(n) = Cyy_C(n) / Cyy_C0
        End Do 


        !--- Integral length scale ---
        Luu_C = 0.0d0
        Do n = 0, (Nx-1)/2, 1
          Luu_C = Luu_C + Cuu_C(n)*dx
          If(Cuu_C(n)<0.0d0) Exit
        End Do


  777   Format(50e25.16e3)
        Return
      End Subroutine

!============================================================================================================================
!============================================================================================================================
      Subroutine Cor_xyz

        Use module_mpi
        Use Definition
        Implicit none 
        Integer i, j, k
        Integer n, i_np
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        Integer fn0, fn1, fn2, fn3, fn4, fn5

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end



        Do n = 0, Nx, 1
          Cuu_xyz(n) = 0.0d0
          Cvv_xyz(n) = 0.0d0
          Cww_xyz(n) = 0.0d0
         ! Cpp_xyz(n) = 0.0d0
         ! Crr_xyz(n) = 0.0d0
         ! Ctt_xyz(n) = 0.0d0
         ! Cyy_xyz(n) = 0.0d0
        End Do

        Do k = ks, ke, 1
        Do j = js, je, 1

          Um_C = 0.0d0 
          Vm_C = 0.0d0 
          Wm_C = 0.0d0 
         ! Pm_C = 0.0d0 
         ! Rm_C = 0.0d0 
         ! Tm_C = 0.0d0 
         ! Ym_C = 0.0d0 


          Do i = 1,Nx-1, 1
            Um_C = Um_C + U (i,j,k) * c1_Nx
            Vm_C = Vm_C + V (i,j,k) * c1_Nx
            Wm_C = Wm_C + W (i,j,k) * c1_Nx
           ! Pm_C = Pm_C + P (i,j,k) * c1_Nx
           ! Rm_C = Rm_C + Rh(i,j,k) * c1_Nx
           ! Tm_C = Tm_C + T (i,j,k) * c1_Nx
           ! Ym_C = Ym_C + Ys(i,j,k) * c1_Nx
          End Do


          Do i = 1,Nx-1, 1
            Do n = 0, (Nx-1)/2, 1
              i_np = i + n
              i_np = max(1   , i_np)
              i_np = min(Nx-1, i_np)

              Cuu_xyz(n) = Cuu_xyz(n)  +  (U (i,j,k)-Um_C)*(U (i_np,j,k)-Um_C)
              Cvv_xyz(n) = Cvv_xyz(n)  +  (V (i,j,k)-Vm_C)*(V (i_np,j,k)-Vm_C)
              Cww_xyz(n) = Cww_xyz(n)  +  (W (i,j,k)-Wm_C)*(W (i_np,j,k)-Wm_C)
             ! Cpp_xyz(n) = Cpp_xyz(n)  +  (P (i,j,k)-Pm_C)*(P (i_np,j,k)-Pm_C)
             ! Crr_xyz(n) = Crr_xyz(n)  +  (Rh(i,j,k)-Rm_C)*(Rh(i_np,j,k)-Rm_C)
             ! Ctt_xyz(n) = Ctt_xyz(n)  +  (T (i,j,k)-Tm_C)*(T (i_np,j,k)-Tm_C)
             ! Cyy_xyz(n) = Cyy_xyz(n)  +  (Ys(i,j,k)-Ym_C)*(Ys(i_np,j,k)-Ym_C)
            End Do
          End Do

        End Do; End Do


        Call mp_allsum_r12( 0,  Nx, Cuu_xyz(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cvv_xyz(0:Nx) )
        Call mp_allsum_r12( 0,  Nx, Cww_xyz(0:Nx) )
       ! Call mp_allsum_r12( 0,  Nx, Cpp_xyz(0:Nx) )
       ! Call mp_allsum_r12( 0,  Nx, Crr_xyz(0:Nx) )
       ! Call mp_allsum_r12( 0,  Nx, Ctt_xyz(0:Nx) )
       ! Call mp_allsum_r12( 0,  Nx, Cyy_xyz(0:Nx) )


        Cuu_xyz0 =  Cuu_xyz(0)
        Cvv_xyz0 =  Cvv_xyz(0)
        Cww_xyz0 =  Cww_xyz(0)
       ! Cpp_xyz0 =  Cpp_xyz(0)
       ! Crr_xyz0 =  Crr_xyz(0)
       ! Ctt_xyz0 =  Ctt_xyz(0)
       ! Cyy_xyz0 =  Cyy_xyz(0)


        Do n = 0, Nx, 1
          Cuu_xyz(n) = Cuu_xyz(n) / Cuu_xyz0
          Cvv_xyz(n) = Cvv_xyz(n) / Cvv_xyz0
          Cww_xyz(n) = Cww_xyz(n) / Cww_xyz0
         ! Cpp_xyz(n) = Cpp_xyz(n) / Cpp_xyz0
         ! Crr_xyz(n) = Crr_xyz(n) / Crr_xyz0
         ! Ctt_xyz(n) = Ctt_xyz(n) / Ctt_xyz0
         ! Cyy_xyz(n) = Cyy_xyz(n) / Cyy_xyz0
        End Do 


        !--- Integral length scale ---
        Luu_xyz = 0.0d0
        Do n = 0, (Nx-1)/2, 1
          Luu_xyz = Luu_xyz + Cuu_xyz(n)*dx
          If(Cuu_xyz(n)<0.0d0) Exit
        End Do


  777   Format(50e25.16e3)
        Return
      End Subroutine



!================================================================================================================================
!*** [ Statistics subroutine      ] *********************************************************************************************
!================================================================================================================================
      Subroutine Sum_xyz_T

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k
        Real*4 Sxx_i, Sxy_i, Sxz_i
        Real*4 Syx_i, Syy_i, Syz_i
        Real*4 Szx_i, Szy_i, Szz_i
        Real*4 Epsi_i
        Real*4 f1,f2,f3
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        UsumT_xyz     = 0.0d0
        VsumT_xyz     = 0.0d0
        WsumT_xyz     = 0.0d0
        UUsumT_xyz    = 0.0d0
        VVsumT_xyz    = 0.0d0
        WWsumT_xyz    = 0.0d0

             !===ïΩãœ===

              RhsumT_xyz     = 0.0d0
               TsumT_xyz     = 0.0d0
               PsumT_xyz     = 0.0d0
               UsumT_xyz     = 0.0d0
               VsumT_xyz     = 0.0d0
               WsumT_xyz     = 0.0d0
              YssumT_xyz     = 0.0d0
              musumT_xyz     = 0.0d0
              nusumT_xyz     = 0.0d0
             UVWsumT_xyz     = 0.0d0

              Rh2sumT_xyz    = 0.0d0
               TTsumT_xyz    = 0.0d0
               PPsumT_xyz    = 0.0d0
               UUsumT_xyz    = 0.0d0
               VVsumT_xyz    = 0.0d0
               WWsumT_xyz    = 0.0d0
              Ys2sumT_xyz    = 0.0d0
             mu2sumT_xyz     = 0.0d0
             UVW2sumT_xyz    = 0.0d0


              Rh3sumT_xyz    = 0.0d0
               T3sumT_xyz    = 0.0d0
               P3sumT_xyz    = 0.0d0
               U3sumT_xyz    = 0.0d0
               V3sumT_xyz    = 0.0d0
               W3sumT_xyz    = 0.0d0
              Ys3sumT_xyz    = 0.0d0

              Rh4sumT_xyz    = 0.0d0
               T4sumT_xyz    = 0.0d0
               P4sumT_xyz    = 0.0d0
               U4sumT_xyz    = 0.0d0
               V4sumT_xyz    = 0.0d0
               W4sumT_xyz    = 0.0d0
              Ys4sumT_xyz    = 0.0d0

              RUUsumT_xyz    = 0.0d0
              RVVsumT_xyz    = 0.0d0
              RWWsumT_xyz    = 0.0d0
             CTKEsumT_xyz    = 0.0d0
               PDsumT_xyz    = 0.0d0
            ForcesumT_xyz    = 0.0d0

              dUdxsumT_xyz   = 0.0d0
              dVdxsumT_xyz   = 0.0d0
              dWdxsumT_xyz   = 0.0d0
              dUdysumT_xyz   = 0.0d0
              dVdysumT_xyz   = 0.0d0
              dWdysumT_xyz   = 0.0d0
              dUdzsumT_xyz   = 0.0d0
              dVdzsumT_xyz   = 0.0d0
              dWdzsumT_xyz   = 0.0d0

              dUdx2sumT_xyz  = 0.0d0
              dVdx2sumT_xyz  = 0.0d0
              dWdx2sumT_xyz  = 0.0d0
              dUdy2sumT_xyz  = 0.0d0
              dVdy2sumT_xyz  = 0.0d0
              dWdy2sumT_xyz  = 0.0d0
              dUdz2sumT_xyz  = 0.0d0
              dVdz2sumT_xyz  = 0.0d0
              dWdz2sumT_xyz  = 0.0d0

              dUdx3sumT_xyz  = 0.0d0
              dVdx3sumT_xyz  = 0.0d0
              dWdx3sumT_xyz  = 0.0d0
              dUdy3sumT_xyz  = 0.0d0
              dVdy3sumT_xyz  = 0.0d0
              dWdy3sumT_xyz  = 0.0d0
              dUdz3sumT_xyz  = 0.0d0
              dVdz3sumT_xyz  = 0.0d0
              dWdz3sumT_xyz  = 0.0d0

              dUdx4sumT_xyz  = 0.0d0
              dVdx4sumT_xyz  = 0.0d0
              dWdx4sumT_xyz  = 0.0d0
              dUdy4sumT_xyz  = 0.0d0
              dVdy4sumT_xyz  = 0.0d0
              dWdy4sumT_xyz  = 0.0d0
              dUdz4sumT_xyz  = 0.0d0
              dVdz4sumT_xyz  = 0.0d0
              dWdz4sumT_xyz  = 0.0d0

              dRhdxsumT_xyz  = 0.0d0
              dRhdysumT_xyz  = 0.0d0
              dRhdzsumT_xyz  = 0.0d0
               dPdxsumT_xyz  = 0.0d0
               dPdysumT_xyz  = 0.0d0
               dPdzsumT_xyz  = 0.0d0
               dTdxsumT_xyz  = 0.0d0
               dTdysumT_xyz  = 0.0d0
               dTdzsumT_xyz  = 0.0d0
              dYsdxsumT_xyz  = 0.0d0
              dYsdysumT_xyz  = 0.0d0
              dYsdzsumT_xyz  = 0.0d0

             dRhdx2sumT_xyz  = 0.0d0
             dRhdy2sumT_xyz  = 0.0d0
             dRhdz2sumT_xyz  = 0.0d0
              dPdx2sumT_xyz  = 0.0d0
              dPdy2sumT_xyz  = 0.0d0
              dPdz2sumT_xyz  = 0.0d0
              dTdx2sumT_xyz  = 0.0d0
              dTdy2sumT_xyz  = 0.0d0
              dTdz2sumT_xyz  = 0.0d0
             dYsdx2sumT_xyz  = 0.0d0
             dYsdy2sumT_xyz  = 0.0d0
             dYsdz2sumT_xyz  = 0.0d0

             dRhdx3sumT_xyz  = 0.0d0
             dRhdy3sumT_xyz  = 0.0d0
             dRhdz3sumT_xyz  = 0.0d0
              dPdx3sumT_xyz  = 0.0d0
              dPdy3sumT_xyz  = 0.0d0
              dPdz3sumT_xyz  = 0.0d0
              dTdx3sumT_xyz  = 0.0d0
              dTdy3sumT_xyz  = 0.0d0
              dTdz3sumT_xyz  = 0.0d0
             dYsdx3sumT_xyz  = 0.0d0
             dYsdy3sumT_xyz  = 0.0d0
             dYsdz3sumT_xyz  = 0.0d0

             dRhdx4sumT_xyz  = 0.0d0
             dRhdy4sumT_xyz  = 0.0d0
             dRhdz4sumT_xyz  = 0.0d0
              dPdx4sumT_xyz  = 0.0d0
              dPdy4sumT_xyz  = 0.0d0
              dPdz4sumT_xyz  = 0.0d0
              dTdx4sumT_xyz  = 0.0d0
              dTdy4sumT_xyz  = 0.0d0
              dTdz4sumT_xyz  = 0.0d0
             dYsdx4sumT_xyz  = 0.0d0
             dYsdy4sumT_xyz  = 0.0d0
             dYsdz4sumT_xyz  = 0.0d0


             UdPdxsumT_xyz   =  0.0d0
             VdPdysumT_xyz   =  0.0d0
             WdPdzsumT_xyz   =  0.0d0


             EpsisumT_xyz = 0.0d0


             rkensumT_xyz = 0.0d0



        Do j = js, je
          Do k = ks, ke
          Do i = 1, Nx-1
             Sxx_i = 0.5d0*( dUdx(i,j,k) + dUdx(i,j,k) )
             Sxy_i = 0.5d0*( dUdy(i,j,k) + dVdx(i,j,k) )
             Sxz_i = 0.5d0*( dUdz(i,j,k) + dWdx(i,j,k) )

             Syx_i = 0.5d0*( dVdx(i,j,k) + dUdy(i,j,k) )
             Syy_i = 0.5d0*( dVdy(i,j,k) + dVdy(i,j,k) )
             Syz_i = 0.5d0*( dVdz(i,j,k) + dWdy(i,j,k) )

             Szx_i = 0.5d0*( dWdx(i,j,k) + dUdz(i,j,k) )
             Szy_i = 0.5d0*( dWdy(i,j,k) + dVdz(i,j,k) )
             Szz_i = 0.5d0*( dWdz(i,j,k) + dWdz(i,j,k) )

             Epsi_i = (   T_xx(i,j,k)*Sxx_i + T_xy(i,j,k)*Sxy_i + T_xz(i,j,k)*Sxz_i        &
                        + T_yx(i,j,k)*Syx_i + T_yy(i,j,k)*Syy_i + T_yz(i,j,k)*Syz_i        &
                        + T_zx(i,j,k)*Szx_i + T_zy(i,j,k)*Szy_i + T_zz(i,j,k)*Szz_i )/ Re

             f1 = + c_d*sqrt(Rh(i,j,k))*Ud(i,j,k)   &
                  + c_s*sqrt(Rh(i,j,k))*Us(i,j,k)    
             f2 = + c_d*sqrt(Rh(i,j,k))*Vd(i,j,k)   &
                  + c_s*sqrt(Rh(i,j,k))*Vs(i,j,k)    
             f3 = + c_d*sqrt(Rh(i,j,k))*Wd(i,j,k)   &
                  + c_s*sqrt(Rh(i,j,k))*Ws(i,j,k)    

            !===ïΩãœ===
            RhsumT_xyz  = RhsumT_xyz  + (   Rh(i,j,k)                                   ) * c1_NxNyNz
            Rh2sumT_xyz = Rh2sumT_xyz + (   Rh(i,j,k)*Rh(i,j,k)                         ) * c1_NxNyNz
            Rh3sumT_xyz = Rh3sumT_xyz + (   Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)               ) * c1_NxNyNz
            Rh4sumT_xyz = Rh4sumT_xyz + (   Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)*Rh(i,j,k)     ) * c1_NxNyNz

             TsumT_xyz =   TsumT_xyz + (    T(i,j,k)               ) * c1_NxNyNz
             PsumT_xyz =   PsumT_xyz + (    P(i,j,k)               ) * c1_NxNyNz
             UsumT_xyz =   UsumT_xyz + (    U(i,j,k)               ) * c1_NxNyNz
             VsumT_xyz =   VsumT_xyz + (    V(i,j,k)               ) * c1_NxNyNz
             WsumT_xyz =   WsumT_xyz + (    W(i,j,k)               ) * c1_NxNyNz
            YssumT_xyz =  YssumT_xyz + (   Ys(i,j,k)               ) * c1_NxNyNz
            musumT_xyz =  musumT_xyz + (   mu(i,j,k)               ) * c1_NxNyNz
            nusumT_xyz =  nusumT_xyz + (   mu(i,j,k)/Rh(i,j,k)     ) * c1_NxNyNz
           UVWsumT_xyz = UVWsumT_xyz + (   U(i,j,k)                ) * c1_NxNyNz
           UVWsumT_xyz = UVWsumT_xyz + (   V(i,j,k)                ) * c1_NxNyNz
           UVWsumT_xyz = UVWsumT_xyz + (   W(i,j,k)                ) * c1_NxNyNz

            TTsumT_xyz =   TTsumT_xyz + (    T(i,j,k) *  T(i,j,k)  ) * c1_NxNyNz
            PPsumT_xyz =   PPsumT_xyz + (    P(i,j,k) *  P(i,j,k)  ) * c1_NxNyNz
            UUsumT_xyz =   UUsumT_xyz + (    U(i,j,k) *  U(i,j,k)  ) * c1_NxNyNz
            VVsumT_xyz =   VVsumT_xyz + (    V(i,j,k) *  V(i,j,k)  ) * c1_NxNyNz
            WWsumT_xyz =   WWsumT_xyz + (    W(i,j,k) *  W(i,j,k)  ) * c1_NxNyNz
           Ys2sumT_xyz =  Ys2sumT_xyz + (   Ys(i,j,k) * Ys(i,j,k)  ) * c1_NxNyNz
           mu2sumT_xyz =  mu2sumT_xyz + (   mu(i,j,k) * mu(i,j,k)  ) * c1_NxNyNz
          UVW2sumT_xyz = UVW2sumT_xyz + (    U(i,j,k) *  U(i,j,k)  ) * c1_NxNyNz
          UVW2sumT_xyz = UVW2sumT_xyz + (    V(i,j,k) *  V(i,j,k)  ) * c1_NxNyNz
          UVW2sumT_xyz = UVW2sumT_xyz + (    W(i,j,k) *  W(i,j,k)  ) * c1_NxNyNz

            T3sumT_xyz =  T3sumT_xyz + (    T(i,j,k) *  T(i,j,k) *  T(i,j,k) ) * c1_NxNyNz
            P3sumT_xyz =  P3sumT_xyz + (    P(i,j,k) *  P(i,j,k) *  P(i,j,k) ) * c1_NxNyNz
            U3sumT_xyz =  U3sumT_xyz + (    U(i,j,k) *  U(i,j,k) *  U(i,j,k) ) * c1_NxNyNz
            V3sumT_xyz =  V3sumT_xyz + (    V(i,j,k) *  V(i,j,k) *  V(i,j,k) ) * c1_NxNyNz
            W3sumT_xyz =  W3sumT_xyz + (    W(i,j,k) *  W(i,j,k) *  W(i,j,k) ) * c1_NxNyNz
           Ys3sumT_xyz = Ys3sumT_xyz + (   Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) ) * c1_NxNyNz

            T4sumT_xyz =  T4sumT_xyz + (    T(i,j,k) *  T(i,j,k) *  T(i,j,k) *  T(i,j,k) ) * c1_NxNyNz
            P4sumT_xyz =  P4sumT_xyz + (    P(i,j,k) *  P(i,j,k) *  P(i,j,k) *  P(i,j,k) ) * c1_NxNyNz
            U4sumT_xyz =  U4sumT_xyz + (    U(i,j,k) *  U(i,j,k) *  U(i,j,k) *  U(i,j,k) ) * c1_NxNyNz
            V4sumT_xyz =  V4sumT_xyz + (    V(i,j,k) *  V(i,j,k) *  V(i,j,k) *  V(i,j,k) ) * c1_NxNyNz
            W4sumT_xyz =  W4sumT_xyz + (    W(i,j,k) *  W(i,j,k) *  W(i,j,k) *  W(i,j,k) ) * c1_NxNyNz
           Ys4sumT_xyz = Ys4sumT_xyz + (   Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) * Ys(i,j,k) ) * c1_NxNyNz

           RUUsumT_xyz = RUUsumT_xyz  + (   Rh(i,j,k) * U(i,j,k) * U(i,j,k)  ) * c1_NxNyNz
           RVVsumT_xyz = RVVsumT_xyz  + (   Rh(i,j,k) * V(i,j,k) * V(i,j,k)  ) * c1_NxNyNz
           RWWsumT_xyz = RWWsumT_xyz  + (   Rh(i,j,k) * W(i,j,k) * W(i,j,k)  ) * c1_NxNyNz
          CTKEsumT_xyz = CTKEsumT_xyz + (   0.5d0*Rh(i,j,k)*( + U(i,j,k)*U(i,j,k)                   &
                                                              + V(i,j,k)*V(i,j,k)                   &
                                                              + W(i,j,k)*W(i,j,k) )  ) * c1_NxNyNz   

            PDsumT_xyz =    PDsumT_xyz + ( + P(i,j,k)*( + dUdx(i,j,k)                  &
                                                        + dVdy(i,j,k)                  &
                                                        + dWdz(i,j,k) )  )* c1_NxNyNz   

         ForcesumT_xyz = ForcesumT_xyz + ( f1*U(i,j,k) + f2*V(i,j,k) + f3*W(i,j,k) ) * c1_NxNyNz   

          dUdxsumT_xyz =   dUdxsumT_xyz + ( dUdx(i,j,k) ) * c1_NxNyNz
          dVdxsumT_xyz =   dVdxsumT_xyz + ( dVdx(i,j,k) ) * c1_NxNyNz
          dWdxsumT_xyz =   dWdxsumT_xyz + ( dWdx(i,j,k) ) * c1_NxNyNz
          dUdysumT_xyz =   dUdysumT_xyz + ( dUdy(i,j,k) ) * c1_NxNyNz
          dVdysumT_xyz =   dVdysumT_xyz + ( dVdy(i,j,k) ) * c1_NxNyNz
          dWdysumT_xyz =   dWdysumT_xyz + ( dWdy(i,j,k) ) * c1_NxNyNz
          dUdzsumT_xyz =   dUdzsumT_xyz + ( dUdz(i,j,k) ) * c1_NxNyNz
          dVdzsumT_xyz =   dVdzsumT_xyz + ( dVdz(i,j,k) ) * c1_NxNyNz
          dWdzsumT_xyz =   dWdzsumT_xyz + ( dWdz(i,j,k) ) * c1_NxNyNz

         dUdx2sumT_xyz =  dUdx2sumT_xyz + ( dUdx(i,j,k)*dUdx(i,j,k) ) * c1_NxNyNz
         dVdx2sumT_xyz =  dVdx2sumT_xyz + ( dVdx(i,j,k)*dVdx(i,j,k) ) * c1_NxNyNz
         dWdx2sumT_xyz =  dWdx2sumT_xyz + ( dWdx(i,j,k)*dWdx(i,j,k) ) * c1_NxNyNz
         dUdy2sumT_xyz =  dUdy2sumT_xyz + ( dUdy(i,j,k)*dUdy(i,j,k) ) * c1_NxNyNz
         dVdy2sumT_xyz =  dVdy2sumT_xyz + ( dVdy(i,j,k)*dVdy(i,j,k) ) * c1_NxNyNz
         dWdy2sumT_xyz =  dWdy2sumT_xyz + ( dWdy(i,j,k)*dWdy(i,j,k) ) * c1_NxNyNz
         dUdz2sumT_xyz =  dUdz2sumT_xyz + ( dUdz(i,j,k)*dUdz(i,j,k) ) * c1_NxNyNz
         dVdz2sumT_xyz =  dVdz2sumT_xyz + ( dVdz(i,j,k)*dVdz(i,j,k) ) * c1_NxNyNz
         dWdz2sumT_xyz =  dWdz2sumT_xyz + ( dWdz(i,j,k)*dWdz(i,j,k) ) * c1_NxNyNz

         dUdx3sumT_xyz =  dUdx3sumT_xyz + ( dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k) ) * c1_NxNyNz
         dVdx3sumT_xyz =  dVdx3sumT_xyz + ( dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k) ) * c1_NxNyNz
         dWdx3sumT_xyz =  dWdx3sumT_xyz + ( dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k) ) * c1_NxNyNz
         dUdy3sumT_xyz =  dUdy3sumT_xyz + ( dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k) ) * c1_NxNyNz
         dVdy3sumT_xyz =  dVdy3sumT_xyz + ( dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k) ) * c1_NxNyNz
         dWdy3sumT_xyz =  dWdy3sumT_xyz + ( dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k) ) * c1_NxNyNz
         dUdz3sumT_xyz =  dUdz3sumT_xyz + ( dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k) ) * c1_NxNyNz
         dVdz3sumT_xyz =  dVdz3sumT_xyz + ( dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k) ) * c1_NxNyNz
         dWdz3sumT_xyz =  dWdz3sumT_xyz + ( dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k) ) * c1_NxNyNz

         dUdx4sumT_xyz =  dUdx4sumT_xyz + ( dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k)*dUdx(i,j,k) ) * c1_NxNyNz
         dVdx4sumT_xyz =  dVdx4sumT_xyz + ( dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k)*dVdx(i,j,k) ) * c1_NxNyNz
         dWdx4sumT_xyz =  dWdx4sumT_xyz + ( dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k)*dWdx(i,j,k) ) * c1_NxNyNz
         dUdy4sumT_xyz =  dUdy4sumT_xyz + ( dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k)*dUdy(i,j,k) ) * c1_NxNyNz
         dVdy4sumT_xyz =  dVdy4sumT_xyz + ( dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k)*dVdy(i,j,k) ) * c1_NxNyNz
         dWdy4sumT_xyz =  dWdy4sumT_xyz + ( dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k)*dWdy(i,j,k) ) * c1_NxNyNz
         dUdz4sumT_xyz =  dUdz4sumT_xyz + ( dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k)*dUdz(i,j,k) ) * c1_NxNyNz
         dVdz4sumT_xyz =  dVdz4sumT_xyz + ( dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k)*dVdz(i,j,k) ) * c1_NxNyNz
         dWdz4sumT_xyz =  dWdz4sumT_xyz + ( dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k)*dWdz(i,j,k) ) * c1_NxNyNz

          dPdxsumT_xyz =   dPdxsumT_xyz + ( dPdx(i,j,k) ) * c1_NxNyNz
          dPdysumT_xyz =   dPdysumT_xyz + ( dPdy(i,j,k) ) * c1_NxNyNz
          dPdzsumT_xyz =   dPdzsumT_xyz + ( dPdz(i,j,k) ) * c1_NxNyNz
          dTdxsumT_xyz =   dTdxsumT_xyz + ( dTdx(i,j,k) ) * c1_NxNyNz
          dTdysumT_xyz =   dTdysumT_xyz + ( dTdy(i,j,k) ) * c1_NxNyNz
          dTdzsumT_xyz =   dTdzsumT_xyz + ( dTdz(i,j,k) ) * c1_NxNyNz

         dPdx2sumT_xyz =   dPdx2sumT_xyz + ( dPdx(i,j,k)*dPdx(i,j,k) ) * c1_NxNyNz
         dPdy2sumT_xyz =   dPdy2sumT_xyz + ( dPdy(i,j,k)*dPdy(i,j,k) ) * c1_NxNyNz
         dPdz2sumT_xyz =   dPdz2sumT_xyz + ( dPdz(i,j,k)*dPdz(i,j,k) ) * c1_NxNyNz
         dTdx2sumT_xyz =   dTdx2sumT_xyz + ( dTdx(i,j,k)*dTdx(i,j,k) ) * c1_NxNyNz
         dTdy2sumT_xyz =   dTdy2sumT_xyz + ( dTdy(i,j,k)*dTdy(i,j,k) ) * c1_NxNyNz
         dTdz2sumT_xyz =   dTdz2sumT_xyz + ( dTdz(i,j,k)*dTdz(i,j,k) ) * c1_NxNyNz

         dPdx3sumT_xyz =   dPdx3sumT_xyz + ( dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k) ) * c1_NxNyNz
         dPdy3sumT_xyz =   dPdy3sumT_xyz + ( dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k) ) * c1_NxNyNz
         dPdz3sumT_xyz =   dPdz3sumT_xyz + ( dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k) ) * c1_NxNyNz
         dTdx3sumT_xyz =   dTdx3sumT_xyz + ( dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k) ) * c1_NxNyNz
         dTdy3sumT_xyz =   dTdy3sumT_xyz + ( dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k) ) * c1_NxNyNz
         dTdz3sumT_xyz =   dTdz3sumT_xyz + ( dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k) ) * c1_NxNyNz

         dPdx4sumT_xyz =   dPdx4sumT_xyz + ( dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k)*dPdx(i,j,k) ) * c1_NxNyNz
         dPdy4sumT_xyz =   dPdy4sumT_xyz + ( dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k)*dPdy(i,j,k) ) * c1_NxNyNz
         dPdz4sumT_xyz =   dPdz4sumT_xyz + ( dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k)*dPdz(i,j,k) ) * c1_NxNyNz
         dTdx4sumT_xyz =   dTdx4sumT_xyz + ( dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k)*dTdx(i,j,k) ) * c1_NxNyNz
         dTdy4sumT_xyz =   dTdy4sumT_xyz + ( dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k)*dTdy(i,j,k) ) * c1_NxNyNz
         dTdz4sumT_xyz =   dTdz4sumT_xyz + ( dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k)*dTdz(i,j,k) ) * c1_NxNyNz

          dRhdxsumT_xyz =   dRhdxsumT_xyz + ( dRhdx(i,j,k) ) * c1_NxNyNz
          dRhdysumT_xyz =   dRhdysumT_xyz + ( dRhdy(i,j,k) ) * c1_NxNyNz
          dRhdzsumT_xyz =   dRhdzsumT_xyz + ( dRhdz(i,j,k) ) * c1_NxNyNz

         dRhdx2sumT_xyz =   dRhdx2sumT_xyz + ( dRhdx(i,j,k)*dRhdx(i,j,k) ) * c1_NxNyNz
         dRhdy2sumT_xyz =   dRhdy2sumT_xyz + ( dRhdy(i,j,k)*dRhdy(i,j,k) ) * c1_NxNyNz
         dRhdz2sumT_xyz =   dRhdz2sumT_xyz + ( dRhdz(i,j,k)*dRhdz(i,j,k) ) * c1_NxNyNz

         dRhdx3sumT_xyz =   dRhdx3sumT_xyz + ( dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k) ) * c1_NxNyNz
         dRhdy3sumT_xyz =   dRhdy3sumT_xyz + ( dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k) ) * c1_NxNyNz
         dRhdz3sumT_xyz =   dRhdz3sumT_xyz + ( dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k) ) * c1_NxNyNz

         dRhdx4sumT_xyz =   dRhdx4sumT_xyz + ( dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k)*dRhdx(i,j,k) ) * c1_NxNyNz
         dRhdy4sumT_xyz =   dRhdy4sumT_xyz + ( dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k)*dRhdy(i,j,k) ) * c1_NxNyNz
         dRhdz4sumT_xyz =   dRhdz4sumT_xyz + ( dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k)*dRhdz(i,j,k) ) * c1_NxNyNz

         dYsdxsumT_xyz =   dYsdxsumT_xyz + ( dYsdx(i,j,k) ) * c1_NxNyNz
         dYsdysumT_xyz =   dYsdysumT_xyz + ( dYsdy(i,j,k) ) * c1_NxNyNz
         dYsdzsumT_xyz =   dYsdzsumT_xyz + ( dYsdz(i,j,k) ) * c1_NxNyNz

        dYsdx2sumT_xyz =   dYsdx2sumT_xyz + ( dYsdx(i,j,k)*dYsdx(i,j,k) ) * c1_NxNyNz
        dYsdy2sumT_xyz =   dYsdy2sumT_xyz + ( dYsdy(i,j,k)*dYsdy(i,j,k) ) * c1_NxNyNz
        dYsdz2sumT_xyz =   dYsdz2sumT_xyz + ( dYsdz(i,j,k)*dYsdz(i,j,k) ) * c1_NxNyNz

        dYsdx3sumT_xyz =   dYsdx3sumT_xyz + ( dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k) ) * c1_NxNyNz
        dYsdy3sumT_xyz =   dYsdy3sumT_xyz + ( dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k) ) * c1_NxNyNz
        dYsdz3sumT_xyz =   dYsdz3sumT_xyz + ( dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k) ) * c1_NxNyNz

        dYsdx4sumT_xyz =   dYsdx4sumT_xyz + ( dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k)*dYsdx(i,j,k) ) * c1_NxNyNz
        dYsdy4sumT_xyz =   dYsdy4sumT_xyz + ( dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k)*dYsdy(i,j,k) ) * c1_NxNyNz
        dYsdz4sumT_xyz =   dYsdz4sumT_xyz + ( dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k)*dYsdz(i,j,k) ) * c1_NxNyNz



          UdPdxsumT_xyz =   UdPdxsumT_xyz + U(i,j,k)*( dPdx(i,j,k) ) * c1_NxNyNz
          VdPdysumT_xyz =   VdPdysumT_xyz + V(i,j,k)*( dPdy(i,j,k) ) * c1_NxNyNz
          WdPdzsumT_xyz =   WdPdzsumT_xyz + W(i,j,k)*( dPdz(i,j,k) ) * c1_NxNyNz


            EpsisumT_xyz = EpsisumT_xyz + (   Epsi_i   ) * c1_NxNyNz

            rkensumT_xyz = rkensumT_xyz + (     Rh(i,j,k)*U(i,j,k)*U(i,j,k)                      &
                                             +  Rh(i,j,k)*V(i,j,k)*V(i,j,k)                      &
                                             +  Rh(i,j,k)*W(i,j,k)*W(i,j,k)   ) * c1_NxNyNz


          End Do; End Do
        End Do


        !---/* MPI */---

        !---------------------------------------------------
        Call mp_allsumr(   RhsumT_xyz  )
        Call mp_allsumr(   Rh2sumT_xyz )
        Call mp_allsumr(   Rh3sumT_xyz )
        Call mp_allsumr(   Rh4sumT_xyz )

        Call mp_allsumr(    TsumT_xyz  )
        Call mp_allsumr(    PsumT_xyz  )
        Call mp_allsumr(    UsumT_xyz  )
        Call mp_allsumr(    VsumT_xyz  )
        Call mp_allsumr(    WsumT_xyz  )
        Call mp_allsumr(   YssumT_xyz  )
        Call mp_allsumr(  UVWsumT_xyz  )

        Call mp_allsumr(   musumT_xyz  )
        Call mp_allsumr(   nusumT_xyz  )
        Call mp_allsumr(  mu2sumT_xyz  )

        Call mp_allsumr(    TTsumT_xyz )
        Call mp_allsumr(    PPsumT_xyz )
        Call mp_allsumr(    UUsumT_xyz )
        Call mp_allsumr(    VVsumT_xyz )
        Call mp_allsumr(    WWsumT_xyz )
        Call mp_allsumr(   Ys2sumT_xyz )
        Call mp_allsumr(  UVW2sumT_xyz )

        Call mp_allsumr(    T3sumT_xyz )
        Call mp_allsumr(    P3sumT_xyz )
        Call mp_allsumr(    U3sumT_xyz )
        Call mp_allsumr(    V3sumT_xyz )
        Call mp_allsumr(    W3sumT_xyz )
        Call mp_allsumr(   Ys3sumT_xyz )

        Call mp_allsumr(    T4sumT_xyz )
        Call mp_allsumr(    P4sumT_xyz )
        Call mp_allsumr(    U4sumT_xyz )
        Call mp_allsumr(    V4sumT_xyz )
        Call mp_allsumr(    W4sumT_xyz )
        Call mp_allsumr(   Ys4sumT_xyz )

        Call mp_allsumr(   RUUsumT_xyz )
        Call mp_allsumr(   RVVsumT_xyz )
        Call mp_allsumr(   RWWsumT_xyz )
        Call mp_allsumr(  CTKEsumT_xyz )

        Call mp_allsumr(    PDsumT_xyz )
        Call mp_allsumr( ForcesumT_xyz )

        Call mp_allsumr(   dUdxsumT_xyz )
        Call mp_allsumr(   dVdxsumT_xyz )
        Call mp_allsumr(   dWdxsumT_xyz )
        Call mp_allsumr(   dUdysumT_xyz )
        Call mp_allsumr(   dVdysumT_xyz )
        Call mp_allsumr(   dWdysumT_xyz )
        Call mp_allsumr(   dUdzsumT_xyz )
        Call mp_allsumr(   dVdzsumT_xyz )
        Call mp_allsumr(   dWdzsumT_xyz )

        Call mp_allsumr(   dUdx2sumT_xyz )
        Call mp_allsumr(   dVdx2sumT_xyz )
        Call mp_allsumr(   dWdx2sumT_xyz )
        Call mp_allsumr(   dUdy2sumT_xyz )
        Call mp_allsumr(   dVdy2sumT_xyz )
        Call mp_allsumr(   dWdy2sumT_xyz )
        Call mp_allsumr(   dUdz2sumT_xyz )
        Call mp_allsumr(   dVdz2sumT_xyz )
        Call mp_allsumr(   dWdz2sumT_xyz )

        Call mp_allsumr(   dUdx3sumT_xyz )
        Call mp_allsumr(   dVdx3sumT_xyz )
        Call mp_allsumr(   dWdx3sumT_xyz )
        Call mp_allsumr(   dUdy3sumT_xyz )
        Call mp_allsumr(   dVdy3sumT_xyz )
        Call mp_allsumr(   dWdy3sumT_xyz )
        Call mp_allsumr(   dUdz3sumT_xyz )
        Call mp_allsumr(   dVdz3sumT_xyz )
        Call mp_allsumr(   dWdz3sumT_xyz )

        Call mp_allsumr(   dUdx4sumT_xyz )
        Call mp_allsumr(   dVdx4sumT_xyz )
        Call mp_allsumr(   dWdx4sumT_xyz )
        Call mp_allsumr(   dUdy4sumT_xyz )
        Call mp_allsumr(   dVdy4sumT_xyz )
        Call mp_allsumr(   dWdy4sumT_xyz )
        Call mp_allsumr(   dUdz4sumT_xyz )
        Call mp_allsumr(   dVdz4sumT_xyz )
        Call mp_allsumr(   dWdz4sumT_xyz )

        Call mp_allsumr(    dPdxsumT_xyz )
        Call mp_allsumr(    dPdysumT_xyz )
        Call mp_allsumr(    dPdzsumT_xyz )
        Call mp_allsumr(    dTdxsumT_xyz )
        Call mp_allsumr(    dTdysumT_xyz )
        Call mp_allsumr(    dTdzsumT_xyz )

        Call mp_allsumr(   dPdx2sumT_xyz )
        Call mp_allsumr(   dPdy2sumT_xyz )
        Call mp_allsumr(   dPdz2sumT_xyz )
        Call mp_allsumr(   dTdx2sumT_xyz )
        Call mp_allsumr(   dTdy2sumT_xyz )
        Call mp_allsumr(   dTdz2sumT_xyz )

        Call mp_allsumr(   dPdx3sumT_xyz )
        Call mp_allsumr(   dPdy3sumT_xyz )
        Call mp_allsumr(   dPdz3sumT_xyz )
        Call mp_allsumr(   dTdx3sumT_xyz )
        Call mp_allsumr(   dTdy3sumT_xyz )
        Call mp_allsumr(   dTdz3sumT_xyz )

        Call mp_allsumr(   dPdx4sumT_xyz )
        Call mp_allsumr(   dPdy4sumT_xyz )
        Call mp_allsumr(   dPdz4sumT_xyz )
        Call mp_allsumr(   dTdx4sumT_xyz )
        Call mp_allsumr(   dTdy4sumT_xyz )
        Call mp_allsumr(   dTdz4sumT_xyz )

        Call mp_allsumr(    dRhdxsumT_xyz )
        Call mp_allsumr(    dRhdysumT_xyz )
        Call mp_allsumr(    dRhdzsumT_xyz )

        Call mp_allsumr(   dRhdx2sumT_xyz )
        Call mp_allsumr(   dRhdy2sumT_xyz )
        Call mp_allsumr(   dRhdz2sumT_xyz )

        Call mp_allsumr(   dRhdx3sumT_xyz )
        Call mp_allsumr(   dRhdy3sumT_xyz )
        Call mp_allsumr(   dRhdz3sumT_xyz )

        Call mp_allsumr(   dRhdx4sumT_xyz )
        Call mp_allsumr(   dRhdy4sumT_xyz )
        Call mp_allsumr(   dRhdz4sumT_xyz )


        Call mp_allsumr(    dYsdxsumT_xyz )
        Call mp_allsumr(    dYsdysumT_xyz )
        Call mp_allsumr(    dYsdzsumT_xyz )

        Call mp_allsumr(   dYsdx2sumT_xyz )
        Call mp_allsumr(   dYsdy2sumT_xyz )
        Call mp_allsumr(   dYsdz2sumT_xyz )

        Call mp_allsumr(   dYsdx3sumT_xyz )
        Call mp_allsumr(   dYsdy3sumT_xyz )
        Call mp_allsumr(   dYsdz3sumT_xyz )

        Call mp_allsumr(   dYsdx4sumT_xyz )
        Call mp_allsumr(   dYsdy4sumT_xyz )
        Call mp_allsumr(   dYsdz4sumT_xyz )

        Call mp_allsumr(   EpsisumT_xyz )

        Call mp_allsumr(  UdPdxsumT_xyz )
        Call mp_allsumr(  VdPdysumT_xyz )
        Call mp_allsumr(  WdPdzsumT_xyz )

        Call mp_allsumr(   rkensumT_xyz )


        Return
      End Subroutine
!================================================================================================================================
!========= [ Subroutine Sumwrite_T ] ============================================================================================
!================================================================================================================================
!      Subroutine Sumwrite_T_xyz
!
!        Use module_mpi
!        Use Definition
!        Implicit none ! à√ñŸÇÃå^êÈåæÇÇµÇ»Ç¢
!        Integer i, j
!        Integer fn0, fn1, fn2, fn3, fn4, fn5
!
!        !=== fname ===
!        fn5 =   nt       / 100000
!        fn4 = ( nt - fn5 * 100000     ) / 10000
!        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
!        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
!        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
!        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1
!
!
!
!        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
!        Open(32, File = 'output/PJET_Sum_time_'                                                                   &
!                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
!                      'bin.', Form = 'unFormatted', Status = 'unknown')
!        Write(32)  nt
!        Write(32)  time
!        Close(32)
!
!        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
!        Open(32, File = 'output/PJET_Sum_F_xyz_'                                                                   &
!                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
!                      'bin.', Form = 'unFormatted', Status = 'unknown')
!
!        Write(32)        TrsumT_xyz
!        Write(32)        PrsumT_xyz
!        Write(32)        UrsumT_xyz
!        Write(32)        VrsumT_xyz
!        Write(32)        WrsumT_xyz
!        Write(32)       YsrsumT_xyz
!
!        Write(32)       TTrsumT_xyz
!        Write(32)       PPrsumT_xyz
!        Write(32)       UUrsumT_xyz
!        Write(32)       VVrsumT_xyz
!        Write(32)       WWrsumT_xyz
!        Write(32)      Ys2rsumT_xyz
!
!        Write(32)       T3rsumT_xyz
!        Write(32)       P3rsumT_xyz
!        Write(32)       U3rsumT_xyz
!        Write(32)       V3rsumT_xyz
!        Write(32)       W3rsumT_xyz
!        Write(32)      Ys3rsumT_xyz
!
!        Write(32)       T4rsumT_xyz
!        Write(32)       P4rsumT_xyz
!        Write(32)       U4rsumT_xyz
!        Write(32)       V4rsumT_xyz
!        Write(32)       W4rsumT_xyz
!        Write(32)      Ys4rsumT_xyz
!
!
!        Write(32)     dUdxrsumT_xyz
!        Write(32)     dVdxrsumT_xyz
!        Write(32)     dWdxrsumT_xyz
!        Write(32)     dUdyrsumT_xyz
!        Write(32)     dVdyrsumT_xyz
!        Write(32)     dWdyrsumT_xyz
!        Write(32)     dUdzrsumT_xyz
!        Write(32)     dVdzrsumT_xyz
!        Write(32)     dWdzrsumT_xyz
!
!        Write(32)    dUdx2rsumT_xyz
!        Write(32)    dVdx2rsumT_xyz
!        Write(32)    dWdx2rsumT_xyz
!        Write(32)    dUdy2rsumT_xyz
!        Write(32)    dVdy2rsumT_xyz
!        Write(32)    dWdy2rsumT_xyz
!        Write(32)    dUdz2rsumT_xyz
!        Write(32)    dVdz2rsumT_xyz
!        Write(32)    dWdz2rsumT_xyz
!
!        Write(32)    dUdx3rsumT_xyz
!        Write(32)    dVdx3rsumT_xyz
!        Write(32)    dWdx3rsumT_xyz
!        Write(32)    dUdy3rsumT_xyz
!        Write(32)    dVdy3rsumT_xyz
!        Write(32)    dWdy3rsumT_xyz
!        Write(32)    dUdz3rsumT_xyz
!        Write(32)    dVdz3rsumT_xyz
!        Write(32)    dWdz3rsumT_xyz
!
!        Write(32)    dUdx4rsumT_xyz
!        Write(32)    dVdx4rsumT_xyz
!        Write(32)    dWdx4rsumT_xyz
!        Write(32)    dUdy4rsumT_xyz
!        Write(32)    dVdy4rsumT_xyz
!        Write(32)    dWdy4rsumT_xyz
!        Write(32)    dUdz4rsumT_xyz
!        Write(32)    dVdz4rsumT_xyz
!        Write(32)    dWdz4rsumT_xyz
!
!        Write(32)     dPdxrsumT_xyz
!        Write(32)     dPdyrsumT_xyz
!        Write(32)     dPdzrsumT_xyz
!        Write(32)     dTdxrsumT_xyz
!        Write(32)     dTdyrsumT_xyz
!        Write(32)     dTdzrsumT_xyz
!
!        Write(32)    dPdx2rsumT_xyz
!        Write(32)    dPdy2rsumT_xyz
!        Write(32)    dPdz2rsumT_xyz
!        Write(32)    dTdx2rsumT_xyz
!        Write(32)    dTdy2rsumT_xyz
!        Write(32)    dTdz2rsumT_xyz
!
!        Write(32)    dPdx3rsumT_xyz
!        Write(32)    dPdy3rsumT_xyz
!        Write(32)    dPdz3rsumT_xyz
!        Write(32)    dTdx3rsumT_xyz
!        Write(32)    dTdy3rsumT_xyz
!        Write(32)    dTdz3rsumT_xyz
!
!        Write(32)    dPdx4rsumT_xyz
!        Write(32)    dPdy4rsumT_xyz
!        Write(32)    dPdz4rsumT_xyz
!        Write(32)    dTdx4rsumT_xyz
!        Write(32)    dTdy4rsumT_xyz
!        Write(32)    dTdz4rsumT_xyz
!
!
!
!        Write(32)    dYsdxrsumT_xyz
!        Write(32)    dYsdyrsumT_xyz
!        Write(32)    dYsdzrsumT_xyz
!
!        Write(32)   dYsdx2rsumT_xyz
!        Write(32)   dYsdy2rsumT_xyz
!        Write(32)   dYsdz2rsumT_xyz
!
!        Write(32)   dYsdx3rsumT_xyz
!        Write(32)   dYsdy3rsumT_xyz
!        Write(32)   dYsdz3rsumT_xyz
!
!        Write(32)   dYsdx4rsumT_xyz
!        Write(32)   dYsdy4rsumT_xyz
!        Write(32)   dYsdz4rsumT_xyz
!
!        Close(32)
!
!
!        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
!        Open(32, File = 'output/PJET_Sum_yz_'                                                                     &
!                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
!                      'bin.', Form = 'unFormatted', Status = 'unknown')
!        Write(32)  RhsumT_xyz
!        Write(32)  Rh2sum_xyz
!        Write(32)  Rh3sum_xyz
!        Write(32)  Rh4sum_xyz
!
!        Write(32)   TsumT_xyz
!        Write(32)   PsumT_xyz
!        Write(32)   UsumT_xyz
!        Write(32)   VsumT_xyz
!        Write(32)   WsumT_xyz
!        Write(32)  YssumT_xyz
!
!        Write(32)   TTsumT_xyz
!        Write(32)   PPsumT_xyz
!        Write(32)   UUsumT_xyz
!        Write(32)   VVsumT_xyz
!        Write(32)   WWsumT_xyz
!        Write(32)  Ys2sumT_xyz
!
!        Write(32)   T3sumT_xyz
!        Write(32)   P3sumT_xyz
!        Write(32)   U3sumT_xyz
!        Write(32)   V3sumT_xyz
!        Write(32)   W3sumT_xyz
!        Write(32)  Ys3sumT_xyz
!
!        Write(32)   T4sumT_xyz
!        Write(32)   P4sumT_xyz
!        Write(32)   U4sumT_xyz
!        Write(32)   V4sumT_xyz
!        Write(32)   W4sumT_xyz
!        Write(32)  Ys4sumT_xyz
!
!        Write(32)  dUdxsumT_xyz
!        Write(32)  dVdxsumT_xyz
!        Write(32)  dWdxsumT_xyz
!        Write(32)  dUdysumT_xyz
!        Write(32)  dVdysumT_xyz
!        Write(32)  dWdysumT_xyz
!        Write(32)  dUdzsumT_xyz
!        Write(32)  dVdzsumT_xyz
!        Write(32)  dWdzsumT_xyz
!
!        Write(32)  dUdx2sumT_xyz
!        Write(32)  dVdx2sumT_xyz
!        Write(32)  dWdx2sumT_xyz
!        Write(32)  dUdy2sumT_xyz
!        Write(32)  dVdy2sumT_xyz
!        Write(32)  dWdy2sumT_xyz
!        Write(32)  dUdz2sumT_xyz
!        Write(32)  dVdz2sumT_xyz
!        Write(32)  dWdz2sumT_xyz
!
!        Write(32)  dUdx3sumT_xyz
!        Write(32)  dVdx3sumT_xyz
!        Write(32)  dWdx3sumT_xyz
!        Write(32)  dUdy3sumT_xyz
!        Write(32)  dVdy3sumT_xyz
!        Write(32)  dWdy3sumT_xyz
!        Write(32)  dUdz3sumT_xyz
!        Write(32)  dVdz3sumT_xyz
!        Write(32)  dWdz3sumT_xyz
!
!        Write(32)  dUdx4sumT_xyz
!        Write(32)  dVdx4sumT_xyz
!        Write(32)  dWdx4sumT_xyz
!        Write(32)  dUdy4sumT_xyz
!        Write(32)  dVdy4sumT_xyz
!        Write(32)  dWdy4sumT_xyz
!        Write(32)  dUdz4sumT_xyz
!        Write(32)  dVdz4sumT_xyz
!        Write(32)  dWdz4sumT_xyz
!
!        Write(32)   dPdxsumT_xyz
!        Write(32)   dPdysumT_xyz
!        Write(32)   dPdzsumT_xyz
!        Write(32)   dTdxsumT_xyz
!        Write(32)   dTdysumT_xyz
!        Write(32)   dTdzsumT_xyz
!
!        Write(32)  dPdx2sumT_xyz
!        Write(32)  dPdy2sumT_xyz
!        Write(32)  dPdz2sumT_xyz
!        Write(32)  dTdx2sumT_xyz
!        Write(32)  dTdy2sumT_xyz
!        Write(32)  dTdz2sumT_xyz
!
!        Write(32)  dPdx3sumT_xyz
!        Write(32)  dPdy3sumT_xyz
!        Write(32)  dPdz3sumT_xyz
!        Write(32)  dTdx3sumT_xyz
!        Write(32)  dTdy3sumT_xyz
!        Write(32)  dTdz3sumT_xyz
!
!        Write(32)  dPdx4sumT_xyz
!        Write(32)  dPdy4sumT_xyz
!        Write(32)  dPdz4sumT_xyz
!        Write(32)  dTdx4sumT_xyz
!        Write(32)  dTdy4sumT_xyz
!        Write(32)  dTdz4sumT_xyz
!
!        Write(32)   dRhdxsumT_xyz
!        Write(32)   dRhdysumT_xyz
!        Write(32)   dRhdzsumT_xyz
!
!        Write(32)  dRhdx2sumT_xyz
!        Write(32)  dRhdy2sumT_xyz
!        Write(32)  dRhdz2sumT_xyz
!
!        Write(32)  dRhdx3sumT_xyz
!        Write(32)  dRhdy3sumT_xyz
!        Write(32)  dRhdz3sumT_xyz
!
!        Write(32)  dRhdx4sumT_xyz
!        Write(32)  dRhdy4sumT_xyz
!        Write(32)  dRhdz4sumT_xyz
!
!
!        Write(32)   dYsdxsumT_xyz
!        Write(32)   dYsdysumT_xyz
!        Write(32)   dYsdzsumT_xyz
!
!        Write(32)  dYsdx2sumT_xyz
!        Write(32)  dYsdy2sumT_xyz
!        Write(32)  dYsdz2sumT_xyz
!
!        Write(32)  dYsdx3sumT_xyz
!        Write(32)  dYsdy3sumT_xyz
!        Write(32)  dYsdz3sumT_xyz
!
!        Write(32)  dYsdx4sumT_xyz
!        Write(32)  dYsdy4sumT_xyz
!        Write(32)  dYsdz4sumT_xyz
!
!        Write(32)  LEpsisumT_xyz
!
!        Close(32)
!
!
!  777   Format(50e25.16e3)
!        Return
!      End Subroutine
!================================================================================================================================
!========= [ Subroutine Sumwrite ] ==============================================================================================
!================================================================================================================================
      Subroutine Sum_Time_development_xyz

        Use module_mpi
        Use Definition
        Implicit none 
        Integer i, j
        Integer fn0, fn1, fn2, fn3, fn4, fn5
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Write(81,777)     &
        Dble(nt),         &
        time+time_0,      &
           UsumT_xyz,     &
           VsumT_xyz,     &
           WsumT_xyz,     &
          RhsumT_xyz,     &
           TsumT_xyz,     &
           PsumT_xyz,     &
        divUsumT_xyz,     &
          musumT_xyz*c1_Re,     &
          nusumT_xyz*c1_Re

        Write(83,777)                                              &!1
        Dble(nt),                                                  &!2
        time+time_0,                                               &!3
         UUsumT_xyz    -     UsumT_xyz*    UsumT_xyz,              &!4
         VVsumT_xyz    -     VsumT_xyz*    VsumT_xyz,              &!5
         WWsumT_xyz    -     WsumT_xyz*    WsumT_xyz,              &!6
         !---Mat1---
          Sqrt(   ( UUsumT_xyz - UsumT_xyz*UsumT_xyz )       &
                + ( VVsumT_xyz - VsumT_xyz*VsumT_xyz )       &
                + ( WWsumT_xyz - WsumT_xyz*WsumT_xyz )  )    &
         /Sqrt(gh*TsumT_xyz*cRT_UU),                         &
         !---Mat2---
          Sqrt(   ( UVW2sumT_xyz - UVWsumT_xyz*UVWsumT_xyz )  )    &
         /Sqrt(     gh*TsumT_xyz*cRT_UU                       ),   &

        Rh2sumT_xyz    -    RhsumT_xyz*   RhsumT_xyz,              &!9
         TTsumT_xyz    -     TsumT_xyz*    TsumT_xyz,              &!10
         PPsumT_xyz    -     PsumT_xyz*    PsumT_xyz,              &!11
         divU2sumT_xyz -  divUsumT_xyz* divUsumT_xyz,              &!12
       ( mu2sumT_xyz    -    musumT_xyz*   musumT_xyz )*c1_Re,     &!13

 ( U3sumT_xyz - 3.0d0*UUsumT_xyz*UsumT_xyz + 2.0d0*UsumT_xyz **3 ) / ( (UUsumT_xyz - UsumT_xyz* UsumT_xyz)**(1.5d0) ),          &!14

 ( U4sumT_xyz - 4.0d0*U3sumT_xyz*UsumT_xyz                                                                                      &!15
              + 6.0d0*UUsumT_xyz*UsumT_xyz*UsumT_xyz - 3.0d0*UsumT_xyz**4 ) / ( (UUsumT_xyz - UsumT_xyz* UsumT_xyz)**(2.0d0) )

        Write(85,777)                                                                                            &
        Dble(nt),                                                                                                &
        time+time_0,                                                                                             &
        Luu_xyz,                                              &
          Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )  /  ( Sqrt( dUdx2sumT_xyz - dUdxsumT_xyz*dUdxsumT_xyz ) ),    &
          Sqrt( VVsumT_xyz - VsumT_xyz*VsumT_xyz )  /  ( Sqrt( dVdy2sumT_xyz - dVdysumT_xyz*dVdysumT_xyz ) ),    &
          Sqrt( WWsumT_xyz - WsumT_xyz*WsumT_xyz )  /  ( Sqrt( dWdz2sumT_xyz - dWdzsumT_xyz*dWdzsumT_xyz ) ),    &

        !---Re_l---
          Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )  /  ( Sqrt( dUdx2sumT_xyz - dUdxsumT_xyz*dUdxsumT_xyz ) )     &
         *Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )                                                               &
          *RhsumT_xyz/musumT_xyz * Re,                                                                           &

        !---eta1---
        ( ( musumT_xyz/RhsumT_xyz*c1_Re)                                                                         &
         *( musumT_xyz/RhsumT_xyz*c1_Re)*                                                                        &
          ( musumT_xyz/RhsumT_xyz*c1_Re)                                                                         &
          / EpsisumT_xyz   )**(1.0d0/4.0d0),                                                                     &

        !---eta2---
        ( (nusumT_xyz*c1_Re)*(nusumT_xyz*c1_Re)*(nusumT_xyz*c1_Re)/ EpsisumT_xyz   )**(1.0d0/4.0d0),             &

        !---Skewness and Flatness---
   ( dUdx3sumT_xyz - 3.0d0*dUdx2sumT_xyz*dUdxsumT_xyz + 2.0d0*dUdxsumT_xyz **3 )                                 &
 / ( (dUdx2sumT_xyz - dUdxsumT_xyz* dUdxsumT_xyz)**(1.5d0) ),                                                    &
   ( dUdx4sumT_xyz - 4.0d0*dUdx3sumT_xyz*dUdxsumT_xyz                                                            &
                   + 6.0d0*dUdx2sumT_xyz*dUdxsumT_xyz*dUdxsumT_xyz - 3.0d0*dUdxsumT_xyz**4 )                     &
    / ( (dUdx2sumT_xyz - dUdxsumT_xyz* dUdxsumT_xyz)**(2.0d0) )

        Write(87,777)                                                                                            &
        Dble(nt),                                                                                                &
        time+time_0,                                                                                             &
        !---Epsi---
           EpsisumT_xyz    ,                                                                                     &
           EpsisumT_xyz*Luu_xyz/(Sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz ))**3.0d0,                               &
         - Epsi_s_sumT_xyz ,                                                                                     &
         - Epsi_c_sumT_xyz ,                                                                                     &
         - Epsi_ih_sumT_xyz,                                                                                     &
          dUdxPsumT_xyz+dVdyPsumT_xyz+dWdzPsumT_xyz

        Write(89,777)                                                                                            &
        Dble(nt),                                                                                                &
        time+time_0,                                                                                             &
        !---Epsi---
            RUUsumT_xyz                                 ,              &!7
            RVVsumT_xyz                                 ,              &!8
            RWWsumT_xyz                                 ,              &!9
           rkensumT_xyz

  777   Format(50e25.16e3)
        Return
      End Subroutine


!================================================================================================================================
!--------------------------------------------------------------------------------------------------------------------------------
!****************[    Inflow_Data_Calc subroutine      ]*************************************************************************
!--------------------------------------------------------------------------------------------------------------------------------
!================================================================================================================================
      Subroutine Inflow_Data_Calc

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Integer n
        Real*4   Count_in(0:Ny)
        Real*4  sum_Uf_in(0:Ny), sum_UfUf_in(0:Ny)
        Real*4  sum_Vf_in(0:Ny), sum_VfVf_in(0:Ny)
        Real*4  sum_Wf_in(0:Ny), sum_WfWf_in(0:Ny)
        Real*4  Uf_in_m(0:Ny)  , UfUf_in_m  (0:Ny)
        Real*4  Vf_in_m(0:Ny)  , VfVf_in_m  (0:Ny)
        Real*4  Wf_in_m(0:Ny)  , WfWf_in_m  (0:Ny)

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Call Inputdata_Inflow  !èâä˙ïœìÆåvéZ

        n_step = int( 0.5d0 + L_in*L_in/(2.0d0*pi*dtime0*D_in) )

      If(my_rank==root) then
        Open( 80,File = 'output/PJET_Log_InitialV.txt', Status = 'unknown' )
        Write(80,*) 'n_step = ', n_step
      End If


        !---[ ägéUï˚íˆéÆéûä‘êiçs ]---
        Do nt_in = 1, n_step ,1
           !=== éûä‘êœï™ ===
           Call Eular_1st

          If(my_rank==root) then
             write(80,*) nt_in
          End If
        End Do

      If(my_rank==root) then
        Close(80)
      End If


        !---[ ìæÇÁÇÍÇΩïœìÆÇï™éUÇ™1Ç∆Ç»ÇÈÇÊÇ§Ç…Ç∑ÇÈ ]---
        !Do j = 1, Ny-1
        !  count_in   (j) = 0.0d0
        !  sum_Uf_in  (j) = 0.0d0
        !  sum_Vf_in  (j) = 0.0d0
        !  sum_Wf_in  (j) = 0.0d0
        !  sum_UfUf_in(j) = 0.0d0
        !  sum_VfVf_in(j) = 0.0d0
        !  sum_WfWf_in(j) = 0.0d0
        !End Do

        !Do j = js, je
        !  Do k = ks, ke
        !  Do i = 1, Nx-1
        !    count_in   (j) = count_in   (j) + 1.0d0
        !    sum_Uf_in  (j) = sum_Uf_in  (j) + Uf_in(i,j,k)                !* c1_NxNz
        !    sum_Vf_in  (j) = sum_Vf_in  (j) + Vf_in(i,j,k)                !* c1_NxNz
        !    sum_Wf_in  (j) = sum_Wf_in  (j) + Wf_in(i,j,k)                !* c1_NxNz
        !    sum_UfUf_in(j) = sum_UfUf_in(j) + Uf_in(i,j,k) * Uf_in(i,j,k) !* c1_NxNz
        !    sum_VfVf_in(j) = sum_VfVf_in(j) + Vf_in(i,j,k) * Vf_in(i,j,k) !* c1_NxNz
        !    sum_WfWf_in(j) = sum_WfWf_in(j) + Wf_in(i,j,k) * Wf_in(i,j,k) !* c1_NxNz
        !  End Do; End Do
        !End Do

        !Call mp_allsum_r12(1,Ny-1,count_in (1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,Uf_in_m  (1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,Vf_in_m  (1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,Wf_in_m  (1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,UfUf_in_m(1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,VfVf_in_m(1:Ny-1))
        !Call mp_allsum_r12(1,Ny-1,WfWf_in_m(1:Ny-1))

        !Do j = 1, Ny-1
        !  Uf_in_m  (j) = sum_Uf_in  (j) / count_in(j)
        !  Vf_in_m  (j) = sum_Vf_in  (j) / count_in(j)
        !  Wf_in_m  (j) = sum_Wf_in  (j) / count_in(j)
        !  UfUf_in_m(j) = sum_UfUf_in(j) / count_in(j)
        !  VfVf_in_m(j) = sum_VfVf_in(j) / count_in(j)
        !  WfWf_in_m(j) = sum_WfWf_in(j) / count_in(j)
        !End Do


        !Do j = js, je
        !  Do k = ks, ke
        !  Do i = 1, Nx-1
        !    Uf_in(i,j,k) = ( Uf_in(i,j,k) - Uf_in_m(j) ) /( Sqrt( UfUf_in_m (j) - Uf_in_m(j)*Uf_in_m(j) ) )
        !    Vf_in(i,j,k) = ( Vf_in(i,j,k) - Vf_in_m(j) ) /( Sqrt( VfVf_in_m (j) - Vf_in_m(j)*Vf_in_m(j) ) )
        !    Wf_in(i,j,k) = ( Wf_in(i,j,k) - Wf_in_m(j) ) /( Sqrt( WfWf_in_m (j) - Wf_in_m(j)*Wf_in_m(j) ) )
        !  End Do; End Do
        !End Do


        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_MPI


        !---ìæÇÁÇÍÇΩinflow_dataÇÃìùåvó ÇåvéZ---
        Do j = 1, Ny-1
          count_in   (j) = 0.0d0
          sum_Uf_in  (j) = 0.0d0
          sum_Vf_in  (j) = 0.0d0
          sum_Wf_in  (j) = 0.0d0
          sum_UfUf_in(j) = 0.0d0
          sum_VfVf_in(j) = 0.0d0
          sum_WfWf_in(j) = 0.0d0
        End Do

        Do j = js, je
          Do k = ks, ke
          Do i = 1, Nx-1
            count_in   (j) = count_in   (j) + 1.0d0
            sum_Uf_in  (j) = sum_Uf_in  (j) + Uf_in(i,j,k)                !* c1_NxNz
            sum_Vf_in  (j) = sum_Vf_in  (j) + Vf_in(i,j,k)                !* c1_NxNz
            sum_Wf_in  (j) = sum_Wf_in  (j) + Wf_in(i,j,k)                !* c1_NxNz
            sum_UfUf_in(j) = sum_UfUf_in(j) + Uf_in(i,j,k) * Uf_in(i,j,k) !* c1_NxNz
            sum_VfVf_in(j) = sum_VfVf_in(j) + Vf_in(i,j,k) * Vf_in(i,j,k) !* c1_NxNz
            sum_WfWf_in(j) = sum_WfWf_in(j) + Wf_in(i,j,k) * Wf_in(i,j,k) !* c1_NxNz
          End Do; End Do
        End Do

        Call mp_allsum_r12(1,Ny-1,count_in (1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,Uf_in_m  (1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,Vf_in_m  (1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,Wf_in_m  (1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,UfUf_in_m(1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,VfVf_in_m(1:Ny-1))
        Call mp_allsum_r12(1,Ny-1,WfWf_in_m(1:Ny-1))

        Do j = 1, Ny-1
          Uf_in_m  (j) = sum_Uf_in  (j) / count_in(j)
          Vf_in_m  (j) = sum_Vf_in  (j) / count_in(j)
          Wf_in_m  (j) = sum_Wf_in  (j) / count_in(j)
          UfUf_in_m(j) = sum_UfUf_in(j) / count_in(j)
          VfVf_in_m(j) = sum_VfVf_in(j) / count_in(j)
          WfWf_in_m(j) = sum_WfWf_in(j) / count_in(j)
        End Do


      If(my_rank==root) then
        Open(80, File = 'output/UVW_in_St.txt', Status = 'unknown' )
        j = (Ny-1)/2
        Write(80,777) Uf_in_m(j),                                  &
                      Vf_in_m(j),                                  &
                      Wf_in_m(j),                                  &
                      Sqrt(UfUf_in_m(j)- Uf_in_m(j)*Uf_in_m(j) ),  &
                      Sqrt(VfVf_in_m(j)- Vf_in_m(j)*Vf_in_m(j) ),  &
                      Sqrt(WfWf_in_m(j)- Wf_in_m(j)*Wf_in_m(j) )
        Close(80)
      End If


      !-------------------------
      !---[ èâä˙íldata_save ]---
      Call Write_UVWf_in





        !==========================================================
        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
        End Do; End Do

        Do k = ks, ke, 1
          If(k==1) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
               U_xy (i,j) =    Uf_in(i,j,k)
            End Do; End Do
          End If

          If(k==(Nz-1)) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
               V_xy (i,j) =    Vf_in(i,j,k)
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then
          Open(201, File = 'output/Uin_xy.txt', Status = 'unknown')
          Do j = 1, Ny-1, 1
          Do i = 1, Nx-1, 1
            Write(201,777)   x(i), y(j),    &
                              U_xy (i,j),   &
                              V_xy (i,j)
          End Do; End Do
          Close(201)
        End If



        !==========================================================
        Do k = 1, Nz-1, 1
        Do i = 1, Nx-1, 1
           U_xz (i,k) = 0.0d0
           V_xz (i,k) = 0.0d0
        End Do; End Do

        Do j = js, je, 1
          If(j==(Ny-1)) then
            Do k = ks, ke, 1
            Do i = 1, Nx-1, 1
               U_xz (i,k) =    Uf_in(i,j,k)
               V_xz (i,k) =    Vf_in(i,j,k)
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    U_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    V_xz(1:Nx-1,1:Nz-1) )

        If(my_rank==root) then
          Open(201, File = 'output/Uin_xz.txt', Status = 'unknown')
          Do k = 1, Nz-1, 1
          Do i = 1, Nx-1, 1
            Write(201,777)   x(i), z(k),    &
                              U_xz (i,k),   &
                              V_xz (i,k)
          End Do; End Do
          Close(201)
        End If

        !==========================================================


   777 Format(50e25.16e3)

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Inputdata_Inflow

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- UÇÃèâä˙íl ---
       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          Call Random_number( Uf_in(i,j,k) )
          Call Random_number( Vf_in(i,j,k) )
          Call Random_number( Wf_in(i,j,k) )
        End Do; End Do; End Do

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1
          Uf_in(i,j,k) = ( Uf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx*dy*dz)
          Vf_in(i,j,k) = ( Vf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx*dy*dz)
          Wf_in(i,j,k) = ( Wf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx*dy*dz)
        End Do; End Do; End Do


        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_MPI

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine BC_Inflow_Data_Calc_Local_MPI

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        Call BC_Periodic_z_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Uf_in )
        Call BC_Periodic_z_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vf_in )
        Call BC_Periodic_z_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wf_in )

        Call BC_Periodic_y_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Uf_in )
        Call BC_Periodic_y_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vf_in )
        Call BC_Periodic_y_dir_MPI_R8(-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wf_in )

        !=== x ===
        Do k = ks, ke, 1
        Do j = js, je, 1
          !--- U ---
          Uf_in(Nx+0,j,k) =  Uf_in(0+1,j,k)
          Uf_in(Nx+1,j,k) =  Uf_in(1+1,j,k)
          Uf_in(Nx+2,j,k) =  Uf_in(2+1,j,k)
          Uf_in(Nx+3,j,k) =  Uf_in(3+1,j,k)
          Uf_in(Nx+4,j,k) =  Uf_in(4+1,j,k)

          Uf_in( 0,j,k)   =  Uf_in(Nx-1-0,j,k)
          Uf_in(-1,j,k)   =  Uf_in(Nx-1-1,j,k)
          Uf_in(-2,j,k)   =  Uf_in(Nx-1-2,j,k)
          Uf_in(-3,j,k)   =  Uf_in(Nx-1-3,j,k)
          Uf_in(-4,j,k)   =  Uf_in(Nx-1-4,j,k)

          !--- V ---
          Vf_in(Nx+0,j,k) =  Vf_in(0+1,j,k)
          Vf_in(Nx+1,j,k) =  Vf_in(1+1,j,k)
          Vf_in(Nx+2,j,k) =  Vf_in(2+1,j,k)
          Vf_in(Nx+3,j,k) =  Vf_in(3+1,j,k)
          Vf_in(Nx+4,j,k) =  Vf_in(4+1,j,k)

          Vf_in( 0,j,k)   =  Vf_in(Nx-1-0,j,k)
          Vf_in(-1,j,k)   =  Vf_in(Nx-1-1,j,k)
          Vf_in(-2,j,k)   =  Vf_in(Nx-1-2,j,k)
          Vf_in(-3,j,k)   =  Vf_in(Nx-1-3,j,k)
          Vf_in(-4,j,k)   =  Vf_in(Nx-1-4,j,k)

          !--- W ---
          Wf_in(Nx+0,j,k) =  Wf_in(0+1,j,k)
          Wf_in(Nx+1,j,k) =  Wf_in(1+1,j,k)
          Wf_in(Nx+2,j,k) =  Wf_in(2+1,j,k)
          Wf_in(Nx+3,j,k) =  Wf_in(3+1,j,k)
          Wf_in(Nx+4,j,k) =  Wf_in(4+1,j,k)

          Wf_in( 0,j,k)   =  Uf_in(Nx-1-0,j,k)
          Wf_in(-1,j,k)   =  Wf_in(Nx-1-1,j,k)
          Wf_in(-2,j,k)   =  Wf_in(Nx-1-2,j,k)
          Wf_in(-3,j,k)   =  Wf_in(Nx-1-3,j,k)
          Wf_in(-4,j,k)   =  Wf_in(Nx-1-4,j,k)
        End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!****************Å@Eular_1stÅ@************************************************************************************************
!=============================================================================================================================
      Subroutine Eular_1st

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_MPI

        Call mp_send_recv_pre_r8(Uf_in,1,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre_r8(Vf_in,1,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call mp_send_recv_pre_r8(Wf_in,1,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

        Call Diff_Eq

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1, 1
          Uf_in(i,j,k) = Uf_in(i,j,k) + dtime0 * dUf(i,j,k) 
          Vf_in(i,j,k) = Vf_in(i,j,k) + dtime0 * dVf(i,j,k) 
          Wf_in(i,j,k) = Wf_in(i,j,k) + dtime0 * dWf(i,j,k) 
        End Do; End Do; End Do

        Return
      End Subroutine
!===========================================================================================================================
!****************Å@Diff_EqÅ@************************************************************************************************
!===========================================================================================================================
      Subroutine Diff_Eq

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

       ! Do k = 1, Nz-1
       ! Do j = 1, Ny-1
        Do k = ks, ke, 1
        Do j = js, je, 1
        Do i = 1, Nx-1

          dUf(i,j,k) = ( - ( - Uf_in(i-1,j  ,k  ) + Uf_in(i  ,j  ,k  ) ) * c1_dx                       &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i+1,j  ,k  ) ) * c1_dx ) * c1_dx   * D_in    &
                     + ( - ( - Uf_in(i  ,j-1,k  ) + Uf_in(i  ,j  ,k  ) ) * c1_dy                       &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i  ,j+1,k  ) ) * c1_dy ) * c1_dy   * D_in    &
                     + ( - ( - Uf_in(i  ,j  ,k-1) + Uf_in(i  ,j  ,k  ) ) * c1_dz                       &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i  ,j  ,k+1) ) * c1_dz ) * c1_dz   * D_in


          dVf(i,j,k) = ( - ( - Vf_in(i-1,j  ,k  ) + Vf_in(i  ,j  ,k  ) ) * c1_dx                       &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i+1,j  ,k  ) ) * c1_dx ) * c1_dx   * D_in    &
                     + ( - ( - Vf_in(i  ,j-1,k  ) + Vf_in(i  ,j  ,k  ) ) * c1_dy                       &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i  ,j+1,k  ) ) * c1_dy ) * c1_dy   * D_in    &
                     + ( - ( - Vf_in(i  ,j  ,k-1) + Vf_in(i  ,j  ,k  ) ) * c1_dz                       &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i  ,j  ,k+1) ) * c1_dz ) * c1_dz   * D_in


          dWf(i,j,k) = ( - ( - Wf_in(i-1,j  ,k  ) + Wf_in(i  ,j  ,k  ) ) * c1_dx                       &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i+1,j  ,k  ) ) * c1_dx ) * c1_dx   * D_in    &
                     + ( - ( - Wf_in(i  ,j-1,k  ) + Wf_in(i  ,j  ,k  ) ) * c1_dy                       &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i  ,j+1,k  ) ) * c1_dy ) * c1_dy   * D_in    &
                     + ( - ( - Wf_in(i  ,j  ,k-1) + Wf_in(i  ,j  ,k  ) ) * c1_dz                       &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i  ,j  ,k+1) ) * c1_dz ) * c1_dz   * D_in



        End Do; End Do; End Do

        Return
      End Subroutine
!===============================================================================================================================
!========= [ Subroutine write_UVWf_in ] ========================================================================================
      Subroutine Write_UVWf_in

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k

        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))

        If(Ityp_data_in .eq. 0) then 

          If( my_rank==root) then
            Open(40, File = 'output/UVW_f_in', Form = 'unFormatted', Status = 'unknown')
          End If

          !--- Uf_in ---
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Uf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) Data_r4
          End If


          !--- Vf_in ---
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Vf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) Data_r4
          End If


          !--- Wf_in ---
          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Wf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) Data_r4
          End If


          If( my_rank==root) then
            Close(40)
          End If

        End If





        If(Ityp_data_in .eq. 1) then 
          !---------------------
          If( my_rank==root) then
            Open(40, File = 'output/Uf_in', Form = 'unFormatted', Status = 'unknown')
          End If

          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Uf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
            Write(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
            Close(40)
          End If


          !---------------------
          If( my_rank==root) then
            Open(40, File = 'output/Vf_in', Form = 'unFormatted', Status = 'unknown')
          End If

          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Vf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
            Write(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
            Close(40)
          End If


          !---------------------
          If( my_rank==root) then
            Open(40, File = 'output/Wf_in', Form = 'unFormatted', Status = 'unknown')
          End If

          Data_r4 = 0.0d0

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Data_r4(i,j,k) = Wf_in(i,j,k)
          End Do; End Do; End Do

          Call mp_sum_r32( -4, Nx+4,  -4, Ny+4,  -4, Nz+4 , Data_r4 )

          If( my_rank==root) then
            Write(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
            Write(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
            Close(40)
          End If


        End If


        Deallocate (Data_r4 )

  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!================================================================================================================================
!========= [ Subroutine Sumwrite_T ] ============================================================================================
!================================================================================================================================
      Subroutine Sumwrite_T

        Use module_mpi
        Use Definition
        Implicit none ! à√ñŸÇÃå^êÈåæÇÇµÇ»Ç¢
        Integer i, j
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
        Open(32, File = 'output/HIT_Sum_yz_'                                                                     &
                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
                      '.bin', Form = 'unFormatted', Status = 'unknown')

        Write(32)       RhsumT
        Write(32)      Rh2sumT
        Write(32)      Rh3sumT
        Write(32)      Rh4sumT

        Write(32)        TsumT
        Write(32)        PsumT
        Write(32)        UsumT
        Write(32)        VsumT
        Write(32)        WsumT
        Write(32)       YssumT

        Write(32)       TTsumT
        Write(32)       PPsumT
        Write(32)       UUsumT
        Write(32)       VVsumT
        Write(32)       WWsumT
        Write(32)      Ys2sumT

        Write(32)       T3sumT
        Write(32)       P3sumT
        Write(32)       U3sumT
        Write(32)       V3sumT
        Write(32)       W3sumT
        Write(32)      Ys3sumT

        Write(32)       T4sumT
        Write(32)       P4sumT
        Write(32)       U4sumT
        Write(32)       V4sumT
        Write(32)       W4sumT
        Write(32)      Ys4sumT

        Write(32)     dUdxsumT
        Write(32)     dVdxsumT
        Write(32)     dWdxsumT
        Write(32)     dUdysumT
        Write(32)     dVdysumT
        Write(32)     dWdysumT
        Write(32)     dUdzsumT
        Write(32)     dVdzsumT
        Write(32)     dWdzsumT

        Write(32)    dUdx2sumT
        Write(32)    dVdx2sumT
        Write(32)    dWdx2sumT
        Write(32)    dUdy2sumT
        Write(32)    dVdy2sumT
        Write(32)    dWdy2sumT
        Write(32)    dUdz2sumT
        Write(32)    dVdz2sumT
        Write(32)    dWdz2sumT

        Write(32)    dUdx3sumT
        Write(32)    dVdx3sumT
        Write(32)    dWdx3sumT
        Write(32)    dUdy3sumT
        Write(32)    dVdy3sumT
        Write(32)    dWdy3sumT
        Write(32)    dUdz3sumT
        Write(32)    dVdz3sumT
        Write(32)    dWdz3sumT

        Write(32)    dUdx4sumT
        Write(32)    dVdx4sumT
        Write(32)    dWdx4sumT
        Write(32)    dUdy4sumT
        Write(32)    dVdy4sumT
        Write(32)    dWdy4sumT
        Write(32)    dUdz4sumT
        Write(32)    dVdz4sumT
        Write(32)    dWdz4sumT

        Write(32)     dPdxsumT
        Write(32)     dPdysumT
        Write(32)     dPdzsumT
        Write(32)     dTdxsumT
        Write(32)     dTdysumT
        Write(32)     dTdzsumT

        Write(32)    dPdx2sumT
        Write(32)    dPdy2sumT
        Write(32)    dPdz2sumT
        Write(32)    dTdx2sumT
        Write(32)    dTdy2sumT
        Write(32)    dTdz2sumT

        Write(32)    dPdx3sumT
        Write(32)    dPdy3sumT
        Write(32)    dPdz3sumT
        Write(32)    dTdx3sumT
        Write(32)    dTdy3sumT
        Write(32)    dTdz3sumT

        Write(32)    dPdx4sumT
        Write(32)    dPdy4sumT
        Write(32)    dPdz4sumT
        Write(32)    dTdx4sumT
        Write(32)    dTdy4sumT
        Write(32)    dTdz4sumT

        Write(32)    dRhdxsumT
        Write(32)    dRhdysumT
        Write(32)    dRhdzsumT

        Write(32)   dRhdx2sumT
        Write(32)   dRhdy2sumT
        Write(32)   dRhdz2sumT

        Write(32)   dRhdx3sumT
        Write(32)   dRhdy3sumT
        Write(32)   dRhdz3sumT

        Write(32)   dRhdx4sumT
        Write(32)   dRhdy4sumT
        Write(32)   dRhdz4sumT


        Write(32)    dYsdxsumT
        Write(32)    dYsdysumT
        Write(32)    dYsdzsumT

        Write(32)   dYsdx2sumT
        Write(32)   dYsdy2sumT
        Write(32)   dYsdz2sumT

        Write(32)   dYsdx3sumT
        Write(32)   dYsdy3sumT
        Write(32)   dYsdz3sumT

        Write(32)   dYsdx4sumT
        Write(32)   dYsdy4sumT
        Write(32)   dYsdz4sumT

        Write(32)     EpsisumT
        Write(32)  Epsi_s_sumT
        Write(32)  Epsi_c_sumT
        Write(32) Epsi_ih_sumT

        Write(32)       UVsumT
        Write(32)       VWsumT
        Write(32)       WUsumT

        Write(32)      UYssumT
        Write(32)      VYssumT
        Write(32)      WYssumT

        Write(32)       UPsumT
        Write(32)       VPsumT
        Write(32)       WPsumT

        Write(32)    dUdxPsumT
        Write(32)    dVdyPsumT
        Write(32)    dWdzPsumT
        Write(32)       nusumT
        Write(32)       musumT
        Write(32)      mu2sumT
        Write(32)     divUsumT
        Write(32)    divU2sumT

        Close(32)


  777   Format(50e25.16e3)
        Return
      End Subroutine

!================================================================================================================================
!========= [ Subroutine Sumwrite_T ] ============================================================================================
!================================================================================================================================
      Subroutine Sumwrite_T_xyz

        Use module_mpi
        Use Definition
        Implicit none ! à√ñŸÇÃå^êÈåæÇÇµÇ»Ç¢
        Integer i, j
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
        Open(32, File = 'output/HIT_Sum_time_'                                                                   &
                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
                      '.bin', Form = 'unFormatted', Status = 'unknown')
        Write(32)  nt
        Write(32)  time
        Close(32)



        !=== ìùåvó ÉfÅ[É^èoóÕópÉtÉ@ÉCÉãÇçÏê¨ ===
        Open(32, File = 'output/HIT_Sum_yz_'                                                                     &
                      //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
                      '.bin', Form = 'unFormatted', Status = 'unknown')
        Write(32)        RhsumT_xyz
        Write(32)       Rh2sumT_xyz
        Write(32)       Rh3sumT_xyz
        Write(32)       Rh4sumT_xyz

        Write(32)         TsumT_xyz
        Write(32)         PsumT_xyz
        Write(32)         UsumT_xyz
        Write(32)         VsumT_xyz
        Write(32)         WsumT_xyz
        Write(32)        YssumT_xyz

        Write(32)        TTsumT_xyz
        Write(32)        PPsumT_xyz
        Write(32)        UUsumT_xyz
        Write(32)        VVsumT_xyz
        Write(32)        WWsumT_xyz
        Write(32)       Ys2sumT_xyz

        Write(32)        T3sumT_xyz
        Write(32)        P3sumT_xyz
        Write(32)        U3sumT_xyz
        Write(32)        V3sumT_xyz
        Write(32)        W3sumT_xyz
        Write(32)       Ys3sumT_xyz

        Write(32)        T4sumT_xyz
        Write(32)        P4sumT_xyz
        Write(32)        U4sumT_xyz
        Write(32)        V4sumT_xyz
        Write(32)        W4sumT_xyz
        Write(32)       Ys4sumT_xyz

        Write(32)      dUdxsumT_xyz
        Write(32)      dVdxsumT_xyz
        Write(32)      dWdxsumT_xyz
        Write(32)      dUdysumT_xyz
        Write(32)      dVdysumT_xyz
        Write(32)      dWdysumT_xyz
        Write(32)      dUdzsumT_xyz
        Write(32)      dVdzsumT_xyz
        Write(32)      dWdzsumT_xyz

        Write(32)     dUdx2sumT_xyz
        Write(32)     dVdx2sumT_xyz
        Write(32)     dWdx2sumT_xyz
        Write(32)     dUdy2sumT_xyz
        Write(32)     dVdy2sumT_xyz
        Write(32)     dWdy2sumT_xyz
        Write(32)     dUdz2sumT_xyz
        Write(32)     dVdz2sumT_xyz
        Write(32)     dWdz2sumT_xyz

        Write(32)     dUdx3sumT_xyz
        Write(32)     dVdx3sumT_xyz
        Write(32)     dWdx3sumT_xyz
        Write(32)     dUdy3sumT_xyz
        Write(32)     dVdy3sumT_xyz
        Write(32)     dWdy3sumT_xyz
        Write(32)     dUdz3sumT_xyz
        Write(32)     dVdz3sumT_xyz
        Write(32)     dWdz3sumT_xyz

        Write(32)     dUdx4sumT_xyz
        Write(32)     dVdx4sumT_xyz
        Write(32)     dWdx4sumT_xyz
        Write(32)     dUdy4sumT_xyz
        Write(32)     dVdy4sumT_xyz
        Write(32)     dWdy4sumT_xyz
        Write(32)     dUdz4sumT_xyz
        Write(32)     dVdz4sumT_xyz
        Write(32)     dWdz4sumT_xyz

        Write(32)      dPdxsumT_xyz
        Write(32)      dPdysumT_xyz
        Write(32)      dPdzsumT_xyz
        Write(32)      dTdxsumT_xyz
        Write(32)      dTdysumT_xyz
        Write(32)      dTdzsumT_xyz

        Write(32)     dPdx2sumT_xyz
        Write(32)     dPdy2sumT_xyz
        Write(32)     dPdz2sumT_xyz
        Write(32)     dTdx2sumT_xyz
        Write(32)     dTdy2sumT_xyz
        Write(32)     dTdz2sumT_xyz

        Write(32)     dPdx3sumT_xyz
        Write(32)     dPdy3sumT_xyz
        Write(32)     dPdz3sumT_xyz
        Write(32)     dTdx3sumT_xyz
        Write(32)     dTdy3sumT_xyz
        Write(32)     dTdz3sumT_xyz

        Write(32)     dPdx4sumT_xyz
        Write(32)     dPdy4sumT_xyz
        Write(32)     dPdz4sumT_xyz
        Write(32)     dTdx4sumT_xyz
        Write(32)     dTdy4sumT_xyz
        Write(32)     dTdz4sumT_xyz

        Write(32)     dRhdxsumT_xyz
        Write(32)     dRhdysumT_xyz
        Write(32)     dRhdzsumT_xyz

        Write(32)    dRhdx2sumT_xyz
        Write(32)    dRhdy2sumT_xyz
        Write(32)    dRhdz2sumT_xyz

        Write(32)    dRhdx3sumT_xyz
        Write(32)    dRhdy3sumT_xyz
        Write(32)    dRhdz3sumT_xyz

        Write(32)    dRhdx4sumT_xyz
        Write(32)    dRhdy4sumT_xyz
        Write(32)    dRhdz4sumT_xyz


        Write(32)     dYsdxsumT_xyz
        Write(32)     dYsdysumT_xyz
        Write(32)     dYsdzsumT_xyz

        Write(32)    dYsdx2sumT_xyz
        Write(32)    dYsdy2sumT_xyz
        Write(32)    dYsdz2sumT_xyz

        Write(32)    dYsdx3sumT_xyz
        Write(32)    dYsdy3sumT_xyz
        Write(32)    dYsdz3sumT_xyz

        Write(32)    dYsdx4sumT_xyz
        Write(32)    dYsdy4sumT_xyz
        Write(32)    dYsdz4sumT_xyz

        Write(32)      EpsisumT_xyz
        Write(32)   Epsi_s_sumT_xyz
        Write(32)   Epsi_c_sumT_xyz
        Write(32)  Epsi_ih_sumT_xyz

        Write(32)        UVsumT_xyz
        Write(32)        VWsumT_xyz
        Write(32)        WUsumT_xyz

        Write(32)       UYssumT_xyz
        Write(32)       VYssumT_xyz
        Write(32)       WYssumT_xyz

        Write(32)        UPsumT_xyz
        Write(32)        VPsumT_xyz
        Write(32)        WPsumT_xyz

        Write(32)     dUdxPsumT_xyz
        Write(32)     dVdyPsumT_xyz
        Write(32)     dWdzPsumT_xyz
        Write(32)        nusumT_xyz
        Write(32)        musumT_xyz
        Write(32)       mu2sumT_xyz
        Write(32)      divUsumT_xyz
        Write(32)     divU2sumT_xyz

        Close(32)


  777   Format(50e25.16e3)
        Return
      End Subroutine

!================================================================================================================================
!========= [ Subroutine Corwrite_T ] ============================================================================================
!================================================================================================================================
      Subroutine Corwrite_T

        Use module_mpi
        Use Definition
        Implicit none ! à√ñŸÇÃå^êÈåæÇÇµÇ»Ç¢
        Integer i, j, n
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        !=== fname ===
        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        !------------------------------------------------------------------------------------------------------------
        Open(32, File = 'output/HIT_Cor_xyz_'                                                          &
                        //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
                        'step.txt', Status = 'unknown' )
        Write(32,*) 'dx Cuu Cvv Cww'
        Do n = 0, (Nx-1)/2, 1
          Write(32,777) Dble(n)*dx,  &
                    Cuu_xyz(n) ,  &
                    Cvv_xyz(n) ,  &
                    Cww_xyz(n) 
        End Do
        Close(32)

       ! !------------------------------------------------------------------------------------------------------------
       ! Open(32, File = 'output/HIT_Cor_xyz_'                                                          &
       !                 //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
       !                 'step.txt', Status = 'unknown' )
       ! Write(32,*) 'dx Cuu Cvv Cww Cpp Crr Ctt Cyy'
       ! Do n = 0, (Nx-1)/2, 1
       !   Write(32,777) Dble(n)*dx,  &
       !             Cuu_xyz(n) ,  &
       !             Cvv_xyz(n) ,  &
       !             Cww_xyz(n) ,  &
       !             Cpp_xyz(n) ,  &
       !             Crr_xyz(n) ,  &
       !             Ctt_xyz(n) ,  &
       !             Cyy_xyz(n)
       ! End Do
       ! Close(32)

       ! !------------------------------------------------------------------------------------------------------------
       ! Open(32, File = 'output/HIT_Cor_C_'                                                          &
       !                 //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//      &
       !                 'step.txt', Status = 'unknown' )
       ! Write(32,*) 'dx Cuu Cvv Cww Cpp Crr Ctt Cyy'
       ! Do n = 0, (Nx-1)/2, 1
       !   Write(32,777) Dble(n)*dx,  &
       !             Cuu_C(n) ,  &
       !             Cvv_C(n) ,  &
       !             Cww_C(n) ,  &
       !             Cpp_C(n) ,  &
       !             Crr_C(n) ,  &
       !             Ctt_C(n) ,  &
       !             Cyy_C(n)
       ! End Do
       ! Close(32)

  777   Format(50e25.16e3)
        Return
      End Subroutine


!===============================================================================================================================
!========= [ Subroutine write_UVWf_in ] ========================================================================================
      Subroutine Read_UVWf_in

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k

        Real*4, dimension(:,:,:),allocatable ::    Data_r4

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------


        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        If(Ityp_data_in .eq. 0) then 

          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))

          !=== Density, Velocity, Temperature ===
          Open(99, File = 'UVW_f_in', Form = 'unFormatted', Status = 'unknown')

          !--- Uf_in ---
          Read(99) Data_r4
          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Uf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do

          !--- Vf_in ---
          Read(99) Data_r4
          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Vf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do

          !--- Wf_in ---
          Read(99) Data_r4
          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Wf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do

        End If


        If(Ityp_data_in .eq. 1) then 
          Allocate (Data_r4 (-4:Nx+4,-4:Ny+4,-4:Nz+4))

          Open(40, File = 'Uf_in', Form = 'unFormatted', Status = 'unknown')
          Read(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
          Read(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
          Close(40)

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Uf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do


          Open(40, File = 'Vf_in', Form = 'unFormatted', Status = 'unknown')
          Read(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
          Read(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
          Close(40)

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Vf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do


          Open(40, File = 'Wf_in', Form = 'unFormatted', Status = 'unknown')
          Read(40) (((Data_r4 (i,j,k),i=-4,(Nx-1)/2       ),j=-4,Ny+4),k=-4,Nz+4)
          Read(40) (((Data_r4 (i,j,k),i=   (Nx-1)/2+1,Nx+4),j=-4,Ny+4),k=-4,Nz+4)
          Close(40)

          Do k = ks,  ke,1
          Do j = js,  je,1
          Do i =  1,Nx-1,1
            Wf_in(i,j,k) = Data_r4(i,j,k)
          End Do; End Do; End Do

        End If


  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!================================================================================================================================
!--------------------------------------------------------------------------------------------------------------------------------
!****************[    Inflow_Data_Calc_CG subroutine      ]**********************************************************************
!--------------------------------------------------------------------------------------------------------------------------------
!================================================================================================================================
      Subroutine Inflow_Data_Calc_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        Integer n
        Integer  Count_in(0:Ny_in)
        Real*8  sum_Uf_in(0:Ny_in), sum_UfUf_in(0:Ny_in)
        Real*8  sum_Vf_in(0:Ny_in), sum_VfVf_in(0:Ny_in)
        Real*8  sum_Wf_in(0:Ny_in), sum_WfWf_in(0:Ny_in)
        Real*8  Uf_in_m(0:Ny_in)  , UfUf_in_m  (0:Ny_in)
        Real*8  Vf_in_m(0:Ny_in)  , VfVf_in_m  (0:Ny_in)
        Real*8  Wf_in_m(0:Ny_in)  , WfWf_in_m  (0:Ny_in)


        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

      If(my_rank==root) then
        If( Ityp_inS == 0) then
          Call Inputdata_Inflow_CG  !èâä˙ïœìÆåvéZ

          n_step = int( 0.5d0 + L_in*L_in/(2.0d0*pi*dtime0*D_in) )

          Open( 80,File = 'output/PJET_Log_InitialV.txt', Status = 'unknown' )
          Write(80,*) 'n_step = ', n_step


          !---[ ägéUï˚íˆéÆéûä‘êiçs ]---
          Do nt_in = 1, n_step ,1
             !=== éûä‘êœï™ ===
             Call Eular_1st_CG
             write(80,*) nt_in
          End Do

          Close(80)
        End If

        If( Ityp_inS == 1 ) Call IFFT_initial


        !---[ ìæÇÁÇÍÇΩïœìÆÇï™éUÇ™1Ç∆Ç»ÇÈÇÊÇ§Ç…Ç∑ÇÈ ]---
        Do j = 1, Ny_in-1
          count_in   (j) = 0
          sum_Uf_in  (j) = 0.0d0
          sum_Vf_in  (j) = 0.0d0
          sum_Wf_in  (j) = 0.0d0
          sum_UfUf_in(j) = 0.0d0
          sum_VfVf_in(j) = 0.0d0
          sum_WfWf_in(j) = 0.0d0
        End Do

        Do j = 1, Ny_in-1
          Do k = 1, Nz_in-1
          Do i = 1, Nx_in-1
            count_in   (j) = count_in   (j) + 1
            sum_Uf_in  (j) = sum_Uf_in  (j) + Uf_in(i,j,k)
            sum_Vf_in  (j) = sum_Vf_in  (j) + Vf_in(i,j,k)
            sum_Wf_in  (j) = sum_Wf_in  (j) + Wf_in(i,j,k)
            sum_UfUf_in(j) = sum_UfUf_in(j) + Uf_in(i,j,k) * Uf_in(i,j,k)
            sum_VfVf_in(j) = sum_VfVf_in(j) + Vf_in(i,j,k) * Vf_in(i,j,k)
            sum_WfWf_in(j) = sum_WfWf_in(j) + Wf_in(i,j,k) * Wf_in(i,j,k)
          End Do; End Do
        End Do

        Do j = 1, Ny_in-1
          Uf_in_m  (j) = sum_Uf_in  (j) / Dble( count_in(j) )
          Vf_in_m  (j) = sum_Vf_in  (j) / Dble( count_in(j) )
          Wf_in_m  (j) = sum_Wf_in  (j) / Dble( count_in(j) )
          UfUf_in_m(j) = sum_UfUf_in(j) / Dble( count_in(j) )
          VfVf_in_m(j) = sum_VfVf_in(j) / Dble( count_in(j) )
          WfWf_in_m(j) = sum_WfWf_in(j) / Dble( count_in(j) )
        End Do

        Do j = 1, Ny_in-1
          Do k = 1, Nz_in-1
          Do i = 1, Nx_in-1
            Uf_in(i,j,k) = ( Uf_in(i,j,k) - Uf_in_m(j) ) /( dSqrt( UfUf_in_m (j) - Uf_in_m(j)*Uf_in_m(j) ) )
            Vf_in(i,j,k) = ( Vf_in(i,j,k) - Vf_in_m(j) ) /( dSqrt( VfVf_in_m (j) - Vf_in_m(j)*Vf_in_m(j) ) )
            Wf_in(i,j,k) = ( Wf_in(i,j,k) - Wf_in_m(j) ) /( dSqrt( WfWf_in_m (j) - Wf_in_m(j)*Wf_in_m(j) ) )
          End Do; End Do
        End Do


        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_CG


        !---ìæÇÁÇÍÇΩinflow_dataÇÃìùåvó ÇåvéZ---
        Do j = 1, Ny_in-1
          count_in   (j) = 0
          sum_Uf_in  (j) = 0.0d0
          sum_Vf_in  (j) = 0.0d0
          sum_Wf_in  (j) = 0.0d0
          sum_UfUf_in(j) = 0.0d0
          sum_VfVf_in(j) = 0.0d0
          sum_WfWf_in(j) = 0.0d0
        End Do

        Do j = 1, Ny_in-1
          Do k = 1, Nz_in-1
          Do i = 1, Nx_in-1
            count_in   (j) = count_in   (j) + 1
            sum_Uf_in  (j) = sum_Uf_in  (j) + Uf_in(i,j,k)
            sum_Vf_in  (j) = sum_Vf_in  (j) + Vf_in(i,j,k)
            sum_Wf_in  (j) = sum_Wf_in  (j) + Wf_in(i,j,k)
            sum_UfUf_in(j) = sum_UfUf_in(j) + Uf_in(i,j,k) * Uf_in(i,j,k)
            sum_VfVf_in(j) = sum_VfVf_in(j) + Vf_in(i,j,k) * Vf_in(i,j,k)
            sum_WfWf_in(j) = sum_WfWf_in(j) + Wf_in(i,j,k) * Wf_in(i,j,k)
          End Do; End Do
        End Do

        Do j = 1, Ny_in-1
          Uf_in_m  (j) = sum_Uf_in  (j) / Dble( count_in(j) )
          Vf_in_m  (j) = sum_Vf_in  (j) / Dble( count_in(j) )
          Wf_in_m  (j) = sum_Wf_in  (j) / Dble( count_in(j) )
          UfUf_in_m(j) = sum_UfUf_in(j) / Dble( count_in(j) )
          VfVf_in_m(j) = sum_VfVf_in(j) / Dble( count_in(j) )
          WfWf_in_m(j) = sum_WfWf_in(j) / Dble( count_in(j) )

        End Do                         

        Open(80, File = 'output/UVW_in_St.txt', Status = 'unknown' )
        j = (Ny_in-1)/2
        Write(80,777) Uf_in_m(j),                                  &
                      Vf_in_m(j),                                  &
                      Wf_in_m(j),                                  &
                      dSqrt(UfUf_in_m(j)- Uf_in_m(j)*Uf_in_m(j) ),  &
                      dSqrt(VfVf_in_m(j)- Vf_in_m(j)*Vf_in_m(j) ),  &
                      dSqrt(WfWf_in_m(j)- Wf_in_m(j)*Wf_in_m(j) )
        Close(80)

        !-------------------------
        !---[ èâä˙íldata_save ]---
        Call Write_UVWf_in_CG


        !==========================================================
        Open(80, File = 'output/Uin_xy_CG.txt', Status = 'unknown')
        Do j = 1, Ny_in-1
        Do i = 1, Nx_in-1
          Write(80,777)                   &
                     x_in(i), y_in(j),    &
                     Uf_in(i,j,1),        &
                     Uf_in(i,j,2)

        End Do; End Do
        Close(80)


        !==========================================================
        Open(80, File = 'output/Uin_xz_CG.txt', Status = 'unknown')
        j = (Ny_in-1)/2
        Do k = 1, Nz_in-1
        Do i = 1, Nx_in-1
          Write(80,777)                   &
                     x_in(i), z_in(k),    &
                     Uf_in(i,j,k),        &
                     Vf_in(i,j,k),        &
                     Wf_in(i,j,k)

        End Do; End Do
        Close(80)

        !==========================================================

      End If

      Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Uf_in )
      Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Vf_in )
      Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Wf_in )

   777 Format(50e25.16e3)

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine Inputdata_Inflow_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------


        !--- UÇÃèâä˙íl ---
        Call Random_number( Uf_in )
        Call Random_number( Vf_in )
        Call Random_number( Wf_in )

        Do k = 1, Nz_in-1
        Do j = 1, Ny_in-1
        Do i = 1, Nx_in-1
          Uf_in(i,j,k) = ( Uf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx_in*dy_in*dz_in)
          Vf_in(i,j,k) = ( Vf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx_in*dy_in*dz_in)
          Wf_in(i,j,k) = ( Wf_in(i,j,k)-0.5d0 ) / (Sqrt(1.0d0/12.0d0)) / Sqrt(dx_in*dy_in*dz_in)
        End Do; End Do; End Do

        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_CG

        Return
      End Subroutine
!=============================================================================================================================
!=============================================================================================================================
      Subroutine BC_Inflow_Data_Calc_Local_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        !=== x ===
        Do k = 1, Nz_in-1, 1
        Do j = 1, Ny_in-1, 1
          !--- U ---
          Uf_in(Nx_in+0,j,k) =  Uf_in(0+1,j,k)
          Uf_in(Nx_in+1,j,k) =  Uf_in(1+1,j,k)
          Uf_in(Nx_in+2,j,k) =  Uf_in(2+1,j,k)
          Uf_in(Nx_in+3,j,k) =  Uf_in(3+1,j,k)
          Uf_in(Nx_in+4,j,k) =  Uf_in(4+1,j,k)

          Uf_in( 0,j,k)   =  Uf_in(Nx_in-1-0,j,k)
          Uf_in(-1,j,k)   =  Uf_in(Nx_in-1-1,j,k)
          Uf_in(-2,j,k)   =  Uf_in(Nx_in-1-2,j,k)
          Uf_in(-3,j,k)   =  Uf_in(Nx_in-1-3,j,k)
          Uf_in(-4,j,k)   =  Uf_in(Nx_in-1-4,j,k)

          !--- V ---
          Vf_in(Nx_in+0,j,k) =  Vf_in(0+1,j,k)
          Vf_in(Nx_in+1,j,k) =  Vf_in(1+1,j,k)
          Vf_in(Nx_in+2,j,k) =  Vf_in(2+1,j,k)
          Vf_in(Nx_in+3,j,k) =  Vf_in(3+1,j,k)
          Vf_in(Nx_in+4,j,k) =  Vf_in(4+1,j,k)

          Vf_in( 0,j,k)   =  Vf_in(Nx_in-1-0,j,k)
          Vf_in(-1,j,k)   =  Vf_in(Nx_in-1-1,j,k)
          Vf_in(-2,j,k)   =  Vf_in(Nx_in-1-2,j,k)
          Vf_in(-3,j,k)   =  Vf_in(Nx_in-1-3,j,k)
          Vf_in(-4,j,k)   =  Vf_in(Nx_in-1-4,j,k)

          !--- W ---
          Wf_in(Nx_in+0,j,k) =  Wf_in(0+1,j,k)
          Wf_in(Nx_in+1,j,k) =  Wf_in(1+1,j,k)
          Wf_in(Nx_in+2,j,k) =  Wf_in(2+1,j,k)
          Wf_in(Nx_in+3,j,k) =  Wf_in(3+1,j,k)
          Wf_in(Nx_in+4,j,k) =  Wf_in(4+1,j,k)

          Wf_in( 0,j,k)   =  Uf_in(Nx_in-1-0,j,k)
          Wf_in(-1,j,k)   =  Wf_in(Nx_in-1-1,j,k)
          Wf_in(-2,j,k)   =  Wf_in(Nx_in-1-2,j,k)
          Wf_in(-3,j,k)   =  Wf_in(Nx_in-1-3,j,k)
          Wf_in(-4,j,k)   =  Wf_in(Nx_in-1-4,j,k)
        End Do; End Do


        !=== y ===
        Do k = 1, Nz_in-1, 1
        Do i = -4, Nx_in+4, 1
          !--- U ---
          Uf_in(i,Ny_in+0,k) =  Uf_in(i,0+1,k)
          Uf_in(i,Ny_in+1,k) =  Uf_in(i,1+1,k)
          Uf_in(i,Ny_in+2,k) =  Uf_in(i,2+1,k)
          Uf_in(i,Ny_in+3,k) =  Uf_in(i,3+1,k)
          Uf_in(i,Ny_in+4,k) =  Uf_in(i,4+1,k)

          Uf_in(i, 0,k)      =  Uf_in(i,Ny_in-1-0,k)
          Uf_in(i,-1,k)      =  Uf_in(i,Ny_in-1-1,k)
          Uf_in(i,-2,k)      =  Uf_in(i,Ny_in-1-2,k)
          Uf_in(i,-3,k)      =  Uf_in(i,Ny_in-1-3,k)
          Uf_in(i,-4,k)      =  Uf_in(i,Ny_in-1-4,k)

          !--- V ---
          Vf_in(i,Ny_in+0,k) =  Vf_in(i,0+1,k)
          Vf_in(i,Ny_in+1,k) =  Vf_in(i,1+1,k)
          Vf_in(i,Ny_in+2,k) =  Vf_in(i,2+1,k)
          Vf_in(i,Ny_in+3,k) =  Vf_in(i,3+1,k)
          Vf_in(i,Ny_in+4,k) =  Vf_in(i,4+1,k)

          Vf_in(i, 0,k)      =  Vf_in(i,Ny_in-1-0,k)
          Vf_in(i,-1,k)      =  Vf_in(i,Ny_in-1-1,k)
          Vf_in(i,-2,k)      =  Vf_in(i,Ny_in-1-2,k)
          Vf_in(i,-3,k)      =  Vf_in(i,Ny_in-1-3,k)
          Vf_in(i,-4,k)      =  Vf_in(i,Ny_in-1-4,k)

          !--- W ---
          Wf_in(i,Ny_in+0,k) =  Wf_in(i,0+1,k)
          Wf_in(i,Ny_in+1,k) =  Wf_in(i,1+1,k)
          Wf_in(i,Ny_in+2,k) =  Wf_in(i,2+1,k)
          Wf_in(i,Ny_in+3,k) =  Wf_in(i,3+1,k)
          Wf_in(i,Ny_in+4,k) =  Wf_in(i,4+1,k)

          Wf_in(i, 0,k)      =  Wf_in(i,Ny_in-1-0,k)
          Wf_in(i,-1,k)      =  Wf_in(i,Ny_in-1-1,k)
          Wf_in(i,-2,k)      =  Wf_in(i,Ny_in-1-2,k)
          Wf_in(i,-3,k)      =  Wf_in(i,Ny_in-1-3,k)
          Wf_in(i,-4,k)      =  Wf_in(i,Ny_in-1-4,k)

        End Do; End Do

        !=== z ===
        Do j = -4, Ny_in+4, 1
        Do i = -4, Nx_in+4, 1
          !--- U ---
          Uf_in(i,j,Nz_in+0) =  Uf_in(i,j,0+1)
          Uf_in(i,j,Nz_in+1) =  Uf_in(i,j,1+1)
          Uf_in(i,j,Nz_in+2) =  Uf_in(i,j,2+1)
          Uf_in(i,j,Nz_in+3) =  Uf_in(i,j,3+1)
          Uf_in(i,j,Nz_in+4) =  Uf_in(i,j,4+1)

          Uf_in(i,j, 0)   =  Uf_in(i,j,Nz_in-1-0)
          Uf_in(i,j,-1)   =  Uf_in(i,j,Nz_in-1-1)
          Uf_in(i,j,-2)   =  Uf_in(i,j,Nz_in-1-2)
          Uf_in(i,j,-3)   =  Uf_in(i,j,Nz_in-1-3)
          Uf_in(i,j,-4)   =  Uf_in(i,j,Nz_in-1-4)

          !--- V ---
          Vf_in(i,j,Nz_in+0) =  Vf_in(i,j,0+1)
          Vf_in(i,j,Nz_in+1) =  Vf_in(i,j,1+1)
          Vf_in(i,j,Nz_in+2) =  Vf_in(i,j,2+1)
          Vf_in(i,j,Nz_in+3) =  Vf_in(i,j,3+1)
          Vf_in(i,j,Nz_in+4) =  Vf_in(i,j,4+1)

          Vf_in(i,j, 0)      =  Vf_in(i,j,Nz_in-1-0)
          Vf_in(i,j,-1)      =  Vf_in(i,j,Nz_in-1-1)
          Vf_in(i,j,-2)      =  Vf_in(i,j,Nz_in-1-2)
          Vf_in(i,j,-3)      =  Vf_in(i,j,Nz_in-1-3)
          Vf_in(i,j,-4)      =  Vf_in(i,j,Nz_in-1-4)

          !--- W ---
          Wf_in(i,j,Nz_in+0) =  Wf_in(i,j,0+1)
          Wf_in(i,j,Nz_in+1) =  Wf_in(i,j,1+1)
          Wf_in(i,j,Nz_in+2) =  Wf_in(i,j,2+1)
          Wf_in(i,j,Nz_in+3) =  Wf_in(i,j,3+1)
          Wf_in(i,j,Nz_in+4) =  Wf_in(i,j,4+1)

          Wf_in(i,j, 0)      =  Wf_in(i,j,Nz_in-1-0)
          Wf_in(i,j,-1)      =  Wf_in(i,j,Nz_in-1-1)
          Wf_in(i,j,-2)      =  Wf_in(i,j,Nz_in-1-2)
          Wf_in(i,j,-3)      =  Wf_in(i,j,Nz_in-1-3)
          Wf_in(i,j,-4)      =  Wf_in(i,j,Nz_in-1-4)

        End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!****************Å@Eular_1stÅ@************************************************************************************************
!=============================================================================================================================
      Subroutine Eular_1st_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- ã´äEèåèÇÃê›íË ---
        Call BC_Inflow_Data_Calc_Local_CG

        Call Diff_Eq_CG

        Do k = 1, Nz_in-1, 1
        Do j = 1, Ny_in-1, 1
        Do i = 1, Nx_in-1, 1
          Uf_in(i,j,k) = Uf_in(i,j,k) + dtime0 * dUf(i,j,k) 
          Vf_in(i,j,k) = Vf_in(i,j,k) + dtime0 * dVf(i,j,k) 
          Wf_in(i,j,k) = Wf_in(i,j,k) + dtime0 * dWf(i,j,k) 
        End Do; End Do; End Do

        Return
      End Subroutine
!===========================================================================================================================
!****************Å@Diff_EqÅ@************************************************************************************************
!===========================================================================================================================
      Subroutine Diff_Eq_CG

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do k = 1, Nz_in-1, 1
        Do j = 1, Ny_in-1, 1
        Do i = 1, Nx_in-1, 1


          dUf(i,j,k) = ( - ( - Uf_in(i-1,j  ,k  ) + Uf_in(i  ,j  ,k  ) ) * c1_dx_in                           &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i+1,j  ,k  ) ) * c1_dx_in ) * c1_dx_in    * D_in    &
                     + ( - ( - Uf_in(i  ,j-1,k  ) + Uf_in(i  ,j  ,k  ) ) * c1_dy_in                           &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i  ,j+1,k  ) ) * c1_dy_in ) * c1_dy_in    * D_in    &
                     + ( - ( - Uf_in(i  ,j  ,k-1) + Uf_in(i  ,j  ,k  ) ) * c1_dz_in                           &
                         + ( - Uf_in(i  ,j  ,k  ) + Uf_in(i  ,j  ,k+1) ) * c1_dz_in ) * c1_dz_in    * D_in


          dVf(i,j,k) = ( - ( - Vf_in(i-1,j  ,k  ) + Vf_in(i  ,j  ,k  ) ) * c1_dx_in                           &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i+1,j  ,k  ) ) * c1_dx_in ) * c1_dx_in    * D_in    &
                     + ( - ( - Vf_in(i  ,j-1,k  ) + Vf_in(i  ,j  ,k  ) ) * c1_dy_in                           &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i  ,j+1,k  ) ) * c1_dy_in ) * c1_dy_in    * D_in    &
                     + ( - ( - Vf_in(i  ,j  ,k-1) + Vf_in(i  ,j  ,k  ) ) * c1_dz_in                           &
                         + ( - Vf_in(i  ,j  ,k  ) + Vf_in(i  ,j  ,k+1) ) * c1_dz_in ) * c1_dz_in    * D_in


          dWf(i,j,k) = ( - ( - Wf_in(i-1,j  ,k  ) + Wf_in(i  ,j  ,k  ) ) * c1_dx_in                           &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i+1,j  ,k  ) ) * c1_dx_in ) * c1_dx_in    * D_in    &
                     + ( - ( - Wf_in(i  ,j-1,k  ) + Wf_in(i  ,j  ,k  ) ) * c1_dy_in                           &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i  ,j+1,k  ) ) * c1_dy_in ) * c1_dy_in    * D_in    &
                     + ( - ( - Wf_in(i  ,j  ,k-1) + Wf_in(i  ,j  ,k  ) ) * c1_dz_in                           &
                         + ( - Wf_in(i  ,j  ,k  ) + Wf_in(i  ,j  ,k+1) ) * c1_dz_in ) * c1_dz_in    * D_in




        End Do; End Do; End Do

        Return
      End Subroutine
!===============================================================================================================================
!========= [ Subroutine write_UVWf_in_CG ] =====================================================================================
      Subroutine Write_UVWf_in_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k


        !--- MPI ---
        Integer js,je,ks,ke
        !-----------


        Open(40, File = 'output/UVW_f_in', Form = 'unFormatted', Status = 'unknown')
        Write(40) Uf_in
        Write(40) Vf_in
        Write(40) Wf_in
        Close(40)

  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!===============================================================================================================================
!========= [ Subroutine write_UVWf_in ] ========================================================================================
      Subroutine Read_UVWf_in_CG

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k


        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        If( my_rank==root ) then
       ! Open(CONVERT='big_endian', File = 'UVW_f_in', Form = 'unFormatted', Status = 'unknown', UNIT=40)
          Open(40, File = 'UVW_f_in', Form = 'unFormatted', Status = 'unknown')
          Read(40) Uf_in
          Read(40) Vf_in
          Read(40) Wf_in
          Close(40)
        End If


        Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Uf_in )
        Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Vf_in )
        Call mp_bcast_r8_32( -4, Nx_in+4,  -4, Ny_in+4,  -4, Nz_in+4 , Wf_in )

  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine
!===============================================================================================================================
!========= [ Subroutine Interpolate_DNS_Inflow ] ===============================================================================
      Subroutine Interpolate_DNS_Inflow

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k ! ïœêî

        Integer x_int, y_int, z_int
        Real*4 ax, by, cz    !ï‚ä‘éÆåWêî

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end


        If( Ityp_Inter==0 ) then

          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1, Nx-1

            U(i,j,k) = Uf_in(i,j,k)
            V(i,j,k) = Vf_in(i,j,k)
            W(i,j,k) = Wf_in(i,j,k)

          End Do; End Do; End Do

        End If

        If( Ityp_Inter==1 ) then

          Do k = ks, ke, 1
          Do j = js, je, 1
          Do i = 1, Nx-1
            x_int = int( x(i) / dx_in + 1.0d0 ) + 1
            y_int = int( y(j) / dy_in + 1.0d0 ) + 1
            z_int = int( z(k) / dz_in + 1.0d0 ) + 1


            ax = ( x_in(x_int)  - (x(i)) ) / dx_in
            by = ( y_in(y_int)  - (y(j)) ) / dy_in
            cz = ( z_in(z_int)  - (z(k)) ) / dz_in


            U(i,j,k) =             cz * (      by  * (            ax  * Uf_in(x_int-1, y_int-1,z_int-1)       &
                                                       + (1.0d0 - ax) * Uf_in(x_int  , y_int-1,z_int-1)  )    &
                                      + (1.0d0-by) * (            ax  * Uf_in(x_int-1, y_int  ,z_int-1)       &
                                                       + (1.0d0 - ax) * Uf_in(x_int  , y_int  ,z_int-1)  )  ) &
                     +    (1.0d0 - cz)* (      by  * (            ax  * Uf_in(x_int-1, y_int-1,z_int  )       &
                                                       + (1.0d0 - ax) * Uf_in(x_int  , y_int-1,z_int  )  )    &
                                      + (1.0d0-by) * (            ax  * Uf_in(x_int-1, y_int  ,z_int  )       &
                                                       + (1.0d0 - ax) * Uf_in(x_int  , y_int  ,z_int  )  )  )


            V(i,j,k) =             cz * (      by  * (            ax  * Vf_in(x_int-1, y_int-1,z_int-1)       &
                                                       + (1.0d0 - ax) * Vf_in(x_int  , y_int-1,z_int-1)  )    &
                                      + (1.0d0-by) * (            ax  * Vf_in(x_int-1, y_int  ,z_int-1)       &
                                                       + (1.0d0 - ax) * Vf_in(x_int  , y_int  ,z_int-1)  )  ) &
                     +    (1.0d0 - cz)* (      by  * (            ax  * Vf_in(x_int-1, y_int-1,z_int  )       &
                                                       + (1.0d0 - ax) * Vf_in(x_int  , y_int-1,z_int  )  )    &
                                      + (1.0d0-by) * (            ax  * Vf_in(x_int-1, y_int  ,z_int  )       &
                                                       + (1.0d0 - ax) * Vf_in(x_int  , y_int  ,z_int  )  )  )


            W(i,j,k) =             cz * (      by  * (            ax  * Wf_in(x_int-1, y_int-1,z_int-1)       &
                                                       + (1.0d0 - ax) * Wf_in(x_int  , y_int-1,z_int-1)  )    &
                                      + (1.0d0-by) * (            ax  * Wf_in(x_int-1, y_int  ,z_int-1)       &
                                                       + (1.0d0 - ax) * Wf_in(x_int  , y_int  ,z_int-1)  )  ) &
                     +    (1.0d0 - cz)* (      by  * (            ax  * Wf_in(x_int-1, y_int-1,z_int  )       &
                                                       + (1.0d0 - ax) * Wf_in(x_int  , y_int-1,z_int  )  )    &
                                      + (1.0d0-by) * (            ax  * Wf_in(x_int-1, y_int  ,z_int  )       &
                                                       + (1.0d0 - ax) * Wf_in(x_int  , y_int  ,z_int  )  )  )

          End Do; End Do; End Do

        End If

    !    If( Ityp_inS==1 ) then
    !
    !      Do k = ks, ke, 1
    !      Do j = js, je, 1
    !      Do i = 1, Nx-1
    !
    !        U(i,j,k) = Uf_in(i,j,k)
    !        V(i,j,k) = Vf_in(i,j,k)
    !        W(i,j,k) = Wf_in(i,j,k)
    !
    !      End Do; End Do; End Do
    !
    !    End If

        !==========================================================
        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
        End Do; End Do

        Do k = ks, ke, 1
          If(k==1) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
               U_xy (i,j) =    U(i,j,k)
            End Do; End Do
          End If

          If(k==(Nz-1)) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
               V_xy (i,j) =    V(i,j,k)
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then
          Open(201, File = 'output/Uin_xy.txt', Status = 'unknown')
          Do j = 1, Ny-1, 1
          Do i = 1, Nx-1, 1
            Write(201,777)   x(i), y(j),    &
                              U_xy (i,j),   &
                              V_xy (i,j)
          End Do; End Do
          Close(201)
        End If



        !==========================================================
        Do k = 1, Nz-1, 1
        Do i = 1, Nx-1, 1
           U_xz (i,k) = 0.0d0
           V_xz (i,k) = 0.0d0
        End Do; End Do

        Do j = js, je, 1
          If(j==(Ny-1)) then
            Do k = ks, ke, 1
            Do i = 1, Nx-1, 1
               U_xz (i,k) =    U(i,j,k)
               V_xz (i,k) =    V(i,j,k)
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    U_xz(1:Nx-1,1:Nz-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Nz-1,    V_xz(1:Nx-1,1:Nz-1) )

        If(my_rank==root) then
          Open(201, File = 'output/Uin_xz.txt', Status = 'unknown')
          Do k = 1, Nz-1, 1
          Do i = 1, Nx-1, 1
            Write(201,777)   x(i), z(k),    &
                              U_xz (i,k),   &
                              V_xz (i,k)
          End Do; End Do
          Close(201)
        End If

        !==========================================================


  777   Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!=============================================================================================================================
!****************Å@Eular_1stÅ@************************************************************************************************
!=============================================================================================================================
      Subroutine Calc_initial_pressure_fluctuations

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k         ! ïœêî
        Integer n_HD            ! ïœêî
        Real*4 :: Rh_sqrt                                            !ñßìx
        Real*4 :: Usd,Vsd,Wsd                                        !ñßìx
        Real*4, dimension(:,:,:),allocatable :: Potential,Source     !ñßìx

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        allocate(    Potential(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(       Source(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- ä÷êîÇÃèÄîı ---
        !---------------------------------------------

        !---É”ÇÃèÄîı---
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Potential (i,j,k) = 0.0d0
         End Do; End Do; End Do

        !---source term(âEï”)ÇÃèÄîı---

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = 0.0d0
         End Do; End Do; End Do

        !--- i=1,j=1 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = U(i,j,k)*U(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=1,j=2 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = U(i,j,k)*V(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=1,j=3 ---
        !---------------

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = U(i,j,k)*W(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------


        !--- i=2,j=1 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = V(i,j,k)*U(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=2,j=2 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = V(i,j,k)*V(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=2,j=3 ---
        !---------------

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = V(i,j,k)*W(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------


        !--- i=3,j=1 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = W(i,j,k)*U(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 1  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=3,j=2 ---
        !---------------
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = W(i,j,k)*V(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 2  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !--- i=3,j=3 ---
        !---------------

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = W(i,j,k)*W(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P    )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P(i,j,k) = dPdx(i,j,k)
         End Do; End Do; End Do
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, P )
         Call Sub_Derivative_MPI( 3  ,   P  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dPdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = Source(i,j,k) + dPdx(i,j,k)
         End Do; End Do; End Do
        !---------------

        !======== Poissonï˚íˆéÆ  Åﬁ^2É” = -g (É”:P, g:source)  ÇÃâ =============================

         Call Poisson_EQ_CG(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Potential )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           P (i,j,k) = P_bg + Potential(i,j,k)
           Tr(i,j,k) = P (i,j,k)
         End Do; End Do; End Do

        deallocate( Potential )
        deallocate( Source    )

        Return
      End Subroutine

!=============================================================================================================================
!****************Å@Eular_1stÅ@************************************************************************************************
!=============================================================================================================================
      Subroutine HD_init

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k         ! ïœêî
        Integer n_HD            ! ïœêî
        Real*4 :: Rh_sqrt                                            !ñßìx
        Real*4 :: Usd,Vsd,Wsd                                        !ñßìx
        Real*4, dimension(:,:,:),allocatable :: Potential,Source     !ñßìx

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        allocate(    Potential(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )
        allocate(       Source(-4:Nx+4,j_sta-5:j_end+5,k_sta-5:k_end+5)    )

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !--- ä÷êîÇÃèÄîı ---
        !---------------------------------------------

        !---É”ÇÃèÄîı---
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Potential (i,j,k) = 0.0d0
         End Do; End Do; End Do

        !---source term(âEï”)ÇÃèÄîı---

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           U(i,j,k) = U(i,j,k)      !U(i,j,k)=Å„Éœ*U
           V(i,j,k) = V(i,j,k)      !V(i,j,k)=Å„Éœ*V
           W(i,j,k) = W(i,j,k)      !W(i,j,k)=Å„Éœ*W
         End Do; End Do; End Do

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, U )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, V )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, W )

         Call Sub_Derivative_MPI( 1  ,   U  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dUdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Call Sub_Derivative_MPI( 2  ,   V  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dVdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )
         Call Sub_Derivative_MPI( 3  ,   W  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5, &
                                        dWdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5  )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Source(i,j,k) = dUdx(i,j,k) + dVdy(i,j,k) + dWdz(i,j,k)
         End Do; End Do; End Do

        !Call Poisson_SOR(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
        Call Poisson_EQ_CG(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

        !--- ë¨ìxÇÃî≠éUê¨ï™åvéZ ---
        !---------------------------------------------

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Potential )

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Ud(i,j,k) = - ( Potential(i+1,j  ,k  ) - Potential(i-1,j  ,k  ) )/(2.0d0*dx)   !Ud(i,j,k)=(Å„Éœ*U)d
           Vd(i,j,k) = - ( Potential(i  ,j+1,k  ) - Potential(i  ,j-1,k  ) )/(2.0d0*dy)   !Vd(i,j,k)=(Å„Éœ*V)d
           Wd(i,j,k) = - ( Potential(i  ,j  ,k+1) - Potential(i  ,j  ,k-1) )/(2.0d0*dz)   !Wd(i,j,k)=(Å„Éœ*W)d
           Us(i,j,k) = U(i,j,k) - Ud(i,j,k)                                               !Us(i,j,k)=(Å„Éœ*U)s
           Vs(i,j,k) = V(i,j,k) - Vd(i,j,k)                                               !Vs(i,j,k)=(Å„Éœ*V)s
           Ws(i,j,k) = W(i,j,k) - Wd(i,j,k)                                               !Ws(i,j,k)=(Å„Éœ*W)s
         End Do; End Do; End Do

        !--- ã´äEÇÃílåvéZ ---
        !---------------------------------------------

          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ud )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Wd )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Us )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vs )
          Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ws )

        !---------------------------------------------

        !---É”sÇÃèÄîı---
         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           Potential (i,j,k) = 0.0d0
         End Do; End Do; End Do

         Do n_HD = 1, 10

           !---source term(âEï”)ÇÃèÄîı---

            Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Us )
            Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Vs )
            Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Ws )

            Call Sub_Derivative_MPI( 1  ,   Us  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                           dUsdx, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            Call Sub_Derivative_MPI( 2  ,   Vs  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                           dVsdy, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )
            Call Sub_Derivative_MPI( 3  ,   Ws  , -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5,  &
                                           dWsdz, -4,Nx+4,j_sta-5,j_end+5, k_sta-5,k_end+5   )

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
              Source    (i,j,k) = dUsdx(i,j,k) + dVsdy(i,j,k) + dWsdz(i,j,k)
            End Do; End Do; End Do

           !Call Poisson_SOR(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)
           Call Poisson_EQ_CG(Potential,Source,-4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5)

           !--- ë¨ìxÇÃî≠éUê¨ï™åvéZ ---
           !---------------------------------------------

            Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, Potential )

            Do k = ks, ke, 1
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
              Usd =  - ( Potential(i+1,j  ,k  ) - Potential(i-1,j  ,k  ) )/(2.0d0*dx)
              Vsd =  - ( Potential(i  ,j+1,k  ) - Potential(i  ,j-1,k  ) )/(2.0d0*dy)
              Wsd =  - ( Potential(i  ,j  ,k+1) - Potential(i  ,j  ,k-1) )/(2.0d0*dz)
              Ud(i,j,k) = Ud(i,j,k) + Usd
              Vd(i,j,k) = Vd(i,j,k) + Vsd
              Wd(i,j,k) = Wd(i,j,k) + Wsd
              Us(i,j,k) = Us(i,j,k) - Usd
              Vs(i,j,k) = Vs(i,j,k) - Vsd
              Ws(i,j,k) = Ws(i,j,k) - Wsd
            End Do; End Do; End Do

         End Do
         !---------------------------------------------

         !call write_xy_init

         Do k = ks, ke, 1
         Do j = js, je, 1
         Do i = 1, Nx-1, 1
           U(i,j,k) = Us(i,j,k)             !U(i,j,k)=U
           V(i,j,k) = Vs(i,j,k)             !V(i,j,k)=V
           W(i,j,k) = Ws(i,j,k)             !W(i,j,k)=W
         End Do; End Do; End Do

         !--- ã´äEÇÃílåvéZ ---
         !---------------------------------------------

         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, U )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, V )
         Call BC_Periodic_xyz( -4,Nx+4,j_sta-5,j_end+5,k_sta-5,k_end+5, W )

         deallocate( Potential )
         deallocate( Source    )

        Return
      End Subroutine

!=============================================================================================================================
!=============================================================================================================================
!****************Å@Eular_1stÅ@************************************************************************************************
!=============================================================================================================================
      Subroutine Make_Dispersion_one

        Use module_mpi
        Use Definition
        Implicit none
        Integer i, j, k         ! ïœêî
        Integer n_HD            ! ïœêî

        !--- MPI ---
        Integer js,je,ks,ke
        !-----------

        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        !===ïΩãœ===

         UsumT_xyz  = 0.0d0
         VsumT_xyz  = 0.0d0
         WsumT_xyz  = 0.0d0
         UUsumT_xyz = 0.0d0
         VVsumT_xyz = 0.0d0
         WWsumT_xyz = 0.0d0

         Do j = js, je
           Do k = ks, ke
           Do i = 1, Nx-1
            UsumT_xyz  = UsumT_xyz  + (   U(i,j,k)            ) * c1_NxNyNz
            VsumT_xyz  = VsumT_xyz  + (   V(i,j,k)            ) * c1_NxNyNz
            WsumT_xyz  = WsumT_xyz  + (   W(i,j,k)            ) * c1_NxNyNz
            UUsumT_xyz = UUsumT_xyz + (   U(i,j,k) * U(i,j,k) ) * c1_NxNyNz
            VVsumT_xyz = VVsumT_xyz + (   V(i,j,k) * V(i,j,k) ) * c1_NxNyNz
            WWsumT_xyz = WWsumT_xyz + (   W(i,j,k) * W(i,j,k) ) * c1_NxNyNz
           End Do; End Do
         End Do

        !---/* MPI */---

        !---------------------------------------------------
        Call mp_allsumr(   UsumT_xyz   )
        Call mp_allsumr(   VsumT_xyz   )
        Call mp_allsumr(   WsumT_xyz   )
        Call mp_allsumr(   UUsumT_xyz  )
        Call mp_allsumr(   VVsumT_xyz  )
        Call mp_allsumr(   WWsumT_xyz  )


        Do k = ks, ke
        Do j = js, je
        Do i = 1, Nx-1
         U(i,j,k) = ( U(i,j,k) - UsumT_xyz )/sqrt( UUsumT_xyz - UsumT_xyz*UsumT_xyz )
         V(i,j,k) = ( V(i,j,k) - VsumT_xyz )/sqrt( VVsumT_xyz - VsumT_xyz*VsumT_xyz )
         W(i,j,k) = ( W(i,j,k) - WsumT_xyz )/sqrt( WWsumT_xyz - WsumT_xyz*WsumT_xyz )
        End Do; End Do; End Do

        Return
      End Subroutine
!=============================================================================================================================
!========= [ Subroutine write_xy ] ===========================================================================================
      Subroutine write_xy_init

        Use module_mpi
        Use Definition
        Implicit none   
        Integer i, j, k
        !--- MPI ---
        Integer js,je,ks,ke
        !-----------
        Integer fn0, fn1, fn2, fn3, fn4, fn5

        Integer i_xyz_min, j_xyz_min, k_xyz_min
        Integer i_xyz_max, j_xyz_max, k_xyz_max
        Integer isize,jsize,ksize

        Integer reclength_inst
        Integer npart
        Real*4 spacex, spacey, spacez, Ox, Oy, Oz

        character(len=80)::binary_form
        character(len=80)::file_description1,file_description2
        character(len=80)::node_id,element_id
        character(len=80)::part,description_part,block
        character(len=80):: desc

        fn5 =   nt       / 100000
        fn4 = ( nt - fn5 * 100000     ) / 10000
        fn3 = ( nt - fn5 * 100000 - fn4 * 10000     ) / 1000
        fn2 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000     ) / 100
        fn1 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100     ) / 10
        fn0 = ( nt - fn5 * 100000 - fn4 * 10000 - fn3 * 1000 - fn2 * 100 - fn1 * 10 ) / 1



        js=j_sta
        je=j_end
        ks=k_sta
        ke=k_end

        Do j = 1, Ny-1, 1
        Do i = 1, Nx-1, 1
          Rh_xy (i,j) = 0.0d0
           U_xy (i,j) = 0.0d0
           V_xy (i,j) = 0.0d0
           W_xy (i,j) = 0.0d0
          Te_xy (i,j) = 0.0d0
           P_xy (i,j) = 0.0d0
        divU_xy (i,j) = 0.0d0
          Ys_xy (i,j) = 0.0d0
        Enst_xy (i,j) = 0.0d0
          Ma_xy (i,j) = 0.0d0

        End Do; End Do

        Do k = ks, ke, 1
          If(k==Nz-1) then
            Do j = js, je, 1
            Do i = 1, Nx-1, 1
              Rh_xy (i,j) =    U(i,j,k)
               U_xy (i,j) =    V(i,j,k)
               V_xy (i,j) =    W(i,j,k)
               W_xy (i,j) =   Us(i,j,k)
              Te_xy (i,j) =   Vs(i,j,k)
               P_xy (i,j) =   Ws(i,j,k)
            divU_xy (i,j) =    U(i,j,k)-Us(i,j,k)
              Ys_xy (i,j) =    V(i,j,k)-Vs(i,j,k)
            Enst_xy (i,j) =    W(i,j,k)-Ws(i,j,k)
              Ma_xy (i,j) =  sqrt(U(i,j,k)**2.0d0 + V(i,j,k)**2.0d0 + W(i,j,k)**2.0d0)/Sqrt(gh*T(i,j,k)*cRT_UU)
            End Do; End Do
          End If
        End Do

        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Rh_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    U_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    V_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    W_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Te_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,    P_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, divU_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ys_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1, Enst_xy(1:Nx-1,1:Ny-1) )
        Call mp_sum_r22( 1, Nx-1, 1, Ny-1,   Ma_xy(1:Nx-1,1:Ny-1) )

        If(my_rank==root) then

          fname_geo = 'output/HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          fname_case = 'output/HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.Case'

          i_xyz_min =  1   ;  i_xyz_max =  Nx-1
          j_xyz_min =  1   ;  j_xyz_max =  Ny-1
          k_xyz_min =  Nz-1;  k_xyz_max =  Nz-1


          !--- Paraview in Ensight Gold Format ---
          binary_form      ='Fortran Binary'
          file_description1='Ensight Model Geometry File Created by '
          file_description2='WriteRectEnsightGeo Routine'
          node_id          ='node id off'
          element_id       ='element id off'
          part             ='part'
          npart            =1
          description_part ='2D DATA'
          block            ='block rectilinear'
          isize=i_xyz_max-i_xyz_min+1
          jsize=j_xyz_max-j_xyz_min+1
          ksize=k_xyz_max-k_xyz_min+1


          !=== Geometry ===
          Open (23, File = fname_geo,form='UNFORMATTED') 
          Write(23) binary_form                     ! 80
          Write(23) file_description1               ! 80
          Write(23) file_description2               ! 80
          Write(23) node_id                         ! 80
          Write(23) element_id                      ! 80
          Write(23) part                            ! 80
          Write(23) npart                           ! 4
          Write(23) description_part                ! 80
          Write(23) block                           ! 80
          Write(23) isize,jsize,ksize               ! 4*3
          Write(23) (Real(x (i)),i=i_xyz_min ,i_xyz_max )   ! 4*isize
          Write(23) (Real(y (j)),j=j_xyz_min ,j_xyz_max )   ! 4*jsize
          Write(23) (Real(z (k)),k=k_xyz_min ,k_xyz_max )   ! 4*ksize
          Close(23)

          fname_geo = 'HIT_xy_'&
          //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.geo'

          !--- Case File ---
          Open (50, File = fname_case, Status = 'unknown' ) 
          Write(50,'(a)') 'FORMAT'
          Write(50,'(a)') 'type:   ensight gold'
          Write(50,'(a)') ' '
          Write(50,'(a)') 'GEOMETRY'
          Write(50,'(a,a)') 'model:  ', trim(fname_geo)
          Write(50,'(a)') ' '
          Write(50,'(a)') 'VARIABLE'


            !--- Inst Data ---
            desc = 'U'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Rh_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'V'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(U_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'W'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(V_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Us'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(W_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Vs'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Te_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Ws'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(P_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Ud'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(divU_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Vd'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart 
            Write(23) block
            Write(23) ((( Real(Ys_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)


            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)

            !--- Inst Data ---
            desc = 'Wd'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) ((( Real(Enst_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)

            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


            !--- Inst Data ---
            desc = 'Ma'
            fname_scl = 'output/HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Open (23, File = fname_scl,form='UNFORMATTED') 
            Write(23) desc
            Write(23) part
            Write(23) npart
            Write(23) block
            Write(23) ((( Real(Ma_xy(i,j)), i=i_xyz_min,i_xyz_max), j=j_xyz_min,j_xyz_max), k=k_xyz_min,k_xyz_max)
            Close(23)

            fname_scl = 'HIT_'//trim(desc)//'_xy_'&
             //char(fn5+48)//char(fn4+48)//char(fn3+48)//char(fn2+48)//char(fn1+48)//char(fn0+48)//'step.bin'
            Write(50,'(a,a,a,a)') 'scalar per node:  ', trim(desc), ' ', trim(fname_scl)


          Close(50)




        End If

    777 Format(E13.6, 200(1x, E13.6)) 

        Return
      End Subroutine

!=================================================================================================================================
!=====[Subroutine IFFT]===========================================================================================================
!=================================================================================================================================
      Subroutine IFFT_initial
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none

         allocate( Fr(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4) )
         allocate( Fi(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4) )
         allocate( a ( 2*(Nx_in+4) ) )
         allocate( theta1(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4) )
         allocate( theta2(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4) )
         allocate( theta3(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4) )

         n1 = Nx_in-1
         n2 = Ny_in-1
         n3 = Nz_in-1

         Call initial_condition_johnsen
         Call initial_random_MT

!========== [U direction]==========================
         call U_IFFT_Johnsen
!========== [V direction]==========================
         call V_IFFT_Johnsen
!========== [W direction]==========================
         call W_IFFT_Johnsen

         deallocate(Fr,Fi)
         deallocate(a)
         deallocate(theta1,theta2,theta3)

         return
      End Subroutine IFFT_initial

!=====[Subroutine initila condition johnsen]===================================================================================
      Subroutine initial_condition_Johnsen
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none

         If(Ityp_Spec == 0) then
           relam0f = Re_ini/L_lambda
           umt0f   = Ma_ini
           urms0f  = rms_in
           lam0f   = L_in/L_lambda

           rel0f   = (relam0f**2)*3.0d0/20.0d0               !R_lambda = (20/3*Re_L)^0.5
           el0f    = lam0f/(dsqrt(10.0d0)*rel0f**(-0.5d0))   !lambda/L = sqrt(10)*Re_L^(-1/2)
           eta0f   = el0f*rel0f**(-0.75d0)                   !eta/L    = Re_L^(-3/4)
           akre0f  = 2.0d0/lam0f

           Asp=coasp*urms0f**2.0d0*16.0d0*dsqrt(2.0d0/pi)/akre0f**5.0d0

           Open (80,File = 'output/HIT_Para_IFFT.txt' , Status = 'unknown' ) 
           Write(80,*) 'ReL    = ', rel0f 
           Write(80,*) 'L      = ', el0f  
           Write(80,*) 'eta0f  = ', eta0f 
           Write(80,*) 'Urms0  = ', urms0f
           Write(80,*) 'k0     = ', akre0f
           Write(80,*) 'A      = ', Asp   
         End If

         If(Ityp_Spec == 1) then

           K_P        = 0.5d0*(3.0d0*rms_in*rms_in)
           Re_lambdaP = Re_ini
           L_P        = L_in
           !Re_LP = Re_ini

           !C_EP     = E_target*L_in/(rms_in*rms_in*rms_in)
           !Re_LP    = Re_ini*Re_ini*C_EP/15.0d0
           !lambda_P = L_in*( sqrt(15.0d0/C_EP)      * Re_LP**(-1.0d0/2.0d0) )
           !eta_P    = L_in*( (C_EP**(-1.0d0/4.0d0)) * Re_LP**(-3.0d0/4.0d0) )
           !C_EP     = E_target*L_in/(rms_in*rms_in*rms_in)

           !Re_lambdaP = sqrt(20.0d0/3.0d0*Re_LP)
           Re_LP      = (3.0d0/20.0d0)*Re_lambdaP*Re_lambdaP
           lambda_P   = L_P*(sqrt(10.0d0) * Re_LP**(-1.0d0/2.0d0) )
           eta_P      = L_P*(               Re_LP**(-3.0d0/4.0d0) )
          ! E_P        = 15.0d0*(rms_in*lambda_P/Re_lambdaP)*(rms_in*rms_in)/(lambda_P*lambda_P)/Re

           E_P        = 15.0d0*(U_ref*rms_in*U_ref*rms_in)/(L_ref*lambda_P*L_ref*lambda_P) &
                         *(P_inf/R_gas/T_inf)                                              &
                          /(  mu_T0 * ( (T_inf/273.0d0 )**(1.5d0)       )                  &
                                 * (273.0d0+110.4d0)/(T_inf+110.4d0) )  
           E_P        = E_P / (U_ref*U_ref*U_ref/L_ref)

           !eta_P      = L_in*eta_P   /L_P
           !lambda_P   = L_in*lambda_P/L_P
           !L_P        = L_in*L_P     /L_P

           Open (80,File = 'output/HIT_Para_IFFT.txt' , Status = 'unknown' ) 
           Write(80,*) 'Kinetic energy = ', K_P        
           Write(80,*) 'Re_lambda      = ', Re_lambdaP 
           Write(80,*) 'lambda         = ', lambda_P   
           Write(80,*) 'eta            = ', eta_P      
           Write(80,*) 'Epsilon        = ', E_P        
         End If

         Close(80)

!////////////////////////////////////////

         return
      End subroutine initial_condition_Johnsen

!=====[Subroutine U_IFFT_johnsen]==================================================================================================
      Subroutine U_IFFT_Johnsen
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none

        Integer i,j,k
        Real*8 L_X,L_Y,L_Z

         L_X = Lx!*L_ini
         L_Y = Ly!*L_ini
         L_Z = Lz!*L_ini

         urms=0.0d0

         do k=0,n3
         do j=0,n2
         do i=0,n1

            if( (i == n1/2) .and. (j == n2/2) .and. (k == n3/2) ) then
              Fr(i,j,k) = 0.0d0
              Fi(i,j,k) = 0.0d0
            else
              nk =   (i-n1/2)*(i-n1/2)/L_X/L_X   &
                   + (j-n2/2)*(j-n2/2)/L_Y/L_Y   &
                   + (k-n3/2)*(k-n3/2)/L_Z/L_Z    
              akk=dble(nk)
              akk=dsqrt(akk)

              if( (i == n1/2) .and. (j == n2/2) ) then
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 =  0.0d0
               ak212 =  1.0d0
              else
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 = dble(i-n1/2)/L_X/ak12
               ak212 = dble(j-n2/2)/L_Y/ak12
              end if

              ak3k  = dble(k-n3/2)/L_Z/akk
              ak12k = ak12/akk
              akre  = akk*2.d0*pi
              !akre  = akre*L_ini
!write(*,*) akre
              If(Ityp_Spec == 0) then
                 eek1 = Asp*(akre**4.d0)*dexp(-2.d0*(akre/akre0f)**2.d0)
              End If
              If(Ityp_Spec == 1) then
                 f_L   = ( akre*L_P/(   (akre*  L_P)**2.0d0 + c_L          )**(1.0d0/2.0d0) )**(5.0d0/3.0d0+p_0)
                 f_eta = exp( -beta*( ( (akre*eta_P)**4.0d0 + c_eta**4.0d0 )**(1.0d0/4.0d0) - c_eta ) )
                 eek1  = K_C*(E_P**(2.0d0/3.0d0))*(akre**(-5.0d0/3.0d0))*f_L*f_eta
              End If
              eek  = dsqrt(2.d0*eek1/(4.d0*pi*akre*akre))
!write(*,*) i,j,k,akre,eek

              are  = eek*dcos(theta1(i,j,k))*dcos(theta3(i,j,k))
              aim  = eek*dsin(theta1(i,j,k))*dcos(theta3(i,j,k))
              bre  = eek*dcos(theta2(i,j,k))*dsin(theta3(i,j,k))
              bim  = eek*dsin(theta2(i,j,k))*dsin(theta3(i,j,k))

              Fr(i,j,k) = ak212*are+ak112*ak3k*bre    !Re(u1(k))
              Fi(i,j,k) = ak212*aim+ak112*ak3k*bim    !Im(u1(k))

              urms=urms+(Fr(i,j,k)**2.d0+Fi(i,j,k)**2.d0)*(2.d0*pi)**3.d0/Lx/Ly/Lz/2.d0

           End If

         End Do; End Do; End Do

         write(*,*)'dumurms=',dsqrt(urms)

         call IFFT_3D

         do k=0,n3
         do j=0,n2
         do i=0,n1
            Uf_in(i,j,k)=Fr(i,j,k)
         end do; end do; end do

        return
     End Subroutine U_IFFT_johnsen


!////////////////////////////////////////////////////

!=====[Subroutine V_IFFT_johnsen]==================================================================================================
     Subroutine V_IFFT_Johnsen
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none
        Integer i,j,k
        Real*8 L_X,L_Y,L_Z

         L_X = Lx!*L_ini
         L_Y = Ly!*L_ini
         L_Z = Lz!*L_ini


         do k=0,n3
         do j=0,n2
         do i=0,n1

            if( (i == n1/2) .and. (j == n2/2) .and. (k == n3/2) ) then
              Fr(i,j,k) = 0.0d0
              Fi(i,j,k) = 0.0d0

            else
              nk =   (i-n1/2)*(i-n1/2)/L_X/L_X   &
                   + (j-n2/2)*(j-n2/2)/L_Y/L_Y   &
                   + (k-n3/2)*(k-n3/2)/L_Z/L_Z    
              akk=dble(nk)
              akk=dsqrt(akk)

              if( (i == n1/2) .and. (j == n2/2) ) then
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 =  0.0d0
               ak212 =  1.0d0
              else
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 = dble(i-n1/2)/L_X/ak12
               ak212 = dble(j-n2/2)/L_Y/ak12
              end if

              ak3k  = dble(k-n3/2)/L_Z/akk
              ak12k = ak12/akk
              akre  = akk*2.d0*pi
              !akre  = akre*L_ini

              If(Ityp_Spec == 0) then
                 eek1 = Asp*(akre**4.d0)*dexp(-2.d0*(akre/akre0f)**2.d0)
              End If
              If(Ityp_Spec == 1) then
                 f_L   = ( akre*L_P/(   (akre*  L_P)**2.0d0 + c_L          )**(1.0d0/2.0d0) )**(5.0d0/3.0d0+p_0)
                 f_eta = exp( -beta*( ( (akre*eta_P)**4.0d0 + c_eta**4.0d0 )**(1.0d0/4.0d0) - c_eta ) )
                 eek1  = K_C*(E_P**(2.0d0/3.0d0))*(akre**(-5.0d0/3.0d0))*f_L*f_eta
              End If
              eek = dsqrt(2.d0*eek1/(4.d0*pi*akre*akre))

              are = eek*dcos(theta1(i,j,k))*dcos(theta3(i,j,k))
              aim = eek*dsin(theta1(i,j,k))*dcos(theta3(i,j,k))
              bre = eek*dcos(theta2(i,j,k))*dsin(theta3(i,j,k))
              bim = eek*dsin(theta2(i,j,k))*dsin(theta3(i,j,k))

              Fr(i,j,k) = ak212*ak3k*bre-ak112*are    !Re(u2(k))
              Fi(i,j,k) = ak212*ak3k*bim-ak112*aim    !Im(u2(k))

           End If

         End Do; End Do; End Do

         call IFFT_3D

         do k=0,n3
         do j=0,n2
         do i=0,n1
            Vf_in(i,j,k)=Fr(i,j,k)
         end do; end do; end do

         return
     End Subroutine V_IFFT_Johnsen


!=====[Subroutine W_IFFT_Johnsen]==================================================================================================
     Subroutine W_IFFT_Johnsen
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none
        Integer i,j,k
        Real*8 L_X,L_Y,L_Z

         L_X = Lx!*L_ini
         L_Y = Ly!*L_ini
         L_Z = Lz!*L_ini


         do k=1,n3
         do j=1,n2
         do i=1,n1

          if( (i == n1/2) .and. (j == n2/2) .and. (k == n3/2) ) then
            Fr(i,j,k) = 0.0d0
            Fi(i,j,k) = 0.0d0

            else
              nk =   (i-n1/2)*(i-n1/2)/L_X/L_X   &
                   + (j-n2/2)*(j-n2/2)/L_Y/L_Y   &
                   + (k-n3/2)*(k-n3/2)/L_Z/L_Z    
              akk=dble(nk)
              akk=dsqrt(akk)

              if( (i == n1/2) .and. (j == n2/2) ) then
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 =  0.0d0
               ak212 =  1.0d0
              else
               nk    = (i-n1/2)*(i-n1/2)/L_X/L_X+(j-n2/2)*(j-n2/2)/L_Y/L_Y
               ak12  = dble(nk)
               ak12  = dsqrt(ak12)
               ak112 = dble(i-n1/2)/L_X/ak12
               ak212 = dble(j-n2/2)/L_Y/ak12
              end if

              ak3k  = dble(k-n3/2)/L_Z/akk
              ak12k = ak12/akk
              akre  = akk*2.d0*pi
              !akre  = akre*L_ini

              If(Ityp_Spec == 0) then
                 eek1 = Asp*(akre**4.d0)*dexp(-2.d0*(akre/akre0f)**2.d0)
              End If
              If(Ityp_Spec == 1) then
                 f_L   = ( akre*L_P/(   (akre*  L_P)**2.0d0 + c_L          )**(1.0d0/2.0d0) )**(5.0d0/3.0d0+p_0)
                 f_eta = exp( -beta*( ( (akre*eta_P)**4.0d0 + c_eta**4.0d0 )**(1.0d0/4.0d0) - c_eta ) )
                 eek1  = K_C*(E_P**(2.0d0/3.0d0))*(akre**(-5.0d0/3.0d0))*f_L*f_eta
              End If
            eek = dsqrt(2.d0*eek1/(4.d0*pi*akre*akre))

            are = eek*dcos(theta1(i,j,k))*dcos(theta3(i,j,k))
            aim = eek*dsin(theta1(i,j,k))*dcos(theta3(i,j,k))
            bre = eek*dcos(theta2(i,j,k))*dsin(theta3(i,j,k))
            bim = eek*dsin(theta2(i,j,k))*dsin(theta3(i,j,k))

            Fr(i,j,k) = -ak12k*bre                  !Re(u3(k))
            Fi(i,j,k) = -ak12k*bim                  !Im(u3(k))

           End If

         End Do; End Do; End Do

         Call IFFT_3D

         do k=0,n3
         do j=0,n2
         do i=0,n1
            Wf_in(i,j,k) = Fr(i,j,k)
         end do; end do; end do

         return
     End Subroutine W_IFFT_Johnsen

!==========[ Subroutine 3D_IFFT ]==================================================================================================
     Subroutine IFFT_3D
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Integer i,j,k,ii,jj,kk
        integer :: loopba,loopre,loopim

!      ///3rd stage(k-direction IFFT)///
       do i=1,n1
        do j=1,n2
         loopba = 0
         do k=n3/2,n3
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do
         do k=1,n3/2-1
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do

         call dfour1(a,n3,1) !FFT=-1,IFFT=1

         do k=1,n3
          kk = 2*k-1
          Fr(i,j,k) = a(kk)
          Fi(i,j,k)   = a(kk+1)
         end do
         Fr(i,j,0) = Fr(i,j,n3)
         Fi(i,j,0) = Fi(i,j,n3)
        end do
       end do
!      ///3rd stage end///

!     ///2nd stage(j-direction IFFT)///
       do i=1,n1
        do k=1,n3
         loopba = 0
         do j=n2/2,n2
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do
         do j=1,n2/2-1
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do
         call dfour1(a,n2,1) !FFT=-1,IFFT=1
         do j=1,n2
          jj = 2*j-1
          Fr(i,j,k) = a(jj)
          Fi(i,j,k)   = a(jj+1)
         end do
         Fr(i,0,k)= Fr(i,n2,k)
         Fi(i,0,k)= Fi(i,n2,k)
        end do
       end do

!      ///2nd stage end//
!      ///1st stage(i-direction IFFT)///
       do k=1,n3
        do j=1,n2
         loopba = 0
         do i=n1/2,n1
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do
         do i=1,n1/2-1
          loopba = loopba+1
          loopre = loopba*2-1
          loopim = loopba*2
          a(loopre)   = Fr(i,j,k)
          a(loopim)   = Fi(i,j,k)
         end do

         call dfour1(a,n1,1) !FFT=-1,IFFT=1

         do i=1,n1
          ii = 2*i-1
          Fr(i,j,k) = a(ii)
          Fi(i,j,k) = a(ii+1)
         end do
         Fr(0,j,k) =Fr(n1,j,k)
         Fi(0,j,k) =Fi(n1,j,k)
        end do
       end do

!     ///1st stage end///

       do k=0,n3
       do j=0,n2
       do i=0,n1
            Fr(i,j,k)=Fr(i,j,k)*((2.0d0*pi)**3.0d0/(Lx*Ly*Lz))**0.5d0
       end do; end do; end do

       return
     End Subroutine IFFT_3D

!=====[Subroutine FFT]=============================================================================================================
!     *********************************************************
!     *         fft subroutine by numerical recipes           *
!     *********************************************************
      Subroutine dfour1(data_IFFT,nn,isign)

      Implicit none

      Integer :: isign,nn
      Integer :: i,istep,j,m,mmax,n
      Real*8  :: data_IFFT(2*nn)
      Real*8  :: tempi,tempr
      Real*8  :: theta,wi,wpi,wpr,wr,wtemp
      Real*8  :: p

      p = dacos(-1.d0)
      n=2*nn
      j=1
      do i=1,n,2
       if(j.gt.i) then
        tempr=data_IFFT(j)
        tempi=data_IFFT(j+1)
        data_IFFT(j)=data_IFFT(i)
        data_IFFT(j+1)=data_IFFT(i+1)
        data_IFFT(i)=tempr
        data_IFFT(i+1)=tempi
       end if
       m=n/2
       do while ((m.ge.2).and.(j.gt.m))
        j=j-m
        m=m/2
       end do
       j=j+m
      end do
      mmax=2
      do while (n.gt.mmax)
       istep=2*mmax
       theta=isign*(2.d0*p/mmax)
       wtemp=dsin(0.5d0*theta)
       wpr=-2.d0*wtemp*wtemp
       wpi=dsin(theta)
       wr=1.d0
       wi=0.d0
       do m=1,mmax,2
        do i=m,n,istep
         j=i+mmax
         tempr=(wr)*data_IFFT(j)-(wi)*data_IFFT(j+1)
         tempi=(wr)*data_IFFT(j+1)+(wi)*data_IFFT(j)
         data_IFFT(j)=data_IFFT(i)-tempr
         data_IFFT(j+1)=data_IFFT(i+1)-tempi
         data_IFFT(i)=data_IFFT(i)+tempr
         data_IFFT(i+1)=data_IFFT(i+1)+tempi
        end do
        wtemp=wr
        wr=wtemp*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
       end do
       mmax=istep
      end do

      return
      end subroutine dfour1

!=====[initial_random]==========================================================================================================
      Subroutine initial_random_MT
        Use module_mpi
        Use Definition
        Use Definition_IFFT
        Implicit none
        Integer :: i,j,k,i2,j2,k2
        Integer :: seed

!     ---ÉÅÉÇÉäÇÃäÑÇËìñÇƒ---
         allocate(  psi1(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)  )
         allocate(  psi2(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)  )
         allocate(  psi3(-4:Nx_in+4,-4:Ny_in+4,-4:Nz_in+4)  )

!     ---(k1,k2,k3) quadrant---
         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(b)
            psi1(i,j,k)  = 2.d0*pi*b
         end do; end do; end do

         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(c)
            psi2(i,j,k)  = 2.d0*pi*c
         end do; end do; end do

         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(d)
            psi3(i,j,k)  = 2.d0*pi*d
         end do; end do; end do
!      ---(-k1,k2,k3) quadrant---
         do k=n3/2  ,n3
         do j=n2/2  ,n2
         do i=n1/2+1,n1
            i2 = n1-i
            psi1(i2,j,k)  = -psi1(i,j,k)
            psi2(i2,j,k)  =  psi2(i,j,k)
            psi3(i2,j,k)  =  psi3(i,j,k)
         end do; end do; end do
!      ---(k1,-k2,k3) quadrant---
         do k=n3/2  ,n3
         do j=n2/2+1,n2
         do i=     1,n1
            j2 = n2-j
            psi1(i,j2,k)  =  psi1(i,j,k)
            psi2(i,j2,k)  = -psi2(i,j,k)
            psi3(i,j2,k)  =  psi3(i,j,k)
         end do; end do; end do
!      ---(k1,k2,-k3) quadrant---
         do k=n3/2+1,n3
         do j=     1,n2
         do i=     1,n1
            k2 = n3-k
            psi1(i,j,k2)  =  psi1(i,j,k)
            psi2(i,j,k2)  =  psi2(i,j,k)
            psi3(i,j,k2)  = -psi3(i,j,k)
         end do; end do; end do
!
         do k=   0,n3
         do j=   0,n2
         do i=n1/2,n1/2
               psi1(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=   0,n3
         do j=n2/2,n2/2
         do i=   0,n1
               psi2(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=n3/2,n3/2
         do j=   0,n2
         do i=   0,n1
               psi3(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=0,n3
         do j=0,n2
         do i=0,n1
            theta1(i,j,k) = psi1(i,j,k)+psi2(i,j,k)+psi3(i,j,k)+0.5d0*pi
         end do; end do; end do
!
!    /// random value for theta2///
!    ---(k1,k2,k3) quadrant---
         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(e)
            psi1(i,j,k)  = 2.d0*pi*e
         end do; end do; end do

         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(f)
            psi2(i,j,k)  = 2.d0*pi*f
         end do; end do; end do

         do k=n3/2,n3
         do j=n2/2,n2
         do i=n1/2,n1
            call MT_random(g)
            psi3(i,j,k)  = 2.d0*pi*g
         end do; end do; end do

!      ---(-k1,k2,k3) quadrant---
         do k=n3/2  ,n3
         do j=n2/2  ,n2
         do i=n1/2+1,n1
            i2 = n1-i
            psi1(i2,j,k)  = -psi1(i,j,k)
            psi2(i2,j,k)  =  psi2(i,j,k)
            psi3(i2,j,k)  =  psi3(i,j,k)
         end do; end do; end do
!      ---(k1,-k2,k3) quadrant---
         do k=n3/2  ,n3
         do j=n2/2+1,n2
         do i=     0,n1
            j2 = n2-j
            psi1(i,j2,k)  =  psi1(i,j,k)
            psi2(i,j2,k)  = -psi2(i,j,k)
            psi3(i,j2,k)  =  psi3(i,j,k)
         end do; end do; end do
!      ---(k1,k2,-k3) quadrant---
         do k=n3/2+1,n3
         do j=     0,n2
         do i=     0,n1
            k2 = n3-k
            psi1(i,j,k2)  =  psi1(i,j,k)
            psi2(i,j,k2)  =  psi2(i,j,k)
            psi3(i,j,k2)  = -psi3(i,j,k)
         end do; end do; end do
!
         do k=   0,n3
         do j=   0,n2
         do i=n1/2,n1/2
               psi1(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=   0,n3
         do j=n2/2,n2/2
         do i=   0,n1
               psi2(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=n3/2,n3/2
         do j=   0,n2
         do i=   0,n1
               psi3(i,j,k)=0.0d0
         end do; end do; end do
!
         do k=0,n3
         do j=0,n2
         do i=0,n1
            theta2(i,j,k) =  psi1(i,j,k)+psi2(i,j,k)+psi3(i,j,k)
         end do; end do; end do

         do k=0,n3
         do j=0,n2
         do i=0,n1
             call MT_random(h)
             theta3(i,j,k) = 2.d0*pi*h
         end do; end do; end do

         deallocate(psi1);deallocate(psi2);deallocate(psi3)

       return
     End subroutine initial_random_MT

!=====[MT_Random]==========================================================================================================
!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

 module mtmod
! Default seed
    integer :: defaultsd,c
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1

! the array for the state vector
    integer, save, dimension(0:N-1) :: mt
    integer, save                   :: mti = N1
    
! Overload procedures for saving and getting mt state
    interface mtsave
      module procedure mtsavef
      module procedure mtsaveu
    end interface
    interface mtget
      module procedure mtgetf
      module procedure mtgetu
    end interface

 contains

!Initialization subroutine
  subroutine sgrnd(seed)
    implicit none
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
    end do
!
    return
  end subroutine sgrnd

!Random number generator
  real(8) function grnd()
    implicit integer(a-z)

! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                            if sgrnd() has not been called,

!        defaultsd = 43570
        call system_clock(count=c)
        defaultsd = c
        call sgrnd( defaultsd )
!                              a default initial seed is used
      end if

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    end if

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      grnd=dble(y)/(2.0d0**32-1.0d0)
    end if

    return
  end function grnd

!State saving subroutines.
! Usage:  call mtsave( file_name, format_character )
!    or   call mtsave( unit_number, format_character )
! where   format_character = 'u' or 'U' will save in unformatted form, otherwise
!         state information will be written in formatted form.
  subroutine mtsavef( fname, forma )

!NOTE: This subroutine APPENDS to the end of the file "fname".

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED', &
            position='APPEND')
       write(10)mti
       write(10)mt

      case default
       open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED', &
            position='APPEND')
       write(10,*)mti
       write(10,*)mt

    end select
    close(10)

    return
  end subroutine mtsavef

  subroutine mtsaveu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       write(unum)mti
       write(unum)mt

      case default
       write(unum,*)mti
       write(unum,*)mt

      end select

    return
  end subroutine mtsaveu

!State getting subroutines.
! Usage:  call mtget( file_name, format_character )
!    or   call mtget( unit_number, format_character )
! where   format_character = 'u' or 'U' will read in unformatted form, otherwise
!         state information will be read in formatted form.
  subroutine mtgetf( fname, forma )

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
       read(10)mti
       read(10)mt

      case default
       open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
       read(10,*)mti
       read(10,*)mt

    end select
    close(10)

    return
  end subroutine mtgetf

  subroutine mtgetu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       read(unum)mti
       read(unum)mt

      case default
       read(unum,*)mti
       read(unum,*)mt

      end select

    return
  end subroutine mtgetu

 end module mtmod

 subroutine MT_random(a)
! this main() outputs first 1000 generated numbers
!
    use mtmod
    implicit none

    double precision :: a
    integer, parameter      :: no=10
    real(8), dimension(0:7) :: r
    integer j,k
!    real(8) grnd
!
!      call sgrnd(4357)
!                         any nonzero integer can be used as a seed

    do j=0,no-1
      r(mod(j,8))=grnd()
      a=r(mod(j,8))
    end do

 end subroutine MT_random


!!=============================================================================================================================
!!=====[ BC_Periodic_z_dir_MPI_R4 ]============================================================================================
!!=============================================================================================================================
!      Subroutine BC_Periodic_z_dir_MPI_R4_init( n11 , n12 , n21 , n22 , n31 , n32 , F1 )
!
!        Use module_mpi
!        Use Definition
!        Implicit none
!        Integer,intent(in) :: n11 , n21 , n31
!        Integer,intent(in) :: n12 , n22 , n32
!        Real*4,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)
!
!        Integer i, j, k ! ïœêî
!
!        !--- MPI ---
!        Integer js,je,ks,ke
!	integer ierror,ierr
!        integer dum_len
!        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
!        !-----------
!
!        js=j_sta
!        je=j_end
!        ks=k_sta
!        ke=k_end
!
!        If(k_sta==1.or.k_end==Nz-1) then
!          Allocate(dumcomxy_r4(1:Nx-1,js:je,5))
!          Allocate(dumcomxy_s4(1:Nx-1,js:je,5))
!          Do icom=0,Ndiv_Ny-1
!            If(j_sta==jjsta(icom)) then
!
!              If(k_end==Nz-1) then
!                Do j = js,je
!                Do i = 1, Nx-1, 1
!                  dumcomxy_s4(i,j,1) = F1(i,j,Nz-5)
!                  dumcomxy_s4(i,j,2) = F1(i,j,Nz-4)
!                  dumcomxy_s4(i,j,3) = F1(i,j,Nz-3)
!                  dumcomxy_s4(i,j,4) = F1(i,j,Nz-2)
!                  dumcomxy_s4(i,j,5) = F1(i,j,Nz-1)
!                End Do; End Do
!                dum_len=(Nx-1-1+1)*(je-js+1)*5
!                Call mpi_isend(dumcomxy_s4(1,js,1),dum_len,MPI_REAL,itable(icom,0),1,MPI_COMM_WORLD,isend(1),ierr)
!                Call mpi_irecv(dumcomxy_r4(1,js,1),dum_len,MPI_REAL,itable(icom,0),1,MPI_COMM_WORLD,isend(2),ierr)
!              End If
!
!              If(k_sta==1) then
!                Do j = js, je
!                Do i = 1, Nx-1, 1
!                  dumcomxy_s4(i,j,1) = F1(i,j,1)
!                  dumcomxy_s4(i,j,2) = F1(i,j,2)
!                  dumcomxy_s4(i,j,3) = F1(i,j,3)
!                  dumcomxy_s4(i,j,4) = F1(i,j,4)
!                  dumcomxy_s4(i,j,5) = F1(i,j,5)
!                End Do; End Do
!                dum_len=(Nx-1-1+1)*(je-js+1)*5
!                Call mpi_isend(dumcomxy_s4(1,js,1),dum_len,MPI_REAL,itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(2),ierr)
!                Call mpi_irecv(dumcomxy_r4(1,js,1),dum_len,MPI_REAL,itable(icom,Ndiv_Nz-1),1,MPI_COMM_WORLD,isend(1),ierr)
!              End if
!
!              Call mpi_wait(isend(1),istatus,ierr)
!              Call mpi_wait(isend(2),istatus,ierr)
!
!            End if
!          End Do
!
!
!          If(k_end==Nz-1) then
!            Do j = js, je
!            Do i = 1, Nx-1, 1
!              F1(i,j,Nz  )=dumcomxy_r4(i,j,1)
!              F1(i,j,Nz+1)=dumcomxy_r4(i,j,2)
!              F1(i,j,Nz+2)=dumcomxy_r4(i,j,3)
!              F1(i,j,Nz+3)=dumcomxy_r4(i,j,4)
!              F1(i,j,Nz+4)=dumcomxy_r4(i,j,5)
!            End Do; End Do
!          End If
!
!          If(k_sta==1) then
!            Do j = js, je
!            Do i = 1, Nx-1, 1
!              F1(i,j, 0)=dumcomxy_r4(i,j,5)
!              F1(i,j,-1)=dumcomxy_r4(i,j,4)
!              F1(i,j,-2)=dumcomxy_r4(i,j,3)
!              F1(i,j,-3)=dumcomxy_r4(i,j,2)
!              F1(i,j,-4)=dumcomxy_r4(i,j,1)
!            End Do; End Do
!          End if
!
!          Deallocate(dumcomxy_r4)
!          Deallocate(dumcomxy_s4)
!        End If
!
!        Return
!      End Subroutine
!!=============================================================================================================================
!!=====[ BC_Periodic_y_dir_MPI_R4 ]============================================================================================
!!=============================================================================================================================
!      Subroutine BC_Periodic_y_dir_MPI_R4_init( n11 , n12      , n21      , n22      , n31      , n32      , F1 )
!            Call BC_Periodic_y_dir_MPI_R4_init( -4  , Nx_old+4 , js_old-5 , je_old+5 , ks_old-5 , ke_old+5 , Data_r4 )
!
!        Use module_mpi
!        Use Definition
!        Implicit none
!        Integer,intent(in) :: n11 , n21 , n31
!        Integer,intent(in) :: n12 , n22 , n32
!        Real*4,intent(inout) :: F1(n11:n12,n21:n22,n31:n32)
!
!        Integer i, j, k ! ïœêî
!
!        !--- MPI ---
!        Integer js,je,ks,ke
!	integer ierror,ierr
!        integer dum_len
!        integer isend(2),irecv(2),istatus(MPI_STATUS_SIZE)
!        !-----------
!
!        js=n21+5
!        je=n22-5
!        ks=n31+5
!        ke=n32-5
!
!        If(js==1.or.je==Ny_old-1) then
!          Allocate(dumcomxy_r4(1:Nx_old-1,5,ks:ke))
!          Allocate(dumcomxy_s4(1:Nx_old-1,5,ks:ke))
!          Do icom=0,Ndiv_Nz-1
!            If(ks==kksta(icom)) then
!
!              If(je==Ny_old-1) then
!                Do k = ks,ke
!                Do i = 1, Nx-1, 1
!                  dumcomxy_s4(i,1,k) = F1(i,Ny_old-5,k)
!                  dumcomxy_s4(i,2,k) = F1(i,Ny_old-4,k)
!                  dumcomxy_s4(i,3,k) = F1(i,Ny_old-3,k)
!                  dumcomxy_s4(i,4,k) = F1(i,Ny_old-2,k)
!                  dumcomxy_s4(i,5,k) = F1(i,Ny_old-1,k)
!                End Do; End Do
!                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
!                Call mpi_isend(dumcomxy_s4(1,1,ks),dum_len,MPI_REAL,itable(0,icom),1,MPI_COMM_WORLD,isend(1),ierr)
!                Call mpi_irecv(dumcomxy_r4(1,1,ks),dum_len,MPI_REAL,itable(0,icom),1,MPI_COMM_WORLD,isend(2),ierr)
!              End If
!
!              If(j_sta==1) then
!                Do k = ks, ke
!                Do i = 1, Nx-1, 1
!                  dumcomxy_s4(i,1,k) = F1(i,1,k)
!                  dumcomxy_s4(i,2,k) = F1(i,2,k)
!                  dumcomxy_s4(i,3,k) = F1(i,3,k)
!                  dumcomxy_s4(i,4,k) = F1(i,4,k)
!                  dumcomxy_s4(i,5,k) = F1(i,5,k)
!                End Do; End Do
!                dum_len=(Nx-1-1+1)*(ke-ks+1)*5
!                Call mpi_isend(dumcomxy_s4(1,1,ks),dum_len,MPI_REAL,itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(2),ierr)
!                Call mpi_irecv(dumcomxy_r4(1,1,ks),dum_len,MPI_REAL,itable(Ndiv_Ny-1,icom),1,MPI_COMM_WORLD,isend(1),ierr)
!              End if
!
!              Call mpi_wait(isend(1),istatus,ierr)
!              Call mpi_wait(isend(2),istatus,ierr)
!
!            End if
!          End Do
!
!
!          If(j_end==Ny-1) then
!            Do k = ks, ke
!            Do i = 1, Nx-1, 1
!              F1(i,Ny  ,k)=dumcomxy_r4(i,1,k)
!              F1(i,Ny+1,k)=dumcomxy_r4(i,2,k)
!              F1(i,Ny+2,k)=dumcomxy_r4(i,3,k)
!              F1(i,Ny+3,k)=dumcomxy_r4(i,4,k)
!              F1(i,Ny+4,k)=dumcomxy_r4(i,5,k)
!            End Do; End Do
!          End If
!
!          If(j_sta==1) then
!            Do k = ks, ke
!            Do i = 1, Nx-1, 1
!              F1(i, 0,k)=dumcomxy_r4(i,5,k)
!              F1(i,-1,k)=dumcomxy_r4(i,4,k)
!              F1(i,-2,k)=dumcomxy_r4(i,3,k)
!              F1(i,-3,k)=dumcomxy_r4(i,2,k)
!              F1(i,-4,k)=dumcomxy_r4(i,1,k)
!            End Do; End Do
!          End if
!
!          Deallocate(dumcomxy_r4)
!          Deallocate(dumcomxy_s4)
!        End If
!
!        Return
!      End Subroutine
!
!
!=============================================================================================================================
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

