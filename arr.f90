program SHOCK_WAVE_S_F
    use mpi
    implicit none
    ! DIFFERENCE UP FLOW ON THREE POINTS
    ! In points X(0),X(1),X(MX-1),X(MX) boundary values of function
    ! Uniform grids (both)
    integer ierr, r, g, id_p, id_z, np, nr, dn, tag, xxx
    integer, parameter:: K=selected_real_kind(14,70)
    real(K),parameter:: S2PI=2.506628274631_K, C23=2._K/3.
    integer,dimension(1):: IM
    integer,dimension(4):: LB  ! boundary point indices
    real(K),dimension(:),allocatable:: D,U,Tt,Tr,T,TET,pXX,fX,omX, D_
    real(K),dimension(4):: DB,UB,TtB,TrB         !  boundary values
    real(K),dimension(:),allocatable:: V,CV, tau, AB,TtR,TrR,Q3,Q4 ! V=abs(Ksi)
    real(K),dimension(:,:),allocatable:: FM5,FP5,FE5,FM7,FP7,FE7,     & ! Fn,Fp,Fom
            &  RM5,RP5,RE5,RM7,RP7,RE7,    & ! (F+)n,p,om
            &  E5,E7                       ! FMn+
    real(K),dimension(:,:),allocatable:: C5,C7, W5,W7, QP,QF,QT  
        ! C - thermal speed. with sign; W- complex H in F+
    real(K) GAM,CM,ttp,s,Z,Tkp,HX,HT,HV,H0,H1,H2, EPS,G53,G53G,G32, Q0 ! G53G=Kappa
    real(K) kP,kF,kTET
    !nzn = 50
    integer ITERMAX,ITER_PRINT,MX,KRT,KHV,MV,MX1,MX2,MX3,M1X,M1,ITER,IDom, I,J
    !-----------------------------------------------------------------------------
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, id_p, ierr )
    call MPI_Comm_size(MPI_COMM_WORLD, np, ierr )
    call OPEN_ALLOC_READ_PREP
    !-----------------------------------------------------------------------------
    !---------------------- Iteration Loop  --------------------------------------
    !-----------------------------------------------------------------------------
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(id_p == 0) then
        write(*,*) 'CONTINUE'
    end if
    MPI_X: if(id_z == 0) then
        xxx = 50
        do i = 1, np - 1
            nr = nr + 1
            dn = i
            tag = nr
            call MPI_Send(W5(nr,1:xxx), xxx, MPI_DOUBLE_PRECISION, dn, tag, MPI_COMM_WORLD, ierr)
        end do
    end if MPI_X
    ITER_CYCLE: do ITER=1,ITERMAX
    !############ FMn+ -function across the field #############################
    TtR=T+G53*kTET*TET; TrR=T-G32*kTET*TET    ! Tt+ & Tr+
    QP=spread(D/S2PI/sqrt(TtR),1,M1);      QT=spread(2.*TtR,1,M1)
    C5=spread(V,2,M1X)-spread(U,1,M1)
    C7=-spread(V,2,M1X)-spread(U,1,M1)
    QF=exp(C5**2/QT);        E5=QP/QF
    QF=exp(C7**2/QT);        E7=QP/QF
    !############ F+ -function across the field ###############################
    QT=spread(TtR,1,M1)
    QP=0.5*spread(kP*pXX/D/TtR,1,M1);   QF=0.2*spread(kF*fX/D/TtR**2,1,M1)
    W5=C5**2/QT;                        W7=C7**2/QT
    RM5=E5*(1.+QP*(W5-1.)+QF*C5*(W5-3.)); RM7=E7*(1.+QP*(W7-1.)+QF*C7*(W7-3.))
    RP5=2.*QT*E5*(1.+QP*(W5-2.)+QF*C5*(W5-1.))
    RP7=2.*QT*E7*(1.+QP*(W7-2.)+QF*C7*(W7-1.))
    select case(IDom); case(2)
        QT=spread(TrR,1,M1); QF=spread(kP*omX/D/TtR/TrR,1,M1)
        RE5=QT*(G53G*RM5+E5*QF*C5);          RE7=QT*(G53G*RM7+E7*QF*C7)
    end select
    !################ F-function in all points of the inner area ##############
    Q3=HT/tau
    CYCLE_X_UP: do I=2,MX2
            AB=1._K+H0*V+Q3(I)
        FM5(:,I)=(FM5(:,I)+Q3(I)*RM5(:,I)+V*(H1*FM5(:,I-1)-H2*FM5(:,I-2)))/AB
        FP5(:,I)=(FP5(:,I)+Q3(I)*RP5(:,I)+V*(H1*FP5(:,I-1)-H2*FP5(:,I-2)))/AB
        select case(IDom); case(2)
            FE5(:,I)=(FE5(:,I)+Q3(I)*RE5(:,I)+V*(H1*FE5(:,I-1)-H2*FE5(:,I-2)))/AB  
        end select
        C_IF1: if(IDom == 1) then
            !FE5(:,I)=(FE*X3);
        end if C_IF1
        call MPI_Bcast(FM5, I, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(FP5, I, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast (FE5, I, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send(RE5, 1, MPI_DOUBLE_PRECISION, 0, RE7, MPI_COMM_WORLD, ierr )
    end do  CYCLE_X_UP
    CYCLE_X_DOWN: do I=MX2,2,-1
            AB=1._K+H0*V+Q3(I)
        FM7(:,I)=(FM7(:,I)+Q3(I)*RM7(:,I)+V*(H1*FM7(:,I+1)-H2*FM7(:,I+2)))/AB
        FP7(:,I)=(FP7(:,I)+Q3(I)*RP7(:,I)+V*(H1*FP7(:,I+1)-H2*FP7(:,I+2)))/AB
        select case(IDom); case(2)
            FE7(:,I)=(FE7(:,I)+Q3(I)*RE7(:,I)+V*(H1*FE7(:,I+1)-H2*FE7(:,I+2)))/AB             
        end select
    end do  CYCLE_X_DOWN
    !############ Macroparameters throughout the field #########################
    D=sum((FM5+FM7)*spread(CV,2,M1X),dim=1);   D(LB)=DB
    U=sum((FM5-FM7)*spread(CV*V,2,M1X),dim=1)/D;  U(LB)=UB
    Tt=sum((C5**2*FM5+FP5+C7**2*FM7+FP7)*spread(CV,2,M1X),dim=1)/3./D;  Tt(LB)=TtB
    pXX=sum((C5**2*FM5+C7**2*FM7)*spread(CV,2,M1X),dim=1)-D*Tt;  pXX(LB)=0.
    fX=0.5*sum((C5**3*FM5+C5*FP5+C7**3*FM7+C7*FP7)*spread(CV,2,M1X),dim=1); fX(LB)=0.
    select case(IDom)
        case(2)
            Tr=sum((FE5+FE7)*spread(CV,2,M1X),dim=1)/G53G/D;   Tr(LB)=TrB 
            omX=sum((C5*FE5+C7*FE7)*spread(CV,2,M1X),dim=1);  omX(LB)=0.
        case(1)
            Tr=Tt; omX=0.
    end select
    tau=ttp*Tt**(s-1.)/D;  TET=Tt-Tr; T=G32*Tt+G53*Tr
    !############ Convergence ##################################################
    EPS=maxval(abs(D-D_)/D)
    D_=D
    !###################### Prompt printing ####################################
    if((ITER/ITER_PRINT)*ITER_PRINT==ITER) then
    IM=maxloc(abs(U(3:MX2)-U(2:MX3)))
    I=IM(1)+1
    write(*,99) ITER,EPS,'  Front=>',I,'   GorbT=>',maxval(Tt(5:MX-5))
    Q3=(D-1.)/(D(MX)-1.);  Q4=0.; Q4(10:MX-10)=(Q3(11:MX-9)-Q3(9:MX-11))/HX/2.
    IM=maxloc(Q4);  J=IM(1)-1
    write(*,998) 'Kn=',1.277*maxval(Q4),'(',J,')'
    998 format(42X,A,F7.4,2X,A,I3,A)

    99 format(1X,I5,ES18.9,A,I4,A,ES15.7)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do ITER_CYCLE
    call WRITE_DEALL_CLOSE
    call MPI_Finalize(ierr)
contains 
    subroutine OPEN_ALLOC_READ_PREP
        !------------------------------------------------------------------------------
        open(unit=1,file='C:\Users\MAXIMUS\source\repos\MPI_Array\MPI_Array\x64\Debug\INIT_SW.IN',status='Old',form='Formatted',action='Read')
        open(unit=2,file='C:\Users\MAXIMUS\source\repos\MPI_Array\MPI_Array\x64\Debug\GRIDS.OUT',status='Old',form='Unformatted',action='Read')
        open(unit=3,file='C:\Users\MAXIMUS\source\repos\MPI_Array\MPI_Array\x64\Debug\FIELDS.OUT',status='Old',form='Unformatted')
        !------------------------------------------------------------------------------
        read(1,1,err=10) ITERMAX,ITER_PRINT, GAM,CM, ttp,s,Z,Tkp, MX,HX,HT, KRT,KHV
        1 format(3/, 2(1X,I7,/),/, 2(1X,F7.6,/),/, 4(1X,F7.6,/),/, 1X,I7,/,2(1X,F7.6,/),&
                                            & /,2(1X,I7,/))
        read(2,err=20) MV
        !--------------------- Service Constants ------------------------------------
        MX1=MX-1;   MX2=MX-2;   MX3=MX-3;   M1X=MX+1;  M1=MV+1
        LB=(/0,1,MX1,MX/)
        !------------------------------------------------------------------------------
        allocate(  D(0:MX),U(0:MX),Tt(0:MX),Tr(0:MX),T(0:MX),TET(0:MX), &
        & pXX(0:MX),fX(0:MX),omX(0:MX),      tau(0:MX), &
        & FM5(0:MV,0:MX),FP5(0:MV,0:MX),FE5(0:MV,0:MX), &
        & FM7(0:MV,0:MX),FP7(0:MV,0:MX),FE7(0:MV,0:MX), &
        & RM5(0:MV,0:MX),RP5(0:MV,0:MX),RE5(0:MV,0:MX), &
        & RM7(0:MV,0:MX),RP7(0:MV,0:MX),RE7(0:MV,0:MX), &
        & E5(0:MV,0:MX),E7(0:MV,0:MX),  D_(0:MX)        )
        allocate( C5(0:MV,0:MX),C7(0:MV,0:MX),  V(0:MV),CV(0:MV),&
        & W5(0:MV,0:MX),  QT(0:MV,0:MX),  QP(0:MV,0:MX),&
        & W7(0:MV,0:MX),                  QF(0:MV,0:MX),&
        & TtR(0:MX),TrR(0:MX),Q3(0:MX),Q4(0:MX),AB(0:MV))
        read(2,err=20) V,CV,H0,H1,H2
        read(3,err=30) FM5,FP5,FE5,FM7,FP7,FE7,  D,U,Tt,Tr,T,TET,pXX,fX,omX
        close(1); close(2)
        !------------------------ Phyisical Contants --------------------------------
        G53=(5._K-3.*GAM)/2.; G53G=G53/(GAM-1.); G32=1.5*(GAM-1.)
        if(ttp<0.001_K)  ttp=(7._K-2.*s)*(5._K-2.*s)/30.
        if(Z<0.001_K) then; write(*,*) 'Z=0 ввводить нельзя'; stop 1000; end if
        kP=1.-ttp; kF=1.-C23*ttp; kTET=1.-1./Z
        if(abs(GAM-1.6666) > 1.E-3) then; IDom=2; else; IDom=1; G53=0.; G32=1.;end if
        tau=ttp*Tt**(s-1.)/D
        DB=(/1._K,1._K,D(MX),D(MX)/); UB=(/U(0),U(0),U(MX),U(MX)/)
        TtB=(/1._K,1._K,Tt(MX),Tt(MX)/);  TrB=TtB
        !------------------------------------------------------------------------------
        D_=D;
        return

        10 write(*,*)'ERROR READ INIT_SW.IN';  stop 1
        20 write(*,*)'ERROR WRITE GRIDS.OUT';   stop 2
        30 write(*,*)'ERROR WRITE FIELDS.OUT';  stop 3
    end subroutine OPEN_ALLOC_READ_PREP
    subroutine WRITE_DEALL_CLOSE
        ! Data output to unformat file
        rewind 3
        write(3,err=40)  FM5,FP5,FE5,FM7,FP7,FE7,  D,U,Tt,Tr,T,TET,pXX,fX,omX
        close(3)  
        return
        40 write(*,*)'ERROR WRITE FIELDS.OUT';  stop 4
    end subroutine WRITE_DEALL_CLOSE
end program SHOCK_WAVE_S_F
