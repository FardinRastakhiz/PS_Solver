    MODULE matrixBuilder
    !USE globalDeclarations

    IMPLICIT NONE
    SAVE

    CONTAINS

    Subroutine soluteCoefMatrix(AI,AJ,BC,nnode,ntube,poreBnd,CONF,Q,POREV,Dmol, VOLUMEV,AREAV, DT, p_origin,Q_tot_node, nnz)
    IMPLICIT NONE
    !utility
    integer:: I,II, J,inlet_edge,outlet_edge,nnz,nnode,ntube
    !input
    integer,dimension(nnode)::poreBnd
    integer,dimension(2,ntube):: conf
    doubleprecision::  Dmol,lDiff,DT
    doubleprecision,dimension(ntube)::Q,VOLUMEV, AREAV
    doubleprecision,dimension(nnode)::Q_tot_node
    doubleprecision,dimension(nnode)::POREV,p_origin

    !internal
    doubleprecision,dimension(nnode)::COEF3
    doubleprecision,dimension(ntube)::COEF1,COEF2,COEF4,COEFB

    integer,dimension(:),allocatable:: ai, aj
    doubleprecision,dimension(:),allocatable:: bc

    lDiff = 0.5d0


    COEFB=0.;COEF1=0.;COEF2=0.;COEF3=0.;COEF4=0.
    !COEFFICIENTS FOR REACTIVE MODEL:
    !WE ADD THE FIRST TWO TERMS OF COEF4 HERE AND ADD THE REST IN THE NEXT LOOP:
    DO I=1,NNODE
        COEF3(I)=POREV(I)/DT + Q_tot_node(I)
    ENDDO

    DO I=1,NTUBE !FIIL TWO COEF. THAT HAVE DIMENTION OF NTUBE
        COEFB(I)=1.0D0+ (Q(I)*DT/VOLUMEV(I))+ &
            + (Dmol*AREAV(I)*2.D0*DT / (VOLUMEV(I)*lDiff))

        COEF1(I)=   ( Q(I)*DT ) &
            /(VOLUMEV(I)*COEFB(I))

        COEF2(I)=(Dmol*AREAV(I)*DT) / (VOLUMEV(I)*lDiff*COEFB(I))
        !AND HERE WE FEEL IN THE THIRD (AND LAST) TERM OF COEF3
        COEF3(CONF(1,I))= COEF3(CONF(1,I))+ &
            Dmol*AREAV(I)/lDiff
        COEF3(CONF(2,I))= COEF3(CONF(2,I))+ &
            Dmol*AREAV(I)/lDiff
    ENDDO

    J=0
    DO I=1,NTube
        if(  poreBnd( CONF(1,I) )  .EQ. 1   .and.   poreBnd( CONF(2,I) ) .EQ. 0  )then
            J=J+1
        elseif(  poreBnd( CONF(1,I) )  .EQ. 0   .and.   poreBnd( CONF(2,I) ) .EQ. 1  )then
            J=J+1
        endif
    ENDDO
    !NOW WE NEED TO KNOW THE NUMBER OF ELEMENTS OF COEF. MATRIX TO ALLOCATE AI AND AJ AND BC MATRICES


    inlet_edge=0
    outlet_edge=0
    DO I=1,NTube
        if(     poreBnd( CONF(1,I) ) .EQ. 1  .or. poreBnd( CONF(2,I) ) .EQ. 1)then
            inlet_edge=inlet_edge+1
        elseif( poreBnd( CONF(1,I) ) .EQ. 2  .or. poreBnd( CONF(2,I) ) .EQ. 2 )then
            outlet_edge=outlet_edge+1
        endif
    ENDDO
    !  Mt=NTube-inlet_edge-outlet_edge




    nnz = NNODE + 2*(NTUBE-inlet_edge) + J
    
    ALLOCATE(AI(nnz),AJ(nnz),BC(nnz))
    AI=0;AJ=0;BC=0.d0

    !NOW FILL THE MATRIX
    J=0
    !FILL THE LEFT BOUNDARY
    DO I=1,NNode
        IF (   poreBnd(I) .EQ. 1 )THEN
            J=J+1
            AI(J)= I
            AJ(J)= I
            BC(J)= 1.0D0
        ELSE
            J=J+1
            AI(J)=I
            AJ(J)=I
            BC(J)= COEF3(I)
        ENDIF
    ENDDO


    DO I=1,NTUBE !SO WE NEED THE UPSTREAM NODE FOR THIS TUBE
        IF( p_origin(CONF(1,I)) .GT. p_origin(CONF(2,I))  )THEN
            II=CONF(2,I) !II=CONF(2,I)= CONSIDERING NODE (DOWNSTREAM NODES)
            IF (   poreBnd( CONF(2,I) ) .NE. 1 )THEN !ACTIVE NODES
                J=J+1
                AI(J)=CONF(2,I)
                AJ(J)=CONF(1,I)
                BC(J)=-Q(I)*(COEF1(I)+COEF2(I))
                !now we add the second sumation which is over all:
                !the second summation is over all nodels. however, we break it into
                !two summation, one over IN one the other one over OUT tubdes:
                !for IN, we need to add in to the obove terms:
                BC(J)=BC(J)-Dmol*AREAV(I)*(COEF1(I)+COEF2(I))/lDiff
            ENDIF
            BC(II) = BC(II) - Q(I)*COEF2(I) !and add exta term to the main diagonal
            BC(II) = BC(II) - Dmol*AREAV(I)*(COEF2(I))/lDiff

            IF (   poreBnd( CONF(1,I) ) .NE. 1  )THEN !ACTIVE NODES
                !for OUT we will have new locations in matrix to fill:
                J=J+1
                AI(J) = CONF(1,I)
                AJ(J) = CONF(2,I)
                BC(J) = -Dmol*AREAV(I)*(COEF1(I)+COEF2(I))/lDiff
                BC(CONF(1,I)) = BC(CONF(1,I)) - Dmol*AREAV(I)*(COEF2(I))/lDiff
            ENDIF
        ELSE
            II=CONF(1,I) !II=CONF(2,I)= CONSIDERING NODE (DOWNSTREAM NODES)
            IF (   poreBnd( CONF(1,I) ) .NE. 1  )THEN !ACTIVE NODES
                J=J+1
                AI(J)=CONF(1,I)
                AJ(J)=CONF(2,I)
                BC(J)=-Q(I)*(COEF1(I)+COEF2(I))
                !now we add the second sumation which is over all:
                BC(J)=BC(J)-Dmol*AREAV(I)*(COEF1(I)+COEF2(I))/lDiff
            ENDIF
            BC(II) = BC(II) -  Q(I)*COEF2(I) !and add exta term to the main diagonal
            BC(II) = BC(II) -  Dmol*AREAV(I)*(COEF2(I))/lDiff

            IF (   poreBnd( CONF(2,I) ) .NE. 1  )THEN !ACTIVE NODES
                !for OUT we will have new locations in matrix to fill:
                J=J+1
                AI(J)=CONF(2,I)
                AJ(J)=CONF(1,I)
                BC(J)=-Dmol*AREAV(I)*(COEF1(I)+COEF2(I))/lDiff
                BC(CONF(2,I))=BC(CONF(2,I))- Dmol*AREAV(I)*(COEF2(I))/lDiff
            ENDIF
        ENDIF
    ENDDO

    IF(SUM(ABS(BC.EQ.0)).GT.0 .AND. Dmol.GT.1D-12)THEN
        PRINT*,"matrix cantains 0 in it!!!"
        print*,"MIN AI",MINVAL(AI)
        print*,"MAX AI",MAXVAL(AI)
        print*,"MIN AJ",MINVAL(AJ)
        print*,"MAX AJ",MAXVAL(AJ)
        PRINT*,"nnz AND BC",nnz,SHAPE(BC);PAUSE
    ENDIF

    End subroutine soluteCoefMatrix


    Subroutine soluteCoefRHS(rhs, nnode,ntube,poreBnd,CONF,Q,POREV,Dmol, VOLUMEV,AREAV, DT, p_origin,conc_node,conc_edje,PASSTO,conc_inlet)
    IMPLICIT NONE
    !utility
    integer:: I,J ,inlet_edge,outlet_edge,nnz,nnode,ntube
    doubleprecision:: conc_inlet
    logical:: PASSTO
    doubleprecision,dimension(nnode)::rhs,conc_node
    doubleprecision,dimension(ntube)::COEFB,conc_edje


    !input
    integer,dimension(nnode)::poreBnd
    integer,dimension(2,ntube):: conf
    doubleprecision::  Dmol,lDiff,DT
    doubleprecision,dimension(ntube)::Q,VOLUMEV, AREAV

    doubleprecision,dimension(nnode)::POREV,p_origin

    !internal
    doubleprecision,dimension(nnode)::Q_tot_node

    lDiff = 0.5d0
    COEFB=0.0d0
    DO I=1,NTUBE !FIIL TWO COEF. THAT HAVE DIMENTION OF NTUBE
        COEFB(I)=1.0d0+ (q(I)*DT/VOLUMEV(I))+ &
            + (Dmol*AREAV(I)*2.D0*DT / (VOLUMEV(I)*lDiff))
    ENDDO

    rhs   = 0.0D0
    rhs(:)=POREV(:)*conc_node(:)/DT !NOW WE HAVE TRANSFERED Ci*V/DT
    DO I=1,NTUBE
        IF (p_origin(CONF(1,I)).GT.p_origin(CONF(2,I))) THEN
            J=CONF(2,I) !J IS DOWNSTREAM NODES
            rhs(J)=rhs(J)+ &
                ( Q(I)*conc_edje(I)/COEFB(I) )
            !DIFFUSIVE FLUX
            rhs(J)=rhs(J)+ &
                ( Dmol*AREAV(I)*conc_edje(I)/ (COEFB(I)*lDiff) )
            J=CONF(1,I) !NOW J IS UPSTREAM FOR DIFFUSIVE
            rhs(J)=rhs(J)+ &
                ( Dmol*AREAV(I)*conc_edje(I)/ (COEFB(I)*lDiff) )
        ELSE
            J=CONF(1,I) !J IS DOWNSTREAM NODES
            rhs(J)=rhs(J)+ &
                ( Q(I)*conc_edje(I)/COEFB(I) )
            !DIFFUSIVE FLUX
            rhs(J)=rhs(J)+ &
                ( Dmol*AREAV(I)*conc_edje(I)/ (COEFB(I)*lDiff) )
            J=CONF(2,I) !NOW J IS UPSTREAM FOR DIFFUSIVE
            rhs(J)=rhs(J)+ &
                ( Dmol*AREAV(I)*conc_edje(I)/ (COEFB(I)*lDiff) )
        ENDIF
    ENDDO

    PASSTO=1
    IF(PASSTO)THEN
        !rhs(1:LNODES,I)=C0_comp(I)
        DO J=1,NNODE
            IF(poreBnd(J).EQ.1 )THEN !.and. J.gt.100
                rhs(J)=conc_inlet
            ENDIF
        ENDDO
    ELSE
        !rhs(1:LNODES,1:NSolCom)=0.00_8
        DO J=1,NNODE
            IF(poreBnd(J).EQ.1 )THEN !.and. J.gt.100
                rhs(J)=0.D0
            ENDIF
        ENDDO
    ENDIF

    End subroutine soluteCoefRHS



    END MODULE matrixBuilder