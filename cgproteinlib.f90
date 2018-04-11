SUBROUTINE ReimageChain(Pos, NAtom, BoxL, ReimagePos)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(OUT), DIMENSION(0:NAtom-1, 0:2) :: ReimagePos
    
    INTEGER :: i
    REAL(8), DIMENSION(0:2) :: Posi, Pos0, Posref, invBoxL
    invBoxl = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Pos0 = 0.d0
    Posref = 0.d0
    ReimagePos = 0.d0

    DO i = 0, NAtom-1
        ! min image i+1 w.r.t. i
        Posi = Pos(i,:)
        Posref = Posi - Pos0
        Posref = Posref - BoxL * DNINT(Posref * invBoxL)
        ReimagePos(i,:) = Posi(:) + Posref
        Pos0 = Posi
    ENDDO
END SUBROUTINE


SUBROUTINE RG_frame(Rg, Pos, NAtom, BoxL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(OUT) :: Rg

    REAL(8), DIMENSION(0:2) :: Posi, Pos0, Posref, PosSum, PosSumsq, invBoxL
    INTEGER :: i,j
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Pos0 = 0.d0
    Posref = 0.d0
    PosSum = 0.d0
    PosSumsq = 0.d0
    Rg = 0.d0

    DO i = 0, NAtom-1
        ! min image i+1 w.r.t i
        Posi = Pos(i,:)
        Posref = Posi - Pos0
        Posref = Posref - BoxL * DNINT(Posref * invBoxL)
        Posi = Pos0 + Posref
        PosSum = PosSum + Posi
        PosSumsq = PosSumsq + Posi * Posi
        Pos0 = Posi
    ENDDO
    Rg = SQRT( SUM( PosSumsq/NAtom - (PosSum*PosSum)/(NAtom*NAtom) ) )
END SUBROUTINE


SUBROUTINE REE_frame(Ree, Pos, NAtom, BoxL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NAtom
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NAtom-1, 0:2) :: Pos
    REAL(8), INTENT(OUT) :: Ree

    REAL(8), DIMENSION(0:2) :: Posi, Pos0, Posref, FirstPos, LastPos, invBoxL
    INTEGER :: i
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Pos0 = 0.d0
    Posref = 0.d0
    FirstPos = 0.d0
    LastPos = 0.d0
    Ree = 0.d0

    DO i = 0, NAtom-1
        ! min image i+1 w.r.t i
        Posi = Pos(i,:)
        Posref = Posi - Pos0
        Posref = Posref - BoxL * DNINT(Posref * invBoxL)
        Posi = Pos0 + Posref
        IF (i==0) FirstPos = Posi
        IF (i==NAtom-1) LastPos = Posi
        Pos0 = Posi
    ENDDO
   Ree = SQRT( SUM( (LastPos - FirstPos) * (LastPos - FirstPos) ) )
END SUBROUTINE


SUBROUTINE ResPairs_frame(ContactMap, ContactDist, ResPos, NRes, &
ResRadius, MinCO, BoxL)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NRes
    REAL(8), INTENT(IN) :: ResRadius, MinCO
    REAL(8), INTENT(IN), DIMENSION(0:2) :: BoxL
    REAL(8), INTENT(IN), DIMENSION(0:NRes-1, 0:2) :: ResPos
    REAL(8), INTENT(OUT), DIMENSION(0:NRes-1, 0:NRes-1) :: ContactMap, ContactDist
    
    REAL(8), DIMENSION(0:2) :: Posi, Posj, invBoxL, rij
    REAL(8) :: rsq
    INTEGER :: i, j
    ContactMap = 0.d0
    ContactDist = 0.d0
    invBoxL = MERGE(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    
    DO i = 0, NRes-2
        Posi = ResPos(i,:)
        DO j = i+1, NRes-1
            Posj = ResPos(j,:)
            rij = Posj - Posi
            rij = rij - BoxL * DNINT(rij * invBoxL) ! min-image
            rsq = SUM(rij * rij)
            ContactDist(i,j) = SQRT(rsq)
            IF ( (rsq <= ResRadius * ResRadius) .AND. (j-i >= MinCO) ) ContactMap(i,j) = 1
        ENDDO
    ENDDO
END SUBROUTINE
