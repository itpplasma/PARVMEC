!  The mgrid interpolation in becoil must reproduce a field that is cubic in
!  R and in Z to rounding error, also next to the grid edges. Bilinear
!  interpolation leaves an O(h**2) error of about 1E-2 here.
PROGRAM test_becoil
   USE stel_kinds, ONLY: rprec
   USE vacmod0, ONLY: nuv3, nv
   USE mgrid_mod, ONLY: nr0b, nz0b, np0b, rminb, rmaxb, zminb, zmaxb,      &
                        delrb, delzb
   USE parallel_vmec_module, ONLY: nuv3min, nuv3max

   IMPLICIT NONE

   INTEGER, PARAMETER :: nr = 9, nz = 11, npts = 400
   REAL(rprec), DIMENSION(nr,nz,1) :: brvac, bpvac, bzvac
   REAL(rprec), DIMENSION(npts) :: rad, zee, br, bp, bz
   REAL(rprec) :: br0, bp0, bz0, err
   INTEGER :: i, j

   nr0b = nr
   nz0b = nz
   np0b = 1
   rminb = 1.0_rprec
   rmaxb = 2.6_rprec
   zminb = -0.9_rprec
   zmaxb = 0.9_rprec
   delrb = (rmaxb - rminb)/(nr - 1)
   delzb = (zmaxb - zminb)/(nz - 1)
   nv = 1
   nuv3 = npts
   nuv3min = 1
   nuv3max = npts

   DO j = 1, nz
      DO i = 1, nr
         CALL field(rminb + (i - 1)*delrb, zminb + (j - 1)*delzb,              &
                    brvac(i,j,1), bpvac(i,j,1), bzvac(i,j,1))
      END DO
   END DO

   !  Points spread over the whole box, including the cells at its edges.
   DO i = 1, npts
      rad(i) = rminb + (rmaxb - rminb)*MOD(0.6180339887_rprec*i, 1.0_rprec)
      zee(i) = zminb + (zmaxb - zminb)*MOD(0.4142135624_rprec*i, 1.0_rprec)
   END DO

   CALL becoil(rad, zee, br, bp, bz, brvac, bpvac, bzvac, .false.)

   err = 0
   DO i = 1, npts
      CALL field(rad(i), zee(i), br0, bp0, bz0)
      err = MAX(err, ABS(br(i) - br0), ABS(bp(i) - bp0), ABS(bz(i) - bz0))
   END DO

   WRITE (*,*) 'max interpolation error', err
   IF (err .gt. 1.0E-12_rprec) THEN
      STOP 1
   END IF

CONTAINS

   SUBROUTINE field(r, z, br, bp, bz)
      REAL(rprec), INTENT(in)  :: r, z
      REAL(rprec), INTENT(out) :: br, bp, bz

      br = r**3 - 2*r*z**2 + z
      bp = 2 - r*z + z**3*r**2
      bz = r**2*z - 0.5_rprec*z**2 + r**3*z**3

   END SUBROUTINE field

END PROGRAM test_becoil
