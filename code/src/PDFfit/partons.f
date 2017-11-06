
      subroutine partons(sq2,sx,fx,nf,isetproton,ippbar)
      implicit none
!
      INTERFACE
!
      subroutine partons_cc(sq2,sx,
     .  fx0,fx1,fx2,fx3,fx4,fx5,fxm1,fxm2,fxm3,fxm4,fxm5,
     .  nf,isetproton,ippbar)
!
      integer, intent(in) :: nf,isetproton,ippbar
      double precision, intent(in) :: sq2,sx
!
      double precision, intent(out) ::
     . fx0,fx1,fx2,fx3,fx4,fx5,fxm1,fxm2,fxm3,fxm4,fxm5
!
      end subroutine
!
      END INTERFACE
!
      integer, intent(in) :: nf,isetproton,ippbar
      double precision, intent(in) :: sq2,sx
!
      double precision, intent(out) :: fx(-5:5)
!
      double precision fx0,fx1,fx2,fx3,fx4,fx5,fxm1,fxm2,fxm3,fxm4,fxm5
      integer ii
!
      call partons_cc(sq2,sx,
     . fx0,fx1,fx2,fx3,fx4,fx5,fxm1,fxm2,fxm3,fxm4,fxm5,
     . nf,isetproton,ippbar)
      
      fx(0)=fx0

      fx(1)=fx1
      fx(2)=fx2
      fx(3)=fx3
      fx(4)=fx4
      fx(5)=fx5

      fx(-1)=fxm1
      fx(-2)=fxm2
      fx(-3)=fxm3
      fx(-4)=fxm4
      fx(-5)=fxm5

      end subroutine partons