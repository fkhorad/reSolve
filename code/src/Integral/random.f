

      subroutine randa(nn,rnd)
!
! Generates N uniformly distributed random no's -- but probably I can find something better
!
      implicit none
!
      integer, intent(in) :: nn
      double precision, intent(out) :: rnd(*)
!
      integer i
!
      do i=1,nn
        rnd(i) = rand()
      enddo
      return
!
      end subroutine randa
!
!
      subroutine randin(seed1)
      implicit none
      integer, intent(in) :: seed1
      call srand(seed1)
      end subroutine randin
