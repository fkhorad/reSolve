c
c PDF FITTING PROGRAM -- NEEDS CERNLIB MINUIT
c
c Original code by D. De Florian -- unpublished as a stand-alone, but
c contained in public codes such as DYres, ...
c
c Extracted and modified by F. Coradeschi, 2017
c
c USAGE:
c
c 1) call fiteador(xtauf,muf2,energy_sector,pdf_label)
c
c where double precision xtauf is the MAXIMUM value of the rate of
c s_partonic / s_hadronic
c to be used in the calculation, double precision muf2 is the squared
c factorization scale, integer energy_sector <= 20 is a label allowing
c for more than one fit to be stored at the same time, and integer
c pdf_label identifies the PDF set to be used.
c
c All of fiteador dummy parameters are INPUTS; the OUTPUT -- the actual
c fit parameters -- is stored in several Common blocks as described below.
c
c 2) For the program to work, the user NEEDS to define, on its own, a
c Fortran function called "partons" which acts as an INTERFACE between the
c the program and the PDF set(s) to be fitted.
c The function template is
c
c partons(sq2,sx,fx,nf,isetproton,ippbar)
c
c with INPUTS:
c double precision sq2, sx
c integer nf, isetproton, ippbar
c
c being the squared factorization scale (sq2), the momentum fraction (sx), the
c number of active flavours (nf), a label for the PDF set to be used (isetproton),
c and a label for the type of hadron (ippbar) respectively. Note that nf and
c ippbar are actually DUMMIES from the point of view of this code: they are
c hardcoded to nf=5 and ippbar=1 (originally meaning a proton). This could
c be modified with some (undocumented) effort on part of the user.
c
c the OUTPUT:
c double precision fx(-nf:nf)
c
c returns an array of x*PDF value for the various flavours, corresponding to
c the given factorization scale, momentum fraction and PDF used.
c
c
c DETAILED DESCRIPTION
c
c This program fits an arbitrary set of PDFs to a standard (though purely
c empirical) analytic function of momentum fraction x, at fixed factorization
c scale mu_f. The ultimate goal is to allow the PDF evaluation in Mellin space.
c The actual functional form (at a given mu_f) is
c
c        f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+a7*x**2
c    .        +a8*X**(aa))
c
c with a1-a8 fit parameters and aa an auxiliary parameter which is set at
c compile time, with aa=2 by default
c
c The core fit subroutine "fiter" uses CernLib Minuit legacy library.
c This MUST pre-installed separately.
c
c The fit is split into integer NFITMAX different kinematical regions
c (according to the partons total pseudorapidity) in order to improve precision.
c NFITMAX is in principle tunable, but changing it to a value different
c from the present NFITMAX=14 requires some (undocumented) handwork on
c part of the user.
c
c Output variable names:
c the resulting fit parameters are saved in the AnXX(m,j) and AnXXp(m,j) variables,
c all of which are contained in CXX1 and CXX2 common blocks: here An stands
c for the par. name a1-a8, XX labels the flavour
c -- Codenames for the flavours are UV=U-Ubar(valence), DV=D-Dbar(val.),
c US=U sea,DS=D sea,SS=S,CH=C,BO=B,GL=g --
c C just stand for Common, m labels the kinematic (pseudorapidity) region,
c and the additional index j allows for fits corresponding to several
c different mu_f values to be stored at once. 1 and 2 refer to beam 1 and
c beam 2 -- an hadron-to-hadron collision is assumed -- and for the variables
c internal to the common blocks, a "p" is appended for beam 2.

c IMPORTANT: The code currently does not distinguish between quarks and antiquarks
c for flavours from strange on and it does not include either the top/antitop
c nor the photon PDFs. This could be modified with some (undocumented) effort
c on part of the user.
c

      subroutine fiteador(xtauf,muf2,energy_sector,pdf_label)
       IMPLICIT NONE
c Dummy parameters
       double precision, intent(in) :: xtauf, muf2
       integer, intent(in) :: energy_sector, pdf_label
c
c Local parameters
       double precision yav(30), xtau1, xtau2, qmin_p, q_temp, etam
       integer kk, jj
c
c Global parameters (commons)
       include 'pdffit.h'

C Maximum rapidity -- limited to avoid reaching the end of phase space
       etam=-0.5*dlog(xtauf)
       etam=IDINT(etam*10)/10d0

c Set auxiliary function parameter aa
       aa=2.5d0

c Store muf2 used for this fit
       if(energy_sector > 20) then
         write(*,*) 'Integer label for the fit exceedes maximum preset'
     .   ,'value of 20'
         stop
       endif
       muf_array(energy_sector) = muf2

c Set additional commons
       NFITMAX=14
       N_ENERGY_SECS = 20
       isetproton = pdf_label
c

       yav(1)=0.5d0
       yav(2)=1.5d0
       yav(3)=2.5d0
       yav(4)=3.25d0
       yav(5)=3.75d0
       yav(6)=4.d0
       yav(7)=4.5d0
       if (etam.lt.1.0001d0) yav(1)=0d0
       if (etam.lt.2.0001d0) yav(2)=yav(1)
       if (etam.lt.3.0001d0) yav(3)=yav(2)
       if (etam.lt.3.50001d0) yav(4)=yav(3)
       if (etam.lt.4.0001d0) yav(5)=yav(4)
       if (etam.lt.4.50001d0) yav(6)=yav(5)
       if (etam.lt.5.0001d0) yav(7)=yav(6)

       yav(8)=-0.5d0
       yav(9)=-1.5d0
       yav(10)=-2.5d0
       yav(11)=-3.25d0
       yav(12)=-3.75d0
       yav(13)=-4.d0
       yav(14)=-4.5d0
       if (etam.le.1.0001d0) yav(8)=0d0
       if (etam.le.2.0001d0) yav(9)=yav(8)
       if (etam.le.3.0001d0) yav(10)=yav(9)
       if (etam.le.3.50001d0) yav(11)=yav(10)
       if (etam.le.4.0001d0) yav(12)=yav(11)
       if (etam.le.4.50001d0) yav(13)=yav(12)
       if (etam.le.5.0001d0) yav(14)=yav(13)

        jj = energy_sector

        write(*,*)'Waiting for PDF fit', jj,'...'

        do kk=1,NFITMAX
c
          xtau1=xtauf**0.5* dexp(yav(kk))

          if (xtau1.gt.1d0) then
            write(6,*)xtauf**0.5,yav(kk),kk,etam
            stop
          endif
c
          call fiter(1,xtau1,muf2,A1UV(kk,jj),A2UV(kk,jj),A3UV(kk,jj),
     .      A4UV(kk,jj),A5UV(kk,jj),A6UV(kk,jj),A7UV(kk,jj),A8UV(kk,jj))
          call fiter(2,xtau1,muf2,A1DV(kk,jj),A2DV(kk,jj),A3DV(kk,jj),
     .      A4DV(kk,jj),A5DV(kk,jj),A6DV(kk,jj),A7DV(kk,jj),A8DV(kk,jj))
          call fiter(3,xtau1,muf2,A1US(kk,jj),A2US(kk,jj),A3US(kk,jj),
     .      A4US(kk,jj),A5US(kk,jj),A6US(kk,jj),A7US(kk,jj),A8US(kk,jj))
          call fiter(4,xtau1,muf2,A1DS(kk,jj),A2DS(kk,jj),A3DS(kk,jj),
     .      A4DS(kk,jj),A5DS(kk,jj),A6DS(kk,jj),A7DS(kk,jj),A8DS(kk,jj))
          call fiter(5,xtau1,muf2,A1SS(kk,jj),A2SS(kk,jj),A3SS(kk,jj),
     .      A4SS(kk,jj),A5SS(kk,jj),A6SS(kk,jj),A7SS(kk,jj),A8SS(kk,jj))
          call fiter(6,xtau1,muf2,A1GL(kk,jj),A2GL(kk,jj),A3GL(kk,jj),
     .      A4GL(kk,jj),A5GL(kk,jj),A6GL(kk,jj),A7GL(kk,jj),A8GL(kk,jj))
          call fiter(7,xtau1,muf2,A1CH(kk,jj),A2CH(kk,jj),A3CH(kk,jj),
     .      A4CH(kk,jj),A5CH(kk,jj),A6CH(kk,jj),A7CH(kk,jj),A8CH(kk,jj))
          call fiter(8,xtau1,muf2,A1BO(kk,jj),A2BO(kk,jj),A3BO(kk,jj),
     .      A4BO(kk,jj),A5BO(kk,jj),A6BO(kk,jj),A7BO(kk,jj),A8BO(kk,jj))
c
          xtau2=xtauf**0.5* dexp(-yav(kk))

          if (xtau2.gt.1d0) then
            write(6,*)xtauf**0.5,-yav(kk),kk,etam
            stop
          endif
c
         call fiter(1,xtau2,muf2,A1UVp(kk,jj),A2UVp(kk,jj),A3UVp(kk,jj),
     . A4UVp(kk,jj),A5UVp(kk,jj),A6UVp(kk,jj),A7UVp(kk,jj),A8UVp(kk,jj))
         call fiter(2,xtau2,muf2,A1DVp(kk,jj),A2DVp(kk,jj),A3DVp(kk,jj),
     . A4DVp(kk,jj),A5DVp(kk,jj),A6DVp(kk,jj),A7DVp(kk,jj),A8DVp(kk,jj))
         call fiter(3,xtau2,muf2,A1USp(kk,jj),A2USp(kk,jj),A3USp(kk,jj),
     . A4USp(kk,jj),A5USp(kk,jj),A6USp(kk,jj),A7USp(kk,jj),A8USp(kk,jj))
         call fiter(4,xtau2,muf2,A1DSp(kk,jj),A2DSp(kk,jj),A3DSp(kk,jj),
     . A4DSp(kk,jj),A5DSp(kk,jj),A6DSp(kk,jj),A7DSp(kk,jj),A8DSp(kk,jj))
         call fiter(5,xtau2,muf2,A1SSp(kk,jj),A2SSp(kk,jj),A3SSp(kk,jj),
     . A4SSp(kk,jj),A5SSp(kk,jj),A6SSp(kk,jj),A7SSp(kk,jj),A8SSp(kk,jj))
         call fiter(6,xtau2,muf2,A1GLp(kk,jj),A2GLp(kk,jj),A3GLp(kk,jj),
     . A4GLp(kk,jj),A5GLp(kk,jj),A6GLp(kk,jj),A7GLp(kk,jj),A8GLp(kk,jj))
         call fiter(7,xtau2,muf2,A1CHp(kk,jj),A2CHp(kk,jj),A3CHp(kk,jj),
     . A4CHp(kk,jj),A5CHp(kk,jj),A6CHp(kk,jj),A7CHp(kk,jj),A8CHp(kk,jj))
         call fiter(8,xtau2,muf2,A1BOp(kk,jj),A2BOp(kk,jj),A3BOp(kk,jj),
     . A4BOp(kk,jj),A5BOp(kk,jj),A6BOp(kk,jj),A7BOp(kk,jj),A8BOp(kk,jj))
c
        enddo
c
        write(*,*)'PDF fit ended'

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine fiter(idist,x,q2,a1,a2,a3,a4,a5,a6,a7,a8)
      implicit double precision (A-H,O-Z)
      external fcng1,chi2
      dimension nprm(8), vstrt(8),stp(8),bl(8),bu(8),arglis(16)
      dimension parsal(8)
      character*10 pnam(8)
      data nprm /1,2,3,4,5,6,7,8/
      data pnam  / 'A1','A2','A3','A4','A5', 'A6', 'A7', 'A8'/
      data vstrt /  1D0, -.1D0, 5D0, 3D0, -1D0, 0D0, 0D0, 0D0/
      data stp   /  1D0, 1D0, 1D0, 10D0, 15D0, 25D0, 25D0, 25D0/
c      data bl    /0D0, -0.8D0, 4D0, -20D0, -40D0, -80D0,-80D0,-80D0/
c      data bu    /0D0, 1.5D0, 19D0, 600D0, 40D0, 80D0, 80D0, 80D0/

c      data bl    /0D0, -0.8D0, 5D0, -5D0, -30D0, -30D0,-30D0,-30D0/
c      data bu    /0D0, 1.2D0, 12D0, 50D0, 30D0, 30D0, 30D0, 30D0/

      data bl    /0D0, -0.8D0, 4D0, -20D0, -60D0, -90D0,-90D0,-90D0/
      data bu    /0D0, 1.5D0, 12D0, 600D0, 60D0, 90D0, 90D0, 90D0/

      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      id=idist
      rx=x
      rq2=q2
c.....initialization :
      call mninit(5,1,1)
c     call mninit(5,6,6)
c.....definitions of the parameters :
      do 11 i = 1, 8
         call mnparm (nprm(i),pnam(i),VSTRT(i),STP(i),BL(i),BU(i),
     .        ierflg,chi2)
         if (ierflg .ne. 0) then
            write (6,*) ' unable to define parameter no.', i
            stop
         end if
 11   continue
c.....output
      arglis(1) = 0.
      call mnexcm(fcng1,'set print',arglis,1,ierflg,chi2)
c.....first call :
      arglis(1) = 1.            !   IFLAG = 1
      call mnexcm(fcng1,'CALL FCN',arglis,1,ierflg,chi2)
c.....simplex fit :

c      goto 100

      arglis(1)=5000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)

      arglis(1)=20000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
      arglis(1)=28000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
 100  continue

c      arglis(1)=500.
c      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
c      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)

cc.....last call :
      arglis(1) = 3             !   iflag = 3
      call mnexcm (fcng1, 'call fcn', arglis, 1, ierflg,chi2)
c.....stop :
      call mnexcm (fcng1,'stop',arglis,1,ierflg,chi2)
      a1=parsal(1)
      a2=parsal(2)
      a3=parsal(3)
      a4=parsal(4)
      a5=parsal(5)
      a6=parsal(6)
      a7=parsal(7)
      a8=parsal(8)
 1200 format(4F8.4)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine fcng1 (npar, g, f, x, iflag ,chi2)
      implicit double precision (a-h, o-z)
      dimension x(*), g(*)
      external chi2
      f = chi2 (x)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function chi2(param)
      implicit double precision (a-h,o-z)
      dimension param(8),xx1(74),xx2(228)
      dimension parsal(8)
      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      common/expo/aa
      data npoints1/74/
      data npoints2/228/

      data xx1/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,9D-5,
     .     1.D-4,2.D-4,3.D-4,4.D-4,5.D-4,6.D-4,8.D-4,
     .     1.D-3,2.D-3,3D-3,4.D-3,5.D-3,6.D-3,7D-3,8d-3,9d-3,
     .     1.D-2,2.D-2,3D-2,4.D-2,5D-2,
     .     6.D-2,6.5d-2,7D-2,7.5d-2,8.D-2,8.5d-2,9d-2,9.5d-2,
     .     .1D0,.11d0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,
     .     .625d0,.65d0,.675d0,.7d0,.725d0,.75d0,.775d0,.8d0,
     .     .825d0,.85d0,.875d0,.9d0,.92d0,.94d0,.96d0,.98d0,1d0/

      data xx2/1.d-5,2.d-5,3.d-5,4.d-5,5.d-5,6.d-5,7.d-5,8.d-5,9D-5,
     . 1.D-4,1.5D-4,2.D-4,2.5D-4,3.D-4,3.5D-4,4.D-4,4.5D-4,5.D-4,
     . 5.5D-4,6.D-4,6.5D-4,7.D-4,7.5D-4,8.D-4,8.5D-4,9.D-4,9.5D-4,
     .0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,
     .0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,
     .0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.050,
     .0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,
     .0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,
     .0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,   !80
     .0.185,0.19,0.195,0.20,0.205,
     .0.21,0.215,0.22,0.225,0.23,0.235,
     .0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,
     .0.30,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,
     .0.36,0.365,0.37,0.375,0.38,0.385,0.39,
     .0.395,0.40,0.405,0.41,0.415,0.42,0.425,0.43,0.435,0.44,.445,0.45,  !134
     .0.455,0.46,0.465,0.47,0.475,0.48,0.485,0.49,0.495,0.50,0.505,
     .0.51,0.515,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,
     .0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.60,0.605,0.61,0.615,0.62,
     .0.625,0.63,0.635,0.64,0.645,0.65,0.655,0.66,0.665,0.67,0.675,0.68,  !180
     .0.685,0.69,0.70,0.705,0.71,0.715,0.72,0.725,0.73,0.735,0.74,
     .0.745,0.75,0.755,0.76,0.765,0.77,0.78,.785,
     .0.79,0.795,0.80,0.805,0.81,
     .0.815,0.82,0.825,0.83,0.835,
     .0.84,0.845,0.85,
     .0.855,0.86,0.865,0.87,0.875,
     .0.88,0.885,0.89,0.895,0.90,
     .0.91,0.93,0.95,0.97,0.99,1.0/! 228

      a1=(param(1))
      a2=(param(2))
      a3=(param(3))
      a4=(param(4))
      a5=(param(5))
      a6=(param(6))
      a7=(param(7))
      a8=(param(8))
      chi2=0.d0

c.....aa can be changed to try to improve the fit
c.....the common changes it at the same time in the subroutine f0moments
      aa=2.5
c      aa=3.5
      npoints=npoints2
      aextra=1
      if (id.eq.6) then
c          npoints=npoints1
c          aextra=30d0
      endif

c.....take less/more points to improve speed/accuracy
      do i=1,npoints,1
         x=rx**( (1-xx2(i)))

c      if (id.eq.6) x=rx**(1-xx1(i))

         if (x.gt.1d0) goto 133

         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))

         call distrit(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
         chi2=chi2+(f2-f)**2*x**1.2*aextra
 133   continue
      enddo

      do i=1,npoints,1
         x=rx**( (1+xx2(i)))

c      if (id.eq.6) x=rx**(1+xx1(i))

         if (x.lt.1d-5) goto 134
         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))
         call distrit(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
         chi2=chi2+(f2-f)**2*x**1.2
 134   continue
        enddo

C.....return the parameters to save the last set
      do i=1,8
         parsal(i)=param(i)
      enddo
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine distrit(x,q,upv,dnv,usea,dsea,str,chm,bot,glu,ippbar)
C.....here I call the new  sets !!!!! (in program prog_pdf.f)
C.....always gives x*distribution!!!!!
      real*8 upv,dnv,usea,dsea,str,chm,bot,glu,x,q
      common/isetproton/isetproton
      real*8 FX(-5:5),sq2,sx
      sq2=(q*q)
      sx=(x)
      call partons(sq2,sx,fx,5,isetproton,ippbar)
      usea=(fx(-1))*x
      dsea=(fx(-2))*x
      str=(fx(-3))*x
      chm=(fx(-4))*x
      bot=(fx(-5))*x
      glu=(fx(0))*x
      upv=((fx(1))*x-usea )
      dnv=((fx(2))*x-dsea )
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



c New routine (2013 - Kora)
      subroutine writepdfout(filename,xx,muf,energy_sector,pdf_label)

      IMPLICIT NONE

      include 'pdffit.h'

      character, intent(in) :: filename*50
      double precision, intent(in) :: xx, muf
      integer, intent(in) :: energy_sector, pdf_label

      integer ii

      open (unit=21, file=filename, status='unknown')

      ii = energy_sector

      write(21,*) 'x_max: ', xx
      write(21,*) 'mu_F: ', muf
      write(21,*) 'pdf_label:', pdf_label
      write(21,*)

      write(21,*) 'Parameters:'
      write(21,*) A1UV(1:NFITMAX,ii)
      write(21,*) A2UV(1:NFITMAX,ii)
      write(21,*) A3UV(1:NFITMAX,ii)
      write(21,*) A4UV(1:NFITMAX,ii)
      write(21,*) A5UV(1:NFITMAX,ii)
      write(21,*) A6UV(1:NFITMAX,ii)
      write(21,*) A7UV(1:NFITMAX,ii)
      write(21,*) A8UV(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1DV(1:NFITMAX,ii)
      write(21,*) A2DV(1:NFITMAX,ii)
      write(21,*) A3DV(1:NFITMAX,ii)
      write(21,*) A4DV(1:NFITMAX,ii)
      write(21,*) A5DV(1:NFITMAX,ii)
      write(21,*) A6DV(1:NFITMAX,ii)
      write(21,*) A7DV(1:NFITMAX,ii)
      write(21,*) A8DV(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1US(1:NFITMAX,ii)
      write(21,*) A2US(1:NFITMAX,ii)
      write(21,*) A3US(1:NFITMAX,ii)
      write(21,*) A4US(1:NFITMAX,ii)
      write(21,*) A5US(1:NFITMAX,ii)
      write(21,*) A6US(1:NFITMAX,ii)
      write(21,*) A7US(1:NFITMAX,ii)
      write(21,*) A8US(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1DS(1:NFITMAX,ii)
      write(21,*) A2DS(1:NFITMAX,ii)
      write(21,*) A3DS(1:NFITMAX,ii)
      write(21,*) A4DS(1:NFITMAX,ii)
      write(21,*) A5DS(1:NFITMAX,ii)
      write(21,*) A6DS(1:NFITMAX,ii)
      write(21,*) A7DS(1:NFITMAX,ii)
      write(21,*) A8DS(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1SS(1:NFITMAX,ii)
      write(21,*) A2SS(1:NFITMAX,ii)
      write(21,*) A3SS(1:NFITMAX,ii)
      write(21,*) A4SS(1:NFITMAX,ii)
      write(21,*) A5SS(1:NFITMAX,ii)
      write(21,*) A6SS(1:NFITMAX,ii)
      write(21,*) A7SS(1:NFITMAX,ii)
      write(21,*) A8SS(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1GL(1:NFITMAX,ii)
      write(21,*) A2GL(1:NFITMAX,ii)
      write(21,*) A3GL(1:NFITMAX,ii)
      write(21,*) A4GL(1:NFITMAX,ii)
      write(21,*) A5GL(1:NFITMAX,ii)
      write(21,*) A6GL(1:NFITMAX,ii)
      write(21,*) A7GL(1:NFITMAX,ii)
      write(21,*) A8GL(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1CH(1:NFITMAX,ii)
      write(21,*) A2CH(1:NFITMAX,ii)
      write(21,*) A3CH(1:NFITMAX,ii)
      write(21,*) A4CH(1:NFITMAX,ii)
      write(21,*) A5CH(1:NFITMAX,ii)
      write(21,*) A6CH(1:NFITMAX,ii)
      write(21,*) A7CH(1:NFITMAX,ii)
      write(21,*) A8CH(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1BO(1:NFITMAX,ii)
      write(21,*) A2BO(1:NFITMAX,ii)
      write(21,*) A3BO(1:NFITMAX,ii)
      write(21,*) A4BO(1:NFITMAX,ii)
      write(21,*) A5BO(1:NFITMAX,ii)
      write(21,*) A6BO(1:NFITMAX,ii)
      write(21,*) A7BO(1:NFITMAX,ii)
      write(21,*) A8BO(1:NFITMAX,ii)
      write(21,*)
      write(21,*)


      write(21,*) A1UVp(1:NFITMAX,ii)
      write(21,*) A2UVp(1:NFITMAX,ii)
      write(21,*) A3UVp(1:NFITMAX,ii)
      write(21,*) A4UVp(1:NFITMAX,ii)
      write(21,*) A5UVp(1:NFITMAX,ii)
      write(21,*) A6UVp(1:NFITMAX,ii)
      write(21,*) A7UVp(1:NFITMAX,ii)
      write(21,*) A8UVp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1DVp(1:NFITMAX,ii)
      write(21,*) A2DVp(1:NFITMAX,ii)
      write(21,*) A3DVp(1:NFITMAX,ii)
      write(21,*) A4DVp(1:NFITMAX,ii)
      write(21,*) A5DVp(1:NFITMAX,ii)
      write(21,*) A6DVp(1:NFITMAX,ii)
      write(21,*) A7DVp(1:NFITMAX,ii)
      write(21,*) A8DVp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1USp(1:NFITMAX,ii)
      write(21,*) A2USp(1:NFITMAX,ii)
      write(21,*) A3USp(1:NFITMAX,ii)
      write(21,*) A4USp(1:NFITMAX,ii)
      write(21,*) A5USp(1:NFITMAX,ii)
      write(21,*) A6USp(1:NFITMAX,ii)
      write(21,*) A7USp(1:NFITMAX,ii)
      write(21,*) A8USp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1DSp(1:NFITMAX,ii)
      write(21,*) A2DSp(1:NFITMAX,ii)
      write(21,*) A3DSp(1:NFITMAX,ii)
      write(21,*) A4DSp(1:NFITMAX,ii)
      write(21,*) A5DSp(1:NFITMAX,ii)
      write(21,*) A6DSp(1:NFITMAX,ii)
      write(21,*) A7DSp(1:NFITMAX,ii)
      write(21,*) A8DSp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1SSp(1:NFITMAX,ii)
      write(21,*) A2SSp(1:NFITMAX,ii)
      write(21,*) A3SSp(1:NFITMAX,ii)
      write(21,*) A4SSp(1:NFITMAX,ii)
      write(21,*) A5SSp(1:NFITMAX,ii)
      write(21,*) A6SSp(1:NFITMAX,ii)
      write(21,*) A7SSp(1:NFITMAX,ii)
      write(21,*) A8SSp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1GLp(1:NFITMAX,ii)
      write(21,*) A2GLp(1:NFITMAX,ii)
      write(21,*) A3GLp(1:NFITMAX,ii)
      write(21,*) A4GLp(1:NFITMAX,ii)
      write(21,*) A5GLp(1:NFITMAX,ii)
      write(21,*) A6GLp(1:NFITMAX,ii)
      write(21,*) A7GLp(1:NFITMAX,ii)
      write(21,*) A8GLp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1CHp(1:NFITMAX,ii)
      write(21,*) A2CHp(1:NFITMAX,ii)
      write(21,*) A3CHp(1:NFITMAX,ii)
      write(21,*) A4CHp(1:NFITMAX,ii)
      write(21,*) A5CHp(1:NFITMAX,ii)
      write(21,*) A6CHp(1:NFITMAX,ii)
      write(21,*) A7CHp(1:NFITMAX,ii)
      write(21,*) A8CHp(1:NFITMAX,ii)
      write(21,*)

      write(21,*) A1BOp(1:NFITMAX,ii)
      write(21,*) A2BOp(1:NFITMAX,ii)
      write(21,*) A3BOp(1:NFITMAX,ii)
      write(21,*) A4BOp(1:NFITMAX,ii)
      write(21,*) A5BOp(1:NFITMAX,ii)
      write(21,*) A6BOp(1:NFITMAX,ii)
      write(21,*) A7BOp(1:NFITMAX,ii)
      write(21,*) A8BOp(1:NFITMAX,ii)

      close(unit=21)

      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c New routine (2013 - Kora)

      subroutine readpdfout(filename,energy_sector)
      IMPLICIT NONE
c
      include 'pdffit.h'
c
      character, intent(in) :: filename*50
      integer, intent(in) :: energy_sector
c
      integer ii
c
c
      aa = 2.5D0
      NFITMAX = 14
c
      ii = energy_sector

      open (unit=21, file=filename, status='old')

      read(21,*)
      read(21,*)
      read(21,*)
      read(21,*)
      read(21,*)

      read(21,*) A1UV(1:NFITMAX,ii)
      read(21,*) A2UV(1:NFITMAX,ii)
      read(21,*) A3UV(1:NFITMAX,ii)
      read(21,*) A4UV(1:NFITMAX,ii)
      read(21,*) A5UV(1:NFITMAX,ii)
      read(21,*) A6UV(1:NFITMAX,ii)
      read(21,*) A7UV(1:NFITMAX,ii)
      read(21,*) A8UV(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1DV(1:NFITMAX,ii)
      read(21,*) A2DV(1:NFITMAX,ii)
      read(21,*) A3DV(1:NFITMAX,ii)
      read(21,*) A4DV(1:NFITMAX,ii)
      read(21,*) A5DV(1:NFITMAX,ii)
      read(21,*) A6DV(1:NFITMAX,ii)
      read(21,*) A7DV(1:NFITMAX,ii)
      read(21,*) A8DV(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1US(1:NFITMAX,ii)
      read(21,*) A2US(1:NFITMAX,ii)
      read(21,*) A3US(1:NFITMAX,ii)
      read(21,*) A4US(1:NFITMAX,ii)
      read(21,*) A5US(1:NFITMAX,ii)
      read(21,*) A6US(1:NFITMAX,ii)
      read(21,*) A7US(1:NFITMAX,ii)
      read(21,*) A8US(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1DS(1:NFITMAX,ii)
      read(21,*) A2DS(1:NFITMAX,ii)
      read(21,*) A3DS(1:NFITMAX,ii)
      read(21,*) A4DS(1:NFITMAX,ii)
      read(21,*) A5DS(1:NFITMAX,ii)
      read(21,*) A6DS(1:NFITMAX,ii)
      read(21,*) A7DS(1:NFITMAX,ii)
      read(21,*) A8DS(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1SS(1:NFITMAX,ii)
      read(21,*) A2SS(1:NFITMAX,ii)
      read(21,*) A3SS(1:NFITMAX,ii)
      read(21,*) A4SS(1:NFITMAX,ii)
      read(21,*) A5SS(1:NFITMAX,ii)
      read(21,*) A6SS(1:NFITMAX,ii)
      read(21,*) A7SS(1:NFITMAX,ii)
      read(21,*) A8SS(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1GL(1:NFITMAX,ii)
      read(21,*) A2GL(1:NFITMAX,ii)
      read(21,*) A3GL(1:NFITMAX,ii)
      read(21,*) A4GL(1:NFITMAX,ii)
      read(21,*) A5GL(1:NFITMAX,ii)
      read(21,*) A6GL(1:NFITMAX,ii)
      read(21,*) A7GL(1:NFITMAX,ii)
      read(21,*) A8GL(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1CH(1:NFITMAX,ii)
      read(21,*) A2CH(1:NFITMAX,ii)
      read(21,*) A3CH(1:NFITMAX,ii)
      read(21,*) A4CH(1:NFITMAX,ii)
      read(21,*) A5CH(1:NFITMAX,ii)
      read(21,*) A6CH(1:NFITMAX,ii)
      read(21,*) A7CH(1:NFITMAX,ii)
      read(21,*) A8CH(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1BO(1:NFITMAX,ii)
      read(21,*) A2BO(1:NFITMAX,ii)
      read(21,*) A3BO(1:NFITMAX,ii)
      read(21,*) A4BO(1:NFITMAX,ii)
      read(21,*) A5BO(1:NFITMAX,ii)
      read(21,*) A6BO(1:NFITMAX,ii)
      read(21,*) A7BO(1:NFITMAX,ii)
      read(21,*) A8BO(1:NFITMAX,ii)
      read(21,*)
      read(21,*)


      read(21,*) A1UVp(1:NFITMAX,ii)
      read(21,*) A2UVp(1:NFITMAX,ii)
      read(21,*) A3UVp(1:NFITMAX,ii)
      read(21,*) A4UVp(1:NFITMAX,ii)
      read(21,*) A5UVp(1:NFITMAX,ii)
      read(21,*) A6UVp(1:NFITMAX,ii)
      read(21,*) A7UVp(1:NFITMAX,ii)
      read(21,*) A8UVp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1DVp(1:NFITMAX,ii)
      read(21,*) A2DVp(1:NFITMAX,ii)
      read(21,*) A3DVp(1:NFITMAX,ii)
      read(21,*) A4DVp(1:NFITMAX,ii)
      read(21,*) A5DVp(1:NFITMAX,ii)
      read(21,*) A6DVp(1:NFITMAX,ii)
      read(21,*) A7DVp(1:NFITMAX,ii)
      read(21,*) A8DVp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1USp(1:NFITMAX,ii)
      read(21,*) A2USp(1:NFITMAX,ii)
      read(21,*) A3USp(1:NFITMAX,ii)
      read(21,*) A4USp(1:NFITMAX,ii)
      read(21,*) A5USp(1:NFITMAX,ii)
      read(21,*) A6USp(1:NFITMAX,ii)
      read(21,*) A7USp(1:NFITMAX,ii)
      read(21,*) A8USp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1DSp(1:NFITMAX,ii)
      read(21,*) A2DSp(1:NFITMAX,ii)
      read(21,*) A3DSp(1:NFITMAX,ii)
      read(21,*) A4DSp(1:NFITMAX,ii)
      read(21,*) A5DSp(1:NFITMAX,ii)
      read(21,*) A6DSp(1:NFITMAX,ii)
      read(21,*) A7DSp(1:NFITMAX,ii)
      read(21,*) A8DSp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1SSp(1:NFITMAX,ii)
      read(21,*) A2SSp(1:NFITMAX,ii)
      read(21,*) A3SSp(1:NFITMAX,ii)
      read(21,*) A4SSp(1:NFITMAX,ii)
      read(21,*) A5SSp(1:NFITMAX,ii)
      read(21,*) A6SSp(1:NFITMAX,ii)
      read(21,*) A7SSp(1:NFITMAX,ii)
      read(21,*) A8SSp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1GLp(1:NFITMAX,ii)
      read(21,*) A2GLp(1:NFITMAX,ii)
      read(21,*) A3GLp(1:NFITMAX,ii)
      read(21,*) A4GLp(1:NFITMAX,ii)
      read(21,*) A5GLp(1:NFITMAX,ii)
      read(21,*) A6GLp(1:NFITMAX,ii)
      read(21,*) A7GLp(1:NFITMAX,ii)
      read(21,*) A8GLp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1CHp(1:NFITMAX,ii)
      read(21,*) A2CHp(1:NFITMAX,ii)
      read(21,*) A3CHp(1:NFITMAX,ii)
      read(21,*) A4CHp(1:NFITMAX,ii)
      read(21,*) A5CHp(1:NFITMAX,ii)
      read(21,*) A6CHp(1:NFITMAX,ii)
      read(21,*) A7CHp(1:NFITMAX,ii)
      read(21,*) A8CHp(1:NFITMAX,ii)
      read(21,*)

      read(21,*) A1BOp(1:NFITMAX,ii)
      read(21,*) A2BOp(1:NFITMAX,ii)
      read(21,*) A3BOp(1:NFITMAX,ii)
      read(21,*) A4BOp(1:NFITMAX,ii)
      read(21,*) A5BOp(1:NFITMAX,ii)
      read(21,*) A6BOp(1:NFITMAX,ii)
      read(21,*) A7BOp(1:NFITMAX,ii)
      read(21,*) A8BOp(1:NFITMAX,ii)

      close(unit=21)

      end
