c
c Common blocks for pdffit
c
c Beam 1 fit parameters
       double precision A1UV(30,20),A2UV(30,20),A3UV(30,20),A4UV(30,20),
     .       A5UV(30,20),A6UV(30,20),A7UV(30,20),A8UV(30,20)
       double precision A1DV(30,20),A2DV(30,20),A3DV(30,20),A4DV(30,20),
     .       A5DV(30,20),A6DV(30,20),A7DV(30,20),A8DV(30,20)
       double precision A1US(30,20),A2US(30,20),A3US(30,20),A4US(30,20),
     .       A5US(30,20),A6US(30,20),A7US(30,20),A8US(30,20)
       double precision A1DS(30,20),A2DS(30,20),A3DS(30,20),A4DS(30,20),
     .       A5DS(30,20),A6DS(30,20),A7DS(30,20),A8DS(30,20)
       double precision A1SS(30,20),A2SS(30,20),A3SS(30,20),A4SS(30,20),
     .       A5SS(30,20),A6SS(30,20),A7SS(30,20),A8SS(30,20)
       double precision A1GL(30,20),A2GL(30,20),A3GL(30,20),A4GL(30,20),
     .       A5GL(30,20),A6GL(30,20),A7GL(30,20),A8GL(30,20)
       double precision A1CH(30,20),A2CH(30,20),A3CH(30,20),A4CH(30,20),
     .       A5CH(30,20),A6CH(30,20),A7CH(30,20),A8CH(30,20)
       double precision A1BO(30,20),A2BO(30,20),A3BO(30,20),A4BO(30,20),
     .       A5BO(30,20),A6BO(30,20),A7BO(30,20),A8BO(30,20)
c
       common/ CUV1/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
       common/ CDV1/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
       common/ CUS1/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
       common/ CDS1/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
       common/ CSS1/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
       common/ CGL1/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
       common/ CCH1/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
       common/ CBO1/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
c
       double precision muf_array
       common/ muf_array/ muf_array(20)
c
c Beam 2 fit parameters
       double precision A1UVp(30,20),A2UVp(30,20),A3UVp(30,20),
     .  A4UVp(30,20),A5UVp(30,20),A6UVp(30,20),A7UVp(30,20),A8UVp(30,20)
       double precision A1DVp(30,20),A2DVp(30,20),A3DVp(30,20),
     .  A4DVp(30,20),A5DVp(30,20),A6DVp(30,20),A7DVp(30,20),A8DVp(30,20)
       double precision A1USp(30,20),A2USp(30,20),A3USp(30,20),
     .  A4USp(30,20),A5USp(30,20),A6USp(30,20),A7USp(30,20),A8USp(30,20)
       double precision A1DSp(30,20),A2DSp(30,20),A3DSp(30,20),
     .  A4DSp(30,20),A5DSp(30,20),A6DSp(30,20),A7DSp(30,20),A8DSp(30,20)
       double precision A1SSp(30,20),A2SSp(30,20),A3SSp(30,20),
     .  A4SSp(30,20),A5SSp(30,20),A6SSp(30,20),A7SSp(30,20),A8SSp(30,20)
       double precision A1GLp(30,20),A2GLp(30,20),A3GLp(30,20),
     .  A4GLp(30,20),A5GLp(30,20),A6GLp(30,20),A7GLp(30,20),A8GLp(30,20)
       double precision A1CHp(30,20),A2CHp(30,20),A3CHp(30,20),
     .  A4CHp(30,20),A5CHp(30,20),A6CHp(30,20),A7CHp(30,20),A8CHp(30,20)
       double precision A1BOp(30,20),A2BOp(30,20),A3BOp(30,20),
     .  A4BOp(30,20),A5BOp(30,20),A6BOp(30,20),A7BOp(30,20),A8BOp(30,20)
c
       common/ CUV2/ A1UVp,A2UVp,A3UVp,A4UVp,A5UVp,A6UVp,A7UVp,A8UVp
       common/ CDV2/ A1DVp,A2DVp,A3DVp,A4DVp,A5DVp,A6DVp,A7DVp,A8DVp
       common/ CUS2/ A1USp,A2USp,A3USp,A4USp,A5USp,A6USp,A7USp,A8USp
       common/ CDS2/ A1DSp,A2DSp,A3DSp,A4DSp,A5DSp,A6DSp,A7DSp,A8DSp
       common/ CSS2/ A1SSp,A2SSp,A3SSp,A4SSp,A5SSp,A6SSp,A7SSp,A8SSp
       common/ CGL2/ A1GLp,A2GLp,A3GLp,A4GLp,A5GLp,A6GLp,A7GLp,A8GLp
       common/ CCH2/ A1CHp,A2CHp,A3CHp,A4CHp,A5CHp,A6CHp,A7CHp,A8CHp
       common/ CBO2/ A1BOp,A2BOp,A3BOp,A4BOp,A5BOp,A6BOp,A7BOp,A8BOp
c
c Additional fit parameters
       double precision aa
       common/expo/aa
       INTEGER NFITMAX, N_ENERGY_SECS
       common/NFITMAX/NFITMAX, N_ENERGY_SECS
       integer isetproton
       common/isetproton/isetproton
