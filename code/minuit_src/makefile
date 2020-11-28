# Compiler
FF=gfortran
FFLAGS = -fno-automatic -fno-f2c -O3
SRC= intrac.F mncuve.F mngrad.F mnmatu.F mnpsdf.F mnstat.F \
minuit.F mnderi.F mnhelp.F mnmigr.F mnrazz.F mnstin.F \
mnamin.F mndxdi.F mnhes1.F mnmnos.F mnread.F mntiny.F \
mnbins.F mneig.F mnhess.F mnmnot.F mnrn15.F mnunpt.F \
mncalf.F mnemat.F mnimpr.F mnparm.F mnrset.F mnvers.F \
mncler.F mnerrs.F mninex.F mnpars.F mnsave.F mnvert.F \
mncntr.F mneval.F mninit.F mnpfit.F mnscan.F mnwarn.F \
mncomd.F mnexcm.F mninpu.F mnpint.F mnseek.F mnwerr.F \
mncont.F mnexin.F mnintr.F mnplot.F mnset.F stand.F \
mncrck.F mnfixp.F mnlims.F mnpout.F mnseti.F \
mncros.F mnfree.F mnline.F mnprin.F mnsimp.F

OBJS= intrac.o mncuve.o mngrad.o mnmatu.o mnpsdf.o mnstat.o \
minuit.o mnderi.o mnhelp.o mnmigr.o mnrazz.o mnstin.o \
mnamin.o mndxdi.o mnhes1.o mnmnos.o mnread.o mntiny.o \
mnbins.o mneig.o mnhess.o mnmnot.o mnrn15.o mnunpt.o \
mncalf.o mnemat.o mnimpr.o mnparm.o mnrset.o mnvers.o \
mncler.o mnerrs.o mninex.o mnpars.o mnsave.o mnvert.o \
mncntr.o mneval.o mninit.o mnpfit.o mnscan.o mnwarn.o \
mncomd.o mnexcm.o mninpu.o mnpint.o mnseek.o mnwerr.o \
mncont.o mnexin.o mnintr.o mnplot.o mnset.o stand.o \
mncrck.o mnfixp.o mnlims.o mnpout.o mnseti.o \
mncros.o mnfree.o mnline.o mnprin.o mnsimp.o



# General compilation instruction for .F sources
%.o: %.F
	$(FF) -c -o $@ $< $(FFLAGS)


minuit.a: $(OBJS)
	ar cr minuit.a $(OBJS)

# Clean
clean:
	rm *.o minuit.a
