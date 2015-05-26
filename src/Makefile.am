

# http://www.jezuk.co.uk/cgi-bin/view/articles/autoconfiscation-part-two
# https://www.lrde.epita.fr/~adl/dl/autotools-handout-4.pdf


AUTOMAKE_OPTIONS = foreign



lib_LTLIBRARIES = libutils.la
lib_LIBRARIES = libutils.a


libutils_la_SOURCES = atomsystem-box.cpp \
atomsystem.cpp \
atomsystem-dist.cpp \
atomsystem-frame.cpp \
atomsystem-nbor.cpp \
bond.cpp \
funcfit-basics.cpp \
mtwister.cpp \
omp-basics.cpp \
param.cpp \
utils-errors.cpp \
utils-rototrans.cpp \
utils-streamio.cpp \
utils-string.cpp \
utils-vecalg.cpp


libutils_a_SOURCES = atomsystem-box.cpp \
atomsystem.cpp \
atomsystem-dist.cpp \
atomsystem-frame.cpp \
atomsystem-nbor.cpp \
bond.cpp \
funcfit-basics.cpp \
mtwister.cpp \
omp-basics.cpp \
param.cpp \
utils-errors.cpp \
utils-rototrans.cpp \
utils-streamio.cpp \
utils-string.cpp \
utils-vecalg.cpp


include_HEADERS = atomsystem.hpp \
bond.hpp \
chisq-basics.hpp \
constants.hpp \
funcfit-basics.hpp \
funcfit-beecolony.hpp \
funcfit-conjgrad.hpp \
funcfit-deriv.hpp \
funcfit-diffevol.hpp \
funcfit-errors.hpp \
funcfit-gravsearch.hpp \
funcfit-gauss-newton.hpp \
funcfit-leve-marq.hpp \
funcfit-powelldogleg.hpp \
funcfit-moldyn.hpp \
funcfit-partswarm.hpp \
funcfit-powell.hpp \
funcfit-simann.hpp \
funcfit-simplexfit.hpp \
mtwister.hpp \
nr-f1dim.hpp \
nr-golden.hpp \
nr-linemethod.hpp \
omp-basics.hpp \
param.hpp \
utils-errors.hpp \
utils.hpp \
utils-math.hpp \
utils-matrix-Choleskydecomp.hpp \
utils-matrix.hpp \
utils-matrixsq3.hpp \
utils-matrix3.hpp \
utils-matrix-LUdecomp.hpp \
utils-matrix-QRdecomp.hpp \
utils-rototrans.hpp \
utils-streamio.hpp \
utils-string.hpp \
utils-vecalg.hpp \
utils-vector3.hpp \
utils-vector.hpp


ma_inc  = -I$(srcdir)
ma_warn = -Wall -Wextra -Wstrict-aliasing
ma_std  = -ansi -pedantic -std=c++98
ma_opt  = -O3
ma_dbg  = -g
ma_omp  = -fopenmp


libutils_la_CPPFLAGS = $(ma_inc) $(ma_warn) $(ma_std) $(ma_opt) $(ma_dbg) $(ma_omp) -fpic
libutils_a_CPPFLAGS  = $(ma_inc) $(ma_warn) $(ma_std) $(ma_opt) $(ma_dbg) $(ma_omp)

libutils_la_LDFLAGS = $(ma_omp) -rpath '$(libdir)' -lrt -lm -shared 
libutils_a_AR = ar rcsvf


# foo_LDADD if foo is a program or a shared library
# foo_LIBADD if foo is a static library




