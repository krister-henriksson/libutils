
.PHONY: all dynamic dynamic2 static static2 all clean


include vars.mk


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Settings:

CC      =  g++ 
WARN    = -Wall -Wextra  #-Wstrict-aliasing
STD     = -ansi -pedantic -std=c++98 
DEBUG   = -g 
OPT     = -O2 #-fstrict-aliasing 
OPENMP  = -fopenmp 
CFLAGS  =  -c $(WARN) $(STD) $(OPT)  $(DEBUG)  # $(INC)

INC     = -I$(INCDIR)

LIB     = -L$(LIBDIR)
LIBFILE = $(LIBDIR)/libutils

LDFLAGS        =  -lrt -lm 
LDFLAGS_STATIC =  $(LDFLAGS) -static 


SRC = atomsystem.cpp exiterrors.cpp param.cpp mtwister.cpp \
	utils-string.cpp utils-streamio.cpp funcfit-basics.cpp \
	omp-basics.cpp

SOURCES = $(addprefix src/,$(SRC))
OBJECTS = $(SOURCES:src/%.cpp=obj/%.o)


#OBJECTS = $(subst src/,obj/,$(TMP))


TARGET_DYNAMIC = $(LIBFILE).so
TARGET_STATIC  = $(LIBFILE).a



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Rules:

all: dirs dynamic2 static2
	@echo "export LIBDIR="$(LIBDIR)"" > bashsettings.text


dynamic: dynamic2
	# strip -s $(LIBFILE).so
	chmod a+rx $(LIBFILE).so

dynamic2: $(OBJECTS)
	$(CC) $(STD) $(WARN) $(DEBUG) $(OPT) $(OPENMP) -shared -o $(LIBFILE).so  $(OBJECTS)  $(LINK)


static: static2
	# strip -s $(LIBFILE).a
	chmod a+rx $(LIBFILE).a

static2: $(OBJECTS)
	ar rcsvf $(LIBFILE).a  $(OBJECTS)


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Pattern rules:

obj/%.o: src/%.cpp
	$(CC) $(STD) $(WARN) $(DEBUG) $(OPT) $(INC) $(OPENMP) -c -fpic $< -o $@


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Phony rules:

dirs:
	-mkdir obj
	-mkdir -p $(INCDIR)
	cp src/*.hpp $(INCDIR)/
	-mkdir -p $(LIBDIR)

clean:
	-rm obj/*.o $(TARGET_DYNAMIC) $(TARGET_STATIC)




# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.


