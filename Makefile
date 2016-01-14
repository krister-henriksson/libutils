

# #######################################################
# Default paths, change values if needed:
# #######################################################
prefix = $(HOME)
INCDIR = $(prefix)/include/libutils
LIBDIR = $(prefix)/lib
# #######################################################
# #######################################################


# #######################################################
# Default settings, change values if needed:
# #######################################################
CC      =  g++ 
WARN    = -Wall -Wextra -Wstrict-aliasing
STD     = -ansi -pedantic -std=c++98 -pg
DEBUG   = -g 
OPT     = -O3 -fstrict-aliasing 
OPENMP  = -fopenmp 
# #######################################################
# #######################################################


CXXFLAGS = -c $(WARN) $(STD) $(OPT)  $(DEBUG)
INC      = -I$(INCDIR)
LIB      = -L$(LIBDIR)
#LIBFILE  = $(LIBDIR)/libutils
LIBFILE  = libutils
LDFLAGS        =  -lrt -lm 
LDFLAGS_STATIC =  $(LDFLAGS) -static 


.PHONY: all default dynamic static clean


# SOURCES = $(addprefix src/,$(SRC))

SRC    := $(wildcard src/*.cpp)
SOURCES = $(SRC)
OBJECTS = $(SOURCES:src/%.cpp=obj/%.o)
DEPS = obj/make.dep

# OBJECTS = $(subst src/,obj/,$(TMP))

TARGET_DYNAMIC = $(LIBFILE).so
TARGET_STATIC  = $(LIBFILE).a
REBUILDABLES   = $(TARGET_DYNAMIC) $(TARGET_STATIC)


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Rules:

default: dirs dynamic static


install: default
	@echo "-----------------------------------------------------------"
	@echo "*** Installing under           :" $(prefix)
	@echo "*** Installing header files in :" $(INCDIR)
	@echo "*** Installing library files in:" $(LIBDIR)
	@echo "-----------------------------------------------------------"
	cp src/*.hpp $(INCDIR)/
	chmod a+rx $(LIBFILE).so
	chmod a+rx $(LIBFILE).a
	cp  $(LIBFILE).so $(LIBFILE).a  $(LIBDIR)/
	@echo "-----------------------------------------------------------"
	@echo '*** If running bash, put these lines into ~/.bashrc :'
	@echo 'export LD_LIBRARY_PATH='$(LIBDIR)':$${LD_LIBRARY_PATH}'
	@echo 'export LD_RUN_PATH='$(LIBDIR)':$${LD_RUN_PATH}'
	@echo "-----------------------------------------------------------"


strip: $(LIBFILE).so $(LIBFILE).a
	strip -g $(LIBFILE).so $(LIBFILE).a
	cp  $(LIBFILE).so $(LIBFILE).a  $(LIBDIR)/


dynamic: $(OBJECTS)
	$(CC) $(STD) $(WARN) $(DEBUG) $(OPT) $(OPENMP) -shared -o $(LIBFILE).so  $(OBJECTS)  $(LINK)

static: $(OBJECTS)
	ar rcsvf $(LIBFILE).a  $(OBJECTS)



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Pattern rules:

obj/%.o: src/%.cpp
	$(CC) $(STD) $(WARN) $(DEBUG) $(OPT) $(INC) $(OPENMP) -c -fpic $< -o $@

$(DEPS): $(SOURCES)
	- mkdir obj
	$(CC) -MM $(INC) $(SOURCES) | sed 's/\(.*\.o\)/obj\/\1/g' > $(DEPS)

include $(DEPS)


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Other rules:

dirs:
	-mkdir -p $(INCDIR)
	-mkdir -p $(LIBDIR)

clean:
	-rm -f obj/*.o $(REBUILDABLES) $(DEPS) $(LIBFILE).so $(LIBFILE).a




# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.


