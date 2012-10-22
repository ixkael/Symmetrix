# ======================================== #

# Directory for DOXYGEN
DOXYGEN=/Applications/Doxygen.app/Contents/Resources/doxygen

# Compiler and options
CC	= gcc
OPT	= -Wall -O3 -g -DSYMTRX_VERSION=\"0.1\" -DSYMTRX_BUILD=\"`git describe`\"
# I MUSTN"T FORGET TO ADD GIT TAGS TO CHANGE THE VERSION!

# ======================================== #

SYMTRXDIR = .
SYMTRXLIB = $(SYMTRXDIR)/lib
SYMTRXINC = $(SYMTRXDIR)/include
SYMTRXBIN = $(SYMTRXDIR)/bin
SYMTRXLIBN= symtrx
SYMTRXSRCMAIN = $(SYMTRXDIR)/src/main/c
SYMTRXSRCTEST = $(SYMTRXDIR)/src/test/c

# ======================================== #

vpath %.c $(SYMTRXSRCMAIN)
vpath %.c $(SYMTRXSRCTEST)
vpath %.h $(SYMTRXINC)

LDFLAGS = -L$(SYMTRXLIB) -l$(SYMTRXLIBN)

FFLAGS  = -I$(SYMTRXINC)

SYMTRXOBJS= $(SYMTRXSRCMAIN)/bisym.o	\
	  $(SYMTRXSRCMAIN)/centrosym.o	\
	  $(SYMTRXSRCMAIN)/miscmath.o	\
	  $(SYMTRXSRCMAIN)/square.o	\
	  $(SYMTRXSRCMAIN)/vector.o

$(SYMTRXSRCMAIN)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(SYMTRXSRCTEST)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

# ======================================== #

.PHONY: default
default: lib test about tidy

.PHONY: all
all: lib doc test about tidy

.PHONY: lib
lib: $(SYMTRXLIB)/lib$(SYMTRXLIBN).a
$(SYMTRXLIB)/lib$(SYMTRXLIBN).a: $(SYMTRXOBJS)
	ar -r $(SYMTRXLIB)/lib$(SYMTRXLIBN).a $(SYMTRXOBJS)

.PHONY: test
test: $(SYMTRXBIN)/symtrx_test
$(SYMTRXBIN)/symtrx_test: $(SYMTRXSRCTEST)/symtrx_test.o $(SYMTRXLIB)/lib$(SYMTRXLIBN).a
	$(CC) $(OPT) $< -o $(SYMTRXBIN)/symtrx_test $(LDFLAGS)
	$(SYMTRXBIN)/symtrx_test

.PHONY: about
about: $(SYMTRXBIN)/about
$(SYMTRXBIN)/about: $(SYMTRXSRCMAIN)/about.o 
	$(CC) $(OPT) $< -o $(SYMTRXBIN)/about
	$(SYMTRXBIN)/about

.PHONY: doc
doc:
	$(DOXYGEN) $(SYMTRXDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(SYMTRXDIR)/doc/html/*

.PHONY: clean
clean:	tidy cleandoc
	rm -f $(SYMTRXLIB)/lib$(SYMTRXLIBN).a
	rm -f $(SYMTRXBIN)/symtrx_test
	rm -f $(SYMTRXBIN)/about

.PHONY: tidy
tidy:
	rm -f $(SYMTRXSRCMAIN)/*.o
	rm -f $(SYMTRXSRCTEST)/*.o
	rm -f *~ 

# ======================================== #