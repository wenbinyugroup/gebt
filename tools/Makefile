
include gebtMake # Change this for specific program

OBJS     = $(OBJSc)  $(OBJS90)  $(OBJSF)

COMPILE  = gfortran
COMPILEc = gcc
LINK     = gfortran
cFLAGS   = -O3 -fbounds-check  -c
cFLAGSc  = -O3    -c
lFLAGS   = -O3   -fbounds-check 
SRCDIR   =../src/

$(PGM):  $(OBJS)   
	$(LINK) $(lFLAGS) -o $(PGM)   $(OBJS)

$(OBJS90): %.o: $(SRCDIR)%.f90
	$(COMPILE) $(cFLAGS) $<

$(OBJSF): %.o: $(SRCDIR)%.f
	$(COMPILE) $(cFLAGS) $<

$(OBJSc): %.o: $(SRCDIR)%.c
	$(COMPILEc) $(cFLAGSc) $<

install: $(PGM)
	cp $< ../gebt

.PHONY: clean
clean :
	rm $(PGM) $(OBJS) 
