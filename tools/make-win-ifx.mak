
TARGET=gebt.exe

SRCDIR=src
SRCF=\
	ddep.f \
	ma28.f \
	mc19.f \
	blas.f \
	lapack.f \
	dnaupd.f \
	dneupd.f
SRCF90=\
	CPUtime.f90 \
	GlobalDataFun.f90 \
	TimeFunction.f90 \
	PrescribedCondition.f90 \
	InternalData.f90 \
	Preprocess.f90 \
	Element.f90 \
	Member.f90 \
	System.f90 \
	Solve.f90 \
	EigenSolve.f90 \
	Analysis.f90 \
	IO.f90 \
	main.f90

OBJF = $(SRCF:.f=.obj)
OBJF90 = $(SRCF90:.f90=.obj)

OBJS = $(OBJF) $(OBJF90)


COMPILE  = ifx
LINK     = ifx
cFLAGS   = /O3 /check:bounds /c
lFLAGS   = /O3 /check:bounds 
INSTALL_DIR=bin

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LINK) $(lFLAGS) $(OBJS) /out:$(TARGET)

# ddep.obj: $(SRCDIR)/ddep.f
# 	$(COMPILE) $(cFLAGS) $(?) /Fo $@

%.obj: $(SRCDIR)/%.f
	$(COMPILE) $(cFLAGS) $(?) /Fo $@

%.obj: $(SRCDIR)/%.f90
	$(COMPILE) $(cFLAGS) $(?) /Fo $@

install: $(TARGET)
	if not exist "$(INSTALL_DIR)" mkdir "$(INSTALL_DIR)"
    copy $(TARGET) "$(INSTALL_DIR)/$(TARGET)"

#	copy $< ..\gebt

clean:
	del /q $(TARGET) $(OBJS) 

print:
    @echo $(SRCF)
    @echo $(SRCF90)
    @echo $(OBJF)
    @echo $(OBJF90)
    @echo $(OBJS)

.PHONY: print all clean install
