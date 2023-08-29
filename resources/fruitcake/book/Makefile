HOMEBOOK = .

DIRPET = pet
DIRSPECT = spect
DIRUTIL = util


SRCDIR = $(HOMEBOOK)/src
SRCPET = $(SRCDIR)/$(DIRPET)
SRCSPECT = $(SRCDIR)/$(DIRSPECT)
SRCUTIL = $(SRCDIR)/$(DIRUTIL)
OBJDIR = $(HOMEBOOK)/objects
INCDIR = $(HOMEBOOK)/include
INCPET = $(INCDIR)/$(DIRPET)
INCSPECT = $(INCDIR)/$(DIRSPECT)
INCUTIL = $(INCDIR)/$(DIRUTIL)
LIBDIR = $(HOMEBOOK)/lib
DEPDIR = $(HOMEBOOK)/depend
DEPFILE = $(HOMEBOOK)/Makefile.d

CC = @gcc 
AR = @ar
MAKE = @make

LIBS = -lm

CFLAGS =  -Wall -fPIC -I$(INCDIR) -g 
LDFLAGS = -Wall  -L$(LIBDIR) $(LIBS)  -g 

vpath %.o   $(OBJDIR)
vpath %.h   $(INCDIR):$(INCPET):$(INCSPECT):$(INCUTIL)
vpath %.c   $(SRCDIR):$(SRCPET):$(SRCSPECT):$(SRCUTIL)
vpath %.cc  $(SRCDIR):$(SRCPET):$(SRCSPECT):$(SRCUTIL)
vpath %.so  $(LIBDIR)
vpath %.a   $(LIBDIR)

objs_pet := $(patsubst $(SRCPET)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCPET)/*.c))
objs_spect := $(patsubst $(SRCSPECT)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCSPECT)/*.c))
objs_util := $(patsubst $(SRCUTIL)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCUTIL)/*.c))

srcs_pet := $(wildcard $(SRCPET)/*.c)
srcs_spect := $(wildcard $(SRCSPECT)/*.c)
srcs_util := $(wildcard $(SRCUTIL)/*.c)

deps 	:= $(patsubst $(SRCPET)/%.c,$(DEPDIR)/%.d,$(srcs_pet)) $(patsubst $(SRCSPECT)/%.c,$(DEPDIR)/%.d,$(srcs_spect)) $(patsubst $(SRCUTIL)/%.c,$(DEPDIR)/%.d,$(srcs_util)) 
headers := $(wildcard $(INCDIR)/*.h)

$(OBJDIR)/%.o : $(SRCPET)/%.c
	@echo Compilando $< ...
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o : $(SRCSPECT)/%.c
	@echo Compilando $< ...
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o : $(SRCUTIL)/%.c
	@echo Compilando $< ...
	$(CC) -c $(CFLAGS) -o $@ $<

all: libubpet.so libubpet.a libubspect.so libubspect.a libubutil.so libubutil.a

libubpet.so: $(objs_pet) 
	@echo Creando $@...	
	$(CC) -shared  $^ -o $(LIBDIR)/libubpet.so

libubpet.a: $(objs_pet)
	@echo Creando $@...	
	$(AR) -rs $(LIBDIR)/libubpet.a $^ 

libubspect.so: $(objs_spect) 
	@echo Creando $@...	
	$(CC) -shared  $^ -o $(LIBDIR)/libubspect.so

libubspect.a: $(objs_spect)
	@echo Creando $@...	
	$(AR) -rs $(LIBDIR)/libubspect.a $^ 

libubutil.so: $(objs_util) 
	@echo Creando $@...	
	$(CC) -shared  $^ -o $(LIBDIR)/libubutil.so

libubutil.a: $(objs_util)
	@echo Creando $@...	
	$(AR) -rs $(LIBDIR)/libubutil.a $^ 

.PHONY: clean
clean: 
	rm -f $(DEPDIR)/*.d
	rm -f $(SRCPET)/*~  $(SRCSPECT)/*~ $(SRCUTIL)/*~ $(OBJDIR)/*.o  $(LIBDIR)/*.so  $(LIBDIR)/*.a *~
	rm -f $(INCPET)/*~  $(INCSPECT)/*~ $(INCUTIL)/*~ 


.PHONY: dep
dep: borrarDepFiles $(deps)
	@echo  Generando fichero $(DEPFILE) ...
	@cat $(deps) > $(DEPFILE)

.PHONY: borrarDepFiles
borrarDepFiles: 
	@echo Borrando dependencias ...
	rm -f $(DEPFILE)
	rm -f $(DEPDIR)/*.d

$(DEPDIR)/%.d: $(SRCPET)/%.c
	@echo Generando dependencias de $< ...
	$(CC) -w -MM -I$(INCDIR) $< -MF $@ -MT $(patsubst $(SRCPET)/%.c,$(OBJDIR)/%.o,./$<)

$(DEPDIR)/%.d: $(SRCSPECT)/%.c
	@echo Generando dependencias de $< ...
	$(CC) -w -MM -I$(INCDIR) $< -MF $@ -MT $(patsubst $(SRCSPECT)/%.c,$(OBJDIR)/%.o,./$<)

$(DEPDIR)/%.d: $(SRCUTIL)/%.c
	@echo Generando dependencias de $< ...
	$(CC) -w -MM -I$(INCDIR) $< -MF $@ -MT $(patsubst $(SRCUTIL)/%.c,$(OBJDIR)/%.o,./$<)

-include $(DEPFILE)
