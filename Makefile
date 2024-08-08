EXEC= MD
CC=g++
CPPFLAGS=-O3 -fopenmp -std=c++11
SUBDIRS= ./analysis ./analysis/asciiMod ./analysis/xyzMod ./analysis/functions
TARGETS: $(EXEC) all

#To make an intel target, MD built with intel is about 10-20% faster
#CC=icc
#CPPFLAGS=-fast -qopenmp

ifeq ($(PREFIX),)
PREFIX=$(HOME)
endif

cuda:
	nvcc -Xcompiler -fopenmp -std=c++14 MDcuda.cu -O3 --gpu-architecture=compute_70 -o MDcuda

.PHONY: all

all: $(TARGETS)
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir; done

.PHONY: generate

generate:
	@echo $(MAKE) -C ./generate all
	

$(TARGET):
	@echo $(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

.PHONY: clean

clean:
	rm -f MDcuda
	rm -f $(EXEC)
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done

install: $(EXEC)
	install -d $(PREFIX)/bin/
	install -m 755 $(EXEC) $(PREFIX)/bin/
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir install; done
