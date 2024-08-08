CC=g++
CPPFLAGS=-O3 -std=c++11
EXEC= abscissas \
 autoCorrelation averageSets collect gaussianKernelClamped \
 gaussianKernel gaussianKernelProj runningAvg \
 skipError splitSets2 splitSets sumSeries sumSets wham2 wham

TARGETS: $(EXEC)

ifeq ($(PREFIX),)
PREFIX=$(HOME)
endif

.PHONY: all

all: $(TARGETS)

$(TARGET):
	@echo $(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

#buildOpenMP: $(TARGET).cpp
#	$(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

.PHONY: clean

clean:
	rm -f $(EXEC)

install: $(EXEC)
	install -d $(PREFIX)/bin/
	install -m 755 $(EXEC) $(PREFIX)/bin/

#test:
#	./test.sh
