CC=g++
CPPFLAGS=-O3 -std=c++11
EXEC= frameSkip framesSeperate move noType oneTypeCount oneType \
oneTypeFrames opticalReduction radiusTube rmFrames \
selectFrames selectRange slice sortFrames sphereMask \
target xyzToBinary

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
