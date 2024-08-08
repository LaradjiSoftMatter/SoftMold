############################################################################
# 'Make all cpp files in this directory'
# Author:  Eric J. Spangler (2021)
# Email: jspngler at memphis dot edu
############################################################################
#  The purpose of this makefile is to compile to executable all cpp files.
#  Makes object files too.
#
###########################################################################
#
#
CC := g++
LDFLAGS := 
CFLAGS_INC := 
CFLAGS := -O3 $(CFLAGS_INC) -std=c++11
#
## ==================- DO NOT CHANGE BELOW THIS LINE ======================
## List of cpp source files.
SRCS := $(wildcard *.cpp)
## List of program targets.
PRGS := $(patsubst %.cpp,%,$(SRCS))
## List of object files.
OBJS := $(patsubst %,%.o,$(PRGS))
##
all : $(PRGS)
##
.SECONDEXPANSION:
OBJ = $(patsubst %,%.o,$@)
	PRG = $@
## Compile the executables
$PRGS : $(OBJS)
	$(CC) $(OBJ)  $(LDFLAGS) -o $(PRG)

## Clean programs and objects
.PHONY: clean
clean:
	$(RM) $(PRGS) $(OBJS)
##
rebuild: clean all
##
## eof Makefile
