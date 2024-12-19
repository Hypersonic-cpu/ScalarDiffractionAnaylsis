AS           = $(CROSS_COMPILE)as 
LD           = $(CROSS_COMPILE)ld
CXX          = $(CROSS_COMPILE)clang++
# CPP          = $(CROSS_COMPILE)cpp
AR           = $(CROSS_COMPILE)ar
NM           = $(CROSS_COMPILE)nm
STRIP        = $(CROSS_COMPILE)strip
OBJCOPY      = $(CROSS_COMPILE)objcopy
OBJDUMP      = $(CROSS_COMPILE)objdump

SRC_DIR := src
BUILD_DIR := build

SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cc,$(BUILD_DIR)/%.o,$(SRC_FILES))
TGT_FILE  := $(BUILD_DIR)/main

CFLAGS 		+= -g -O2
CFLAGS 		+= -Wall
CFLAGS 		+= -I/usr/local/Cellar/llvm/19.1.3/include
CFLAGS 		+= -I/usr/local/include
CFLAGS 		+= -I../Symmetric/src
CFLAGS 		+= -std=c++23
LDFLAGS		+= -L/usr/local/lib -lgsl -lgslcblas -lm 

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CFLAGS) -c -o $@ $<

# $(BUILD_DIR)/%.elf: $(BUILD_DIR)/%.o
# 	$(CXX) -o $@ $< $(LDFLAGS) 

$(TGT_FILE): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean pre_build all
.PRECIOUS: $(TGT_FILE)
clean:
	rm -rf $(OBJ_FILES) $(TGT_FILE) 

pre_build:
	mkdir -p build

all: pre_build $(TGT_FILE)
