# ============================================================
# Project root (Makefile lives in scripts/)
# ============================================================

ROOT_DIR := .
SRC_DIR  := $(ROOT_DIR)/src
INC_DIR  := $(ROOT_DIR)/include
BUILD_DIR := $(ROOT_DIR)/build
OBJ_DIR  := $(BUILD_DIR)/obj
BIN_DIR  := $(BUILD_DIR)/bin

PROJECT := stars

# ============================================================
# OS detection
# ============================================================

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CXX := clang++
	CC  := clang
	EXE :=
endif

ifeq ($(UNAME_S),Linux)
	CXX := g++
	CC  := gcc
	EXE :=
endif

ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
	CXX := g++
	CC  := gcc
	EXE := .exe
endif

TARGET := $(BIN_DIR)/$(PROJECT)$(EXE)

# ============================================================
# Build configuration
# ============================================================

CONFIG ?= release

CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -Wno-writable-strings -Wno-unused-but-set-variable -Wno-comment -I$(INC_DIR)
CFLAGS_COMMON   := -Wall -Wextra -I$(INC_DIR)

ifeq ($(CONFIG),debug)
	CXXFLAGS := $(CXXFLAGS_COMMON) -O0 -g
	CFLAGS   := $(CFLAGS_COMMON) -O0 -g
else
	CXXFLAGS := $(CXXFLAGS_COMMON) -O3
	CFLAGS   := $(CFLAGS_COMMON) -O3
endif

LDFLAGS :=

# ============================================================
# Source files
# ============================================================

CPP_SOURCES := \
	aim_functions.cpp \
	aim_modules.cpp \
	aircraft_functions.cpp \
	aircraft_modules.cpp \
	class_functions.cpp \
	dsp_utils.cpp \
	execution.cpp \
	flat3_modules.cpp \
	global_functions.cpp \
	utility_functions.cpp

C_SOURCES := \
	kiss_fft.c

CPP_SRCS := $(addprefix $(SRC_DIR)/,$(CPP_SOURCES))
C_SRCS   := $(addprefix $(SRC_DIR)/,$(C_SOURCES))

OBJ_CPP := $(CPP_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJ_C   := $(C_SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

OBJS := $(OBJ_CPP) $(OBJ_C)

# ============================================================
# Build rules
# ============================================================

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean

