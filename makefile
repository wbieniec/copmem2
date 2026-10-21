CXX := g++

TARGET := copmem2
MEMFILL := memoryfill

SRC_COPMEM := \
    CopMEM2.src/CopMEM2.cpp \
    common/StringUtils.cpp \
    common/StopWatch.cpp \
    common/city.cpp \
    common/metrohash64.cpp \
    fmt/format.cc

SRC_MEMFILL := \
    MemoryFill.src/MemoryFill.cpp

INCLUDES := \
    -Icommon \
    -Ifmt

COMMON_FLAGS := \
    -std=c++17 \
    -pthread \
    $(INCLUDES)

RELEASE_FLAGS := \
    -march=native \
    -O3 \
    -funroll-loops \
    -DNDEBUG

DEBUG_FLAGS := \
    -O0 \
    -g3 \
    -ggdb3 \
    -fno-omit-frame-pointer

LDFLAGS := -pthread

###########################################################################
# Release configuration
###########################################################################

RELEASE_DIR := x64/Release
RELEASE_OBJ := $(RELEASE_DIR)/obj

RELEASE_COPMEM_OBJ := \
    $(RELEASE_OBJ)/CopMEM2.src/CopMEM2.o \
    $(RELEASE_OBJ)/common/StringUtils.o \
    $(RELEASE_OBJ)/common/StopWatch.o \
    $(RELEASE_OBJ)/common/city.o \
    $(RELEASE_OBJ)/common/metrohash64.o \
    $(RELEASE_OBJ)/fmt/format.o

RELEASE_MEMFILL_OBJ := \
    $(RELEASE_OBJ)/MemoryFill.src/MemoryFill.o

###########################################################################
# Debug configuration
###########################################################################

DEBUG_DIR := x64/Debug
DEBUG_OBJ := $(DEBUG_DIR)/obj

DEBUG_COPMEM_OBJ := \
    $(DEBUG_OBJ)/CopMEM2.src/CopMEM2.o \
    $(DEBUG_OBJ)/common/StringUtils.o \
    $(DEBUG_OBJ)/common/StopWatch.o \
    $(DEBUG_OBJ)/common/city.o \
    $(DEBUG_OBJ)/common/metrohash64.o \
    $(DEBUG_OBJ)/fmt/format.o

DEBUG_MEMFILL_OBJ := \
    $(DEBUG_OBJ)/MemoryFill.src/MemoryFill.o

###########################################################################
# Default target
###########################################################################

all: release

###########################################################################
# Release
###########################################################################

release: CXXFLAGS := $(COMMON_FLAGS) $(RELEASE_FLAGS)
release: $(RELEASE_DIR)/$(TARGET) $(RELEASE_DIR)/$(MEMFILL)

$(RELEASE_DIR)/$(TARGET): $(RELEASE_COPMEM_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(RELEASE_DIR)/$(MEMFILL): $(RELEASE_MEMFILL_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

###########################################################################
# Debug
###########################################################################

debug: CXXFLAGS := $(COMMON_FLAGS) $(DEBUG_FLAGS)
debug: $(DEBUG_DIR)/$(TARGET) $(DEBUG_DIR)/$(MEMFILL)

$(DEBUG_DIR)/$(TARGET): $(DEBUG_COPMEM_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(DEBUG_DIR)/$(MEMFILL): $(DEBUG_MEMFILL_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $^ -o $@ $(LDFLAGS)

###########################################################################
# Release object rules
###########################################################################

$(RELEASE_OBJ)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(RELEASE_OBJ)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

###########################################################################
# Debug object rules
###########################################################################

$(DEBUG_OBJ)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DEBUG_OBJ)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

###########################################################################
# Clean
###########################################################################

clean:
	rm -rf x64

.PHONY: all release debug clean