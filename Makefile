include config.mk

# define
# CXX = g++
CXXFLAGS = -Wall -std=c++14
INCLUDES = -I./include

# if USE_OPENMP=1，add -fopenmp
ifeq ($(USE_OPENMP),1)
    CXXFLAGS += -fopenmp
endif

MAIN_DIR = ./src/main
SETUP_DIR = ./src/setup
OBJ_DIR = ./obj
BIN_DIR = ./bin
TARGET = $(BIN_DIR)/xeno.sph

# obtain source files (.cpp files)
SRCS = $(MAIN_DIR)/main.cpp \
       $(MAIN_DIR)/config.cpp \
       $(MAIN_DIR)/density.cpp \
       $(MAIN_DIR)/evolve.cpp \
       $(MAIN_DIR)/force.cpp \
       $(MAIN_DIR)/kdtree.cpp \
       $(MAIN_DIR)/particles.cpp \
       $(MAIN_DIR)/random.cpp \
       $(MAIN_DIR)/utils.cpp \
       $(MAIN_DIR)/eos.cpp \
       $(SETUP_DIR)/$(SETUP)

# create object files (.o 文件)
OBJS = $(SRCS:$(MAIN_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJS := $(OBJS:$(SETUP_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# target
all: $(TARGET)

# compile rules
$(OBJ_DIR)/%.o: $(MAIN_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SETUP_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# chain rule
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# clean 
clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET)

# create obj directory (if not exist)
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# create bin directory (if not exist)
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: all clean