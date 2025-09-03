include config.mk

# define
CXX = g++
CXXFLAGS = -Wall -std=c++14
INCLUDES = -I./include

SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
TARGET = $(BIN_DIR)/xeno.sph

# obtain source files (.cpp files)
SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/config.cpp $(SRC_DIR)/density.cpp $(SRC_DIR)/evolve.cpp $(SRC_DIR)/force.cpp $(SRC_DIR)/kdtree.cpp $(SRC_DIR)/particles.cpp $(SRC_DIR)/random.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/$(SETUP) $(SRC_DIR)/eos.cpp

# create object files (.o 文件)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# target
all: $(TARGET)

# compile rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# cbain rule
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