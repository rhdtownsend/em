# File     : Makefile
# Purpose  : top-level makefile

# Build frontend executables
FRONTENDS ?= yes

# Enable debugging (with a performance penalty)
DEBUG ?= no

# Build & link against shared libraries (with a possible performance
# penalty)
SHARED ?= no

# Enable OpenMP parallelization
OMP ?= yes

# Enable FPE checks
FPE ?= yes

# Enable correctly rounded math functions
CRMATH ?= yes

# Use IEEE fortran features
IEEE ?= yes

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# Export options

export FRONTENDS
export DEBUG
export SHARED
export OMP
export FPE
export CRMATH
export IEEE

# General make settings

export SH = /bin/bash
export MAKEFLAGS += --no-print-directory

# Paths

export BIN_DIR ?= $(CURDIR)/bin
export LIB_DIR ?= $(CURDIR)/lib
export PKG_DIR ?= $(LIB_DIR)/pkgconfig
export INC_DIR ?= $(CURDIR)/include
export SCT_DIR ?= $(CURDIR)/scripts

export SRC_DIR := $(CURDIR)/src
export SRC_DIRS := $(SRC_DIR)/lib $(SRC_DIR)/mesa $(SRC_DIR)/gyre $(SRC_DIR)/freq \
                   $(SRC_DIR)/math $(SRC_DIR)/include $(SRC_DIR)/frontend $(SRC_DIR)/amp

# Rules

install : build | $(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR)
	@$(MAKE) -C build $@

build : $(PKG_DIR)/mesa.pc $(PKG_DIR)/gyre.pc
	@$(MAKE) -C build $@

$(PKG_DIR)/mesa.pc : $(SRC_DIR)/pkgconfig/mesa.pc.in | $(PKG_DIR)
	@$(SCT_DIR)/configure_mesa.sh $< $@

$(PKG_DIR)/gyre.pc : $(SRC_DIR)/pkgconfig/gyre.pc.in | $(PKG_DIR)
	@$(SCT_DIR)/configure_gyre.sh $< $@

clean :
	@$(MAKE) -C build $@
	@rm -rf $(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR)

check_versions :
	@

check_src :
	@$(MAKE) -C build $@

.PHONY: install build clean check_src

$(BIN_DIR) $(LIB_DIR) $(PKG_DIR) $(INC_DIR) :
	@mkdir -p $@
