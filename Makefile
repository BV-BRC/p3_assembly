TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime

SRC_PYTHON = $(wildcard scripts/*.py)

all: bin 

bin: $(BIN_PYTHON)

deploy: deploy-client 
deploy-all: deploy-client 
deploy-client: deploy-scripts 

include $(TOP_DIR)/tools/Makefile.common.rules
