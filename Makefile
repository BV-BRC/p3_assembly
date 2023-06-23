TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime

WRAP_PYTHON_TOOL = wrap_python3

APP_SERVICE = app_service

all: bin 

bin: $(BIN_PYTHON) $(BIN_PERL) $(BIN_SERVICE_PERL)

deploy: deploy-client 
deploy-all: deploy-client 
deploy-client: deploy-scripts 
deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

include $(TOP_DIR)/tools/Makefile.common.rules
