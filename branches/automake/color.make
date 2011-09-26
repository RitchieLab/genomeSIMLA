
ifeq ($(TERM), xterm-color)
TPUT=$(shell which tput)
TERM_CYAN=--@$(TPUT) setaf 4
TERM_PINK=--@$(TPUT) setaf 5
TERM_GREEN=--@$(TPUT) setaf 2
TERM_RESET=--@$(TPUT) sgr0
endif