UNAME := $(shell uname)
current_dir = $(shell pwd)


ifeq ($(UNAME), Linux)
TARGET =  Rfunc.so
endif
ifeq ($(UNAME), Darwin)
TARGET = Rfunc_mac.so
endif

$(TARGET): Sim.h Sim.cpp Tree.h Tree.cpp Rfunc.cpp
	g++ -m64 -fpic  -O3 -pipe -I/$(current_dir) -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -shared $^ -o $@
