UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
TARGET =  Rfunc.so
endif
ifeq ($(UNAME), Darwin)
TARGET = Rfunc_mac.so
endif

$(TARGET): Sim.h Sim.cpp Tree.h Tree.cpp Rfunc.cpp
	g++ -m64 -I/Users/opennopasswordrequired/Downloads/boost_1_55_0 -I/home/bop13mc/include/boost_1_55_0 -I/usr/include/R -DNDEBUG  -I/usr/local/include -fpic  -O3 -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -shared $^ -o $@
