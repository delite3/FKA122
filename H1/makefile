CC = gcc
LOCAL=${HOME}/.local
LIBDIR=lib

CFLAGS = \
	 -g \
	 -pedantic \
	 -Werror \
	 -Wall \
	 -Iinclude \
	 -Wl,-rpath,${LOCAL} \
	 -I${LOCAL}/include \
	 -Wl,-rpath,${LIBDIR} \
	 -L${LIBDIR} \
	 -L${LOCAL}/lib \
	 -Isrc/nepc/

CFLAGS_OPT = \
	     -O0 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas \
      -lnep

OBJ = \
      obj/tools.o \
      obj/lattice.o \
      obj/potential.o \
      obj/run.o

OBJ_LIB = \
	obj_nep/nep.o \
	obj_nep/interface.o

CXXFLAGS = \
	   -DUSE_TABLE_FOR_RADIAL_FUNCTIONS \
	   -O3 \
	   -fPIC

MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

H1: $(MAIN) $(OBJ) lib
	$(CC) $(CFLAGS) $(MAIN) $(OBJ) $(LIB) -o $@

obj/%.o: src/%.c | obj
	$(CC) -MMD -c $(CFLAGS) $< -o $@

obj_nep/%.o: src/nepc/%.cpp | obj
	$(CXX) -MMD -c $(CXXFLAGS) $< -o $@

lib: $(OBJ_LIB)
	$(CXX) -shared $^ -o ${LIBDIR}/libnep.so

obj:
	mkdir -p obj
	mkdir -p obj_nep
	mkdir -p ${LIBDIR}

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	find -iname "*.Identifier" -exec rm {} \;
	
	rm -f H1 run-test
	rm -rf obj obj_test obj_nep ${LIBDIR}
	rm -f energy_time.txt

.PHONY: clean
