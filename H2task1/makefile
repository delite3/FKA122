CC = gcc

LOCAL=${HOME}/.local
CFLAGS = \
	 -pedantic \
	 -Werror \
	 -Wall \
	 -g \
	 -Iinclude \
	 -I${LOCAL}/include

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LDFLAGS = \
         -Wl,-rpath,${LOCAL}/lib \
         -L${LOCAL}/lib

LIB = \
      -lm

OBJ = \
      obj/tools.o \
      obj/run.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

H2: obj metropolis_simulation

metropolis_simulation: obj/metropolis_simulation.o | obj
	$(CC) $(LDFLAGS) $^ -o metropolis_simulation $(LIB)

obj/%.o: src/%.c | obj
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find . -name "*.o" -exec rm {} \;
	find . -name "*.d" -exec rm {} \;
	find . -name "*Identifier" -exec rm {} \;
	rm -f H2 metropolis_simulation
	rm -rf obj obj_test

.PHONY: clean