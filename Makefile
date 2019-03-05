CXXFLAGS:=-O3 -Wall -std=c++17 -flto -g -march=native -I. -fno-exceptions -fopenmp
LDFLAGS:=-flto -lpthread -fopenmp
CXX:=g++

ALL_SRCS:=$(wildcard **/*.cc) $(wildcard *.cc)
BIN_SRCS:=$(wildcard main/**/*.cc) $(wildcard main/*.cc)
TEST_SRCS:=$(wildcard **/*_test.cc) $(wildcard *_test.cc)
BINS:=$(BIN_SRCS:main/%.cc=bin/%)
TESTS:=$(TEST_SRCS:%.cc=build/%)
RUNTEST:=$(TESTS:%=.test_outputs/%)
SRCS:=$(filter-out ${BIN_SRCS}, ${ALL_SRCS})
SRCS:=$(filter-out ${TEST_SRCS}, ${SRCS})
ALL_OBJS:=$(ALL_SRCS:%.cc=build/%.o)
OBJS:=$(SRCS:%.cc=build/%.o)
DEPS:=$(ALL_SRCS:%.cc=.deps/%.d)
DIRS:=$(dir ${ALL_OBJS}) $(dir ${DEPS}) $(dir ${BINS}) $(dir ${TESTS}) \
	$(dir ${RUNTEST}) build

$(shell mkdir -p $(DIRS))

all: ${BINS}

.deps/%.d: %.cc Makefile
	${CXX} $< -M -MM -MP -MT $(patsubst .deps/%.d,build/%.o,$@) -o $@ ${CXXFLAGS}

build/%_test.o: %_test.cc .deps/%_test.d
	${CXX} $< -c -o $@ ${CXXFLAGS} $(shell pkg-config --cflags gmock gtest)

build/%.o: %.cc .deps/%.d
	${CXX} $< -c -o $@ ${CXXFLAGS}

build/%_test: build/%_test.o ${OBJS}
	${CXX} $^ -o $@ ${CXXFLAGS} ${LDFLAGS} $(shell pkg-config --libs gmock gtest_main)

bin/%: build/main/%.o ${OBJS}
	${CXX} $^ -o $@ ${CXXFLAGS} ${LDFLAGS}

test: ${RUNTEST}

.test_outputs/%: %
	./$^ &> $@ || ( cat $@ && exit 1 )

.PHONY: clean

clean:
	rm -rf bin/ build/ .deps/ .test_outputs/

.PRECIOUS: ${DEPS}

-include ${DEPS}

