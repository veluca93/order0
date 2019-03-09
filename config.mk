CXXFLAGS:=-O3 -Wall -std=c++17 -flto -g -march=native -I. -fno-exceptions \
	-fopenmp -march=native
LDFLAGS:=-flto -lpthread
CXX:=g++

${TGT}bin/optimizer: ${TGT}build/optimizers/annealing.o.tar \
	${TGT}build/optimizers/line_search.o.tar \
	${TGT}build/optimizers/simulated_gradient.o.tar \
	${TGT}build/optimizers/cma_es.o.tar \
	${TGT}build/optimizers/nelder_mead.o.tar
