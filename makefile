all:
	mpic++ -O2 -std=c++17 main.cpp -o boruvka_mpi.out -lstdc++fs

debug:
	mpic++ -g -std=c++17 main.cpp -o boruvka_mpi.out -lstdc++fs