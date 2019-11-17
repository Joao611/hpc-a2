all:
	mpic++ -O2 main.cpp -o boruvka_mpi.out

debug:
	mpic++ -g main.cpp -o boruvka_mpi.out