#include <iostream>
#include <fstream>
#include <sstream>
#include <limits> // std::numeric_limits
#include "mpi.h"

using namespace std;

static int proc = -1, numProcs = -1;

struct Edge {
    int from;
    int to;
    double weight;
};

struct Forest {
    
};

struct Graph {
    int nVerts = -1;
    int nEdges = -1;
    Edge *edges = nullptr;

    Graph(const string &filename) {
        ifstream ifs(filename);
        if (!ifs.is_open()) {
            cerr << "Could not open file " << filename << "\n";
            MPI_Finalize();
            exit(0);
        }
        // Bypass comments/header.
        while (ifs.peek() == '%') {
            string line;
            getline(ifs, line);
        }
        // ignore rows since it's a square matrix (rows = cols = number of Vertices)
        ifs.ignore(numeric_limits<streamsize>::max(), ' ');
        ifs >> nVerts >> nEdges;
        string line;
        while (getline(ifs, line)) {
            istringstream iss(line);
            Edge e;
            iss >> e.from >> e.to >> e.weight;
        }
        ifs.close();
    }
};

void handleUsage(int argc) {
    if (argc != 2) {
        if (proc == 0) {
            printf("Usage:\n\tmpirun [-np X] boruvka_mpi.out file\n");
        }
        MPI_Finalize();
        exit(0);
    }
}

void distributeEdges() {

}

int main(int argc, char *argv[]) {
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Get_processor_name(processor_name, &namelen);
    printf("Process %d on %s\n", proc, processor_name);
    
    handleUsage(argc);

    double startTime = -1.0, readEndTime = -1.0, endTime = -1.0;
    if (proc == 0) {
        startTime = MPI_Wtime();
    }
    Graph g(argv[1]);
    if (proc == 0) {
        readEndTime = MPI_Wtime();
        printf("Graph loaded with %d vertices and %d edges in %f seconds\n",
                g.nVerts, g.nEdges, readEndTime - startTime);
    }

    bool edgesFoundBefore = true;
    while (edgesFoundBefore) {
        edgesFoundBefore = false;
        distributeEdges();
    }

    if (proc == 0) {
        endTime = MPI_Wtime();
        printf("Result computed in %f seconds\n", endTime - readEndTime);
    }
    
    MPI_Finalize();

    return 0;
}