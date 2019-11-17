#include <iostream>
#include <fstream>
#include <sstream>
#include <limits> // std::numeric_limits
#include <memory>
#include <cstddef> // offsetof
#include <vector>

#include "mpi.h"

using namespace std;

static int proc = -1, numProcs = -1;
static MPI_Datatype mpiEdgeType;

struct Edge {
    int from;
    int to;
    double weight;
};

// TODO: Might need MPI_Type_create_resized() to be able to send multiple structs
void createMpiEdgeType() {
    const int nItems = 3;
    int blockLengths[nItems] = { 1, 1, 1 };
    MPI_Aint offsets[nItems] = { offsetof(Edge, from), offsetof(Edge, to), offsetof(Edge, weight) };
    MPI_Datatype types[nItems] = { MPI_INT, MPI_INT, MPI_DOUBLE };

    MPI_Type_create_struct(nItems, blockLengths, offsets, types, &mpiEdgeType);
    MPI_Type_commit(&mpiEdgeType);
}

/* Forest of vertices represented by a disjoint-set structure.
 * Each index is a vertice, and the value is its parent.
 * A vertice with itself as value means it's its forest's root.
 * Rank is an optional property for each vertice to optimize performance.
 */
struct Forest {
    unique_ptr<int[]> parents;
    vector<int> ranks;

    Forest(int nVertices) {
        parents = make_unique<int[]>(nVertices);
        ranks = vector<int>(nVertices, 0);
    }

    int find(int v) {
        if (parents[v] != v) {
            parents[v] = find(parents[v]);
        }
        return parents[v];
    }

    void merge(int v1, int v2) {
        int root1 = find(v1);
        int root2 = find(v2);

        // merge smaller rank root into the other root
        if (ranks[root1] < ranks[root2]) {
            parents[root1] = root2;
            if (ranks[root2] == ranks[root1]) {
                ranks[root2]++;
            }
        } else if (ranks[root2] < ranks[root1]) {
            parents[root2] = root1;
            if (ranks[root1] == ranks[root2]) {
                ranks[root1]++;
            }
        }
    }
};

struct Graph {
    int nVerts = 0;
    int nEdges = 0;
    unique_ptr<Edge[]> edges;

    void read(const string &filename) {
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
        edges = make_unique<Edge[]>(nEdges);
        string line;
        for (int i = 0; getline(ifs, line); i++) {
            istringstream iss(line);
            Edge e;
            iss >> e.from >> e.to >> e.weight;
            edges[i] = e;
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

void distributeEdges(const Graph &g, Edge *localEdges, int nLocalEdges) {
    MPI_Scatter(&(g.edges[0]), nLocalEdges, mpiEdgeType,
            localEdges, nLocalEdges, mpiEdgeType,
            0, MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Get_processor_name(processor_name, &namelen);
    printf("Process %d on %s\n", proc, processor_name);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    createMpiEdgeType();
    
    handleUsage(argc);

    double startTime = -1.0, readEndTime = -1.0, endTime = -1.0;
    if (proc == 0) {
        startTime = MPI_Wtime();
    }
    Graph g;
    Graph msf; // minimum spanning forest
    if (proc == 0) {
        g.read(argv[1]);
        readEndTime = MPI_Wtime();
        printf("Graph loaded with %d vertices and %d edges in %f seconds\n",
                g.nVerts, g.nEdges, readEndTime - startTime);
    }
    MPI_Bcast(&g.nVerts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g.nEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Forest forest(g.nVerts);

    int nLocalEdges = (g.nEdges + numProcs - 1) / numProcs;
    if (proc == numProcs - 1) {
        nLocalEdges = g.nEdges % nLocalEdges;
    }
    Edge *localEdges = (Edge *) malloc(nLocalEdges * sizeof(Edge));
    distributeEdges(g, localEdges, nLocalEdges);

    bool edgesFoundBefore = true;
    while (edgesFoundBefore) {
        edgesFoundBefore = false;
    }

    if (proc == 0) {
        endTime = MPI_Wtime();
        printf("Result computed in %f seconds\n", endTime - readEndTime);
    }
    
    MPI_Type_free(&mpiEdgeType);
    MPI_Finalize();

    return 0;
}