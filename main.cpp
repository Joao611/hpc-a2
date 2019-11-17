#include <iostream>
#include <fstream>
#include <sstream>
#include <limits> // std::numeric_limits
#include <memory>
#include <cstddef> // offsetof
#include <vector>

#include "mpi.h"

using namespace std;

#define SELF_ROOT -1

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
 * A vertice with value SELF_ROOT means it's its forest's root.
 * Rank is an optional property for each vertice to optimize performance.
 */
struct Forest {
    vector<int> parents;
    vector<int> ranks;

    Forest(int nVertices) {
        parents = vector<int>(nVertices, SELF_ROOT);
        ranks = vector<int>(nVertices, 0);
    }

    int find(int v) {
        if (parents[v] != SELF_ROOT) {
            parents[v] = find(parents[v]);
            return parents[v];
        } else {
            return v;
        }
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

    Graph() {}

    Graph(int nVerts) : nVerts(nVerts) {}

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

void distributeEdges(const Graph &g, Edge *localEdges, int *nLocalEdges) {
    printf("P%d: Scattering %d edges in %d chunks\n", proc, g.nEdges, *nLocalEdges);
    MPI_Scatter(&(g.edges[0]), *nLocalEdges, mpiEdgeType,
            localEdges, *nLocalEdges, mpiEdgeType,
            0, MPI_COMM_WORLD);
    // If |E| isn't a multiple of numProcs, correct the number of edges for the last process.
    if (proc == numProcs - 1 && g.nEdges % *nLocalEdges != 0) {
        *nLocalEdges = g.nEdges % *nLocalEdges;
    }
}

void syncClosestEdges(const Graph &g, unique_ptr<Edge[]> &closestEdges) {
    Edge *closestEdgesRecvd = new Edge[g.nVerts];
    // iteratively get all local edges to proc 0
    for (int iter = 1; iter < numProcs; iter *= 2) {
        for (int currProc = 0; currProc < numProcs; currProc++) {
            if (currProc % iter == 0) {
                int dest = currProc - iter;
                MPI_Send(&(closestEdges[0]), g.nVerts, mpiEdgeType, dest, 0, MPI_COMM_WORLD);
            } else if (currProc % (2 * iter) == 0) {
                int sender = currProc + iter;
                if (sender < numProcs) {
                    MPI_Recv(closestEdgesRecvd, g.nVerts, mpiEdgeType, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // combine results with own closest edges
                    for (int v = 0; v < g.nVerts; v++) {
                        if (closestEdgesRecvd[v].weight < closestEdges[v].weight) {
                            closestEdges[v] = closestEdgesRecvd[v];
                        }
                    }
                }
            }
        }
    }
    // send global closest edges to all processes
    MPI_Bcast(&(closestEdges[0]), g.nVerts, mpiEdgeType, 0, MPI_COMM_WORLD);
    delete[] closestEdgesRecvd;
}

void printResults(const Graph &msf, const double readEndTime, const double endTime) {
    printf("P%d: Result computed in %f seconds\n", proc, endTime - readEndTime);
    printf("P%d: Minimum Spanning Forest (MSF) has %d vertices and %d edges\n", proc, msf.nVerts, msf.nEdges);
    double weight = 0.0;
    for (int e = 0; e < msf.nEdges; e++) {
        weight += msf.edges[e].weight;
    }
    printf("P%d: MSF weight: %f\n", proc, weight);
}

int main(int argc, char *argv[]) {
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Get_processor_name(processor_name, &namelen);
    printf("Process %d on %s\n", proc, processor_name);

    createMpiEdgeType();
    
    handleUsage(argc);

    double startTime = -1.0, readEndTime = -1.0, endTime = -1.0;
    if (proc == 0) {
        startTime = MPI_Wtime();
    }
    Graph g;
    Forest forest(g.nVerts);

    if (proc == 0) {
        g.read(argv[1]);
        readEndTime = MPI_Wtime();
        printf("P%d: Graph loaded with %d vertices and %d edges in %f seconds\n",
                proc, g.nVerts, g.nEdges, readEndTime - startTime);
    }
    MPI_Bcast(&g.nVerts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g.nEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
    Graph msf(g.nVerts); // minimum spanning forest

    int nLocalEdges = (g.nEdges + numProcs - 1) / numProcs;
    Edge *localEdges = new Edge[nLocalEdges];
    distributeEdges(g, localEdges, &nLocalEdges);

    auto closestEdges = make_unique<Edge[]>(g.nVerts);
    bool edgesFound = true;

    while (msf.nEdges < g.nVerts - 1 && edgesFound) {
        edgesFound = false;
        for (int v = 0; v < g.nVerts; v++) {
            closestEdges[v].weight = numeric_limits<double>::max();
        }
        for (int e = 0; e < nLocalEdges; e++) {
            int roots[2] = { forest.find(localEdges[e].from),
                            forest.find(localEdges[e].to) };
            if (roots[0] != roots[1]) { // different tree
                for (int r = 0; r < 2; r++) {
                    if (localEdges[roots[r]].weight < closestEdges[roots[r]].weight) {
                        closestEdges[roots[r]] = localEdges[roots[r]];
                    }
                }
            }
        }

        syncClosestEdges(g, closestEdges);

        for (int v = 0; v < g.nVerts; v++) {
            const auto e = closestEdges[v];
            // if edge exists and it unites different sets
            if (e.weight < numeric_limits<double>::max()
                    && forest.find(e.from) != forest.find(e.to)) {
                edgesFound = true; // continue loop as we may find more
                if (proc == 0) {
                    msf.edges[msf.nEdges] = e;
                }
                forest.merge(e.from, e.to);
            }
        }
    }

    if (proc == 0) {
        endTime = MPI_Wtime();
        printResults(msf, readEndTime, endTime);
    }
    
    delete[] localEdges;
    MPI_Type_free(&mpiEdgeType);
    MPI_Finalize();

    return 0;
}