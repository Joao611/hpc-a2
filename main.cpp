#include <iostream>
#include <fstream>
#include <sstream>
#include <limits> // std::numeric_limits
#include <memory>
#include <cstddef> // offsetof
#include <vector>
/* DAS-4 has GCC 6.3.0, which doesn't support std::filesystem. */
#include <experimental/filesystem>

#include "mpi.h"

using namespace std;
namespace fs = std::experimental::filesystem;

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
        if (root1 == root2) {
            return;
        }

        if (ranks[root1] < ranks[root2]) {
            int temp = root1;
            root1 = root2;
            root2 = temp;
        }

        parents[root2] = root1;
        if (ranks[root1] == ranks[root2]) {
            ranks[root1]++;
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
            exit(1);
        }
        // Bypass comments/header.
        while (ifs.peek() == '%') {
            string line;
            getline(ifs, line);
        }
        string line;
        getline(ifs, line);
        istringstream iss(line);
        // ignore rows since it's a square matrix (rows = cols = number of Vertices)
        iss.ignore(numeric_limits<streamsize>::max(), ' ');
        iss >> nVerts >> nEdges;
        edges = make_unique<Edge[]>(nEdges);
        for (int i = 0; getline(ifs, line); i++) {
            istringstream iss(line);
            Edge e;
            iss >> e.from >> e.to >> e.weight;
            // correct file's 1-indexing
            e.from--; e.to--;
            edges[i] = e;
        }
        ifs.close();
    }

    void printEdges() const {
        for (int e = 0; e < nEdges; e++) {
            cout << edges[e].from << '\t' << edges[e].to << '\t' << edges[e].weight << '\n';
        }
    }
};

// Minimum Spanning Forest.
struct MSF : Graph {
    MSF(int nVerts) : Graph(nVerts) {
        edges = make_unique<Edge[]>(nVerts - 1); // max edges for MSF (fully connected case)
    }
};

void handleUsage(int argc) {
    if (argc != 2) {
        if (proc == 0) {
            printf("Usage:\n\tmpirun [-np X] boruvka_mpi.out file\n");
        }
        MPI_Finalize();
        exit(1);
    }
}

void distributeEdges(const Graph &g, Edge *localEdges, int *nLocalEdges) {
    printf("P%d: Scattering %d edges in %d chunks\n", proc, g.nEdges, *nLocalEdges);
    MPI_Scatter(g.edges.get(), *nLocalEdges, mpiEdgeType,
            localEdges, *nLocalEdges, mpiEdgeType,
            0, MPI_COMM_WORLD);
    // If |E| isn't a multiple of numProcs, correct the number of edges for the last process.
    if (proc == numProcs - 1 && g.nEdges % *nLocalEdges != 0) {
        *nLocalEdges = g.nEdges % *nLocalEdges;
    }
}

/**
 * Gets all processes with the same closest edges,
 * prevailing the ones with the lowest weight when in conflict.
 * Processes send to each other the closest edges they obtained.
 * All data sent to P0 in log2(n) steps.
 */
void syncClosestEdges(const Graph &g, unique_ptr<Edge[]> &closestEdges) {
    printf("P%d: Syncing edges\n", proc);
    Edge *closestEdgesRecvd = new Edge[g.nVerts];
    // iteratively get all local edges to proc 0
    for (int iter = 1; iter < numProcs; iter *= 2) {
        printf("P%d: Step %d\n", proc, iter);
        if (proc % (2 * iter) == 0) {
            int sender = proc + iter;
            if (sender < numProcs) {
                printf("P%d: Receiving from %d\n", proc, sender);
                MPI_Recv(closestEdgesRecvd, g.nVerts, mpiEdgeType, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // combine results with own closest edges
                for (int v = 0; v < g.nVerts; v++) {
                    if (closestEdgesRecvd[v].weight < closestEdges[v].weight) {
                        closestEdges[v] = closestEdgesRecvd[v];
                    }
                }
            }
        } else if (proc % iter == 0) {
            int dest = proc - iter;
            printf("P%d: Sending to %d\n", proc, dest);
            MPI_Send(closestEdges.get(), g.nVerts, mpiEdgeType, dest, 0, MPI_COMM_WORLD);
        }
    }
    // send global closest edges to all processes
    printf("P%d: Broadcasting edges\n", proc);
    MPI_Bcast(closestEdges.get(), g.nVerts, mpiEdgeType, 0, MPI_COMM_WORLD);
    delete[] closestEdgesRecvd;
}

void printResults(const Graph &msf, const double readEndTime, const double endTime) {
    printf("-------------------------------------\n");
    printf("P%d: Result computed in %f seconds\n", proc, endTime - readEndTime);
    printf("P%d: Minimum Spanning Forest (MSF) has %d vertices and %d edges\n", proc, msf.nVerts, msf.nEdges);
    // msf.printEdges();
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
    createMpiEdgeType();
    handleUsage(argc);

    if (proc == 0) {
        fs::path p = string(argv[1]);
        double fileSize = fs::file_size(p) / 1000; // KB
        cout << "Name: " << processor_name << "\n";
        cout << "Number of processors: " << numProcs << "\n";
        cout << "File: " << argv[1] << " - ";
        if (fileSize > 1000) {
            cout << fileSize / 1000 << " MB\n";
        } else {
            cout << fileSize << " KB\n";
        }
        cout << "---------------------------------------\n";
    }

    double startTime = -1.0, readEndTime = -1.0, endTime = -1.0;
    if (proc == 0) {
        startTime = MPI_Wtime();
    }
    Graph g;

    if (proc == 0) {
        g.read(argv[1]);
        readEndTime = MPI_Wtime();
        printf("P%d: Graph loaded with %d vertices and %d edges in %f seconds\n",
                proc, g.nVerts, g.nEdges, readEndTime - startTime);
    }
    MPI_Bcast(&g.nVerts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g.nEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MSF msf(g.nVerts); // minimum spanning forest, only meant for proc 0
    Forest forest(g.nVerts);

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
                    if (localEdges[e].weight < closestEdges[roots[r]].weight) {
                        closestEdges[roots[r]] = localEdges[e];
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
                    msf.nEdges++;
                }
                forest.merge(e.from, e.to);
            }
        }
        if (proc == 0) {
            printf("P%d: MSF: %d edges\n", proc, msf.nEdges);
        }
    }

    if (proc == 0) {
        endTime = MPI_Wtime();
        printResults(msf, readEndTime, endTime);
    }
    
    delete[] localEdges;
    MPI_Type_free(&mpiEdgeType);
    
    if (proc == 0) {
        cout << "--------------------------------------------------\n";
        cout << "You may see an MPI_ABORT error below. Not to fret, as with "
            << "certain kinds of graphs the processes with rank != 0 may not notice "
            << "a solution has been found and are forcefully terminated this way, "
            << "as all required output has already been provided above.\n";
    }
    // the only main loop condition applying to non-P0's processes
    if (edgesFound) {
        // non-P0 processes may have not detected a solution was found and are
        // blocked in a MPI_Send() to P0, so all processes are forcefully terminated
        MPI_Abort(MPI_COMM_WORLD, 0);
    } else {
        MPI_Finalize();
    }

    return 0;
}