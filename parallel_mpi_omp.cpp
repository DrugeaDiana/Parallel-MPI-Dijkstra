#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <utility>
#include <set>
#include <climits>
#include <pthread.h>
#include <chrono>
#include <cmath>
#include "constants.h"
#include <mpi.h>
#include <omp.h>

using namespace std;

#define ROOT 0

// Global variables are initialized with 0
int source;
unsigned short int dist[NUM_NODES + 1][NUM_NODES + 1];
ofstream fout ("out_mpi_omp");
int dij[NUM_NODES + 1];
bool sptSet[NUM_NODES + 1];

#ifdef DEBUG
void print_matrix() {
    for (int i = 1; i <= NUM_NODES; i++) {
        for (int j = 1; j <= NUM_NODES; j++) {
            cout << dist[i][j] << " ";
        }
        cout << endl;
    }
}
#endif

void generate_matrix()
{
    // Set the seed of the random number generator
    srand(SEED);
    // fill the array with infinity
    fill(&dist[0][0], &dist[0][0] + sizeof(dist) / sizeof(dist[0][0]), INF_GEN);
    // Set the diagonal to 0
    for (int i = 1; i <= NUM_NODES; i++) {
        dist[i][i] = 0;
    }
    int x, y, i, cost;
    const int edge_per_line = NUM_EDGES / NUM_NODES + 1;
    for (x = 1; x <= NUM_NODES; x++) {
        for (i = 0; i < edge_per_line; i++) {
            do {
                y = rand() % NUM_NODES + 1;
            } while (x == y || dist[x][y] != INF_GEN);
            cost = rand() % MAX_COST + 1;
            dist[x][y] = cost;
        }
    }
    // Choose a random source
    source = rand() % NUM_NODES + 1;
}

void printSolution()
{
    for (int i = 1; i <= NUM_NODES; i++)
        if (dij[i] == INF) {
            fout << "NaN" << " ";
        } else {
            fout << dij[i] << " ";
        }
}

void dijkstra(int proc_no, int rank)
{
    if (rank == ROOT) {
        fill(&dij[0], &dij[0] + sizeof(dij) / sizeof(dij[0]), INF);
        fill(&sptSet[0], &sptSet[0] + sizeof(sptSet) / sizeof(sptSet[0]), false);
        dij[source] = 0;
    }
    
    const int local_size = NUM_NODES / proc_no;
    
    int dij_local[local_size] = {0};
    bool sptSet_local[local_size] = {false};
    for (int count = 1; count <= NUM_NODES; count++) {
        // Scatter dij to the workers
        if (rank == ROOT) {
            MPI_Scatter(&dij[1], local_size, MPI_INT, dij_local, local_size, MPI_INT, ROOT, MPI_COMM_WORLD);
        } else {
            MPI_Scatter(NULL, 0, MPI_INT, dij_local, local_size, MPI_INT, ROOT, MPI_COMM_WORLD);
        }

        // calculate min distance from the portion of array the proc has
        int index = 0;
        int min_dist = INF;
        for (int v = 0; v < local_size; v++)
            if (sptSet_local[v] == false && dij_local[v] < min_dist)
                index = local_size * rank + v + 1, min_dist = dij_local[v];

        // gather all index and min dist
        int index_array[proc_no];
        int min_dist_array[proc_no];
        if (rank == ROOT) {
            MPI_Gather(&index, 1, MPI_INT, index_array, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&min_dist, 1, MPI_INT, min_dist_array, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        } else {
            MPI_Gather(&index, 1, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&min_dist, 1, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
        }

        // root determine the min distance and index
        if (rank == ROOT) {
            for (int i = 1; i < proc_no; i++) {
                if (min_dist_array[i] < min_dist) {
                    min_dist = min_dist_array[i];
                    index = index_array[i];
                }
            }
            if (min_dist == INF) {
                break;
            }
            sptSet[index] = true;
        }

        // broadcast the index and min_dist
        MPI_Bcast(&index, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        MPI_Bcast(&min_dist, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

        // Worker who has index has to update sptSet as well
        if (index >= local_size * rank + 1 && index <= local_size * (rank + 1)) {
            sptSet_local[index - local_size * rank - 1] = true;
        }

        // Scatter the current line of dist matrix aka dist[idx]
        unsigned short int dist_local[local_size];
        if (rank == ROOT) {
            MPI_Scatter(&dist[index][1], local_size, MPI_UNSIGNED_SHORT, dist_local, local_size, MPI_UNSIGNED_SHORT, ROOT, MPI_COMM_WORLD);
        } else {
            MPI_Scatter(NULL, 0, MPI_UNSIGNED_SHORT, dist_local, local_size, MPI_UNSIGNED_SHORT, ROOT, MPI_COMM_WORLD);
        }

        // Update shortest paths
        #pragma omp parallel for schedule(static, 2)
        for (int v = 0; v <= local_size; v++) {
            if (!sptSet_local[v] && dist_local[v] != INF_GEN)
                if (min_dist + dist_local[v] < dij_local[v])
                    dij_local[v] = min_dist + dist_local[v];
        }

        // Gather the results
        if (rank == ROOT) {
            MPI_Gather(dij_local, local_size, MPI_INT, &dij[1], local_size, MPI_INT, ROOT, MPI_COMM_WORLD);
        } else {
            MPI_Gather(dij_local, local_size, MPI_INT, NULL, 0, MPI_INT, ROOT, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, proc_no;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_no);
    omp_set_num_threads(NUM_THREADS);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    // Generate matrix onl by master to make sure it's consistent
    if (rank == ROOT) {
        start = std::chrono::system_clock::now();
        generate_matrix();
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        cout << "Generate time: " << elapsed_seconds.count() << "s\n";
    }

#ifdef DEBUG
    if (rank == ROOT)
        print_matrix();
#endif

    if (rank == ROOT) {
        start = std::chrono::system_clock::now();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    dijkstra(proc_no, rank);
    MPI_Finalize();

    if (rank == ROOT) {
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        cout << "Dijkstra time: " << elapsed_seconds.count() << "s\n";
        printSolution();
	    fout.close();
    }

    return 0;
}