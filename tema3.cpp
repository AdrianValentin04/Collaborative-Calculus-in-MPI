#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>
#include <fstream>

using namespace std;

struct boundaries {
    int start;
    int end;
};


int minimum(int x, int y) {
    
    if (x < y) {
        return x;
    } else {
        return y;
    }

}

int maximum(int x, int y) {

    if (x > y) {
        return x;
    } else {
        return y;
    }

}

int isCoord(int rank) {

    if (rank < 3) {
        return 1;
    } else {
        return 0;
    }

}

boundaries getLimits(int index, int k, int numberWorkers) {
    boundaries limits;

    limits.start = index * (double)k / numberWorkers;
    limits.end = minimum((index + 1) * k / numberWorkers, k);

    return limits; 
}

void readInput(int workers[][100], int *noWorkers, int rank) {
    string rankString = to_string(rank);
    string file = "cluster" + rankString + ".txt";
    ifstream fin(file);

    fin >> noWorkers[rank];

    for (int i = 1; i <= noWorkers[rank]; i ++) {
        fin >> workers[rank][i - 1];
    }
    fin.close();
}

void switchVector(int *x, int *y, int size) {

    for (int i = 0; i < size; i++) {
        x[i] = y[i];
    }

}

int main(int argc, char *argv[]) {

    int nProcesses, rank;
    string result;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // Declare the variables for topology - TASK 1
    int numberCoord = 3; //number of coordinators
    int maxWorkers = 100; //max number of workers for every coordinator
    int workers[numberCoord][100];

    int noWorkers[numberCoord]; //number of workers for every coordinator

    // For workers
    int rankPersonalCoord; //Rank of coordinator
    int defection = atoi(argv[2]);
        

    if (defection == 0){

        if (isCoord(rank) == 1) {
            // Read fron file
            
            readInput(workers, noWorkers, rank);

            // Send the topology to the other coordinators

            // First coordinator - send number of workers
            result = "M(" + to_string(rank) + ',' + to_string((rank + 1) % 3) + ")";
            cout << result << "\n";
            MPI_Send(&noWorkers[rank], 1, MPI_INT, ((rank + 1) % 3), 0, MPI_COMM_WORLD);
            

            // Second coordinator - send number of workers
            result = "M(" + to_string(rank) + ',' + to_string((rank + 2) % 3) + ")";
            cout << result << "\n";
            MPI_Send(&noWorkers[rank], 1, MPI_INT, (rank + 2) % 3, 0, MPI_COMM_WORLD);
            

            // Receive the number of workers from others
            for (int i = 1; i < 3; i ++) {
                MPI_Recv(&noWorkers[(rank + i) % 3], 1, MPI_INT, (rank + i) % 3, 0, MPI_COMM_WORLD, &status);
            }

            // First coordinator - send workers
            result = "M(" + to_string(rank) + ',' + to_string((rank + 1) % 3) + ")";
            MPI_Send(workers[rank], noWorkers[rank], MPI_INT, (rank + 1) % 3, 0, MPI_COMM_WORLD);
            cout << result << "\n";

            // Second coordinator - send workers
            result = "M(" + to_string(rank) + ',' + to_string((rank + 2) % 3) + ")";
            MPI_Send(workers[rank], noWorkers[rank], MPI_INT, (rank + 2) % 3, 0, MPI_COMM_WORLD);
            cout << result << "\n";

            // Receive the workers from the others
            for (int i = 1; i < 3; i ++) {
                MPI_Recv(workers[(rank + i) % 3], noWorkers[(rank + i) % 3], MPI_INT, (rank + i) % 3, 0, MPI_COMM_WORLD, &status);
            }

            // Send the topology to personal workers
            for (int i = 0; i < noWorkers[rank]; i++) {
                
                // Send info to every personal worker
                result = "M(" + to_string(rank) + ',' + to_string(workers[rank][i]) + ")";
                MPI_Send(&noWorkers, numberCoord, MPI_INT, workers[rank][i], 0, MPI_COMM_WORLD);
                cout << result << "\n";

                for (int j = 0; j < numberCoord; ++j) {
                    
                    // Every coordinator's workers
                    result = "M(" + to_string(rank) + ',' + to_string(workers[rank][i]) + ")";
                    MPI_Send(workers[j], noWorkers[j], MPI_INT, workers[rank][i], j, MPI_COMM_WORLD);
                    cout << result << "\n";

                }
            }

        } else {

            // Receive number of workers
            MPI_Recv(noWorkers, numberCoord, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            rankPersonalCoord = status.MPI_SOURCE;

            // Receive the topology
            for (int i = 0; i < numberCoord; i++) {
                MPI_Recv(workers[i], noWorkers[i], MPI_INT, rankPersonalCoord, i, MPI_COMM_WORLD, &status);
            }
        }

        
        // Every process needs to print the topology
        string topology = to_string(rank) + " ->";

        for (int i = 0; i < numberCoord; i++) {
            topology = topology + ' ' + to_string(i) + ':';

            for (int j = 0; j < noWorkers[i]; ++j) {
                topology = topology + to_string(workers[i][j]) + ",";
            }

            // Eliminate the last ","
            topology.pop_back();

            cout << topology << "\n";
        }

        // TASK 2

        // Declare variables for task 2

        int v[10000];
        int k = atoi(argv[1]);

        if(isCoord(rank) == 1) {
            
            if(rank == 0) {
                
                // Initialise the vector
                for (int i = 0; i < k; i++) {
                    v[i] = i;
                }

                // Send the vector to the other coordinators
                result = "M(" + to_string(rank) + ',' + to_string((rank + 1) % 3) + ")";
                MPI_Send(v, k, MPI_INT, (rank + 1) % 3, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string((rank + 2) % 3) + ")";
                MPI_Send(v, k, MPI_INT, (rank + 2) % 3, 0, MPI_COMM_WORLD);
                cout << result << "\n";

            } else {

                // The other coordinators receive the vector
                MPI_Recv(v, k, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

            }
        }

        boundaries limit;
        limit.start = 100000;
        limit.end = -1;
        
        if (isCoord(rank) == 1) {

            for (int i = 0; i < noWorkers[rank]; ++i) {
                int processRank = workers[rank][i];

                int index = i;
                for (int j = 1; j <= rank; j++) {
                    index = index + noWorkers[j - 1];
                }

                int noWorkersTotal = nProcesses - numberCoord;
                boundaries vectorLimit = getLimits(index, k, noWorkersTotal);
                

                // Send the vector and the indices
                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(v, k, MPI_INT, processRank, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(&vectorLimit.start, 1, MPI_INT, processRank, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(&vectorLimit.end, 1, MPI_INT, processRank, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                // Change the limits
                limit.start = minimum(vectorLimit.start, limit.start);
                limit.end = maximum(limit.end, vectorLimit.end);

                // Receive the calculus
                MPI_Recv(v + vectorLimit.start, vectorLimit.end - vectorLimit.start, MPI_INT, processRank, 0, MPI_COMM_WORLD, &status);

            }

            if (rank != 0){

                // Send to coordinator 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.end, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(v + limit.start, limit.end - limit.start, MPI_INT, 0, 0, MPI_COMM_WORLD);
                cout << result << "\n";

            } else {
                
                // Receive from coordinator 1
                MPI_Recv(&limit.start, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&limit.end, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(v + limit.start, limit.end - limit.start, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            
                // Receive from coordinator 2
                MPI_Recv(&limit.start, 1, MPI_INT, 2, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&limit.end, 1, MPI_INT, 2, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(v + limit.start, limit.end - limit.start, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

            }

        } else {

            // Receive the elements
            boundaries lim;
            MPI_Recv(v, k, MPI_INT, rankPersonalCoord, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&lim.start, 1, MPI_INT, rankPersonalCoord, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lim.end, 1, MPI_INT, rankPersonalCoord, 2, MPI_COMM_WORLD, &status);

            // Calculate the double of each element
            for (int i = lim.start; i < lim.end; i++) {
                v[i] = v[i] * 2;
            }

            // Send back
            result = "M(" + to_string(rank) + ',' + to_string(rankPersonalCoord) + ")";
            MPI_Send(v + lim.start, lim.end - lim.start, MPI_INT, rankPersonalCoord, 0, MPI_COMM_WORLD);
            cout << result << "\n";

        }    

        // Print the final vector
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            result = "Rezultat:";

            for (int i = 0; i < k; i++) {
                result += ' ' + to_string(v[i]);
            }
            
            cout << result << "\n";
        }
    

    } else {

        // The case when we don't have the 0 -> 1 comunication
        int to1Tag = 5;
        int to0Tag = 6;

        if (isCoord(rank) == 1) {

            readInput(workers, noWorkers, rank);

            if (rank == 0) {

                // Send both messages to coord. 2
                
                // This is for coord 2
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[rank], 1, MPI_INT, 2, 0, MPI_COMM_WORLD);

                // This is for coord 2
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(workers[rank], noWorkers[rank], MPI_INT, 2, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                // Receive from coord 2 - information for 2 
                MPI_Recv(&noWorkers[2], 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Receive from coord 2 - information for 1 
                MPI_Recv(&noWorkers[1], 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                // Receive from coord 2 - information for 2 
                MPI_Recv(workers[2], noWorkers[2], MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                // Receive from coord 2 - information for 1
                MPI_Recv(workers[1], noWorkers[1], MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


            } else if (rank == 1) {

                // Send both messages to coord. 2
                // This is for coord 2
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[rank], 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
            
                // This is for coord 2
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(workers[rank], noWorkers[rank], MPI_INT, 2, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                // Receive from coord 2 - information for 2 
                MPI_Recv(&noWorkers[2], 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Receive from coord 2 - information for 0 
                MPI_Recv(&noWorkers[0], 1, MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                // Receive from coord 2 - information for 2 
                MPI_Recv(workers[2], noWorkers[2], MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                
                // Receive from coord 2 - information for 0 
                MPI_Recv(workers[0], noWorkers[0], MPI_INT, 2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


            } else if(rank == 2) {

                // Send information to 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[rank], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

                // Send information to 1
                result = "M(" + to_string(rank) + ',' + to_string(1) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[rank], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

                // Receive info from 0
                MPI_Recv(&noWorkers[0], 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Send information for 0 to 1
                result = "M(" + to_string(rank) + ',' + to_string(1) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);


                // Receive info from 1
                MPI_Recv(&noWorkers[1], 1, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Send information for 1 to 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                cout << result << "\n";
                MPI_Send(&noWorkers[1], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);


                //Send information to 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(workers[rank], noWorkers[rank], MPI_INT, 0, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                // Send information to 1
                result = "M(" + to_string(rank) + ',' + to_string(1) + ")";
                MPI_Send(workers[rank], noWorkers[rank], MPI_INT, 1, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                // Receive from 0
                MPI_Recv(workers[0], noWorkers[0], MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Send information for 0 to 1
                result = "M(" + to_string(rank) + ',' + to_string(1) + ")";
                MPI_Send(workers[0], noWorkers[0], MPI_INT, 1, 0, MPI_COMM_WORLD);
                cout << result << "\n";


                // Receive from 1
                MPI_Recv(workers[1], noWorkers[1], MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                // Send information for 1 to 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(workers[1], noWorkers[1], MPI_INT, 0, 0, MPI_COMM_WORLD);
                cout << result << "\n";

            }

            // Send the topology to personal workers
            for (int i = 0; i < noWorkers[rank]; i++) {
                
                // Send info to every personal worker
                result = "M(" + to_string(rank) + ',' + to_string(workers[rank][i]) + ")";
                MPI_Send(&noWorkers, numberCoord, MPI_INT, workers[rank][i], 0, MPI_COMM_WORLD);
                cout << result << "\n";

                for (int j = 0; j < numberCoord; ++j) {
                    
                    // Every coordinator's workers
                    result = "M(" + to_string(rank) + ',' + to_string(workers[rank][i]) + ")";
                    MPI_Send(workers[j], noWorkers[j], MPI_INT, workers[rank][i], j, MPI_COMM_WORLD);
                    cout << result << "\n";

                }
            }

        } else {

             // Receive number of workers
            MPI_Recv(noWorkers, numberCoord, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            rankPersonalCoord = status.MPI_SOURCE;

            // Receive the topology
            for (int i = 0; i < numberCoord; i++) {
                MPI_Recv(workers[i], noWorkers[i], MPI_INT, rankPersonalCoord, i, MPI_COMM_WORLD, &status);
            }
        }

        // Every process needs to print the topology
        string topology = to_string(rank) + " ->";

        for (int i = 0; i < numberCoord; i++) {
            topology = topology + ' ' + to_string(i) + ':';

            for (int j = 0; j < noWorkers[i]; ++j) {
                topology = topology + to_string(workers[i][j]) + ",";
            }

            // Eliminate the last ","
            topology.pop_back();

            cout << topology << "\n";
        }



        // TASK 2
        // Declare variables for task 2

        int v[10000];
        int k = atoi(argv[1]);

        if(isCoord(rank) == 1) {
            
            if(rank == 0) {
                
                // Initialise the vector
                for (int i = 0; i < k; i++) {
                    v[i] = i;
                }

                // Send the vector to coordinator 2
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(v, k, MPI_INT, 2, 0, MPI_COMM_WORLD);
                cout << result << "\n";


            } else if (rank == 1) {

                // Receive from 2
                MPI_Recv(v, k, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

            } else {

                //Receive from coordinator 0
                MPI_Recv(v, k, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

                // Send to coordinator 1
                result = "M(" + to_string(rank) + ',' + to_string(1) + ")";
                MPI_Send(v, k, MPI_INT, 1, 0, MPI_COMM_WORLD);
                cout << result << "\n";

            }
        }

        boundaries limit;
        limit.start = 100000;
        limit.end = -1;
        
        if (isCoord(rank) == 1) {

            for (int i = 0; i < noWorkers[rank]; ++i) {
                int processRank = workers[rank][i];

                int index = i;
                for (int j = 1; j <= rank; j++) {
                    index = index + noWorkers[j - 1];
                }

                int noWorkersTotal = nProcesses - numberCoord;
                boundaries vectorLimit = getLimits(index, k, noWorkersTotal);
                

                // Send the vector and the indices
                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(v, k, MPI_INT, processRank, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(&vectorLimit.start, 1, MPI_INT, processRank, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(processRank) + ")";
                MPI_Send(&vectorLimit.end, 1, MPI_INT, processRank, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                // Change the limits
                limit.start = minimum(vectorLimit.start, limit.start);
                limit.end = maximum(limit.end, vectorLimit.end);

                // Receive the calculus
                MPI_Recv(v + vectorLimit.start, vectorLimit.end - vectorLimit.start, MPI_INT, processRank, 0, MPI_COMM_WORLD, &status);

            }

            if (rank == 0) {

                // Receive vector from cluster 2
                MPI_Recv(&limit.start, 1, MPI_INT, 2, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&limit.end, 1, MPI_INT, 2, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(v + limit.start, limit.end - limit.start, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

                // Receive vector from cluster 1
                MPI_Recv(&limit.start, 1, MPI_INT, 2, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&limit.end, 1, MPI_INT, 2, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(v + limit.start, limit.end - limit.start, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

            } else if (rank == 1) {

                // Send to coordinator 0
                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(&limit.start, 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(&limit.end, 1, MPI_INT, 2, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(2) + ")";
                MPI_Send(v + limit.start, limit.end - limit.start, MPI_INT, 2, 0, MPI_COMM_WORLD);
                cout << result << "\n";

            } else if (rank == 2) {

                // Send to coordinator 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.end, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(v + limit.start, limit.end - limit.start, MPI_INT, 0, 0, MPI_COMM_WORLD);
                cout << result << "\n";

                // Receive from 1
                MPI_Recv(&limit.start, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&limit.end, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(v + limit.start, limit.end - limit.start, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);

                // Send info from 1 to 0
                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(&limit.end, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                cout << result << "\n";

                result = "M(" + to_string(rank) + ',' + to_string(0) + ")";
                MPI_Send(v + limit.start, limit.end - limit.start, MPI_INT, 0, 0, MPI_COMM_WORLD);
                cout << result << "\n";


            }

        } else {

            // Receive the elements
            boundaries lim;
            MPI_Recv(v, k, MPI_INT, rankPersonalCoord, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&lim.start, 1, MPI_INT, rankPersonalCoord, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lim.end, 1, MPI_INT, rankPersonalCoord, 2, MPI_COMM_WORLD, &status);

            // Calculate the double of each element
            for (int i = lim.start; i < lim.end; i++) {
                v[i] = v[i] * 2;
            }

            // Send back
            result = "M(" + to_string(rank) + ',' + to_string(rankPersonalCoord) + ")";
            MPI_Send(v + lim.start, lim.end - lim.start, MPI_INT, rankPersonalCoord, 0, MPI_COMM_WORLD);
            cout << result << "\n";

        }    

        // Print the final vector
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            result = "Rezultat:";

            for (int i = 0; i < k; i++) {
                result += ' ' + to_string(v[i]);
            }
            
            cout << result << "\n";
        }

    }

    MPI_Finalize();
    return 0;
        
}