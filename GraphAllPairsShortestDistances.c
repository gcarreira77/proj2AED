//
// Algoritmos e Estruturas de Dados --- 2024/2025
//
// Joaquim Madeira - Dec 2024
//
// GraphAllPairsShortestDistances
//

// Student Name :
// Student Number :
// Student Name :
// Student Number :

/*** COMPLETE THE GraphAllPairsShortestDistancesExecute FUNCTION ***/

#include "GraphAllPairsShortestDistances.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "GraphBellmanFordAlg.h"

struct _GraphAllPairsShortestDistances {
  int** distance;  // The 2D matrix storing the all-pairs shortest distances
                   // It is stored as an array of pointers to 1D rows
                   // Idea: an INDEFINITE distance value is stored as -1
  Graph* graph;
};

// Allocate memory and initialize the distance matrix
// Compute the distances between vertices by running the Bellman-Ford algorithm
GraphAllPairsShortestDistances* GraphAllPairsShortestDistancesExecute(Graph* g) {
  assert(g != NULL);

  // Get the number of vertices in the graph
  unsigned int numVertices = GraphGetNumVertices(g);

  // Allocate memory for the result structure
  GraphAllPairsShortestDistances* result =
      (GraphAllPairsShortestDistances*)malloc(sizeof(struct _GraphAllPairsShortestDistances));
  assert(result != NULL);

  result->graph = g;

  // Allocate memory for the distance matrix
  result->distance = (int**)malloc(numVertices * sizeof(int*));
  assert(result->distance != NULL);

  // Allocate memory for each row in the distance matrix
  for (unsigned int i = 0; i < numVertices; i++) {
    result->distance[i] = (int*)malloc(numVertices * sizeof(int));
    assert(result->distance[i] != NULL);
  }

  // Initialize the distance matrix
  for (unsigned int i = 0; i < numVertices; i++) {
    for (unsigned int j = 0; j < numVertices; j++) {
      if (i == j) {
        result->distance[i][j] = 0;  // Distance to itself is 0
      } else {
        result->distance[i][j] = -1;  // No path initially
      }
    }
  }

  // Compute the shortest paths for each vertex using Bellman-Ford
  for (unsigned int startVertex = 0; startVertex < numVertices; startVertex++) {
    GraphBellmanFordAlg* bfResult = GraphBellmanFordAlgExecute(g, startVertex);
    
    // Copy the results from Bellman-Ford algorithm into the distance matrix
    for (unsigned int v = 0; v < numVertices; v++) {
      int dist = GraphBellmanFordAlgDistance(bfResult, v);
      result->distance[startVertex][v] = dist;
    }

    // Clean up the Bellman-Ford result
    GraphBellmanFordAlgDestroy(&bfResult);
  }

  return result;
}


void GraphAllPairsShortestDistancesDestroy(GraphAllPairsShortestDistances** p) {
  assert(*p != NULL);

  GraphAllPairsShortestDistances* aux = *p;
  unsigned int numVertices = GraphGetNumVertices(aux->graph);

  for (unsigned int i = 0; i < numVertices; i++) {
    free(aux->distance[i]);
  }

  free(aux->distance);

  free(*p);
  *p = NULL;
}

// Getting the result

int GraphGetDistanceVW(const GraphAllPairsShortestDistances* p, unsigned int v,
                       unsigned int w) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));
  assert(w < GraphGetNumVertices(p->graph));

  return p->distance[v][w];
}

// DISPLAYING on the console

void GraphAllPairsShortestDistancesPrint(
    const GraphAllPairsShortestDistances* p) {
  assert(p != NULL);

  unsigned int numVertices = GraphGetNumVertices(p->graph);
  printf("Graph distance matrix - %u vertices\n", numVertices);

  for (unsigned int i = 0; i < numVertices; i++) {
    for (unsigned int j = 0; j < numVertices; j++) {
      int distanceIJ = p->distance[i][j];
      if (distanceIJ == -1) {
        // INFINITY - j was not reached from i
        printf(" INF");
      } else {
        printf(" %3d", distanceIJ);
      }
    }
    printf("\n");
  }
}
