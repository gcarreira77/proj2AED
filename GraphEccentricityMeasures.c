//
// Algoritmos e Estruturas de Dados --- 2024/2025
//
// Joaquim Madeira - Dec 2024
//
// GraphEccentricityMeasures
//

// Student Name : Guilherme Carreira
// Student Number : 120159
// Student Name : Daniel Couto
// Student Number : 50138

/*** COMPLETE THE GraphEccentricityMeasuresCompute FUNCTION ***/
/*** COMPLETE THE GraphGetCentralVertices FUNCTION ***/
/*** COMPLETE THE GraphEccentricityMeasuresPrint FUNCTION ***/

#include "GraphEccentricityMeasures.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "GraphAllPairsShortestDistances.h"

struct _GraphEccentricityMeasures {
  unsigned int*
      centralVertices;  // centralVertices[0] = number of central vertices
                        // array size is (number of central vertices + 1)
  int* eccentricity;    // the eccentricity value of each vertex
  Graph* graph;         // the graph
  int graphRadius;      // the graph radius
  int graphDiameter;    // the graph diameter
};


static void computeEccentricities(Graph* g, GraphAllPairsShortestDistances* apsd, 
                                  int* eccentricity, int* radius, int* diameter);

static unsigned int* computeCentralVertices(const int* eccentricity, unsigned int numVertices, int radius);                                 

// Allocate memory
// Compute the vertex eccentricity values
// Compute graph radius and graph diameter
// Compute the set of central vertices
GraphEccentricityMeasures* GraphEccentricityMeasuresCompute(Graph* g) {
  assert(g != NULL);

  // COMPLETE THE CODE
  // CREATE AUXILIARY (static) FUNCTIONS, IF USEFUL
  // Graph radius --- the smallest vertex eccentricity value
  // Graph diameter --- the largest vertex eccentricity value
  // Do not forget that -1 represents an IDEFINITE value

  // Computing the set of central vertices
  // Allocate the central vertices array : number of central vertices + 1
  // Fill in the central vertices array
  
  unsigned int numVertices = GraphGetNumVertices(g);

  // Allocate the structure
  GraphEccentricityMeasures* measures = malloc(sizeof(GraphEccentricityMeasures));
  measures->eccentricity = calloc(numVertices, sizeof(int));
  measures->graph = g;
  measures->graphRadius = -1;
  measures->graphDiameter = -1;
  
  // Compute all-pairs shortest distances
  GraphAllPairsShortestDistances* allPairsShortestDistances = GraphAllPairsShortestDistancesExecute(g);
  
  // Compute eccentricities, radius, and diameter
  computeEccentricities(g, allPairsShortestDistances, measures->eccentricity, &measures->graphRadius, &measures->graphDiameter);
  
  // Compute central vertices
  measures->centralVertices = computeCentralVertices(measures->eccentricity, numVertices, measures->graphRadius);
  
  // Destroy the APSD structure
  GraphAllPairsShortestDistancesDestroy(&allPairsShortestDistances);

  return measures;
}

void GraphEccentricityMeasuresDestroy(GraphEccentricityMeasures** p) {
  assert(*p != NULL);

  GraphEccentricityMeasures* aux = *p;

  free(aux->centralVertices);
  free(aux->eccentricity);

  free(*p);
  *p = NULL;
}

// Getting the computed measures

int GraphGetRadius(const GraphEccentricityMeasures* p) {
  assert(p != NULL);

  return p->graphRadius;
}

int GraphGetDiameter(const GraphEccentricityMeasures* p) {
  assert(p != NULL);

  return p->graphDiameter;
}

int GraphGetVertexEccentricity(const GraphEccentricityMeasures* p,
                               unsigned int v) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));
  assert(p->eccentricity != NULL);

  return p->eccentricity[v];
}

// Getting a copy of the set of central vertices
// centralVertices[0] = number of central vertices in the set
unsigned int* GraphGetCentralVertices(const GraphEccentricityMeasures* p) {
  assert(p != NULL);
  assert(p->centralVertices != NULL);

  // COMPLETE THE CODE
  unsigned int size = p->centralVertices[0];
  unsigned int* copy = malloc((size + 1) * sizeof(unsigned int));
  for (unsigned int i = 0; i <= size; i++) {
      copy[i] = p->centralVertices[i];
  }
  return copy;
}

// Print the graph radius and diameter
// Print the vertex eccentricity values
// Print the set of central vertices
void GraphEccentricityMeasuresPrint(const GraphEccentricityMeasures* p) {
  // COMPLETE THE CODE
  assert(p != NULL);

  printf("Graph Radius: %d\n", p->graphRadius);
  printf("Graph Diameter: %d\n", p->graphDiameter);
  printf("Vertex Eccentricities:\n");
  for (unsigned int v = 0; v < GraphGetNumVertices(p->graph); v++) {
      printf("Vertex %u: %d\n", v, p->eccentricity[v]);
  }
  printf("Central Vertices (Eccentricity = %d):\n", p->graphRadius);
  for (unsigned int i = 1; i <= p->centralVertices[0]; i++) {
      printf("Vertex %u\n", p->centralVertices[i]);
  }
}

// Compute eccentricity, radius, and diameter
static void computeEccentricities(Graph* g, GraphAllPairsShortestDistances* allPairsShortestDistances, 
                                  int* eccentricity, int* radius, int* diameter) {
    unsigned int numVertices = GraphGetNumVertices(g);

    *radius = -1;
    *diameter = -1;

    for (unsigned int v = 0; v < numVertices; v++) {
        int maxDistance = -1;

        for (unsigned int w = 0; w < numVertices; w++) {
            int distance = GraphGetDistanceVW(allPairsShortestDistances, v, w);
            if (distance > maxDistance) {
                maxDistance = distance;
            }
        }

        eccentricity[v] = maxDistance;

        if (*radius == -1 || maxDistance < *radius) {
            *radius = maxDistance;
        }

        if (maxDistance > *diameter) {
            *diameter = maxDistance;
        }
    }
}

// Compute the set of central vertices
static unsigned int* computeCentralVertices(const int* eccentricity, unsigned int numVertices, int radius) {
    unsigned int centralCount = 0;

    // Count the number of central vertices
    for (unsigned int v = 0; v < numVertices; v++) {
        if (eccentricity[v] == radius) {
            centralCount++;
        }
    }

    // Allocate and populate the central vertices array
    unsigned int* centralVertices = malloc((centralCount + 1) * sizeof(unsigned int));
    centralVertices[0] = centralCount;

    unsigned int idx = 1;
    for (unsigned int v = 0; v < numVertices; v++) {
        if (eccentricity[v] == radius) {
            centralVertices[idx++] = v;
        }
    }

    return centralVertices;
}
