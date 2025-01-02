//
// Algoritmos e Estruturas de Dados --- 2024/2025
//
// Joaquim Madeira - Dec 2024
//
// GraphBellmanFord - Bellman-Ford Algorithm
//

// Student Name : Guilherme Carreira
// Student Number : 120159
// Student Name : Daniel Couto
// Student Number : 50138

/*** COMPLETE THE GraphBellmanFordAlgExecute FUNCTION ***/

#include "GraphBellmanFordAlg.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "Graph.h"
#include "IntegersStack.h"
#include "instrumentation.h"

struct _GraphBellmanFordAlg {
  unsigned int* marked;  // To mark vertices when reached for the first time
  int* distance;  // The number of edges on the path from the start vertex
                  // distance[i]=-1, if no path found from the start vertex to i
  int* predecessor;  // The predecessor vertex in the shortest path
                     // predecessor[i]=-1, if no predecessor exists
  Graph* graph;
  unsigned int startVertex;  // The root of the shortest-paths tree
};

/* Inicial unchanged function
GraphBellmanFordAlg* GraphBellmanFordAlgExecute(Graph* g,
                                                unsigned int startVertex) {
  assert(g != NULL);
  assert(startVertex < GraphGetNumVertices(g));
  assert(GraphIsWeighted(g) == 0);

  GraphBellmanFordAlg* result =
      (GraphBellmanFordAlg*)malloc(sizeof(struct _GraphBellmanFordAlg));
  assert(result != NULL);

  // Given graph and start vertex for the shortest-paths
  result->graph = g;
  result->startVertex = startVertex;

  unsigned int numVertices = GraphGetNumVertices(g);

  //
  // TO BE COMPLETED !!
  //
  // CREATE AND INITIALIZE
  // result->marked
  // result->distance
  // result->predecessor
  //

  // Mark all vertices as not yet visited, i.e., ZERO
  
  // No vertex has (yet) a (valid) predecessor
  
  // No vertex has (yet) a (valid) distance to the start vertex
  
  // THE ALGORTIHM TO BUILD THE SHORTEST-PATHS TREE

  return NULL;
}

*/

GraphBellmanFordAlg* GraphBellmanFordAlgExecute(Graph* g,
                                                unsigned int startVertex) {
  assert(g != NULL);
  assert(startVertex < GraphGetNumVertices(g));
  assert(GraphIsWeighted(g) == 0); // Ensure the graph is unweighted

  GraphBellmanFordAlg* result =
      (GraphBellmanFordAlg*)malloc(sizeof(struct _GraphBellmanFordAlg));
  assert(result != NULL);

  // Initialize attributes
  unsigned int numVertices = GraphGetNumVertices(g);
  result->graph = g;
  result->startVertex = startVertex;

  // Allocate memory for marked, distance, and predecessor arrays
  result->marked = (unsigned int*)calloc(numVertices, sizeof(unsigned int));
  result->distance = (int*)malloc(numVertices * sizeof(int));
  result->predecessor = (int*)malloc(numVertices * sizeof(int));

  assert(result->marked != NULL);
  assert(result->distance != NULL);
  assert(result->predecessor != NULL);

  // Initialize all vertices as unvisited and distances as "infinity" (-1)
  for (unsigned int i = 0; i < numVertices; i++) {
    result->marked[i] = 0;
    result->distance[i] = -1;  // -1 represents infinity
    result->predecessor[i] = -1;
  }

  // Initialize the start vertex
  result->marked[startVertex] = 1;
  result->distance[startVertex] = 0;

  int counterEx = 0; //counter exterior loop
  int counterMid = 0; //counter middle loop
  int counterInt= 0; //counter interior loop

  // Perform the Bellman-Ford algorithm
  for (unsigned int i = 0; i < numVertices - 1; i++) {
    counterEx += 1;
    for (unsigned int j = 0; j < numVertices; j++) {
      counterMid += 1;

      unsigned int* adjacents = GraphGetAdjacentsTo(g, j);   //adjacent vertexes
      
      unsigned int num_adj = adjacents[0];        //Number of adjacent vertexes

      for (unsigned int k=1; k<=num_adj; k++) {
        counterInt += 1;
        unsigned int v = adjacents[k];

        //If distance to vertex is smaller than j, or if j doesnt have a path yet
        if ((result->distance[v] == -1 || result->distance[j]+1 < result->distance[v]) && result->distance[j] != -1) { 
          result->marked[v] = 1;   //Marks the path to v
          result->distance[v] = result->distance[j] + 1;//Saves distance to the vertex
          result->predecessor[v] = j;  
        }
      }
      free(adjacents);    
    }
  }

  printf("Iteraçoes Exteriores: %d; Iterações Medio: %d; Iterações Interiores: %d;\n", counterEx, counterMid, counterInt);

  return result;
}



void GraphBellmanFordAlgDestroy(GraphBellmanFordAlg** p) {
  assert(*p != NULL);

  GraphBellmanFordAlg* aux = *p;

  free(aux->marked);
  free(aux->predecessor);
  free(aux->distance);

  free(*p);
  *p = NULL;
}

// Getting the paths information

int GraphBellmanFordAlgReached(const GraphBellmanFordAlg* p, unsigned int v) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));

  return p->marked[v];
}

int GraphBellmanFordAlgDistance(const GraphBellmanFordAlg* p, unsigned int v) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));

  return p->distance[v];
}
Stack* GraphBellmanFordAlgPathTo(const GraphBellmanFordAlg* p, unsigned int v) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));

  Stack* s = StackCreate(GraphGetNumVertices(p->graph));

  if (p->marked[v] == 0) {
    return s;
  }

  // Store the path
  for (unsigned int current = v; current != p->startVertex;
       current = p->predecessor[current]) {
    StackPush(s, current);
  }

  StackPush(s, p->startVertex);

  return s;
}

// DISPLAYING on the console

void GraphBellmanFordAlgShowPath(const GraphBellmanFordAlg* p, unsigned int v) {
  assert(p != NULL);
  assert(v < GraphGetNumVertices(p->graph));

  Stack* s = GraphBellmanFordAlgPathTo(p, v);

  while (StackIsEmpty(s) == 0) {
    printf("%d ", StackPop(s));
  }

  StackDestroy(&s);
}

// Display the Shortest-Paths Tree in DOT format
void GraphBellmanFordAlgDisplayDOT(const GraphBellmanFordAlg* p) {
  assert(p != NULL);

  Graph* original_graph = p->graph;
  unsigned int num_vertices = GraphGetNumVertices(original_graph);

  // The paths tree is a digraph, with no edge weights
  Graph* paths_tree = GraphCreate(num_vertices, 1, 0);

  // Use the predecessors array to add the tree edges
  for (unsigned int w = 0; w < num_vertices; w++) {
    // Vertex w has a predecessor vertex v?
    int v = p->predecessor[w];
    if (v != -1) {
      GraphAddEdge(paths_tree, (unsigned int)v, w);
    }
  }

  // Display the tree in the DOT format
  GraphDisplayDOT(paths_tree);

  // Housekeeping
  GraphDestroy(&paths_tree);
}
