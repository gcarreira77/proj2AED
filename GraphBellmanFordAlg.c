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

GraphBellmanFordAlg* GraphBellmanFordAlgExecute(Graph* g, unsigned int startVertex) {
  assert(g != NULL);  // Ensure the graph is not NULL
  assert(startVertex < GraphGetNumVertices(g));  // Ensure the start vertex is valid
  assert(GraphIsWeighted(g) == 0); // Assuming the graph is unweighted

  // Allocate memory for the result structure
  GraphBellmanFordAlg* result = (GraphBellmanFordAlg*)malloc(sizeof(struct _GraphBellmanFordAlg));
  assert(result != NULL);  // Ensure memory allocation was successful

  unsigned int numVertices = GraphGetNumVertices(g);  // Get the number of vertices
  result->graph = g;  // Store the graph in the result structure
  result->startVertex = startVertex;  // Store the start vertex

  // Initialize arrays for marked, distance, and predecessor
  result->marked = (unsigned int*)malloc(numVertices * sizeof(unsigned int));
  result->distance = (int*)malloc(numVertices * sizeof(int));
  result->predecessor = (int*)malloc(numVertices * sizeof(int));

  // Initialize arrays: all vertices are unvisited, no path initially, no predecessors
  for (unsigned int i = 0; i < numVertices; i++) {
    result->marked[i] = 0;  // Not visited
    result->distance[i] = -1;  // No path initially
    result->predecessor[i] = -1;  // No predecessor initially
  }

  result->distance[startVertex] = 0;  // The distance to the start vertex is 0

  // Bellman-Ford algorithm: Relax edges up to numVertices-1 times
  for (unsigned int i = 0; i < numVertices - 1; i++) {
    for (unsigned int v = 0; v < numVertices; v++) {
      unsigned int* neighbors = GraphGetAdjacentsTo(g, v);  // Get the neighbors of vertex v
      
      unsigned int numNeighbors;

      // Check if the graph is directed or not
      if (GraphIsDigraph(g)) {
        // For directed graphs, use the out-degree (number of outgoing edges)
        numNeighbors = GraphGetVertexOutDegree(g, v);
      } else {
        // For undirected graphs, use the vertex degree (number of edges connected)
        numNeighbors = GraphGetVertexDegree(g, v);
      }

      // Relax the edges (v, w)
      for (unsigned int j = 0; j < numNeighbors; j++) {
        unsigned int w = neighbors[j];  // Get the neighbor w of vertex v
        // Relaxation: if there's a shorter path to w through v
        if (result->distance[v] != -1 && (result->distance[w] == -1 || result->distance[w] > result->distance[v] + 1)) {
          result->distance[w] = result->distance[v] + 1;  // Update the distance to w
          result->predecessor[w] = v;  // Update the predecessor of w
        }
      }
    }
  }

  return result;  // Return the Bellman-Ford result
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
