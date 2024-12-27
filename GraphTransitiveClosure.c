//
// Algoritmos e Estruturas de Dados --- 2024/2025
//
// Joaquim Madeira - Dec 2024
//
// GraphTransitiveClosure - Transitive Closure of a directed graph
//

// Student Name : Guilherme Carreira
// Student Number : 120159
// Student Name : Daniel Couto
// Student Number : 50138

/*** COMPLETE THE GraphComputeTransitiveClosure FUNCTION ***/

#include "GraphTransitiveClosure.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "GraphBellmanFordAlg.h"
#include "instrumentation.h"

// Compute the transitive closure of a directed graph
// Return the computed transitive closure as a directed graph
// Use the Bellman-Ford algorithm
Graph* GraphComputeTransitiveClosure(Graph* g) {
  assert(g != NULL);             
  assert(GraphIsDigraph(g));           
  assert(GraphIsWeighted(g) == 0);       

  unsigned int numVertices = GraphGetNumVertices(g); // numero de vertices do grafo
  Graph* transitiveClosure = GraphCreate(numVertices, 1, 0); //cria um grafo direcionado vaziozz

  // Itera por cada vertice
  for (unsigned int u = 0; u < numVertices; u++) {
    //Utiliza o algoritmo de BellmanFord para obter os vertices alcançaveis
    GraphBellmanFordAlg* result = GraphBellmanFordAlgExecute(g, u);
    assert(result != NULL);

    //verifica se o vertice é alcançavel
    for (unsigned int v = 0; v < numVertices; v++) {
      if (u != v && GraphBellmanFordAlgReached(result, v)) {
        // em caso positivo, adiciona uma aresta direta
        GraphAddEdge(transitiveClosure, u, v);
      }
    }

    //Liberta memória da variavel result
    GraphBellmanFordAlgDestroy(&result);
  }
  return transitiveClosure;
}