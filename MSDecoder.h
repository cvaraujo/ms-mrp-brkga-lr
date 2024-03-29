/*
 * Created by Carlos - 07/09/2020
 */

#ifndef MSDECODER_H
#define MSDECODER_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/config.hpp>
#include "edmonds_optimum_branching.hpp"

#include <vector>
#include <iostream>
#include <iomanip>
#include <bits/ios_base.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

using namespace boost;
using namespace std;

struct Arc
{
  Arc(int dest, int del, int jit) : j(dest), delay(del), jitter(jit) {}
  // Destine vertex, delay and jitter of the edge
  int j, delay, jitter;
};

class MSDecoder
{

public:
  typedef adjacency_list<vecS, vecS, directedS, property<vertex_index_t, int>, property<edge_weight_t, double>> BoostGraph;
  typedef graph_traits<BoostGraph>::edge_descriptor Edge;
  typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;
  property_map<BoostGraph, edge_weight_t>::type weightmap;
  property_map<BoostGraph, vertex_index_t>::type indexmap;

  int n, m, root, paramDelay, paramJitter, paramVariation, paramBandwidth;
  int incumbent;
  // BoostGraph graph;
  vector<vector<Arc>> arcs;
  vector<vector<pair<int, int>>> costs;
  vector<int> terminals = vector<int>(), nonTerminals = vector<int>(), DuS = vector<int>();
  vector<bool> removeNodes;
  vector<vector<bool>> removeArcs;
  vector<double> mapping;
  map<string, double> visitedArbs;
  float localSearchChance = 0.1;

  MSDecoder();

  ~MSDecoder();

  void loadInstance(string instance, string param, string preprocessing);

  int getN() const;

  int getM() const;

  int getIncumbent() const;

  double decode(const std::vector<double> &chromosome);

  int decodeFinal(const std::vector<double> &chromosome);

  void loadBias(string instance, bool ls);

  double getBias(int i) const;

  double already_computed();

  bool checkSolution(vector<int> &predecessors, vector<int> &delayPaths, vector<int> &jitterPaths, vector<bool> &notAttended, int notAttendedTerm);
};

#endif
