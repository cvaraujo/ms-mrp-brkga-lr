/*
 * Created by Carlos - 07/09/2020
 */

#include "MSDecoder.h"
#include <string>

MSDecoder::MSDecoder() {}

MSDecoder::~MSDecoder() {}

void MSDecoder::loadInstance(string instance, string param, string preprocessing)
{
  int u, v;
  double delay, jitter, bandwidth, ldp, paramDelayToken,
      paramJitterToken, paramVariationToken, paramBandwidthToken;
  int delayInt, jitterInt, variationInt;
  string token;
  ifstream fileBoostGraph, fileParam, filePrep;

  // Read the QoS parameters
  fileParam.open(param, fstream::in);
  while (!fileParam.eof())
  {
    fileParam >> token;
    if (token == "Delay")
    {
      fileParam >> token;
      if (token == "variation")
      {
        fileParam >> token >> paramVariationToken;
        paramVariation = int(1e5 * paramVariationToken);
      }
      else
      {
        fileParam >> paramDelayToken;
        paramDelay = int(1e5 * paramDelayToken);
      }
    }
    if (token == "Jitter")
    {
      fileParam >> token >> paramJitterToken;
      paramJitter = int(1e6 * paramJitterToken);
    }
    if (token == "Bandwidth")
    {
      fileParam >> token >> paramBandwidthToken;
      paramBandwidth = int(paramBandwidthToken);
    }
  }

  // Read the graph
  fileBoostGraph.open(instance, fstream::in);
  while (!fileBoostGraph.eof())
  {
    fileBoostGraph >> token;
    if (token == "Nodes")
    {
      fileBoostGraph >> n, n++;
      arcs = vector<vector<Arc>>(n, vector<Arc>());
      costs = vector<vector<pair<int, int>>>(n, vector<pair<int, int>>(n));
      removeNodes = vector<bool>(n);
      removeArcs = vector<vector<bool>>(n, vector<bool>(n));

      // Read preprocessing file
      filePrep.open(preprocessing, fstream::in);
      while (!filePrep.eof())
      {
        filePrep >> token;
        if (token == "MVE")
        {
          filePrep >> u;
          removeNodes[u] = true;
        }
        if (token == "MAE")
        {
          filePrep >> u >> v;
          removeArcs[u][v] = true;
        }
      }
    }

    if (token == "Edges")
    {
      fileBoostGraph >> m;
      m *= 2;
    }
    if (token == "E")
    {
      fileBoostGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
      if (removeNodes[u] || removeNodes[v] || bandwidth < paramBandwidth)
      {
        m--;
        continue;
      }
      else
      {
        delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter);
        if (!removeArcs[u][v])
        {
          arcs[u].push_back(Arc(v, delayInt, jitterInt));
          costs[u][v] = make_pair(delayInt, jitterInt);
        }
        if (!removeArcs[v][u])
        {
          arcs[v].push_back(Arc(u, delayInt, jitterInt));
          costs[v][u] = make_pair(delayInt, jitterInt);
        }
      }
    }
    if (token == "Root")
    {
      fileBoostGraph >> root;
      // Update artificial node
      arcs[root].push_back(Arc(0, 1, 1));
      costs[root][0] = make_pair(1, 1);
    }
    if (token == "T")
      fileBoostGraph >> u, terminals.push_back(u), DuS.push_back(u);
  }

  for (int i = 1; i < n; ++i)
  {
    if (i != root)
    {
      // Update artificial arcs
      arcs[0].push_back(Arc(i, paramDelay, paramJitter));
      costs[0][i] = make_pair(paramDelay, paramJitter);
      // Setting the non-terminal nodes
      if (find(terminals.begin(), terminals.end(), i) == terminals.end())
      {
        nonTerminals.push_back(i), DuS.push_back(i);
      }
    }
  }
  m += (DuS.size() + 1);
  incumbent = terminals.size();
  cout << "Load graph successfully" << endl;
}

double MSDecoder::decode(const std::vector<double> &chromosome)
{
  bool found;
  Edge ed;
  int i, j, edNum = 0, actual, jitter, delay, bestCand, obj = 0, count = 0;
  vector<int> delayPaths = vector<int>(n), jitterPaths = vector<int>(n);
  vector<bool> notAttended = vector<bool>(n);
  double alpha;

  // Update the arcs cost
  BoostGraph graph(n);
  for (int i = 0; i < n; i++)
    for (auto e : arcs[i])
      add_edge(i, e.j, -chromosome[edNum++], graph);

  Vertex roots[2];
  roots[0] = roots[1] = root;
  vector<Edge> branching = vector<Edge>();
  // Compute the MST
  edmonds_optimum_branching<false, true, true>(graph,
                                               indexmap,
                                               weightmap,
                                               roots,
                                               roots + 1,
                                               back_inserter(branching));

  // Get the paths form the arborescence
  vector<int> predecessors = vector<int>(n);
  for (std::vector<Edge>::iterator ei = branching.begin(); ei != branching.end(); ++ei)
  {
    i = source(*ei, graph), j = target(*ei, graph);
    predecessors[j] = i;
  }

  // Return the computed arborescence
  string strArc, strHash = "";
  for (i = 0; i < n; i++)
  {
    if (i != root)
    {
      stringstream ss;
      ss << predecessors[i] << i;
      ss >> strArc;
      strHash += strArc;
    }
  }

  if (visitedArbs.find(strHash) != visitedArbs.end())
  {
    return visitedArbs[strHash];
  }

  // Set the serviced terminals
  for (auto k : DuS)
  {
    actual = k, jitter = 0, delay = 0;
    while (actual != root)
    {
      jitter += costs[predecessors[actual]][actual].second,
          delay += costs[predecessors[actual]][actual].first;
      actual = predecessors[actual];
    }
    delayPaths[k] = delay, jitterPaths[k] = jitter;
    if (delay > paramDelay || jitter > paramJitter)
      notAttended[k] = true, obj++;
  }

  if (float(rand() / RAND_MAX) < localSearchChance)
  {
    // Select the pivot terminal
    random_shuffle(terminals.begin(), terminals.end());
    int diffDelay, diffJitter, losts, minSelectedDelay = numeric_limits<int>::max(), selected = -1;
    bool canMove;

    for (auto k : terminals)
      if (!notAttended[k] && delayPaths[k] < minSelectedDelay)
        selected = k, minSelectedDelay = delayPaths[k];

    // Local Search procedure
    if (selected != -1)
    {
      vector<int> changed, delayPathAux = vector<int>(n),
                           jitterPathAux = vector<int>(n),
                           predecessorsAux = vector<int>(n);
      vector<bool> notAttendedAux = vector<bool>(n);
      vector<vector<int>> sub;

      // get the values of this path
      for (auto k : terminals)
      {
        if (k != selected)
        {
          bestCand = -1, delay = paramDelay + 1, jitter = paramJitter + 1;

          if (delayPaths[k] < (delayPaths[selected] - paramVariation) ||
              delayPaths[k] > (delayPaths[selected] + paramVariation))
          {

            for (auto arc : arcs[k])
            {
              i = arc.j;
              // Get the first (lowest delay) candidate to move
              if (i != predecessors[k] && k != predecessors[i])
              {
                if (delayPaths[i] + arc.delay >= (delayPaths[selected] - paramVariation) &&
                    delayPaths[i] + arc.delay <= (delayPaths[selected] + paramVariation))
                {
                  bestCand = i;
                  delay = delayPaths[i] + arc.delay, jitter = jitterPaths[i] + arc.jitter;
                }
              }
            }

            if (bestCand != -1)
            {
              // create the temporary vectors
              losts = 0;
              notAttended[k] = false, canMove = true;
              sub = vector<vector<int>>(n, vector<int>());
              changed = vector<int>();

              for (i = 0; i < n; i++)
              {
                delayPathAux[i] = delayPaths[i], jitterPathAux[i] = jitterPaths[i];
                predecessorsAux[i] = predecessors[i], notAttendedAux[i] = notAttended[i];
                sub[predecessors[i]].push_back(i);
              }

              // Evaluate the move
              diffDelay = delay - delayPathAux[k], diffJitter = jitter - jitterPathAux[k];
              delayPathAux[k] = delay, jitterPathAux[k] = jitter;
              predecessorsAux[k] = bestCand;

              // Propagates the QoS changes
              changed.push_back(k);
              while (!changed.empty())
              {
                actual = changed.back();
                changed.pop_back();

                for (auto j : sub[actual])
                {
                  delayPathAux[j] += diffDelay, jitterPathAux[j] += diffJitter;
                  if (delayPathAux[j] > paramDelay || jitterPathAux[j] > paramJitter)
                    notAttendedAux[j] = true;
                  else
                    notAttendedAux[j] = false;
                  changed.push_back(j);
                }
              }

              for (auto terminal : terminals)
                if (!notAttended[terminal] && notAttendedAux[terminal])
                  losts += 1;
                else if (notAttended[terminal] && !notAttendedAux[terminal])
                  losts -= 1;

              if (losts < 0)
                for (i = 0; i < n; i++)
                {
                  delayPaths[i] = delayPathAux[i], jitterPaths[i] = jitterPathAux[i],
                  predecessors[i] = predecessorsAux[i], notAttended[i] = notAttendedAux[i];
                }
            }
          }
        }
      }
    }
  }

  // Computing the difference in the metrics to make terminals serviced
  vector<int> delays = vector<int>();
  vector<double> minMetric = vector<double>(3, numeric_limits<int>::max());
  double auxMetric;

  for (auto k : terminals)
  {
    if (delayPaths[k] > paramDelay)
    {
      notAttended[k] = true;
      auxMetric = double(delayPaths[k] - paramDelay) / paramDelay;
      if (auxMetric < minMetric[0])
        minMetric[0] = auxMetric;
    }
    else if (jitterPaths[k] <= paramJitter)
    {
      notAttended[k] = false;
      delays.push_back(delayPaths[k]);
    }

    if (jitterPaths[k] > paramJitter)
    {
      notAttended[k] = true;
      auxMetric = double(jitterPaths[k] - paramJitter) / paramJitter;
      if (auxMetric < minMetric[1])
        minMetric[1] = auxMetric;
    }
  }

  // Applying the selection using the delay variation constraints
  sort(delays.begin(), delays.end());
  int actMax = 0, p2 = 0, maxi = 0, fd, larggest = 0, smallest = 0;
  for (i = 0; i < delays.size(); i++)
  {
    fd = delays[i];

    while (delays[p2] <= fd + paramVariation && p2 < delays.size())
      p2++;

    if ((p2 - i) > maxi)
    {
      maxi = (p2 - i);
      larggest = delays[p2];
      smallest = delays[i];
    }
  }

  // Delay variation penalty
  int penalty = 0, toUse = 0, maxDiff = 0;
  for (auto k : terminals)
  {
    toUse = -1;
    if (delayPaths[k] < smallest)
      toUse = larggest - delayPaths[k];
    else if (delayPaths[k] > larggest)
      toUse = delayPaths[k] - smallest;

    if (toUse > paramVariation)
    {
      auxMetric = double(toUse - paramVariation) / paramVariation;
      if (auxMetric < minMetric[2])
        minMetric[2] = auxMetric;
    }
  }

  count = terminals.size() - maxi;
  // alpha = *min_element(minMetric.begin(), minMetric.end());
  alpha = (minMetric[0] + minMetric[1] + minMetric[2]) / 3;

  if (count < incumbent)
  {
    // cout << "Delay: " << minMetric[0] << ", Jitter: " << minMetric[1] << ", Variation: " << minMetric[2] << endl;
    incumbent = count;
    // if (!checkSolution(predecessors, delayPaths, jitterPaths, notAttended, count))
    // {
    //   cout << "O PORRA" << endl;
    // }
  }

  // Save this arborescence in the map
  strHash = "";
  for (i = 0; i < n; i++)
  {
    if (i != root)
    {
      stringstream ss;
      ss << predecessors[i] << i;
      ss >> strArc;
      strHash += strArc;
    }
  }
  visitedArbs[strHash] = count * terminals.size() + alpha;
  // Return the value
  return (count * terminals.size() + alpha);
}

int MSDecoder::getM() const
{
  return m;
}

int MSDecoder::getN() const
{
  return n;
}

int MSDecoder::getIncumbent() const
{
  return incumbent;
}

void MSDecoder::loadBias(string instance, bool ls)
{
  ifstream file;
  file.open(instance, fstream::in);
  string token;
  int i, j, freq;
  vector<vector<int>> edges = vector<vector<int>>(n, vector<int>(n));
  mapping = vector<double>(m);

  int big = -1, les = 10000;
  while (!file.eof())
  {
    file >> token;

    if (ls)
    {
      if (token == "FU")
      {

        file >> i >> j >> freq;

        if (freq < les)
          les = freq;
        else if (freq > big)
          big = freq;

        edges[i][j] = freq;
      }
    }
    else
    {
      if (token == "FL")
      {
        file >> i >> j >> freq;
        if (freq < les)
          les = freq;
        else if (freq > big)
          big = freq;

        edges[i][j] = freq;
      }
    }
  }

  int ed = 0;
  for (i = 0; i < n; i++)
  {
    for (auto arc : arcs[i])
    {
      j = arc.j;
      mapping[ed++] = (edges[i][j] - les) / (big - les);
    }
  }
}

double MSDecoder::getBias(int i) const
{
  return mapping[i];
}

bool MSDecoder::checkSolution(vector<int> &predecessors, vector<int> &delayPaths, vector<int> &jitterPaths, vector<bool> &notAttended, int notAttendedTerm)
{
  int actual, jitter, delay;
  for (auto k : DuS)
  {
    actual = k, jitter = 0, delay = 0;
    while (actual != root)
    {
      jitter += costs[predecessors[actual]][actual].second,
          delay += costs[predecessors[actual]][actual].first;
      actual = predecessors[actual];
    }
    if ((delayPaths[k] != delay || jitterPaths[k] != jitter) ||
        ((delay > paramDelay || jitter > paramJitter) && (!notAttended[k])))
    {
      return false;
    }
  }

  // Computing the difference in the metrics to make terminals serviced
  vector<int> delays = vector<int>();
  for (auto k : terminals)
    if (!notAttended[k])
      delays.push_back(delayPaths[k]);

  // Applying the selection using the delay variation constraints
  sort(delays.begin(), delays.end());
  int actMax = 0, p2 = 0, maxi = 0, fd, larggest = 0, smallest = 0;
  for (int i = 0; i < delays.size(); i++)
  {
    fd = delays[i];

    while (delays[p2] <= fd + paramVariation && p2 < delays.size())
      p2++;

    if ((p2 - i) > maxi)
    {
      maxi = (p2 - i);
      larggest = delays[p2];
      smallest = delays[i];
    }
  }

  return notAttendedTerm != (terminals.size() - maxi) ? false : true;
}