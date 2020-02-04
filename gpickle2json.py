#!/usr/bin/env python
#-*- coding:utf-8 -*-

# Converts gpickle or other NetworkX-supported file format to .json file
# 
# This code is compatible with both networkx 1.9.1 and 2.x versions, but the version used to save the file should match!
#
# input argument 1: input file path with extension, e.g. /Users/username/data/MyGraph.gpickle
# 
# input argument 2: output file path with extension, e.g. /Users/username/data/MyGraph_converted.json
#
# Renaud Marquis @ FBMlab, March 2019

import os
import sys
import json
import networkx as nx
from networkx.readwrite import json_graph

G = nx.read_gpickle(sys.argv[1])
J = json_graph.node_link_data(G)
S = json.dumps(J,indent=4)

OutFile = open(sys.argv[2],'w')
OutFile.write(S)
OutFile.close()
