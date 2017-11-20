#!/usr/bin/python

import sys, math, random, itertools
from fastkml import kml
from fastkml.styles import Style, LineStyle, StyleMap
from shapely.geometry import LineString
import pandas as pd
from gurobipy import *
from random import choice

random.seed(1)

# Code extracted from Google TSP tutorial @ developers.google.com/optimization/routing/tsp/tsp
# Haversine function
def haversine(angle):
    return math.sin(angle / 2) ** 2
# Calculate distances between two geographical points
def distance(lat1, long1, lat2, long2):
    # Note: The formula used in this function is not exact, as it assumes
    # the Earth is a perfect sphere.
    # Mean radius of Earth in miles
    radius_earth = 3959
    km_per_mile = 1.60934
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi / 180.0
    phi1 = lat1 * degrees_to_radians
    phi2 = lat2 * degrees_to_radians
    lambda1 = long1 * degrees_to_radians
    lambda2 = long2 * degrees_to_radians
    dphi = phi2 - phi1
    dlambda = lambda2 - lambda1
    # Calculate distance
    a = haversine(dphi) + math.cos(phi1) * math.cos(phi2) * haversine(dlambda)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius_earth * c
    return d * km_per_mile

# Callback - use lazy constraints to eliminate sub-tours
def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._varsx)
        selected = tuplelist((i,j) for i,j in model._varsx.keys() if vals[i,j] > 0.5)
        tour = subtour(selected)
        if len(tour) < n:
            not_tour = [node for node in list(visit_costs.keys()) if node not in tour]
            for jnu in not_tour:
                iters = itertools.combinations(tour, 2)
                model.cbLazy(quicksum(model._varsx[i,j] for i,j in iters)
                             <= len(tour) - 1 + big_M *
                             (len(tour) + 1 - sum(model._varsy[i] for i in iters) - model._varsy[jnu]))


# Given a tuplelist of edges, find the shortest subtour
def subtour(edges):
    unvisited = list(range(n))
    cycle = range(n+1) # initial length has 1 more city
    while unvisited: # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i,j in edges.select(current,'*') if j in unvisited]
        if len(cycle) > len(thiscycle):
            cycle = thiscycle
    return cycle


## Import geo data from KML files
# Open file with KML data points and fit into features dictionary
with open("Nodos de Venta.kml", 'rt', encoding='utf-8') as mf:
    doc = mf.read().encode('utf-8')
    kml_doc = kml.KML()
    kml_doc.from_string(doc)
# Fit KML object into dictionaries with info
features = {feat.name: [feat._geometry.geometry.x, feat._geometry.geometry.y] for feat in kml_doc._features[0]._features}
points = {feat.name: feat._geometry.geometry for feat in kml_doc._features[0]._features}

## Import parameter values from CSV files
# CSV file names
par = "parameters.csv"
par_cts = "parameters_constants.csv"
# Read CSVs into DFs
df_pars = pd.read_csv(par)
df_pars_cts = pd.read_csv(par_cts)
# Obtain constants
alphaval = df_pars_cts.alpha[0]
betaval = df_pars_cts.beta[0]
K = df_pars_cts.bk[0]
K_min = df_pars_cts.sk[0]
tau = df_pars_cts.tau[0]
# Obtain parameters
pi = list(df_pars.pi)
di = list(df_pars.di)
Di = list(df_pars.Di)
ci2 = list(df_pars.ci2)
epi = [round(val, 2) for val in list(df_pars.epi)]
efi = [round(val, 2) for val in list(df_pars.efi)]
evi = [round(val, 2) for val in list(df_pars.evi)]
en = [0 for _ in range(len(efi))]
# Create list of event modifiers throughout week
ei = [[random.choice([epi[i], efi[i], en[i], evi[i], en[i]]) for j in range(5)] for i in range(len(en))]
# Check feasibility within parameters
if sum(di) > K:
    print("Product availability is less than minimum total demand. Problem is infeasible!")
    exit(-2)


# Definition of parameters
# K = (sum(di) + sum(Di))/2
big_M = 5 * K

## Prepare navigation variables
# Calculate distances in between nodes. Nodes are alphabetically ordered by name. Matrix must
# be symmetrical and have null diagonal elements (distance from node i to node i)
distances = []
for k1 in sorted(features):
    k1dists = []
    for k2 in sorted(features):
        k1dists.append(distance(features[k1][1], features[k1][0], features[k2][1], features[k2][0]))
    distances.append(k1dists)
# Check parameters and nodes
n = len(distances)
if n != len(en):
    print("Mismatch between parameter and nodes dimensionality!")
    exit(-1)

# Distance costs - c_{ij}^1
dist = {(i, j): distances[i][j]*alphaval*betaval for i in range(n) for j in range(i)}
# Selling benefits - Accompany z_i
cost_values = [pi[i]*(1+choice(ei[i])) for i in range(n)]
costs = {i: cost_values[i] for i in range(n)}
# Visiting costs - c_{i}^2
visit_costs = {i: ci2[i]*alphaval for i in range(n)}

# Create a model object and create variables
m = Model()
sales_vars = m.addVars(costs.keys(), lb=0, ub=Di, obj=costs, vtype=GRB.INTEGER, name='z')
arc_vars = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='x')
visit_vars = m.addVars(visit_costs.keys(), obj=visit_costs, vtype=GRB.BINARY, name='y')
# Correction: Edges in opposite direction are equal
for i,j in arc_vars.keys():
    arc_vars[j, i] = arc_vars[i, j]

## Problem restrictions
# R1 Edge limit for every visited node
for i in range(n):
    n_noi = [x for x in range(n) if x != i]
    m.addConstr(sum(arc_vars[i, j] for j in n_noi) == visit_vars[i])
    m.addConstr(sum(arc_vars[j, i] for j in n_noi) == visit_vars[i])

m.addConstrs(big_M*visit_vars[i] >= sales_vars[i] for i in range(n))
m.addConstr(visit_vars[1] >= 1, name='entry')
m.addConstrs(sales_vars[i] >= di[i] * visit_vars[i] for i in range(n))
m.addConstr(sum(sales_vars[i] for i in range(n)) <= K, name='max_stock')
m.addConstr(sum(sales_vars[i] for i in range(n)) >= K_min, name='min_stock')
# R2 Subtour elimination: Occurs in callback while optimizing
# R3 Add stock constraints to sales

# R5 Add lower bound upon activation constraint

# R6 Total sale time constraint
# m.addConstr(betaval * quicksum(arc_vars[i, j] * distances[i][j] for i, j in dist) +
#             quicksum(visit_vars[i] * visit_costs[i] for i in range(n)) <= tau, name='max_time')
# R7 Start at node 0

# R8 Activation constraints


# Optimize model
m._varsx = arc_vars
m._varsy = visit_vars
m.Params.lazyConstraints = 1
m.setObjective(sales_vars.prod(costs) - arc_vars.prod(dist) - visit_vars.prod(visit_costs), GRB.MAXIMIZE)

m.optimize(subtourelim)

vals = m.getAttr('x', arc_vars)
selected = tuplelist((i,j) for i,j in vals.keys() if vals[i,j] > 0.5)

tour = subtour(selected)
# assert len(tour) == n

print('')
print('Optimal tour: {}'.format(str(tour)))
print('Optimal cost: {}'.format(m.objVal))
print('')
