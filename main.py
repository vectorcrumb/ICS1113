import math, random, itertools
from gurobipy import *
from fastkml import kml
from shapely.geometry import LineString
from fastkml.styles import Style, LineStyle, StyleMap
import pandas as pd


# Set random seed
random.seed(1)

## Code extracted from Gurobi TSP tutorial @ gurobi.com/documentation/7.5/examples/tsp_py.html
## and @ examples.gurobi.com/traveling-salesman-problem/
# Lazy constraint to eliminate subtours
def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:
        # Create list of currently selected edges
        vals = model.cbGetSolution(model._vars)
        # print(model.cbGetSolution(model._vars).values())
        selected = tuplelist((i, j) for i, j in model._vars.keys() if vals[i, j] > 0.5)
        # Find shortest cycle in selected edges list
        tour = subtour(selected)
        if tour == []:
            active_nodes = list(set([n[0] for n in selected]))
            for i in active_nodes:
                model.cbLazy(model._vars.sum(i, '*') == 2)
        else:
            for i in tour:
                model.cbLazy(model._vars.sum(i, '*') == 2)
        if len(tour) < n:
            model.cbLazy(quicksum(model._vars[i, j] for i, j in itertools.combinations(tour, 2)) <= len(tour) - 1)
# Find shortest subtour given a tuplelist of edges
def subtour(edges):
    unvisited = list(range(n))
    cycle = range(n + 1)
    while unvisited:
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*') if j in unvisited]
        if len(cycle) > len(thiscycle):
            cycle = thiscycle
    return cycle

## Code extracted from Google TSP tutorial @ developers.google.com/optimization/routing/tsp/tsp
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
k = df_pars_cts.k[0]
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
if sum(di) > k:
    print("Product availability is less than minimum total demand. Problem is infeasible!")
    exit(-2)

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
# Take subdiagonal triangular matrix and fit into dictionary
# according to the Gurobi implementation of TSP.
dist = {(i, j): distances[i][j]
        for i in range(n) for j in range(i)}
# Definition of parameters
# K = k
# For nice simulation purposes
K = (sum(di) + sum(Di))/2
velocity = 3
T = tau * velocity
cost_to_dist = velocity/alphaval
# Typical cost values
# cost_values = [pi[i]*(1+ei[i][0]) for i in range(n)]
# Test with transforming costs to distance units (km)
cost_values = [pi[i]*(1+ei[i][0])*cost_to_dist for i in range(n)]
costs = {i: cost_values[i] for i in range(n)}

# Create models and variables
m = Model()
sales_vars = m.addVars(costs.keys(), lb=di, ub=Di, obj=costs, vtype=GRB.INTEGER, name='p')
variables = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='e')
for i, j in variables.keys():
    variables[j, i] = variables[i, j]  # Edge in opposite direction is the same variable

# Add 2nd degree constraint to each node
# m.addConstrs(variables.sum(i, '*') == 2 for i in range(n))
m.addConstr(variables.sum(1,'*') == 2, name='entry')
# Add total sale constraints to sales
# m.addConstr(sales_vars.sum(i) <= K for i in range(n))
m.addConstr(quicksum(sales_vars[i] for i in range(n)) <= K, name='max_stock')
# Add time constraint
m.addConstr(quicksum(variables[i,j]*distances[i][j] for i,j in dist) <= T, name='max_time')

# Optimize model
m._vars = variables
m.Params.lazyConstraints = 1
m.setObjective(sales_vars.prod(costs), GRB.MAXIMIZE)

m.optimize(subtourelim)

print("Difference from upper bound:", sum(list((list(m.getAttr('x', sales_vars).values())[i] - Di[i] for i in range(n)))))
print("Difference from lower bound:", sum(list((list(m.getAttr('x', sales_vars).values())[i] - di[i] for i in range(n)))))
print("Total sold:", sum(m.getAttr('x', sales_vars).values()))

vals = m.getAttr('x', variables)
selected = tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

tour = subtour(selected)
# assert len(tour) == n
print('Optimal tour: {}'.format(str(tour)))
print('Optimal cost: {}'.format(m.objVal))

ord_features = [sorted(features)[i] for i in tour]
print('Order of nodes:', *ord_features, sep='\n - ')


## Export created optimal path to new KML file
ns = '{http://www.opengis.net/kml/2.2}'
# Create styles for normal and highlighted lines
style_line_normal = Style(styles = [LineStyle(ns, color='8f781414', width=5)], id='line_normal')
style_line_highlight = Style(styles = [LineStyle(ns, color='8f781414', width=8)], id='line_highlight')
# Place both styles inside a StyleMap with an URL reference
styles_map = StyleMap(ns, id='line_styles', normal=style_line_normal, highlight=style_line_highlight)
# Create a document to store placemarks in and append it to the KML object
document = kml.Document(ns, 'opt_paths', 'Optimal Path', 'Draws the optimal TSP solution', styles=[styles_map])
kml_doc.append(document)
# Create a styled placemark to join every corresponding adjacent node and append to document
for i in range(len(ord_features)):
    place = kml.Placemark(ns, 'p_{}'.format(i), 'path_{}'.format(i))
    place.geometry = LineString([points[ord_features[i-1]], points[ord_features[i]]])
    place.styleUrl = '#line_styles'
    document.append(place)
# Display generated KML for debugging
# print(kml_doc.to_string(prettyprint=True))
# with open('Camino Optimo.kml', 'w') as file:
#     file.write(kml_doc.to_string(prettyprint=True))
print("\n\nFinished!")