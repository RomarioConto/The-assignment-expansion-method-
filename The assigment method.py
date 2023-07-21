!pip install -q pyomo
from pyomo.environ import *
from pyomo.opt import *

# Set the solver
solver = solvers.SolverFactory('gurobi')

vec=[-1,1]     #first factor levels
vec2=[-1,0,1]  #Second factor levels
k = 2          #number of factors
beta = 5       #maximum number of level changes allowed
  
import itertools
import numpy as np

# factors 
factors = {i for i in range(k)}

# Corridas
runs = itertools.product(vec, vec2)
runs = [i for i in runs]
runs = {j+1:runs[j] for j in range(len(runs))}
positions = {i+1 for i in range(len(runs))}
runs

response = {i:i for i in positions}
distances = {}
for k1,v1 in runs.items():
  for k2,v2 in runs.items():
    distances[(k1,k2)] = [a == b for a, b in zip(v1, v2)].count(False)
arcos = set(distances.keys())

model = ConcreteModel(name="AssignmentPriority")

# Sets
model.POSITIONS = Set(ordered = False, initialize=positions)                     
model.RUNS = Set(ordered = False, initialize=positions)                          
model.FACTORS = Set(ordered = False, initialize=factors)   
model.ARCOS = Set(ordered = False, initialize=arcos)

# define parámetros
model.runs = Param(model.RUNS, initialize=runs, within=Any)
model.response = Param(model.POSITIONS, initialize=response, within=Any)
model.distance = Param(model.ARCOS, initialize=distances, within=Any)
model.beta = Param(initialize=beta, within=Any)

# Declarión variables de decisión
model.x = Var(model.POSITIONS, model.RUNS, domain=Binary)
model.s = Var(model.FACTORS, domain=NonNegativeReals)
model.s_max = Var(domain=NonNegativeReals)
model.p = Var(model.POSITIONS, model.RUNS, model.RUNS, domain=Binary)


# Objetive Funtion
def obj_rule(model):
    return model.s_max 
model.objetive = Objective(sense=minimize, rule=obj_rule)

# linearizing maximum
def s_max_c(model, k):
    return  model.s[k] <= model.s_max
model.s_max_c = Constraint(model.FACTORS,rule=s_max_c)

# linearizing absolute value
def linear1(model, k):
    return sum(model.response[i]*model.runs[j][k]*model.x[i,j]for i in model.POSITIONS for j in model.RUNS) <= model.s[k]
model.linear1 = Constraint(model.FACTORS,rule=linear1)
def linear2(model, k):
    return -sum(model.response[i]*model.runs[j][k]*model.x[i,j]for i in model.POSITIONS for j in model.RUNS) <= model.s[k]
model.linear2 = Constraint(model.FACTORS,rule=linear2)


# one run per position
def run_1(model,i):
    return sum(model.x[i,j] for j in model.RUNS) == 1
model.run_1 = Constraint(model.POSITIONS,rule=run_1)

# one position per run
def position_1(model,j):
    return sum(model.x[i,j] for i in model.POSITIONS) == 1
model.posicion_1 = Constraint(model.RUNS,rule=position_1)

# shorten distance
def distance_c(model):
    return sum(model.distance[(k,j)]*model.p[i,k,j] for i in model.POSITIONS for k in model.RUNS for j in model.RUNS) <= model.beta
model.distancia_c = Constraint(rule=distance_c)

# liniealizing distance
def product1(model,i,k,j):
    return model.p[i,k,j] <= model.x[i,k]
model.product1 = Constraint(model.POSITIONS, model.RUNS, model.RUNS, rule=product1)

def product2(model,i,k,j):
  if i < max(positions):
    return model.p[i,k,j] <= model.x[i+1,j]
  else:
    return Constraint.Skip
model.product2 = Constraint(model.POSITIONS, model.RUNS, model.RUNS, rule=product2)

def product3(model,i,k,j):
  if i < max(positions):
    return model.p[i,k,j] >= model.x[i,k] + model.x[i+1,j] -1
  else:
    return Constraint.Skip
model.product3 = Constraint(model.POSITIONS, model.RUNS, model.RUNS, rule=product3)

# Cuts 
# def cuts1(model,i,k,j):
#   if i < max(posiciones) and model.distancia[(k,j)]>1:
#     return model.x[i,k] + model.x[i+1,j] <= 1
#   else:
#     return Constraint.Skip
# model.cuts1 = Constraint(model.POSICIONES, model.CORRIDAS, model.CORRIDAS, rule=cuts1)

results = solver.solve(model, tee=True)
term_cond = results.solver.termination_condition
print ("Termination condition={}".format(term_cond))
if term_cond== TerminationCondition.optimal:
# Get funtion objetive
  obj_val = model.objetive.expr()
  print("Objective function value: ", obj_val)
  
def get_results(model):
  results_dict = {}
  # get order  
  res = model.x.get_values()
  order = []
  for key, value in res.items():
    if value != None and value > 0:
      order.append(key)
  results_dict['solution'] = order
  run_order = [runs[val[1]] for val in order]
  results_dict['extended_solution'] = run_order
  # biases
  biases = [0 for f in factors]
  for factor in range(len(factors)):
    for i in range(len(run_order)):
      biases[factor]+=(i+1)*run_order[i][factor]
  results_dict['biases'] = biases
  results_dict['max_bias'] = max(list(map(abs, biases)))
  # level changes between runs
  distanceT = 0
  for i in range(len(order)-1):
    distanceT += distances[(order[i][1], order[i+1][1])]
  distanceT
  results_dict['distance'] = distanceT
  
  return results_dict

results = get_results(model)
print("Optimal run order:")
print(results['extended_solution'])
print("factor biases: ", results['biases'])
print("Maximum Bias: ", results['max_bias'])
print("Maximum distance: ", results['distance'])

