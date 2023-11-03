import pyautocad
from shapely.geometry import Polygon
from rtree import index
import ezdxf
import time
from cvxopt import solvers, matrix
import cvxopt
import pyautocad
from shapely.geometry import Polygon
from rtree import index
import ezdxf
import time
from cvxopt import solvers, matrix
import numpy as np
import pygame
from shapely.geometry import Point
import shapely.geometry.polygon 
from matplotlib.patches import Polygon
from sympy import linear_eq_to_matrix, symbols
import math
import matplotlib.pyplot as plt
from colormath.color_objects import LabColor, XYZColor
from colormath.color_conversions import convert_color
import matplotlib.transforms as mtransforms
from matplotlib.colors import Normalize
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from scipy.optimize import minimize, LinearConstraint, linprog, least_squares, root,fsolve
from scipy.optimize import least_squares
from sympy import Matrix

from pyomo.environ import ConcreteModel, RangeSet, Var, Objective, Constraint
from pyomo.opt import SolverFactory

import json

from data_entry import data_entry
from autocad_data_extraction import autocad_data_extraction
from loading_extraction import loading_extraction
from problem_building_up import problem_building_up
from static_approach_cvxopt import static_approach_cvxopt
from kinematic_approach_cvxopt import kinematic_approach_cvxopt
from iterative_approach_cvxopt import iterative_approach_cvxopt
from mix_data import mix_data
from stat_ipopt import stat_ipopt
from ipopt_data import ipopt_data
from kin_ipopt import kin_ipopt
from iterative_ipopt import iterative_ipopt 

friction_ang=24
fname= ["","C:\\Users\\mmoussa\\Desktop\\backup_before_osserain_trip\\Python-projects\\mechanic_data_pont.txt"]
dxf_file ="C:\\Users\\mmoussa\\Desktop\\backup_before_osserain_trip\\Python-projects\\pont (6).dxf" 
cohesion=0
max_nl=9999
min_nl=0
max_nl_1 =9999
max_nl_2=9999
dilantancy=24

data = data_entry( cohesion, friction_ang , dilantancy, max_nl , min_nl,max_nl_1 ,max_nl_2,fname,dxf_file  )

data = autocad_data_extraction(data,pont = True,culee_bloc = False,arch = False)
data = loading_extraction(data,remplissage = True)
data = problem_building_up(data,case=3,arch =False, pont = True, case_pont = 1)


#Static approach

#data = static_approach_cvxopt(data,graph_rep = True ,scale_factor=10)      

#Kinematic approach

#data = kinematic_approach_cvxopt(data,graph_rep = False,scale_factor=6) 

data = iterative_approach_cvxopt(data, alpha=0.3, beta=0.6, tolerance=0.001 ,graph_rep = True,scale_factor =20, rescale =False )


data  = mix_data(data)

if False:
    def custom_encoder(obj):
        if isinstance(obj, np.ndarray):
            return np.array(obj).tolist()  # Convert NumPy array to a list
        if isinstance(obj, cvxopt.base.matrix):
            obj = np.array(obj)
            return obj.tolist()  # Convert cvxopt matrix to a list
        #Add more custom encoding logic for other data types if needed
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


    # Save the dictionary to a JSON file using the custom encoder
    with open("data.json", "w") as file:
        json.dump(data, file, default=custom_encoder)

 
A_df, b_df, c_df, len_df = ipopt_data(data["Amix"],data["bmix"],data["cmix"],data["edges_length"])

#statsol,x,nb,nc,obj,result = stat_ipopt(data,A_df, b_df, c_df, len_df,graph_rep = True)

#statsol,x,nb,nc,obj,result = kin_ipopt(data,A_df, b_df, c_df, len_df)

#iterative_ipopt(data,graph_rep = True)

plt.show(block = True)