import os
import pandas as pd  # for pandas.DataFrame objects (https://pandas.pydata.org/)
import numpy as np  # for numpy.matrix objects (https://numpy.org/)
from amplpy import AMPL
import time
import math
from ipopt_data import ipopt_data
from stat_ipopt import stat_ipopt

def iterative_ipopt(data,graph_rep):

    A_df, b_df, c_df, len_df = ipopt_data(data["Amix"],data["bmix"],data["cmix"],data["edges_length"])
    statsol,x,nb,nc,obj,result,ampl = stat_ipopt(data,A_df, b_df, c_df, len_df,graph_rep = False)
    
    Amix = data["Amix"].copy()
    bmix = data["bmix"].copy()
    length = data["edges_length"].copy()
    cmix = data["cmix"].copy()         
            
      
            
    strength = data["strength"]

    activeEdgesInd = data["activeEdgesInd"]
    tol =1
    itr = 0
    alpha=0.2 
    beta=0.2 
    #tolerance=0.001
    tolerance=0.001

    #Cstat = np.array(Cstat)
    for element in range(0,nc,3):
        Amix[element*2+3*nb][1+element] = abs(Amix[element*2+3*nb][1+element]*alpha)

        Amix[element*2+1+3*nb][1+element] = abs(Amix[element*2+1+3*nb][1+element]*alpha)

        start = time.time()


        
    while tol >tolerance: 

        # Check if the solver found an optimal solution
        if result  == "solved":
            
            if itr>1:
                normalvalues_old = normalvalues
            else:
                normalvalues_old = np.zeros((int((nc)/3),1))

            statsolver_old = obj    
            shearvalues = np.empty((0,1))
            normalvalues = np.empty((0,1))
            momentvalues = np.empty((0,1))
            rup_test = np.empty((0,1))
            for element in range(0,len(statsol)-2,3):
                shearvalues = np.append(shearvalues,statsol[element])
                normalvalues = np.append(normalvalues,statsol[element+1])
                momentvalues = np.append(momentvalues,statsol[element+2])


        for element in range(0,int((nc)/3),1):
            ind = np.where( strength[:,0]  == str(int(activeEdgesInd[element])) )[0]

            if itr >1:
                bmix[3*nb + element*6] = 0.00001*max(normalvalues) + data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element])

                bmix[3*nb + element*6+1] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element])
            else:

                bmix[3*nb + element*6] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*normalvalues[element] 

                bmix[3*nb + element*6+1] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*normalvalues[element]


        #Cstat = matrix(Cstat)

        statsol_old = statsol
        A_df, b_df, c_df, len_df = ipopt_data(Amix,bmix,cmix,length)
        statsol,x,nb,nc,obj,result,ampl  = stat_ipopt(data,A_df, b_df, c_df, len_df ,graph_rep = False)
        
        end = time.time()

        if result == "solved":

            if itr >1:
                tol = abs(obj - statsolver_old)/abs(obj )
                print("tolerance", tol)
        else:
            statsol = statsol_old
            
        if alpha <= 0.0001:
            alpha = 0.0001
        else:
            alpha = alpha*0.7
       

        for element in range(0,nc,3):
            Amix[element*2+3*nb][1+element] = abs(Amix[element*2+3*nb][1+element]*0.7)

            Amix[element*2+3*nb][1+element] = abs(Amix[element*2+3*nb][1+element]*0.7)


        itr = itr+1
    
    #x = ampl.get_variable("x")
    
    
    statsol = np.empty((0,nc*3+1));
    print("hereeeee")
    # Print it
    for i in range(0,3*nc+1):
        statsol = np.append(statsol,x[i+1].value()) 
    
    
    constraint1 = ampl.get_constraint("equalities1")
    constraint1_values = constraint1.getValues()
    ub = np.empty(3*nb)
    
    ind = 0
    for i in constraint1_values:
        ub[ind] = i[1]
        ind = ind+1
    
    crushing_top = ampl.get_constraint("crushing_top")
    crushing_top_values = crushing_top.getValues()
    
    gamma1 = np.empty(nc)
    ind = 0;
    for i in crushing_top_values:
        gamma1[ind] = i[1]
        ind = ind+1
    
    
    
    crushing_bot = ampl.get_constraint("crushing_bot")
    crushing_bot_values = crushing_bot.getValues()
    
    ind = 0
    
    gamma2 = np.empty(nc)
    for i in crushing_bot_values:
        
        gamma2[ind] = i[1]
        ind = ind +1
    
    
    inequalities2 = ampl.get_constraint("inequalities2")
    inequalities2_values = inequalities2.getValues()
    #print(inequalities2_values)
    d = np.empty(6*nc)
    ind =0
    for i in inequalities2_values:
        d[ind] = i[1]
        ind = ind+1
    
    # Create an empty list to store the vectors
    gamma = np.array([])

    # Add your vectors to the list
    gamma= np.append(gamma,gamma1)
    gamma= np.append(gamma,gamma2)
    

    print("I am here")
    if graph_rep == True:
        from graphical_rep import graphical_rep
        graphical_rep(data,statsol,ub,d,gamma, 20,graph_title = "Static approach")