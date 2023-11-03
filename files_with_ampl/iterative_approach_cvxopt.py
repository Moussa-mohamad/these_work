def iterative_approach_cvxopt(data, alpha, beta, tolerance ,graph_rep,scale_factor,rescale):
    
    import time
    from cvxopt import solvers, matrix
    import numpy as np
    import math
    
    title = "Iterative approach"
    print("=" * 50)
    print(title.center(50))
    print("=" * 50)
    
    Astat = data["Astat"]  
    bstat = data["bstat"] 
    Cstat = data["Gstat"] 
    cohstat = np.array(data["hstat"]).copy() 
    cohstat = matrix(cohstat)
    strength = data["strength"]
    
    
    edges = data["edges"]
    allpolygons = data["allpolygons"]
    blocks = data["blocks"]
    activeEdgesInd = data["activeEdgesInd"]
    zstat = data["cstat"]
    Dead_loadsvect = data["dead_loadvect"]
    Live_loadsvect = data["live_loadvect"]
    B = data["equilibrium_matrix"] 
    coh = data["hstat"]
    
    Astat =data["Astat"]
    bstat = data["bstat"] 
    Cstat = data["Gstat"] 
    cohstat = data["hstat"] 
    zstat = data["cstat"] 

    
    statsolver = solvers.lp(zstat,Cstat,cohstat,Astat,-bstat)
    statsol = statsolver['x']
    
    if not rescale:
        tol =1
        itr = 0
        
        Cstat = np.array(Cstat)
        for element in range(0,Astat.size[1]-1,3):
            Cstat[element*2][1+element] = abs(Cstat[element*2][1+element]*alpha)

            Cstat[element*2+1][1+element] = abs(Cstat[element*2+1][1+element]*alpha)
        
        start = time.time()
        
        while tol >tolerance: 

            
            if statsolver['x']  != None:
                if itr>1:
                    normalvalues_old = normalvalues
                else:
                    normalvalues_old = np.zeros((int((len(statsol)-1)/3),1))
                                                
                statsolver_old = statsolver['primal objective']    
                shearvalues = np.empty((0,1))
                normalvalues = np.empty((0,1))
                momentvalues = np.empty((0,1))
                rup_test = np.empty((0,1))
                for element in range(0,len(statsol)-2,3):
                    shearvalues = np.append(shearvalues,statsol[element])
                    normalvalues = np.append(normalvalues,statsol[element+1])
                    momentvalues = np.append(momentvalues,statsol[element+2])


            for element in range(0,int((Astat.size[1]-1)/3),1):
                ind = np.where( strength[:,0]  == str(int(activeEdgesInd[element])) )[0]

                if itr >1:
                    cohstat[element*6] = 0.00001*max(normalvalues) + data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element])

                    cohstat[element*6+1] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element])
                else:

                    cohstat[element*6] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*normalvalues[element] 

                    cohstat[element*6+1] = 0.00001*max(normalvalues) +data["hstat"][element*6] + (1+alpha)*math.tan(math.radians(float(strength[ind[0]][2])))*normalvalues[element]


            Cstat = matrix(Cstat)
            

            statsolver = solvers.lp(zstat,Cstat,cohstat,Astat,-bstat)
            end = time.time()
            
            if statsolver['x'] != None:
                statsol = statsolver['x']
            #print(Astat.size,bstat.size,Cstat.size,cohstat.size,zstat.size)
                print(statsolver['primal objective'])
                if itr >1:
                    tol = abs(statsolver['primal objective']-statsolver_old)/abs(statsolver['primal objective'])
                    print("tolerance", tol)
            
            if alpha <= 0.0001:
                alpha = 0.0001
            else:
                alpha = alpha/2
            Cstat = np.array(Cstat)
            for element in range(0,Astat.size[1]-1,3):
                Cstat[element*2][1+element] = abs(Cstat[element*2][1+element]/2)

                Cstat[element*2+1][1+element] = abs(Cstat[element*2+1][1+element]/2)
            
            
            itr = itr+1
            
            
        lagrange_stateq = np.array(statsolver['y'])
        lagrange_statineq = np.array(statsolver['z'])
        fourth = np.append(lagrange_stateq,lagrange_statineq)
        data["statvect_itr"] = statsolver['x']
        data["kinvect_itr"] = fourth
    else:
        statsolver['x'] = data["statvect_itr"]
        fourth = data["kinvect_itr"]
    
    lagrange_stateq = []
    for i in statsolver['y']:
        lagrange_stateq.append(i)
        
    statvect = statsolver['x']
    lagrange_statineq = np.array(statsolver['z'])
   
    d = np.empty(6*(Astat.size[1]-1))
    ind =0
    for i in statsolver['z']:
        d[ind] = i
        ind = ind+1
        
    
    fourth = np.append(lagrange_stateq,lagrange_statineq)
    
    if graph_rep:
        from graphical_rep import graphical_rep
        graphical_rep(data,statvect,lagrange_stateq,d,[],scale_factor,graph_title = "Static approach")
    
    end = time.time()
    print("Iterative method:", end-start, "seconds")
    effcoh_energy = 0
    for element in range(0,len(lagrange_statineq),6):
        effcoh_energy = effcoh_energy + 0.00001*max(normalvalues)*abs(lagrange_statineq[element] - lagrange_statineq[element+1])
    print("Energy generated by C0 ", effcoh_energy)
    
    Dead_energie = 0
    Dead_energie = sum([Dead_energie  + lagrange_stateq[element]*Dead_loadsvect[element] for element in range(0,len(lagrange_stateq))])
    print("Dead loads energie ", + Dead_energie)
    
    Live_energie = 0
    Live_energie = sum([Live_energie  + lagrange_stateq[element]*Live_loadsvect[element] for element in range(0,len(lagrange_stateq))])
    
    print("live loads energie ", + Live_energie*statsolver['x'][-1:])
    
    print("Internal loads energy",+ np.dot(np.dot(B.T,lagrange_stateq).T,statsolver['x'][0:Astat.size[1]-1]))
    print("Prm",+ np.dot(lagrange_statineq.T,coh))
    return data
