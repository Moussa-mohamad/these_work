def problem_building_up(data,case,arch,pont,case_pont):
    
    import numpy as np
    import time 
    import matplotlib
    
    import matplotlib.pyplot as plt
    #plt.ion()
    from shapely.geometry import Point
    import shapely.geometry.polygon 
    import math
    from cvxopt import  matrix
    from sympy import linear_eq_to_matrix, symbols
    Live_loadsvect = np.empty((0,1))
    Dead_loadsvect = np.empty((0,1))
    count_edges = 0;
    
    strength = data["strength"]
    dilantancy  = data["dilantancy"]         # dilatancy angle
    max_nl = data["max_normalload"]           # Maximum allowable normal load 
    min_nl = data["min_normalload"]           # minimum allowable normal load 
    blocks = data["blocks"]
    edges = data["edges"]
    nodes = data["nodes"]
    Live_loads = data["Live_loads"]
    Dead_loads = data["Dead_loads"]
    support = data["support"]
    
    start_problem =  time.time()
    allequations =[];                # define a variable for system of equilibrium equations 
    allvar = [];                     #define a variable to collect all interfaces variables
    blockscenters = np.array([])
    activeEdgesInd = np.empty((0,1)) # define a variable to collect all active edges indices 
    i=0                              #block indexing variable
    testvect = np.empty((0,1))       # to check whether an edge variable is equal to zero
    strength = [[float(strength[element1][element2]) for element2 in range(0,strength.shape[1]) ]for element1 in range(0,strength.shape[0])]
    strength = np.array(strength) 
    C = np.array([[]]);                 #Matrix to define the rupture criteria matrix
    N_mat = np.array([[]]);             #Matrix to define the displacement directions 
    coh = np.empty((0,1))               #cohesion vector
    allpolygons = dict()                # collect all blocks as polygons
    edgesIndCoor = dict()               # collect all edges  (with repetition)
    blockedgesInd = dict()             #collect block index for each edge in the variable edgesIndCoor
    
    fig, ax = plt.subplots(figsize=(12, 6))
    # set the backend to inline
    edgeslength = np.array([])
    all_edges_length = np.array([])
    
    while i < len(blocks):              #enumerate all blocks
        
        
        
        Live_loadsvect = np.append(Live_loadsvect, np.array([[0],[0],[0]]),axis=0) #add lines for the block j live loads
        Dead_loadsvect = np.append(Dead_loadsvect, np.array([[0],[0],[0]]),axis=0) #add lines for the block j dead loads
        j = int(blocks[i][0] )          #current block index

        c=np.where( blocks[:,0] ==str(j ))    # collect the coordinates of vertices belonging to the block 

        x = blocks[c[0]]
        y = [[   int(edges[int(x[element][1])-1][1]), int(edges[int(x[element][1])-1][2])  ]  for element in range(len(x)) ]
        z1 = [[   float(nodes[int(y[element][0])-1][1]), float(nodes[int(y[element][0])-1][2])  ]  for element in range(len(y)) ]
        z2 = [[   float(nodes[int(y[element][1])-1][1]), float(nodes[int(y[element][1])-1][2])  ]  for element in range(len(y)) ]
        blockjnodes = z1+z2
        blockjnodes = np.unique(blockjnodes,axis=0)
        blockjnodes = list(zip(*blockjnodes))

        blockjnodesX = np.array(blockjnodes[0])
        blockjnodesY = np.array(blockjnodes[1])

        A = [blockjnodesX,blockjnodesY] 
        A = A - np.mean(A,1)[:, None]
        A = np.argsort(np.arctan2(A[1,:],A[0,:]))
        blockjnodesXsort = blockjnodesX[A]           # x coordinates of all vertices of block j  
        blockjnodesYsort = blockjnodesY[A]           # y coordinates of all vertices of block j
        x = np.column_stack((blockjnodesXsort,blockjnodesYsort))
        polygon = shapely.geometry.polygon.Polygon(x)       #convert x to an area
        allpolygons[str(j)] = x
        area = polygon.area                       # calculate the area of the block j
        blockjcenter = np.array([polygon.centroid.xy[0][0],polygon.centroid.xy[1][0]]) # center of gravity of the block j
        blockscenters = np.concatenate((blockscenters,np.array([polygon.centroid.xy[0][0],polygon.centroid.xy[1][0]])))
        
        (x, y) = (blockjcenter[0],blockjcenter[1])
        
        #ax.text(x, y, str(j), ha='center', va='center', fontsize=5, fontweight='normal')
        
        #Add blocks dead and live loads

        #check volumetric live loads

        blockjLoadcheck = np.where(Live_loads['Blocks']['volume'][:,0] ==str(j))[0]
        blockjliveX = 0
        blockjliveY = 0
        blockjmoment = 0

        # Calculate live loads vector for block j
        #step 1: volumetric live loads

        if not len(blockjLoadcheck) ==0:
            for element in blockjLoadcheck:

                blockjliveX = blockjliveX + float(Live_loads['Blocks']['volume'][int(element)][1])*area
                blockjliveY = blockjliveY + float(Live_loads['Blocks']['volume'][int(element)][2])*area

        # step 2: concentrated live loads

        blockjLoadcheck = np.where(Live_loads['Blocks']['concentrated'][:,0] ==str(j))[0]
        if not len(blockjLoadcheck) ==0:
            for element in blockjLoadcheck:

                refpointX = float(nodes[int(Live_loads['Blocks']['concentrated'][element][4])-1][1])
                refpointY = float(nodes[int(Live_loads['Blocks']['concentrated'][element][4])-1][2])

                pointloadX = refpointX + float(Live_loads['Blocks']['concentrated'][element][5])
                pointloadY =  refpointY + float(Live_loads['Blocks']['concentrated'][element][6]) #are the coor sorted !
                loadX = float(Live_loads['Blocks']['concentrated'][element][1])
                loadY = float(Live_loads['Blocks']['concentrated'][element][2])
                moment = float(Live_loads['Blocks']['concentrated'][element][3])
                blockjliveX = blockjliveX + loadX  
                blockjliveY = blockjliveY + loadY
                blockjmoment = blockjmoment + moment + \
                               np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]

        Live_loadsvect[-3:] = np.array([[blockjliveX],[blockjliveY],[blockjmoment]]) # add block j live loads to the live loads vector

        # Calculate dead loads vector for block j

        blockjLoadcheck = np.where(Dead_loads['Blocks']['volume'][:,0] ==str(j))[0]
        blockjdeadX = 0
        blockjdeadY = 0
        blockjmoment = 0

        #step 1: volumetric dead loads

        if not len(blockjLoadcheck) ==0:
            for element in blockjLoadcheck:
                
                
                if pont:
                    blockjdeadX = blockjdeadX + float(Dead_loads['Blocks']['volume'][int(element)][1])*area
                    blockjdeadY = blockjdeadY + float(Dead_loads['Blocks']['volume'][int(element)][2])*area*2.6/0.142857142
                else:
                    blockjdeadX = blockjdeadX + float(Dead_loads['Blocks']['volume'][int(element)][1])*area
                    blockjdeadY = blockjdeadY + float(Dead_loads['Blocks']['volume'][int(element)][2])*area
                
                
                    
                
        # step 2: concentrated dead loads

        blockjLoadcheck = np.where(Dead_loads['Blocks']['concentrated'][:,0] ==str(j))[0]
        if not len(blockjLoadcheck) ==0:
            for element in blockjLoadcheck:

                refpointX = float(nodes[int(Dead_loads['Blocks']['concentrated'][element][4])-1][1])
                refpointY = float(nodes[int(Dead_loads['Blocks']['concentrated'][element][4])-1][2])

                pointloadX = refpointX + float(Dead_loads['Blocks']['concentrated'][element][5])
                pointloadY =  refpointY + float(Dead_loads['Blocks']['concentrated'][element][6]) #are the coor sorted !


                loadX = float(Dead_loads['Blocks']['concentrated'][element][1])
                loadY = float(Dead_loads['Blocks']['concentrated'][element][2])
                moment = float(Dead_loads['Blocks']['concentrated'][element][3])

                blockjdeadX = blockjdeadX + loadX  
                blockjdeadY = blockjdeadY + loadY
                blockjmoment = blockjmoment + moment + \
                               np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]

        Dead_loadsvect[-3:] = np.array([[blockjdeadX],[blockjdeadY],[blockjmoment]]) # add block j dead loads to the live loads vector
        blockjvector = []; # edges vectors of the block j
        
        while int(blocks[i][0]) == j:          # extracte all edges of the block j

            vertex1 = nodes[int(edges[int(blocks[i][1])-1][1])-1] #extract the data of the first vertex of an edge
            vertex2 = nodes[int(edges[int(blocks[i][1])-1][2])-1] #extract the data of the second vertex of an edge
            vertex1coor = np.array([float(vertex1[element]) for element in range(1,len(vertex1))  ])
            vertex2coor = np.array([float(vertex2[element]) for element in range(1,len(vertex2)) ])
            edgeindex = edges[int(blocks[i][1])-1][0] # current edge index
            edgesIndCoor[str(i+1)] = np.array([[edgeindex,vertex1coor[0],vertex2coor[0],vertex1coor[1],vertex2coor[1]]]) #find edges indices with repetition
            blockedgesInd[str(i+1)] = j #find blocks the current edge belongs to
            plt.plot([vertex1coor[0],vertex2coor[0]], [vertex1coor[1],vertex2coor[1]],'k' ) #mohamad
            
            edgelength = math.sqrt((vertex1coor[0]-vertex2coor[0])**2 + (vertex1coor[1]-vertex2coor[1])**2 )          

            #Add edge dead and live loads

            edgeLoadcheck = np.where( Live_loads['Edges']['continuous'][:,0] == float(edgeindex))[0]
            edgeliveX = 0
            edgeliveY = 0
            edgemoment = 0

            #step 1: continuous Live loads

            if not len(edgeLoadcheck) ==0:

                for element in edgeLoadcheck:
                    refpointX = float(nodes[int(Live_loads['Edges']['continuous'][element][3])-1][1])
                    refpointY = float(nodes[int(Live_loads['Edges']['continuous'][element][3])-1][2])

                    startpointloadX = refpointX + float(Live_loads['Edges']['continuous'][element][4])
                    startpointloadY =  refpointY + float(Live_loads['Edges']['continuous'][element][5]) #are the coor sorted !

                    endpointloadX = refpointX + float(Live_loads['Edges']['continuous'][element][6])
                    endpointloadY =  refpointY + float(Live_loads['Edges']['continuous'][element][7])

                    pointloadX = (startpointloadX +  endpointloadX)/2
                    pointloadY = (startpointloadY +  endpointloadY)/2

                    Loadlength = math.sqrt((startpointloadX-endpointloadX)**2 + (startpointloadY-endpointloadY)**2 ) 

                    loadX = float(Live_loads['Edges']['continuous'][element][1])*Loadlength
                    loadY = float(Live_loads['Edges']['continuous'][element][2])*Loadlength

                    edgeliveX = edgeliveX + loadX  
                    edgeliveY = edgeliveY  + loadY
                    edgemoment = edgemoment + \
                                   np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]


            #step 2: concentrated live loads

            edgeLoadcheck = np.where( Live_loads['Edges']['concentrated'][:,0] == float(edgeindex))[0]
            if not len(edgeLoadcheck) ==0:
                
                for element in edgeLoadcheck:
                    refpointX = float(nodes[int(Live_loads['Edges']['concentrated'][element][4])-1][1])
                    refpointY = float(nodes[int(Live_loads['Edges']['concentrated'][element][4])-1][2])

                    pointloadX = refpointX + float(Live_loads['Edges']['concentrated'][element][5])
                    pointloadY =  refpointY + float(Live_loads['Edges']['concentrated'][element][6])

                    loadX = float(Live_loads['Edges']['concentrated'][element][1])
                    loadY = float(Live_loads['Edges']['concentrated'][element][2])
                    moment = float(Live_loads['Edges']['concentrated'][element][3])

                    edgeliveX = edgeliveX + loadX  
                    edgeliveY = edgeliveY  + loadY
                    edgemoment = edgemoment + moment + \
                                   np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]
            
            Live_loadsvect[-3:] = Live_loadsvect[-3:] + \
            np.array([[edgeliveX],[edgeliveY],[edgemoment]]) # add edges live loads to the dead loads vector  



            #step 3: continuous dead loads

            edgeLoadcheck = np.where( Dead_loads['Edges']['continuous'][:,0] == float(edgeindex))[0]
            edgedeadX = 0
            edgedeadY = 0
            edgemoment = 0

            if not len(edgeLoadcheck) ==0:
                for element in edgeLoadcheck:
                    refpointX = float(nodes[int(Dead_loads['Edges']['continuous'][element][3])-1][1])
                    refpointY = float(nodes[int(Dead_loads['Edges']['continuous'][element][3])-1][2])

                    startpointloadX = refpointX + float(Dead_loads['Edges']['continuous'][element][4])
                    startpointloadY =  refpointY + float(Dead_loads['Edges']['continuous'][element][5]) #are the coor sorted !

                    endpointloadX = refpointX + float(Dead_loads['Edges']['continuous'][element][6])
                    endpointloadY =  refpointY + float(Dead_loads['Edges']['continuous'][element][7])

                    pointloadX = (startpointloadX +  endpointloadX)/2
                    pointloadY = (startpointloadY +  endpointloadY)/2

                    Loadlength = math.sqrt((startpointloadX-endpointloadX)**2 + (startpointloadY-endpointloadY)**2 ) 

                    loadX = float(Dead_loads['Edges']['continuous'][element][1])*Loadlength
                    loadY = float(Dead_loads['Edges']['continuous'][element][2])*Loadlength

                    edgedeadX = edgedeadX + loadX  
                    edgedeadY = edgedeadY  + loadY
                    edgemoment = edgemoment + \
                                   np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]


            #step 4: concentrated dead loads

            edgeLoadcheck = np.where( Dead_loads['Edges']['concentrated'][:,0] == float(edgeindex))[0]
          
            if not len(edgeLoadcheck) ==0:
                count_edges = count_edges+1;
                for element in edgeLoadcheck:
                    refpointX = float(nodes[int(Dead_loads['Edges']['concentrated'][element][4])-1][1])
                    refpointY = float(nodes[int(Dead_loads['Edges']['concentrated'][element][4])-1][2])

                    pointloadX = refpointX + float(Dead_loads['Edges']['concentrated'][element][5])
                    pointloadY =  refpointY + float(Dead_loads['Edges']['concentrated'][element][6])

                    loadX = float(Dead_loads['Edges']['concentrated'][element][1])
                    loadY = float(Dead_loads['Edges']['concentrated'][element][2])
                    moment = float(Dead_loads['Edges']['concentrated'][element][3])

                    edgedeadX = edgedeadX + loadX  
                    edgedeadY = edgedeadY  + loadY
                    edgemoment = edgemoment + moment + \
                                   np.cross(-np.array(blockjcenter)+np.array([[pointloadX,pointloadY]]),np.array([[loadX,loadY]]))[0]

            Dead_loadsvect[-3:] = Dead_loadsvect[-3:] + \
            np.array([[edgedeadX],[edgedeadY],[edgemoment]]) # add edges dead loads to the dead loads vector    


            if not int(edges[int(blocks[i][1])-1][3]) == 0:

                edgevector = vertex1coor - vertex2coor  # find the direction of the edge
                edgevector = edgevector / np.linalg.norm(edgevector) # find a unit vector 

                normalvector = np.array([-edgevector[1], edgevector[0]]) #find the normal to the edge
                         #check whether yhe normal is inward or outward ?
                edgecenter = ((vertex1coor[0]+vertex2coor[0]) /2, (vertex1coor[1]+vertex2coor[1]) /2) # find the center of the edge
                testpoint = edgecenter + 0.001*normalvector # find a new point very close to the edge center following the direction 'normal vector'  
                #A = [blockjnodesX,blockjnodesY] 
                #A = A - np.mean(A,1)[:, None]
                #A = np.argsort(np.arctan2(A[1,:],A[0,:]))
                #blockjnodesX = blockjnodesX[A]    # x coordinates of all vertices of block j  
                #blockjnodesY = blockjnodesY[A]     # y coordinates of all vertices of block j
                #x = np.column_stack((blockjnodesX,blockjnodesY))
                #polygon = Polygon(x)  #convert x to an area
                #blockjcenter = np.array([polygon.centroid.xy[0][0],polygon.centroid.xy[1][0]]) # center of gravity of the block j
                #print( Point(testpoint[0],testpoint[1]).within(polygon))
                if  not Point(testpoint[0],testpoint[1]).within(polygon):
                    normalvector = -normalvector
                    edgevector = -edgevector 
                #print(normalvector,edgevector)
                blockjvector.append(edgevector) 
                blockjvector.append(normalvector) 
                blockjvector.append(vertex1) 
                blockjvector.append(vertex2) 
                blockjvector.append(int(edges[int(blocks[i][1])-1][3])) # 0 1 or 2
                blockjvector.append(int(edges[int(blocks[i][1])-1][0])) # edge index
                #print(blockjvector[0])
            i+=1
            if i == len(blocks):
                break
        blockjvar =[]; # all edge variables (t,n &m) for the block j
        blockjeq1 = 0 # block equilibrium equation in x direction
        blockjeq2 = 0 # block equilibrium equation in y direction
        blockjeq3 = 0 # block equilibrium equation (moment)  
        
        
        for v in range(0,len(blockjvector),6):
            testx = 1
            testy = 1
            testm = 1
            currentEdgeInd = blockjvector[v+5] 
            minInd = min(blockjvector[v+2][0],blockjvector[v+3][0]) #find the min index of the vertex of an edge
            maxInd = max(blockjvector[v+2][0],blockjvector[v+3][0]) #find the max index of the vertex of an edge 
            nsym= symbols('n_'+minInd+'_'+maxInd) # normal reaction at an edge, symbolic representation 
            tsym= symbols('t'+minInd+'_'+maxInd) # shear at an edge, symbolic representation 
            msym = symbols('m'+minInd+'_'+maxInd) # moment at an edge, symbolic representation 
            if blockjvector[v+4] == 2:

                testx = (support[(np.where(support[:,0] == str( blockjvector[v+5])))[0][0]][1] == '0')*1 #test if the support prevents disp in x direction
                testy = (support[(np.where(support[:,0] == str( blockjvector[v+5])))[0][0]][2] == '0')*1 #test if the support prevents disp in y direction
                testm = (support[(np.where(support[:,0] == str( blockjvector[v+5])))[0][0]][3] == '0')*1  #test if the support prevents rotation around z


            n= blockjvector[v+1]*nsym # normal vector, symbolic representation
            t= blockjvector[v]*tsym # shear vector, symbolic representation
            
            edgenode1 = [float(blockjvector[v+2][1]),float(blockjvector[v+2][2]) ] #first node coordinates
            edgenode2 = [float(blockjvector[v+3][1]),float(blockjvector[v+3][2]) ] # second node coordinates
            edgecenter = [np.mean([float(blockjvector[v+2][1]),float(blockjvector[v+3][1])]), np.mean([float(blockjvector[v+2][2]),float(blockjvector[v+3][2])])]
            edgelength = np.linalg.norm(np.array(edgenode1)-np.array(edgenode2))
            
            
            
            if len(np.where(blockjvar == nsym)[0]) == 0:
                 #Add the equations corresponding to a new edge to the problem of the block j 
                blockjeq1 = blockjeq1 + (n[0]+t[0])
                blockjeq2 = blockjeq2 + (n[1]+t[1])
                blockjeq3 = blockjeq3 + np.cross(-np.array(blockjcenter)+np.array(edgecenter),np.array([[n[0] +t[0] ,0]])) + \
                np.cross(-np.array(blockjcenter)+np.array(edgecenter),np.array([[0,n[1] +t[1] ]]))

                if blockjvector[v+1][1] ==0:
                    if blockjvector[v+1][0] > 0:
                        blockjeq3 = blockjeq3 - msym*edgelength/2
                    else:
                        blockjeq3 = blockjeq3 + msym*edgelength/2
                else:
                    if blockjvector[v+1][1] > 0:
                        blockjeq3 = blockjeq3 + msym*edgelength/2
                    else:
                        blockjeq3 = blockjeq3 - msym*edgelength/2

                #Add the new edge variables to the previous variables        

                if len(np.where(allvar == nsym)[0]) == 0:
                    all_edges_length = np.append(all_edges_length,edgelength)
                    activeEdgesInd = np.append(activeEdgesInd,np.array([currentEdgeInd]))
                    testvect = np.append(testvect,np.array([testx,testy,testm]))
                    blockjvar = np.append(blockjvar,np.array([tsym,nsym,msym ]))
                    allvar = np.append(allvar,np.array([tsym,nsym,msym ]))
                    # Create the rupture criteria matrix C
                    ind = np.where( strength[:,0] == currentEdgeInd)[0] #find the index of an edge in the strength array
                    if not len(ind) ==0:
                     
                        #if int(data["edges"][int(currentEdgeInd)-1][3]) == 2 and float(data["nodes"][int(data["edges"][int(currentEdgeInd)-1][1])-1][1]) >35 and float(data["nodes"][int(data["edges"][int(currentEdgeInd)-1][1])-1][1]) > 39:
                         #   coh_add = np.array([999999999,9999999999,99999999,999999999,99999999,9999999])
                        #else:
                         #   coh_add = np.array([strength[ind[0]][1]*edgelength,strength[ind[0]][1]*edgelength,0,0,strength[ind[0]][4]*edgelength,strength[ind[0]][3]*edgelength])
                            

                        
                        C_add =np.array([[1,-math.tan(math.radians(strength[ind[0]][2])),0],[-1,-math.tan(math.radians(strength[ind[0]][2])),0],
                           [0,-1,1],[0,-1,-1],[0,-1,0],[0,1,0]]) #find new block of the matrix C
                        C = np.block([[C,np.zeros((C.shape[0],3))],[np.zeros((6,C.shape[1])),C_add]]) #Add a new block to the matrix C 

                        C_add  =np.array([[1,-math.tan(math.radians(dilantancy)),0],[-1,-math.tan(math.radians(dilantancy)),0],
                           [0,-1,1],[0,-1,-1],[0,-1,0],[0,1,0]]) #find new block of the matrix C_non ass

                        N_mat =   np.block([[N_mat,np.zeros((N_mat.shape[0],3))],[np.zeros((6,N_mat.shape[1])),C_add]]) #Add a new block to the matrix C non ass                 

                        coh_add = np.array([strength[ind[0]][1]*edgelength,strength[ind[0]][1]*edgelength,0,0,strength[ind[0]][4]*edgelength,strength[ind[0]][3]*edgelength])

                        coh = np.append (coh,coh_add)
                        edgeslength = np.append(edgeslength,edgelength)

        #Add the problems of the block j to the global problem
        allequations = np.append(allequations,np.array(blockjeq1))
        allequations = np.append(allequations,np.array(blockjeq2))
        allequations = np.append(allequations,np.array(blockjeq3))
    #Write the linear equations problem in matrix form
    data["count_edges"] = count_edges
    
    #plt.xlim(0,5)
    nodes = data["nodes"]
    xlimmin = float(nodes[0,1])
    xlimmax=  float(nodes[0,1])
    ylimmin=  float(nodes[0,2])
    ylimmax = float(nodes[0,2])   
    for element in range(nodes.shape[0]):
        xlimmin = min([xlimmin,float(nodes[element,1])])
        xlimmax = max([xlimmax,float(nodes[element,1])])
        ylimmin = min([ylimmin,float(nodes[element,2])])
        ylimmax = max([ylimmax,float(nodes[element,2])])

    plt.xlim(xlimmin, xlimmax)
    plt.ylim(ylimmin, ylimmax)
    
    ax.set_xlim(xlimmin, xlimmax)
    ax.set_ylim(ylimmin, ylimmax)
    ax.set_axis_off()
    
    plt.savefig('structure_shape.jpg', dpi=300, format='jpg')
    plt.show(block = False)
    
    
    
    #display(plt.gcf())
    B,b = linear_eq_to_matrix(allequations,allvar)
    B = np.array(B)

    #Astat=B
    
    C = np.delete(C,0,axis=0)
    N_mat = np.delete(N_mat,0,axis=0)
    Cstat = np.append(C,np.zeros((C.shape[0],1)),axis=1)

    B = B.astype(float)
    testrankmatrix = np.dot(B,B.T)
    
    #print("B.TxB", testrankmatrix)
    #testrankmatrix = B.T
    #print(np.linalg.eigvals(testrankmatrix))
    #print("Matrix (BB^t)^-1: dimension ", np.linalg.inv(testrankmatrix)) 
    print("Matrix B.BT: dimension " + str(testrankmatrix.shape[0])+"x" + str(testrankmatrix.shape[1]) + \
          " Rank " + str(np.linalg.matrix_rank(testrankmatrix)))



    end_problem = time.time()
    print("Problem building-up",end_problem-start_problem , "seconds")
    
    #Dead_loadsvect = Dead_loadsvect/5*5
    if arch:
        
        for element in range(0,len(Dead_loadsvect),3):    
            Live_loadsvect[element] =  Live_loadsvect[element]*20/0.142857142
            Dead_loadsvect[element+1] = Dead_loadsvect[element+1]*20/0.142857142
    

    
    
    
    
    if arch:
        if case ==3:
            Live_loadsvect = Live_loadsvect*0
            Live_loadsvect[data["layer_loading"]*3 + 1] = -1
            Live_loadsvect[data["layer_loading"]*3 ]=0.0001*0
    
        
    
    
    
    
    if arch:
        if case == 2:
            data["dead_loadvect"] = Dead_loadsvect*0
            data["live_loadvect"] = Dead_loadsvect
    
    
    
    
    
    if case_pont == 2:
        data["dead_loadvect"] = Dead_loadsvect*0
        data["live_loadvect"] = Dead_loadsvect
        
   
    
    
    Astat = B
    Cstat = np.append(C,np.zeros((C.shape[0],1)),axis=1)
    bstat = Dead_loadsvect
    
    for element in np.where(testvect == 0)[0]: # add boundary conditions to A
        new_line = np.zeros(Astat.shape[1]) 
        new_line[element] = 1        # create a line corresponding to a boundary condition
        Astat = np.append(Astat,[new_line],axis=0) #add the created line to A
        Live_loadsvect = np.append(Live_loadsvect,np.array([[0]]),axis=0)
        bstat = np.append(bstat,np.array([[0]]),axis=0)
        #cohkin = np.append(cohkin,np.array([[0]]))
    
    Astat = np.append(Astat,Live_loadsvect,axis=1)
    
    count = len(bstat)
    while not count == Astat.shape[0]:
        bstat = np.append(bstat,np.array([[0]]))
        count = count+1
   
 

    zstat = matrix(np.array([*np.zeros(Astat.shape[1] -1),-1.0])) 
    Astat = Astat.astype(np.double)
    bstat = Dead_loadsvect 
 
    Cstat = matrix(Cstat)
    cohstat = matrix(coh)
    Astat = matrix(Astat)
    bstat = matrix(bstat)
    
    data["Astat"] = Astat
    data["bstat"] = bstat
    data["Gstat"] = Cstat
    data["hstat"] = cohstat
    data["cstat"] = zstat 
         
     

    B = B.astype(np.double)
    cohkin = coh
    
    #Akin = Astat.T
    
    zkin = -matrix(np.append(Dead_loadsvect,-cohkin))
    Akin = np.append(B.T,C.T,axis=1)

    for element in np.where(testvect == 0)[0]: # add boundary conditions to A
        new_line = np.zeros((Akin.shape[0],1)) 
        new_line[element] = 1        # create a line corresponding to a boundary condition
        Akin = np.append(Akin,np.array(new_line),axis=1) #add the created line to A
        cohkin = np.append(cohkin,np.array([[0]]))
        
    Akin = np.array( [np.append(Live_loadsvect,np.zeros(len(zkin)-len(Live_loadsvect))), *Akin ])
    bkin = np.array([1,*np.zeros(Akin.shape[0]-1)]) # put the 1 back for pext
    Gkin = np.block([np.zeros((len(zkin)-len(Dead_loadsvect),len(Dead_loadsvect))),-np.eye((len(zkin)-len(Dead_loadsvect)))])
    hkin = np.zeros(len(zkin)-len(Dead_loadsvect))
    
    Gkin = matrix(Gkin)
    hkin = matrix(hkin)
    Akin = matrix(Akin)
    bkin=matrix(bkin)
     
     
     
    data["Akin"] = Akin
    data["bkin"] = bkin
    data["Gkin"] = Gkin
    data["hkin"] = hkin
    data["zkin"] = zkin
    
    
    
    
    
    data["testvect"] = testvect
    data["allpolygons"] = allpolygons
    data["blockscenters"] = blockscenters
    data["edgesIndCoor"] = edgesIndCoor 
    data["activeEdgesInd"] = activeEdgesInd
    data["blockedgesInd"] = blockedgesInd
    data["edgeslength"] = edgeslength
    
    
    
    
    data["criteria_matrix"] = C
    data["disp_direction"] =  N_mat
    data["equilibrium_matrix"] =  B
    data["cohesion_vector"] =  coh
    data["dead_loadvect"] = Dead_loadsvect
    data["edges_length"] = all_edges_length 
    data["live_loadvect"] = Live_loadsvect
    return data