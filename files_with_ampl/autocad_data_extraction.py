def autocad_data_extraction(data,pont,culee_bloc,arch):
    
    import pyautocad
    from rtree import index
    import ezdxf
    import time

    import pyautocad
    from shapely.geometry import Polygon
    import ezdxf
    import time
    import numpy as np
    #import pygame
    #from shapely.geometry import Point
    import shapely.geometry.polygon 
 
    
    
    friction_ang = data["friction_angle"]   #friction angle
    cohesion = data["cohesion"]               #cohesion
    dilantancy  = data["dilantancy"]         # dilatancy angle
    max_nl = data["max_normalload"]           # Maximum allowable normal load 
    min_nl = data["min_normalload"]           # minimum allowable normal load 
    fname = data["notepad_files"]             # text file names (the file in the first position is no longer used in this code so no need to change it while the second one contain the loading data)
    dxf_file = data["autocad_file"]           #DXF file location
    cohesion_reg = data["cohesion_region"] 
    friction_reg = data["friction_ang_region"]  
    
    max_nl_1 = data["max_nl_1"]
    max_nl_2 = data["max_nl_2"]
    
    blocks = []                            #define an empty vector for block references
    start_cad = time.time()                #Begin timing the execution of this function
    doc = ezdxf.readfile(dxf_file)         # import the dxf file
    msp = doc.modelspace()                 #read data in the file

    block_references = msp.query("INSERT") #import block references from the dxf file
    yy= block_references
    support_edges = np.empty((0,2))
    edges = []                             # Initialize empty list to store edges
    ind =0                                 #initialize block index 
    allnodes = set()                       #define a set to collect all nodes from the drawing
    acad = pyautocad.Autocad()
    
    print(len(block_references)) 
    # Extract all the block references data
    for block_reference in block_references:
        
        points = block_reference.explode() #explode blocks references into polylines
        block_polygon = []                 #define a variable to temporarily store current block references in form of polygon 
        
        print(points[0])
        
        vertices = points[0].get_points()  # Get polygon vertices
        for element in range(0,len(vertices)): #enumerate all block references vertices
            node_coordinates = tuple([round(vertices[element][0],3),round(vertices[element][1],3)]) #extract nodes coordinates of the current vertex
            allnodes.add(node_coordinates)    #add the coordinate of the current point to allnodes set
            block_polygon.append( [round(vertices[element][0],3),round(vertices[element][1],3) ]) #store the current vertex coordinates in block poygon 

        block_polygon  = shapely.geometry.polygon.Polygon(block_polygon)   # convert the current block to a "python" polygon


        blocks.append({ 'id': ind,'polygon': block_polygon}) #add current block's data to the blocks dictionary which collects all blocks data
        ind = ind+1 # increase blocks indexing variable by 1

    # Extract all polylines data
    # Open DXF file
    doc = ezdxf.readfile(dxf_file)  # import the dxf file
    msp = doc.modelspace() #read data in the file
    polygons = msp.query('LWPOLYLINE') # Extract all polygon entities
    polygons.extend(msp.query('3DFACE'))

    for poly in polygons: # Iterate over polygons
        
        layer_name = poly.dxf.layer
    
        if layer_name == 'Layer1':
            data["layer_loading"] = ind
        
        block_polygon = []           #define an empty vector for block polygons
        vertices = poly.get_points() # Get polygon vertices
        
        for element in range(0,len(vertices)): #enumerate all polylines vertices

            node_coordinates = tuple([round(vertices[element][0],3),round(vertices[element][1],3)]) #extract nodes coordinates of the current vertex
            allnodes.add(node_coordinates)


            block_polygon.append( [round(vertices[element][0],3),round(vertices[element][1],3) ]) 

        block_polygon  = shapely.geometry.polygon.Polygon(block_polygon) # convert the current block to a "python" polygon  
        #add current block's data to the blocks dictionary which collects all blocks data
        blocks.append({'id': ind,'polygon': block_polygon})
        ind = ind+1 # increase blocks indexing variable by 1
    
    
    
    #reorganize nodes data

    ind = 1 # nodes indexing variable
    nodes = np.empty((0,3)) #create an array to reorganize nodes data previously extracted
    for element in allnodes: #enumerate all nodes 
        nodes= np.append(nodes,np.array([[ind,element[0],element[1]]]),axis=0) # organize data for the current node (index,xcoor,ycoor)
        ind += 1   # increase node indexing variable by 1  


    #Extract attributes

    doc = ezdxf.readfile(dxf_file)               #import the dxf file
    msp = doc.modelspace()                       #read data in the file
    block_references = msp.query("INSERT")       #extract references
    attrib_points =[] 
    for block_reference in block_references:
        entity = block_reference                 # process all entities
        attrib_points = np.empty((0,2))          # define attribute points vector
        
        for attrib in entity.attribs:
            
            attrib_point = ( attrib.dxf.insert[0],
                             attrib.dxf.insert[1]) # calculate position of an attribute
            attrib_points = np.append(attrib_points,np.array([[round(attrib_point[0],3),round(attrib_point[1],3)]]),axis=0)
        
        #B = attrib_points - np.mean(attrib_points,1)[:, None]
        #B = np.argsort(np.arctan2(B[1,:],B[0,:]))
        #attrib_points = attrib_points.T
        #attrib_points = attrib_points[B].T
        
        for element in range(0,len(attrib_points)-1):
            condition1= (nodes[:,1] == attrib_points[element,0]) & (nodes[:,2] == attrib_points[element,1]) #find the position of the first point

            condition2 = (nodes[:,1] == attrib_points[element+1,0]) & (nodes[:,2] == attrib_points[element+1,1]) #find the position of the first point

            if (len(nodes[condition1,:]) !=0) & (len(nodes[condition2,:]) !=0): #if no intersection between the current edge and is found 

                node1_ind = nodes[condition1,:][0][0]   # index of the first node of the attribute
                node2_ind = nodes[condition2,:][0][0] # index of the second node of the attribute

                support_edges = np.append(support_edges,np.array([[min(node1_ind,node2_ind),max(node1_ind,node2_ind)]]),axis=0) #add indices to the global matrix
        attrib_points =[]


    # Create an R-Tree index to find neighboring blocks
    p = index.Property()
    ind = index.Index(properties=p)            # extract R-Tree tools
    for i, block in enumerate(blocks):
        ind.insert(i, block['polygon'].bounds) #add blocks boundaries 

    # Iterate over the blocks and check for connections
    edges = np.empty((0,4)) #define edges array
    blocks_data = np.empty((0,2)) #define new block
    edgeind = 1
    unique_nodes = set() #find nodes without repetition


    #Get the final shape of the vectors "blocks", "nodes" and "edges"
    print(len(blocks))
    
    
    for i in range(len(blocks)):
        print(blocks[i])
        
        block = blocks[i]
        blockunique_nodes = set() #all nodes for the current block

        nearby_blocks = ind.intersection(block['polygon'].bounds) # Find all nearby blocks for the current block using the R-Tree index

        for nearby_block in nearby_blocks:

            # Check if the nearby block is actually connected to the current block
            if block['polygon'].intersects(blocks[nearby_block]['polygon']) or block['polygon'].touches(blocks[nearby_block]['polygon']):

                intersection = block['polygon'].intersection(blocks[nearby_block]['polygon'])

                if isinstance(intersection, shapely.geometry.polygon.Polygon): #if the intersection is a polygon

                    xx, yy = intersection.exterior.coords.xy
                    for i in range(0,len(xx)):
                        node_coordinates = tuple([xx[i],yy[i]])
                        unique_nodes.add(node_coordinates)
                        blockunique_nodes.add(node_coordinates)

                else: #if the intersection is a line

                    for i in range(0,len(intersection.xy[0])):
                        node_coordinates = tuple([intersection.xy[0][i],intersection.xy[1][i]]) 
                        unique_nodes.add(node_coordinates) # add the current node to the nodes list
                        blockunique_nodes.add(node_coordinates) #add the current node to all current block nodes


        block_nodes = np.empty((0,2))  #variable to reorganize the block's data      
        for element in blockunique_nodes:
            block_nodes= np.append(block_nodes,np.array([[element[0],element[1]]]),axis=0) # collect all the current block nodes in a specific array form

        block_nodes = block_nodes.T
        blockjnodesX = np.array(block_nodes[0])
        blockjnodesY = np.array(block_nodes[1])

        A = [blockjnodesX,blockjnodesY] 
        A = A - np.mean(A,1)[:, None]
        A = np.argsort(np.arctan2(A[1,:],A[0,:]))
        blockjnodesX = blockjnodesX[A]      
        blockjnodesY = blockjnodesY[A]
        blockjnodesX = np.append(blockjnodesX,np.array([blockjnodesX[0]]))
        blockjnodesY = np.append(blockjnodesY,np.array([blockjnodesY[0]]))
        for element in range(0,len(blockjnodesX)-1):

            condition = (nodes[:, 1] == blockjnodesX[element]) & (nodes[:, 2] == blockjnodesY[element]) #look for the first vertex
            vertex1 = nodes[condition, :] #extract first vertex data
            condition = (nodes[:, 1] == blockjnodesX[element+1]) & (nodes[:, 2] == blockjnodesY[element+1]) #look for the first vertex
            vertex2 = nodes[condition, :] #extract first vertex data
            vertex1 = vertex1[0]
            print(vertex1)
            vertex2 = vertex2[0]
            min_ind = int(min(vertex1[0],vertex2[0]))
            max_ind = int(max(vertex1[0],vertex2[0]))
            condition = (edges[:,1] ==min_ind) & (edges[:,2] ==max_ind) # search the current edge in the egdes list
            if len(edges[condition, :]) ==0:                            #if the edge does not already exist in the list
                edges = np.append(edges,np.array([[int(edgeind),int(min_ind),int(max_ind),int(0)]]),axis =0)
                blocks_data = np.append(blocks_data,np.array([[block['id']+1,edgeind]]),axis =0)
                edgeind +=1
            else:                                                       #if the edge already exists in the list
                rep_edgeind = edges[condition, :][0][0] #find the repeated edge index in the list to add it to the block edges column

                edges[int(rep_edgeind)-1][3]=1

                blocks_data = np.append(blocks_data,np.array([[block['id']+1,rep_edgeind]]),axis =0)
   
    #support_edges[1][0]=3 add this for the special example
 
    # Give the value 2 to edges on supports 
    for element in support_edges:
        edge_elements  = edges[(edges[:,1] == element[0]) & (edges[:,2] == element[1]) ,:][0] # extract current edge elements
        edges[int(edge_elements[0])-1][3] = 2
        # Print the connections


    support_ind = np.empty((0,1)) 
    contact_ind = np.empty((0,1))
    for element in edges:
        if element[3] == 1:
            contact_ind = np.append(contact_ind,np.array([[str(int(element[0]))]]))
        if element[3] == 2:
            support_ind = np.append(support_ind,np.array([[str(int(element[0]))]]))

    blocks  =  blocks_data #create an  array for blocks

    strength = np.zeros((len(contact_ind) + len(support_ind),5)) #create an empty array for edges characteristics
    strength = strength.astype(str)
    strength[:,0] = [*contact_ind,*support_ind]
    strength[:,1] = (len(contact_ind) + len(support_ind))*[str(cohesion)]
    strength[:,2] = (len(contact_ind) + len(support_ind)) * [str(friction_ang)]
    strength[:,3] = (len(contact_ind) + len(support_ind)) * [str(max_nl)]
    strength[:,4] = (len(contact_ind) + len(support_ind)) * [str(min_nl)]
    
    
    
    if len(support_ind): 
        if arch:
            for element in support_ind:
                ind = np.where(strength[:,0] == element)[0]  
                if nodes[int(edges[int(element)-1][1])-1][1] < 19:

                    strength[int(ind),3] =  str(max_nl_1)
                else:
                    strength[int(ind),3] =  str(max_nl_2)
        if pont and culee_bloc:
            for element in support_ind:
                ind = np.where(strength[:,0] == element)[0]  
                if nodes[int(edges[int(element)-1][1])-1][1] > 51:

                    strength[int(ind),1] =  str(9999)
        
            

        
        support =  np.zeros((len(support_ind),4)) #create an empty array for the support data
        support = support.astype(str)
        support[:,0] = support_ind
        support[:,1] = (len(support_ind))*[str(0)]
        support[:,2] = (len(support_ind))*[str(0)]
        support[:,3] = (len(support_ind))*[str(0)]

        strength[-1:][0][1:3] = [cohesion_reg,friction_reg]

    edges_perm = np.empty((0,4))
    for element in range(0,len(edges)):
        edges_perm  = np.append(edges_perm ,np.array([[str(int(edges[element][0])),str(int(edges[element][1])),str(int(edges[element][2])),str(int(edges[element][3]))]]),axis=0)
    edges = edges_perm 

    blocks_perm = np.empty((0,2))
    for element in range(0,len(blocks)):
        blocks_perm  = np.append(blocks_perm ,np.array([[str(int(blocks[element][0])),str(int(blocks[element][1]))]]),axis=0)
    blocks = blocks_perm    

    nodes_perm = np.empty((0,3))
    for element in range(0,len(nodes)):
        nodes_perm  = np.append(nodes_perm ,np.array([[str(int(nodes[element][0])),str(nodes[element][1]),str(nodes[element][2])]]),axis=0)
    nodes = nodes_perm   

    end_cad = time.time()
    print("CAD: " ,end_cad-start_cad, "seconds")
    blocks = blocks_perm
    nodes = nodes_perm
    edges = edges_perm
    
    data["support_ind"] = support_ind
    data["strength"] = strength
    data["support"] = support
    data["blocks"] = blocks
    data["nodes"] = nodes
    data["edges"] = edges

    return data