def loading_extraction(data,remplissage):
    
    import pyautocad
    #from rtree import index
    import ezdxf
    import time

    import pyautocad
    from shapely.geometry import Polygon
    import ezdxf
    import time
    import numpy as np
    #import pygame
    #from shapely.geometry import Point
    #import shapely.geometry.polygon 
    
    
    fname = data["notepad_files"]
    dxf_file = data["autocad_file"]
    start_loading_ext = time.time()
    
    
    Dead_loads = dict()
    Dead_loads['Edges'] = dict()
    Dead_loads['Blocks'] = dict()
    Dead_loads['Edges']['concentrated']= np.empty((0,7))
    Dead_loads['Edges']['continuous']= np.empty((0,8))
    Dead_loads['Blocks']['concentrated']= np.empty((0,7))
    Dead_loads['Blocks']['volume']= np.empty((0,3))
    Live_loads = dict()
    Live_loads['Edges'] = dict()
    Live_loads['Blocks'] = dict()
    Live_loads['Edges']['concentrated']= np.empty((0,7))
    Live_loads['Edges']['continuous']= np.empty((0,8))
    Live_loads['Blocks']['concentrated']= np.empty((0,7))
    Live_loads['Blocks']['volume']= np.empty((0,3))
    
    

    # Read the mechanic file 
    try:
        file = open(fname[1],'r') #check if the file exists
    except: 
        print('File can not be opened:',fname) # if the file is not found
        quit()
    Textlines = file.readlines() 
    file.close() 
    #extract strength data

    #extract Loads data
    for line in range(4,len(Textlines)):
        if not Textlines[line].startswith('block'):

            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip() 
                line_split = line_content.split("\t")
                Dead_loads['Edges']['concentrated'] = np.append(Dead_loads['Edges']['concentrated'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6]]]),axis =0) 

        else:
            line_block_Dccload = line
            break

    for line in range(line_block_Dccload+1,len(Textlines)):
        if not Textlines[line].startswith('Continuous'):
            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip() 

                line_split = line_content.split("\t")

                Dead_loads['Blocks']['concentrated'] = np.append(Dead_loads['Blocks']['concentrated'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6]]]),axis =0)  
        else:
            line_Edge_Dcnload = line
            break
    for line in range(line_Edge_Dcnload+2,len(Textlines)):
        if not Textlines[line].startswith('Block'):
            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip() 
                line_split = line_content.split("\t")
                Dead_loads['Edges']['continuous'] = np.append(Dead_loads['Edges']['continuous'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6],line_split[7]]]),axis =0)  
        else:
            line_block_Dcnload = line
            break
    for line in range(line_block_Dcnload+1,len(Textlines)):
        if not Textlines[line].startswith('Live'):
            if  not Textlines[line]=='\n':
                line_content1 = Textlines[line].rstrip() 
                line_content = Textlines[line_block_Dcnload+1].rstrip()
                line_split = line_content.split("\t")
                line_split1 = line_content1.split("\t")

                Dead_loads['Blocks']['volume'] = np.append(Dead_loads['Blocks']['volume'],np.array([[line_split1[0],line_split[1],line_split[2]]]),axis =0)  
        else:
            line_edge_Lccload = line
            break

    for line in range(line_edge_Lccload+3,len(Textlines)):
        if not Textlines[line].startswith('block'):
            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip() 
                line_split = line_content.split("\t")


                Live_loads['Edges']['concentrated'] = np.append(Live_loads['Edges']['concentrated'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6]]]),axis =0)  
        else:
            line_block_Lccload = line
            break
    for line in range(line_block_Lccload+1,len(Textlines)):
        if not Textlines[line].startswith('Continuous'):
            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip()

                line_split = line_content.split("\t")


                Live_loads['Blocks']['concentrated'] = np.append(Live_loads['Blocks']['concentrated'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6]]]),axis =0)  
        else:
            line_Edge_Lcnload = line
            break

    for line in range(line_Edge_Lcnload+2,len(Textlines)):
        if not Textlines[line].startswith('Block'):
            if  not Textlines[line]=='\n':
                line_content = Textlines[line].rstrip() 
                line_split = line_content.split("\t")
                Live_loads['Edges']['continuous'] = np.append(Live_loads['Edges']['continuous'],np.array([[line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6],line_split[7]]]),axis =0)  
        else:
            line_block_Lcnload = line
            break
    for line in range(line_block_Lcnload+1,len(Textlines)):

        if  not Textlines[line]=='\n':
            line_content1 = Textlines[line].rstrip()
            line_content = Textlines[line_block_Lcnload+1].rstrip()

            line_split = line_content.split("\t")
            line_split1 = line_content1.split("\t")

            Live_loads['Blocks']['volume'] = np.append(Live_loads['Blocks']['volume'],np.array([[line_split1[0],line_split[1],line_split[2]]]),axis =0)  
    
    #Add fills and pavement weight
    
    doc = ezdxf.readfile(dxf_file)  # Import the DXF file
    msp = doc.modelspace()  # Read data in the file
    lines = msp.query('LINE')  # Extract all LINE entities

    # Create a list to store all vertices
    # Create a list to store all vertices
    all_vertices = []
    all_vertices_ind = []
    all_edges = []
    nodes = data["nodes"]
    nodes = np.array([[round(float(value), 3) for value in row] for row in nodes])

    for line in lines:  # Iterate over lines
        # Collect start and end points
        start_point = (round(line.dxf.start.x,3), round(line.dxf.start.y,3), round(line.dxf.start.z,3))

        end_point = (round(line.dxf.end.x,3), round(line.dxf.end.y,3), round(line.dxf.end.z,3))

        # Add start and end points to the list
        all_vertices.append(start_point)
        all_vertices.append(end_point)

    # Sort the vertices by x-coordinate
    all_vertices_sorted = sorted(all_vertices, key=lambda vertex: vertex[0])

    for vertex in range(0,len(all_vertices_sorted),2):

        condition = np.all(
        np.column_stack((
            np.abs(nodes[:,1] - all_vertices_sorted[vertex][0]) < 0.01,
            np.abs(nodes[:,2] - all_vertices_sorted[vertex][1]) < 0.01
        )),
        axis=1
        )

        # Use the condition to index nodes or perform other operations
        filtered_nodes_1 = nodes[condition]


        condition = np.all(
        np.column_stack((
            np.abs(nodes[:,1] - all_vertices_sorted[vertex+1][0]) < 0.01,
            np.abs(nodes[:,2] - all_vertices_sorted[vertex+1][1]) < 0.01
        )),
        axis=1
        )

        # Use the condition to index nodes or perform other operations
        filtered_nodes_2 = nodes[condition]

        all_vertices_ind.append( filtered_nodes_1[0][0])
        all_vertices_ind.append( filtered_nodes_2[0][0])

        edges = data["edges"].astype(float)


        condition = np.all(
        np.column_stack((
            edges[:,1] == min([filtered_nodes_1[0][0],filtered_nodes_2[0][0]]) ,
            edges[:,2] == max([filtered_nodes_1[0][0],filtered_nodes_2[0][0]])
        )),
        axis=1
        )

        edge = edges[condition][0][0] 
        all_edges.append(edge)
    
    
    def get_unit_anticlockwise_normal(direction):
    # Calculate the clockwise normal vector
        dx, dy = direction
        normal = (dy, -dx)

        # Calculate the length of the normal vector
        length = np.sqrt(normal[0]**2 + normal[1]**2)

        # Calculate the unit normal vector
        unit_normal = (normal[0] / length, normal[1] / length)

        return unit_normal
    

    for i in range(0,len(all_edges)): 
        
        if remplissage:
            fact = 1
        else:
            fact =0;
            
        direction = (all_vertices[2*i+1][0]-all_vertices[2*i][0], all_vertices[2*i+1][1]-all_vertices[2*i][1])  # direction de l'arret
        
        unit_normal = get_unit_anticlockwise_normal(direction)
        
        load_value = fact*2.6*abs((16.03-0.5*(all_vertices[2*i][1] + all_vertices[2*i+1][1]))*(all_vertices[2*i][0]- all_vertices[2*i+1][0]))
        
        Dead_loads['Edges']['concentrated'] = np.append(Dead_loads['Edges']['concentrated'], \
            np.array([[all_edges[i],unit_normal[0]*load_value,unit_normal[1]*load_value ,0,int(all_vertices_ind[2*i]),0.5*(all_vertices[2*i+1][0]-all_vertices[2*i][0]), 0.5*(all_vertices[2*i+1][1]-all_vertices[2*i][1]) ]]),axis =0)  
        
        if True:
            if i > 121:
                print(i)
                load_value = 6.87*abs((16.03-0.5*(all_vertices[2*i][1] + all_vertices[2*i+1][1]))*(all_vertices[2*i][0]- all_vertices[2*i+1][0]))

                Live_loads['Edges']['concentrated'] = np.append(Live_loads['Edges']['concentrated'], \
            np.array([[all_edges[i],unit_normal[0]*load_value,unit_normal[1]*load_value ,0,int(all_vertices_ind[2*i]),0.5*(all_vertices[2*i+1][0]-all_vertices[2*i][0]), 0.5*(all_vertices[2*i+1][1]-all_vertices[2*i][1]) ]]),axis =0)  
                
    
    
    end_loading_ext = time.time() # stop recording time
    print("Loading extraction",end_loading_ext-start_loading_ext, "seconds")
    
    #print(Dead_loads['Blocks']['volume'][:,2])
    #Dead_loads['Blocks']['volume'][:,2] = Dead_loads['Blocks']['volume'][:,2]*5
   
    
    data["Dead_loads"]= Dead_loads
    data["Live_loads"]= Live_loads
    
    
    return data 