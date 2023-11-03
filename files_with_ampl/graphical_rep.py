def graphical_rep(data,effort_vect,disp_vect,lag_1,lag_2, scale_factor,graph_title): # A ne pas oublier: il faut tenir compte des zones ayant une cohesion et un angle de frotement différents dans la représentation graphique "basculement, glissement, aucun critère !!"
    
    
    from shapely.geometry import Polygon
   
    import time
   

    from shapely.geometry import Polygon
    import math
    import time
    import numpy as np
    from shapely.geometry import Point
    import shapely.geometry.polygon 
    from matplotlib.patches import Polygon
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.colors import Normalize
    from matplotlib.colors import BoundaryNorm
    from matplotlib.cm import ScalarMappable
    import matplotlib.pyplot as plt
   

    
    
    blockedgesInd = data["blockedgesInd"]
    activeEdgesInd = data["activeEdgesInd"]
    edgesIndCoor = data["edgesIndCoor"]
    cohesion = data["cohesion"] 
    friction_ang = data["friction_angle"] 
    nodes = data["nodes"]
    blockscenters = data["blockscenters"]
    allpolygons = data["allpolygons"] 
    strength = data["strength"]
    kinsol = disp_vect*scale_factor
    edgeslength = data["edgeslength"]  
    
    if max(max( lag_1),- min(lag_1)) >0:
        lag_1 = lag_1/ max(max( lag_1),- min(lag_1))
    if len(lag_2) != 0:
        if max(max( lag_2),- min(lag_2)) >0:
        
            lag_2 = lag_2/ max(max( lag_2),- min(lag_2))
    print(kinsol)
    xlimmin = float(nodes[0,1])
    xlimmax=  float(nodes[0,1])
    ylimmin=  float(nodes[0,2])
    ylimmax = float(nodes[0,2])   
    for element in range(nodes.shape[0]):
        xlimmin = min([xlimmin,float(nodes[element,1])])
        xlimmax = max([xlimmax,float(nodes[element,1])])
        ylimmin = min([ylimmin,float(nodes[element,2])])
        ylimmax = max([ylimmax,float(nodes[element,2])])


    statsol = effort_vect

    shearvalues = np.empty((0,1))
    normalvalues = np.empty((0,1))
    momentvalues = np.empty((0,1))
    rup_test = np.empty((0,1))
    normalstress = np.empty((0,1))
    
    for element in range(0,len(statsol)-2,3):
        
        
        ind_strength = int(np.where(strength[:,0] == str(int(activeEdgesInd[int(element/3)])))[0])
        

        shearvalues = np.append(shearvalues,statsol[element])
        normalvalues = np.append(normalvalues,statsol[element+1])
        normalstress = np.append(normalstress,statsol[element+1]/edgeslength[int(element/3)])
        momentvalues = np.append(momentvalues,statsol[element+2])
        if len(lag_2) !=0:
        
            rup_test = np.append(rup_test,1*(round(abs(lag_1[int((element/3)*6)]),4) > 0 or  round(abs(lag_1[int((element/3)*6)+1]),4) > 0 ) + \
                                  2*(round(abs(lag_2[int((element/3))]),4) > 0 or round(abs(lag_2[int(len(lag_2)/2) + int(element/3) ]),4) > 0 ) + 4*(round(abs(lag_1[int(element/3)*6+4]),4) > 0) + \
                                  8*(round(abs(lag_1[int(element/3)*6+5]),4) > 0))
        else:
            rup_test = np.append(rup_test,1*(round(abs(lag_1[int((element/3)*6)]),4) > 0 or  round(abs(lag_1[int((element/3)*6)+1]),4) > 0 ) + \
                                  2*(round(abs(lag_1[int((element/3)*6)+2]),4) > 0 or round(abs(lag_1[int((element/3)*6)+3]),4) > 0 ) + 4*(round(abs(lag_1[int(element/3)*6+4]),4) > 0) + \
                                  8*(round(abs(lag_1[int(element/3)*6+5]),4) > 0))
            
            
            #rup_test =  np.append(rup_test,1*(round(abs(statsol[element]),3) >=round(statsol[element+1]*math.tan(math.radians(friction_ang)) + cohesion,3)) + \
             #                     2*(round(abs(statsol[element+2]),3)>=round(statsol[element+1],3) ) + 4*(round(statsol[element+1],3) >= round(float(strength[ind_strength][3])*edgeslength[int(element/3)],3)) + \
              #                    8*(round(statsol[element+1],3) <= round(float(strength[ind_strength][4])*edgeslength[int(element/3)],3) ))     
    data["rup_test"] = rup_test
    
    if True:
        # Normalize the data
        #norm = Normalize(vmin=shearvalues.min(), vmax=shearvalues.max())
        fig, ax = plt.subplots(figsize=(12, 6))
        # Normalize the data
        boundaries = np.linspace(shearvalues.min(), shearvalues.max(), num=5)

        # Create a BoundaryNorm object that maps data values to colors based on the boundaries
        norm = BoundaryNorm(boundaries, ncolors=256)

        # Create a ScalarMappable object
        sm = ScalarMappable(norm=norm, cmap='rainbow')
        # get the color of each element data
        datacolors = sm.to_rgba(shearvalues)

        plt.xlim(0,2)
        plt.ylim(0, 2)
        for  element in edgesIndCoor:
            if not len(np.where(activeEdgesInd == int(edgesIndCoor[element][0][0]))[0]) == 0:
                ind = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0])) #find edge index in activeedges array
                edgecolor = datacolors[ind][0]
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],color = (edgecolor[0],edgecolor[1],edgecolor[2]) )
            else:
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],'k' )

        # Add a colorbar
        ax.set_xlim(xlimmin-0.1, xlimmax+0.1)
        ax.set_ylim(ylimmin-0.1, ylimmax+0.1)


        #plt.title(graph_title)
        # Create a colorbar
        cbar = plt.colorbar(sm)

        # Add a label to the colorbar
        cbar.set_label("Effort tranchant aux joints (KN)")
        ax.set_axis_off()
        plt.savefig(graph_title+'_shear.jpg', dpi=300)

        plt.show(block = False)

        fig, ax = plt.subplots(figsize=(12, 6))

        # Normalize the data
        #norm = Normalize(vmin=normalvalues.min(), vmax=normalvalues.max())

        # Normalize the data
        boundaries = np.linspace(normalvalues.min(), normalvalues.max(), num=5)

        # Create a BoundaryNorm object that maps data values to colors based on the boundaries
        norm = BoundaryNorm(boundaries, ncolors=256)


        # Create a ScalarMappable object
        sm = ScalarMappable(norm=norm, cmap='rainbow')
        # get the color of each element data
        datacolors = sm.to_rgba(normalvalues)

        plt.xlim(xlimmin-0.1, xlimmax+0.1)
        plt.ylim(ylimmin-0.1, ylimmax+0.1)
        for  element in edgesIndCoor:
            if not len(np.where(activeEdgesInd == int(edgesIndCoor[element][0][0]))[0]) == 0:
                ind = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0])) #find edge index in activeedges array
                edgecolor = datacolors[ind][0]
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],color = (edgecolor[0],edgecolor[1],edgecolor[2]) )
            else:
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],'k' )

        # Add a colorbar
        ax.set_xlim(xlimmin-0.1, xlimmax+0.1)
        ax.set_ylim(ylimmin-0.1, ylimmax+0.1)
        #plt.title(graph_title)
        # Create a colorbar
        cbar = plt.colorbar(sm)

        # Add a label to the colorbar
        cbar.set_label("Effort normal aux joints (KN)")
        ax.set_axis_off()
        plt.savefig(graph_title+'_normalload.jpg', dpi=300, format='jpg')

        plt.show(block = False)


        #Normal stress 

        fig, ax = plt.subplots(figsize=(12, 6))

        # Normalize the data
        #norm = Normalize(vmin=normalvalues.min(), vmax=normalvalues.max())

        # Normalize the data
        boundaries = np.linspace(normalstress.min(), normalstress.max(), num=5)

        # Create a BoundaryNorm object that maps data values to colors based on the boundaries
        norm = BoundaryNorm(boundaries, ncolors=256)


        # Create a ScalarMappable object
        sm = ScalarMappable(norm=norm, cmap='rainbow')
        # get the color of each element data
        datacolors = sm.to_rgba(normalstress)

        plt.xlim(xlimmin-0.1, xlimmax+0.1)
        plt.ylim(ylimmin-0.1, ylimmax+0.1)
        for  element in edgesIndCoor:
            if not len(np.where(activeEdgesInd == int(edgesIndCoor[element][0][0]))[0]) == 0:
                ind = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0])) #find edge index in activeedges array
                edgecolor = datacolors[ind][0]
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],color = (edgecolor[0],edgecolor[1],edgecolor[2]) )
            else:
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],'k' )

        # Add a colorbar
        ax.set_xlim(xlimmin-0.1, xlimmax+0.1)
        ax.set_ylim(ylimmin-0.1, ylimmax+0.1)
        #plt.title(graph_title)
        # Create a colorbar
        cbar = plt.colorbar(sm)

        # Add a label to the colorbar
        cbar.set_label("Contrainte aux joints (KPa.m)")
        ax.set_axis_off()
        plt.savefig(graph_title+'_normalstress.jpg', dpi=300, format='jpg')

        plt.show(block = False)

        #Moment 


        fig, ax = plt.subplots(figsize=(12, 6))

        # Normalize the data
        #norm = Normalize(vmin=momentvalues.min(), vmax=momentvalues.max())

        # Normalize the data
        boundaries = np.linspace(momentvalues.min(), momentvalues.max(), num=5)

        # Create a BoundaryNorm object that maps data values to colors based on the boundaries
        norm = BoundaryNorm(boundaries, ncolors=256)


        # Create a ScalarMappable object
        sm = ScalarMappable(norm=norm, cmap='rainbow')

        # get the color of each element data
        datacolors = sm.to_rgba(momentvalues)

        plt.xlim(0,2)
        plt.ylim(0, 2)
        for  element in edgesIndCoor:
            if not len(np.where(activeEdgesInd == int(edgesIndCoor[element][0][0]))[0]) == 0:
                ind = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0])) #find edge index in activeedges array
                edgecolor = datacolors[ind][0]
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],color = (edgecolor[0],edgecolor[1],edgecolor[2]) )
            else:
                plt.plot([float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][2])],  [float(edgesIndCoor[element][0][3]),float(edgesIndCoor[element][0][4])],'k' )

        # Add a colorbar
        ax.set_xlim(xlimmin-0.1, xlimmax+0.1)
        ax.set_ylim(ylimmin-0.1, ylimmax+0.1)
        #plt.title(graph_title)
        # Create a colorbar
        cbar = plt.colorbar(sm)

        # Add a label to the colorbar
        cbar.set_label("Moment aux joints (KN.m/m)")
        ax.set_axis_off()
        plt.savefig(graph_title+'_moment.jpg', dpi=300, format='jpg')

        plt.show(block = False)

    fig, ax = plt.subplots(figsize=(12, 6))

    # Normalize the data
    #norm = Normalize(vmin=rup_test.min(), vmax=rup_test.max())

    # Create a ScalarMappable object
    #sm = ScalarMappable(norm=norm, cmap='rainbow')

    # get the color of each element data
    #datacolors = sm.to_rgba(rup_test)
    
    
     
    labelcolors = ['k', 'r', 'g', 'b', 'c', 'm', 'pink', 'indigo', 'purple' , 'y', 'orange','gray', 'teal', 'olive','lime' , 'navy']
    #leg = plt.legend(['Aucun critère','Glissement(G)','Basculement (B)' ,'G + B ','tassement (T)','T+G','T+B','T+G+B','Soulevement(S)','G+S','S+B', 'S+G+B','T+S','T+G+S','T+B+S','T+B+S+G'],loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
  

    
    def rotate_point(point, center, angle):
        #angle = np.radians(angle)
        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                    [np.sin(angle), np.cos(angle)]])
        translated_point = np.array([point[0] - center[0], point[1] - center[1]])
        rotated_point = np.dot(rotation_matrix, translated_point)
        return rotated_point + center
    
    xlimmin = np.inf
    ylimmin = np.inf
    
    xlimmax = -np.inf
    ylimmax = -np.inf
    
    for  element in edgesIndCoor:

        if not len(np.where(activeEdgesInd == int(edgesIndCoor[element][0][0]))[0]) == 0:

            #edgeind = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0]))[0][0] #find edge index in activeedges array
            ind = blockedgesInd[element]-1 #block index 
            indcol = np.where (activeEdgesInd == int(edgesIndCoor[element][0][0])) #find edge index in activeedges array
            (x, y) = (blockscenters[2*ind],blockscenters[2*ind+1])

            point1 = [kinsol[3*ind],kinsol[3*ind+1]]+ rotate_point((float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][3])), (x,y), kinsol[3*ind+2])
            point2 = [kinsol[3*ind],kinsol[3*ind+1]] + rotate_point((float(edgesIndCoor[element][0][2]),float(edgesIndCoor[element][0][4])), (x,y), kinsol[3*ind+2])
            plt.plot([point1[0],point2[0]], \
                     [point1[1],point2[1]],labelcolors[int(rup_test[indcol])] )
        else:
            ind = blockedgesInd[element]-1 #block index
            (x, y) = (blockscenters[2*ind],blockscenters[2*ind+1])

            point1 = [kinsol[3*ind],kinsol[3*ind+1]] + rotate_point((float(edgesIndCoor[element][0][1]),float(edgesIndCoor[element][0][3])), (x,y), kinsol[3*ind+2])
            point2 = [kinsol[3*ind],kinsol[3*ind+1]]  + rotate_point((float(edgesIndCoor[element][0][2]),float(edgesIndCoor[element][0][4])), (x,y), kinsol[3*ind+2])
            
            xlimmin = min([xlimmin,min([point1[0],point2[0]])])
            xlimmax = max([xlimmax,max([point1[0],point2[0]])])
            
            ylimmin = min([ylimmin,min([point1[1],point2[1]])])
            
            ylimmax = max([ylimmax,max([point1[1],point2[1]])])
            plt.plot([point1[0],point2[0]], \
                     [point1[1],point2[1]],'k' )
    


    # Add a label to the colorbar
    #cbar.set_label("Moment aux joints")
    
    labelcolors = ['k', 'r', 'g', 'b', 'c', 'm', 'pink', 'indigo', 'purple' , 'y', 'orange','gray', 'teal', 'olive','lime' , 'navy']
    leg = plt.legend(['Aucun critère','Glissement(G)','Basculement (B)' ,'G + B ','tassement (T)','T+G','T+B','T+G+B','Soulevement(S)','G+S','S+B', 'S+G+B','T+S','T+G+S','T+B+S','T+B+S+G'],loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
  
    for i, j in enumerate(leg.legendHandles):
        
        print(i,j)
        j.set_color(labelcolors[i])

    # Add a colorbar
    ax.set_xlim(xlimmin-0.1, xlimmax+0.1)
    ax.set_ylim(ylimmin-0.1, ylimmax+0.1)
    #plt.title(graph_title)
    # Create a colorbar
    ax.set_axis_off()
    plt.savefig(graph_title+'_edgesdisp.jpg', dpi=300, format='jpg')
    
    plt.show(block = False)
    
    
    
    
    fig, ax = plt.subplots(figsize=(12, 6))
    

    for ind, element in enumerate(allpolygons):
        polygonplot = Polygon(allpolygons[element], color='w',linestyle=':' ,lw=1, ec='black', fill='false') #use polygon function of mathplotlib to put the polygon on the suitable form for plotting    
        #ax.add_patch(polygonplot) # plot the blockj
        
       

    for ind, element in enumerate(allpolygons):
        polygonplot = Polygon(allpolygons[element], lw=1, ec='black', fill='false',facecolor="none") #use polygon function of mathplotlib to put the polygon on the suitable form for plotting    
        polygonplotdef = polygonplot
        
        #ax.add_patch(polygonplot) # plot the blockj
        
        (x, y) = (blockscenters[2*ind],blockscenters[2*ind+1])

        # Create the transformation
        transform = mtransforms.Affine2D().translate(-x, -y)+ mtransforms.Affine2D().translate(kinsol[3*ind], kinsol[3*ind+1]) + \
        mtransforms.Affine2D().rotate_deg(kinsol[3*ind+2]*180/3.14) + mtransforms.Affine2D().translate(x, y) + ax.transData

        # Apply the transformation to the polygon
        polygonplotdef.set_transform(transform)
        #transform =  mtransforms.Affine2D().translate(kinsol[3*ind], kinsol[3*ind+1]) + ax.transData
        #polygonplot.set_transform(transform)
        ax.add_patch(polygonplotdef) # plot the blockj after displacement

    
    plt.xlim(xlimmin-0.5, xlimmax+0.5)
    plt.ylim(ylimmin-0.5, ylimmax+0.5)
    #plt.title(graph_title)
    ax.set_axis_off()
    plt.savefig(graph_title+'_blocksdisp.jpg', dpi=300, format='jpg')
    
    plt.show(block = False)