
import rhinoscriptsyntax as rs
import gmsh
import numpy as np
import importlib
import plotly.graph_objects as go
import meshio
from cvxopt import matrix, solvers, spmatrix
from scipy.sparse import coo_matrix
import sys
import time
#import pandas as pd
import os
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad, VtkVertex
import mesh


def calculate_moment(point_coords, force_coords, applied_point_coords):
    """
    Calculate the moment of a force around a point.

    Parameters:
    - point_coords: Coordinates of the point where the moment is calculated (numpy array).
    - force_coords: Coordinates of the force vector (numpy array).
    - applied_point_coords: Coordinates of the point where the force is applied (numpy array).

    Returns:
    - moment: Moment of the force around the point (numpy array).
    """
    r = - point_coords + applied_point_coords
    F = force_coords

    # Calculate the cross product to get the moment
    moment = np.cross(r, F)

    return moment
def rotate_x(theta):
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])

def rotate_y(theta):
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

def rotate_z(theta):
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
def transform_point(point, displacement, rotations):
    

    # Apply rotations
    for rotation in rotations:
        axis, angle = rotation[0], rotation[1]
        if axis == 'x':
            rotation_matrix = rotate_x(angle)
        elif axis == 'y':
            rotation_matrix = rotate_y(angle)
        elif axis == 'z':
            rotation_matrix = rotate_z(angle)
        else:
            raise ValueError("Invalid rotation axis. Use 'x', 'y', or 'z'.")

        point = np.dot(rotation_matrix, point)
        # Apply displacement
    point = np.array(point) + np.array(displacement)
    return point.tolist()

def check_outward_normal(block_centroid,face_centroid,face_normal,local_ref):

    block_centroid = [ block_centroid[0], block_centroid[1], block_centroid[2]]
   
    if np.dot(face_normal, np.array(block_centroid )- np.array(face_centroid)) > 0:
        face_normal = [-1*element for element in face_normal]
        
        
    if face_normal[1] != 0 or face_normal[2] != 0:
        x_local = [face_normal[1]**2 + face_normal[2]**2, -face_normal[0]*face_normal[1] , -face_normal[0]*face_normal[2] ]
        y_local = [0, face_normal[2] , -face_normal[1] ]
    else:
        if face_normal[0] > 0:
            x_local = [0,1,0]
            y_local = [0,0,1]
        else:
            x_local = [0,-1,0]
            y_local = [0,0,1]
  
        
    x_local /= np.linalg.norm(x_local)
    y_local /= np.linalg.norm(y_local)

    x_local = [round(element,4) for element in x_local ]
    y_local = [round(element,4) for element in y_local ]
    

    local_ref = np.append(local_ref, np.array([[ face_ind, face_normal[0], face_normal[1],face_normal[2], \
                                              x_local[0],x_local[1],x_local[2], y_local[0],y_local[1],y_local[2]  ]]), axis =0)
        
  
    return local_ref

def find_element_in_matrix(matrix, element):
    for i, row in enumerate(matrix):
        if element in row:
            return i  # Return the index of the row where the element is found
    return -1  # Return -1 if the element is not found in any row

def stat_parav(all_unique_points,noc_triangles_coor,noc_triangles,triangles_num,all_triangles,Ntriangles_num,stat_sol):
    FILE_PATH = "C:\\Users\\mmoussa\\Desktop\\rhino_test\\stat_stress"
    print("Running unstructured...")
    
    normal_load = [ round(stat_sol[ind+2],2)   for ind in range(0,len(stat_sol[0:-1:]),3) ]
    shear_loadx = [ round(stat_sol[ind],2)   for ind in range(0,len(stat_sol[0:-1:]),3) ]
    shear_loady = [ round(stat_sol[ind+1],2)   for ind in range(0,len(stat_sol[0:-1:]),3) ]
    
    
    
    all_unique_points = np.array(all_unique_points)

    # Define vertices
    x = [float(element) for element in all_unique_points[:, 0]]
    y = [float(element) for element in all_unique_points[:, 1]]
    z = [float(element) for element in all_unique_points[:, 2]]

    
    x = np.array(x, dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z = np.array(z,dtype = np.float64)
    

    # Define connectivity or vertices that belong to each element
    conn = np.zeros(3*triangles_num, dtype = np.int32)
    # Define offset of the last vertex of each element
    offset = np.zeros(triangles_num , dtype = np.int32)
    # Define cell types
    ctype = np.zeros(triangles_num)
   
    
    tri_ind = 0
 
    for face_triangles in all_triangles:
      
        for triangle in face_triangles:
        
            conn[3*tri_ind], conn[3*tri_ind + 1], conn[3*tri_ind +2] = int(triangle[0]-1), int(triangle[1]-1), int(triangle[2]-1)  # first triangle
            offset[tri_ind] = 3*tri_ind +3
            ctype[tri_ind] = VtkTriangle.tid
            tri_ind += 1

    
    #cd = np.random.rand(2)
    cellData = {"pressure": []}

    # Define displacement components
    normal_stress = np.array(normal_load)
    shear_stress1 = np.array(shear_loadx)
    shear_stress2 = np.array(shear_loady)
    

    # Combine displacement components into separate arrays
    pointData = {"normal_stress":normal_stress, "shear_stress1": shear_stress1, "shear_stress2": shear_stress2  }

  
    
    # Add combined displacement to pointData
    #pointData["displacement"] = displacement
    
    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype, pointData=pointData)
    
    
    FILE_PATH = "C:\\Users\\mmoussa\\Desktop\\rhino_test\\stat_body"
    
    
  
    # Define vertices
    
    
    noc_triangles_coor = np.array(noc_triangles_coor)

    # Define vertices
    x = [float(element) for element in noc_triangles_coor[:, 0]]
    y = [float(element) for element in noc_triangles_coor[:, 1]]
    z = [float(element) for element in noc_triangles_coor[:, 2]]

   
    x = np.array(x, dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z = np.array(z,dtype = np.float64)
    

    # Define connectivity or vertices that belong to each element
    conn = np.zeros(3*Ntriangles_num, dtype = np.int32)
    # Define offset of the last vertex of each element
    offset = np.zeros(Ntriangles_num , dtype = np.int32)
    # Define cell types
    ctype = np.zeros(Ntriangles_num)
   
    
    tri_ind = 0
 
    for face_triangles in noc_triangles:
      
        for triangle in face_triangles:
        
            conn[3*tri_ind], conn[3*tri_ind + 1], conn[3*tri_ind +2] = int(triangle[0]-1), int(triangle[1]-1), int(triangle[2]-1)  # first triangle
            offset[tri_ind] = 3*tri_ind +3
            ctype[tri_ind] = VtkTriangle.tid
            tri_ind += 1
    
    pointData = {'d': np.ones(x.shape[0])}

    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype, pointData=pointData)


def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()



def block_disp_plot(manifold_data,blocks,faces,points,ub):
    ub = [ 1*element for element in ub ]
 
    
    fig = go.Figure()
    new_points = []
    initial_points = np.empty((0,3))
    points_disp = np.empty((0,3))


    FILE_PATH = "C:\\Users\\mmoussa\\Desktop\\rhino_test\\disp"
  
    # Define connectivity or vertices that belong to each element
    conn = []
    # Define offset of the last vertex of each element
    offset = []
    # Define cell types
    ctype = []
  
    for ind, block_centroid in enumerate(blocks_centroid):
      
        faces_ind = blocks[np.where( blocks[:,0] == ind+1)][:,1]
        for face in faces_ind:
            points_ind = faces[ np.where(faces[:,0] == face)][:,2]
            points_ind = [ int(element) -1 for element in points_ind] 
            
            points_coor = points[points_ind][:,-3:]
            points_coor.astype(float)
            points_coor = np.array([[float(element) for element in row] for row in points_coor])
            

            all_new = []
            reference_point = block_centroid[-3:]
       
            
            displacement_vector = [ ub[6*ind], ub[6*ind +1], ub[6*ind +2] ]
      
            
            rotations = [('x', ub[6*ind +3]), ('y',ub[6*ind +4]),('z',  ub[6*ind +5])]
            
            for element in points_coor:
                #conn.extend( initial_points.shape[0])
                #initial_points.append(element)
                initial_points = np.append(initial_points, np.array([element]), axis = 0)
                

                conn.extend( [initial_points.shape[0] - 1 ])
                
                point_new_ref = [element[element1] - reference_point[element1] for element1 in range(3) ]
               
                
                new_point = transform_point(point_new_ref , displacement_vector, rotations)
                
               
                new_point = [new_point[element] + reference_point[element] for element in range(3) ]
                
                all_new.append(new_point)

                points_disp = np.append(points_disp, np.array( [new_point[:] - element[:] ] ), axis = 0)
                #points_disp.append( new_point[:] - element[:] ) 
                
            ctype.extend ([VtkQuad.tid])
            offset.extend( [initial_points.shape[0] ])
            
            new_points.append(all_new[0]) 
            
            all_new = np.array(all_new)   
            
            

            fig.add_trace( go.Scatter3d(x= all_new[:,0] ,
                y= all_new[:,1] ,
                z= all_new[:,2] ,
                mode='lines',
               
                line=dict(color='blue', width=2),
              showlegend=False  ))

            
            # Add the second mesh3d surface to the existing figure
            
            fig.add_trace(go.Mesh3d(x=all_new[:,0],
                    y=all_new[:,1],
                    z=all_new[:,2],

                    i = [element for element in range(points_coor.shape[0]-2)],
                    j = [element for element in range(1,points_coor.shape[0]-1)],
                    k = [points_coor.shape[0]-1 for element in range(1,points_coor.shape[0]-1) ],
                    opacity=0.9,
                    color='grey',
                    flatshading=True,
                    alphahull=-1   )  )
            
    
    points = new_points
   
    point = np.array([[float(element) for element in row] for row in points])
    max_coor = np.max( point,axis = 0 )
    min_coor = np.min( point,axis = 0 )
    max_dist = max([max_coor[0] -min_coor[0], max_coor[1] -min_coor[1], max_coor[2] -min_coor[2]  ])
  
    center = np.mean( np.array([max_coor,min_coor]) ,axis = 0)
    fig.update_layout(scene=dict(xaxis=dict(range=[center[0] - max_dist/2, center[0] + max_dist/2]),
                                 yaxis=dict(range=[center[1] - max_dist/2, center[1] + max_dist/2]),
                                 zaxis=dict(range=[center[2] - max_dist/2, center[2] + max_dist/2])))
    fig.show()
      # Save the figure as an HTML file
    fig.write_html("C:/Users/mmoussa/Desktop/rhino_test/figure_disp.html")  

  #paraview part
    
    # Define vertices
    x = [float(element) for element in initial_points[:, 0]]
    y = [float(element) for element in initial_points[:, 1]]
    z = [float(element) for element in initial_points[:, 2]]

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # Define displacement components
    dx = np.array(points_disp[:,0])
    dy = np.array(points_disp[:,1])
    dz = np.array(points_disp[:,2])

    pointData = {'dx': dx, 'dy': dy, 'dz': dz } 
    conn = np.array(conn)
    offset = np.array(offset)
    ctype = np.array(ctype)
    

    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype, pointData=pointData)


def blocks_plot(face_points,fig):
    



    fig.add_trace(go.Scatter3d(x=face_points[:,0],
            y=face_points[:,1],
            z=face_points[:,2],
            mode='lines',
           
            line=dict(color='blue', width=2),
          showlegend=False   ))

    
    fig.add_trace( go.Mesh3d( x=face_points[:,0],
            y=face_points[:,1],
            z=face_points[:,2],

            i = [0],
            j = [1],
            k = [2 ],
            opacity=0.9,
            color='grey',
            flatshading=True,
            alphahull=-1  ) )
    
           
                   
    return fig

# Reload the RhinoScriptSyntax module
importlib.reload(rs)
objs = rs.AllObjects()

faces_rep = []
pts_coor = np.empty((0,3))

block_ind = 1
face_ind = 1

blocks = np.empty((0, 7), dtype=np.int32)
faces = np.empty((0,3))
blocks_centroid = np.empty((0,4))
blocks_volume = np.empty((0,2))
blocks_att = []

unsorted_faces_rep = []
plans_data = np.empty((0,4))
plan_ind = 1
f_points = [] #faces points
keys_of_interest = ["fx", "fy", "fz", "mx", "my", "mz"]
local_ref = np.empty((0,10))
blocks_brep = [] 
bars = []
supports_type = [] 

contacts_ind = []
contacts_nodes = []
contacts_FE = [0]

faces_ind = []
faces_nodes = []
faces_FE = [0]

phi_type = []
coh_type = []
lc_type = []
for obj in objs:
    obj_type = rs.ObjectType(obj)
    if obj_type == 4:
        bars.append(obj) 

    if obj_type == 8: # Check if the object is a surface
      
        points = rs.SurfacePoints(obj)
     

        for i, point in enumerate(points):
            print(point)
            plans_data = np.append(plans_data, np.array([[ plan_ind, np.round(point[0],5) , np.round( point[1],5), np.round(point[2],5) ]]), axis = 0)
        
        plan_ind += 1
        att_lcase = [item.lower() for item in  rs.GetUserText(obj)]
       

        if 'type' in att_lcase:
            index = att_lcase.index('type')
            if rs.GetUserText(obj, rs.GetUserText(obj)[index]) == '3':
                supports_type.append(3)
            elif rs.GetUserText(obj, rs.GetUserText(obj)[index]) == '2':
                supports_type.append(2)
        else:
            supports_type.append(-1)

        
        if 'c' in att_lcase:
            index = att_lcase.index('c')
            coh_type.append(float(rs.GetUserText(obj, rs.GetUserText(obj)[index])))
            
        else:
            coh_type.append(-1)
        
        if 'phi' in att_lcase:
            index = att_lcase.index('phi')
            phi_type.append(float(rs.GetUserText(obj, rs.GetUserText(obj)[index])))
            
        else:
            phi_type.append(-1)
        
        if 'lc' in att_lcase:
            index = att_lcase.index('lc')
            lc_type.append(float(rs.GetUserText(obj, rs.GetUserText(obj)[index])))
            
        else:
            lc_type.append(-1)
            

   
    if obj_type == 16:  # Check if the object is a solid
        block_name = rs.ObjectName(obj)
        
        blocks_brep.append(obj)
        # Explode the solid into its faces
        faces_data = rs.ExplodePolysurfaces(obj)
        
        if faces_data:
            for face in faces_data:

                boundary_curve_id = rs.DuplicateSurfaceBorder(face)
                if not boundary_curve_id:
                    print("Failed to extract surface boundary.")
                    
    
                # Get points along the boundary curve
                points = rs.CurvePoints(boundary_curve_id)
                points = points[0:-1:]
                # Get the control points of each face
                #points = rs.SurfacePoints(face)
                
                face_pts = np.empty((0,3))
                # Loop through the control points
                for i, point in enumerate(points):
                  
                    # Print the coordinates of each control point
                    face_pts = np.append(face_pts, np.array([[ np.round(point[0],5), np.round(point[1],5), np.round(point[2],5) ]]), axis = 0)
                 
              
                
                b_centroid = rs.SurfaceVolumeCentroid(obj)[0]
                
                face_centroid = np.mean(face_pts, axis = 0)
                
                face_normal = np.cross(face_centroid - face_pts[0], face_pts[1] - face_pts[0])
            
                face_normal /= np.linalg.norm(face_normal)
                
                
                
                
                
                pts_coor = np.append(pts_coor, face_pts, axis=0)

                # Get unique points and their indices
                pts_coor , unique_indices = np.unique(pts_coor, axis=0, return_index=True)

                # Reorder the unique points based on their original order
                pts_coor = pts_coor[np.argsort(unique_indices)]

                added_positions = []

                for point in face_pts:
                    pt_position = np.where(np.all(pts_coor == point, axis=1))[0][0]
                    added_positions.append(pt_position +1 )
                
                

                #added_positions_copy = [element +1 for element in added_positions_copy]

                unsorted_face_rep = '#'.join(map(str, added_positions.copy() ))

                added_positions_copy = added_positions.copy()

                sorted_point_indices = sorted(set(added_positions_copy))
                #sorted_point_indices = [element + 1 for element in sorted_point_indices]

                face_rep = '#'.join(map(str, sorted_point_indices ))
               
                if (not face_rep in faces_rep) or block_ind == 1:
                    faces_rep.append(face_rep)
                    f_points.append(face_pts)
                    unsorted_faces_rep.append(unsorted_face_rep)

                    faces_ind += [faces_rep.index(face_rep) + 1 ]
                    faces_nodes += added_positions
                    faces_FE.append(faces_FE[-1:][0] + len(added_positions))


                    local_ref = check_outward_normal(b_centroid,face_centroid,face_normal,local_ref)
                    local_ref[-1:, 0] = len(faces_rep)
                
                    block_name = -99
                    blocks = np.append(blocks, np.array([[block_ind, len(faces_rep), block_name, face_ind  , 0,-1, 11 ]]), axis = 0)
                    for ind in added_positions:
                        faces = np.append(faces, np.array([[len(faces_rep), face_ind, ind  ]]),axis = 0) #ind or ind+1
                else:
                    contacts_ind += [faces_rep.index(face_rep) + 1 ]
                    contacts_nodes += added_positions
                  
                    contacts_FE.append(contacts_FE[-1:][0] + len(added_positions))

                    block_name = -99
                    blocks = np.append(blocks, np.array([[block_ind, faces_rep.index(face_rep) + 1, block_name, face_ind, 0 , -1, 22 ]]), axis = 0)
  
                # Do something with each face (e.g., print face area)
                face_area = rs.SurfaceArea(face)[1]
                
                
                
                face_ind += 1

                # Optionally, you can delete the exploded face
                rs.DeleteObject(face)
                rs.DeleteObject(boundary_curve_id)
                
        
       
        b_centroid = rs.SurfaceVolumeCentroid(obj)[0]
        b_volume = rs.SurfaceVolume(obj)[0]

        blocks_volume = np.append(blocks_volume, np.array([[block_ind, b_volume ]]), axis = 0)

        blocks_centroid = np.append(blocks_centroid, np.array([[block_ind, b_centroid[0], b_centroid[1], b_centroid[2] ]]), axis = 0)
        
        atb_keys = rs.GetUserText(obj)
        
        if atb_keys:
            
            b_att = [block_ind]

            for atb_key in atb_keys:
                if atb_key in keys_of_interest:
                    key_ind = keys_of_interest.index(atb_key)
                    

                    b_att.extend([key_ind +1, rs.GetUserText(obj, atb_key) ])
            
            blocks_att.append(b_att)

                    

   
        
        block_ind += 1
       
# Get the number of rows in the array
num_rows = pts_coor.shape[0]

# Create an array of 1-based indices
indices_column = (np.arange(num_rows) + 1).reshape(-1, 1)

# Add the indices column to the original array
points = np.hstack((indices_column, pts_coor))

plans_rep = []
supports_ind = []

if plans_data.shape[0] != 0:
    for ind in range(plan_ind-1):
     
        plan_elements = plans_data[ np.where(plans_data[:,0] == ind +1 )[0] ]
  
        plan_pts = []
        for  element in plan_elements:
            pt_coor = element[-3:]
            pt_ind = np.where( np.all(points[:,-3:] == pt_coor, axis=1)  )[0]
            print(pt_ind)
            plan_pts.extend(pt_ind)
        
        sorted_plan_pts = sorted(set(plan_pts))
        sorted_plan_pts = [element + 1 for element in sorted_plan_pts]

        plan_rep = '#'.join(map(str, sorted_plan_pts ))
        
        plans_rep.append(plan_rep)
        
    
for ind in contacts_ind:    blocks[np.where(blocks[:,1] == ind)[0], 4] = 1

print("hereee",contacts_ind)
print(supports_type)

supports_rep = set(plans_rep) &  set(faces_rep)
print(supports_type)
print(plans_rep)
for ind in supports_rep :
    if supports_type[plans_rep.index(ind)] == 2:
        
        blocks[np.where(blocks[:,1] == faces_rep.index(ind) +1)[0] , 4] = 2
        
        contacts_ind += [faces_rep.index(ind) +1]
        
        
        unsorted_rep = unsorted_faces_rep[faces_rep.index(ind)] 
      
        nodes_ind = [int(num) for num in unsorted_rep.split('#')]
        
        contacts_nodes += nodes_ind
        
        contacts_FE += [ contacts_FE[-1:][0] + len(nodes_ind)]
    
    if supports_type[plans_rep.index(ind)] == 3:

        supports_ind += [faces_rep.index(ind) +1]
        blocks[np.where(blocks[:,1] == faces_rep.index(ind) +1)[0] , 4 ] = 3

        contacts_ind += [faces_rep.index(ind) +1]
       
        unsorted_rep = unsorted_faces_rep[faces_rep.index(ind)] 
        
        nodes_ind = [int(num) for num in unsorted_rep.split('#')]
        
        contacts_nodes += nodes_ind
        
        contacts_FE += [ contacts_FE[-1:][0] + len(nodes_ind)]
    if coh_type[plans_rep.index(ind)] != -1:
        blocks[np.where(blocks[:,1] == faces_rep.index(ind) +1)[0] , 5 ] = float(coh_type[plans_rep.index(ind)])

print("hereee",contacts_ind)
# Select the curve(s) you want to split
#curve_ids = rs.GetObjects("Select curve(s) to split", rs.filter.curve)

contacts_ind_sorted = sorted(contacts_ind)


supports_pos = []
for ind in supports_ind:
    supports_pos += [contacts_ind_sorted.index(ind)]



curve_ids = bars

if not curve_ids:
    print("No curve selected. Exiting...")
    #quit()
rs.EnableRedraw(False)


bpts_inter = np.empty((0,3))
bars_data = np.empty((0,7))

if True:
    for block_ind, solid_id in enumerate(blocks_brep):
        pts_pos = set()
        # Split the curves with the solid
        for curve_id in curve_ids:
            added_positions = []
            att_num = 0

            atb_keys = rs.GetUserText(curve_id)   
            if atb_keys:
                
                for atb_key in atb_keys:
                    if atb_key == 'A':
                        bar_area = rs.GetUserText(curve_id, 'A')
                        att_num += 1
                    if atb_key == 'fy':
                        bar_fy = rs.GetUserText(curve_id, 'fy')
                        att_num += 1
                        
            if att_num == 2:
                exploded_curve = rs.ExplodeCurves(curve_id)
            
                for line in exploded_curve:
                    start_pt = rs.CurvePoints(line)[0]
                    end_pt = rs.CurvePoints(line)[1]
                
                    bar_dir = np.array([start_pt[0]- end_pt[0]  , start_pt[1]- end_pt[1], start_pt[2]- end_pt[2]  ] )
                    
                    bar_dir /= np.linalg.norm(bar_dir)  
                    bar_dir = [np.round(element, 5)  for element in bar_dir]
                    # Check if the curve intersects with the solid
                    if rs.CurveBrepIntersect(line, solid_id):
            
                        
                        inter1_coor = rs.PointCoordinates((rs.CurveBrepIntersect(line, solid_id)[1][0]))
                        inter1_coor = [inter1_coor[0], inter1_coor[1], inter1_coor[2] ]
                        part_to_add = np.array( [inter1_coor])

                        if len(rs.CurveBrepIntersect(line, solid_id)[1]) == 2:
                            inter2_coor = rs.PointCoordinates((rs.CurveBrepIntersect(line, solid_id)[1][1]))
                            inter2_coor = [inter2_coor[0], inter2_coor[1], inter2_coor[2]]
                            part_to_add = np.array( [inter1_coor, inter2_coor])
                    
                        part_to_add = np.round(part_to_add, decimals=5)

                        bpts_inter = np.append(bpts_inter, part_to_add , axis = 0)
                        # Get unique points and their indice
                        bpts_inter , unique_indices = np.unique(bpts_inter, axis=0, return_index=True)
                
                        # Reorder the unique points based on their original order
                        bpts_inter = bpts_inter[np.argsort(unique_indices)]
                    
                        # Use broadcasting to compare all points in points_array with each point in set_of_points
                        matches = np.all(np.expand_dims( bpts_inter, axis=1) == part_to_add, axis=2)

                        # Find row indices where any row in set_of_points matches points_array
                        added_positions = np.where(matches.any(axis=1))[0]
                    
                
                        
                        if len(added_positions) != 0:
                            
                            for  pt_pos in added_positions:
                                
                                pt_coor = bpts_inter[pt_pos]
                                if np.dot(bar_dir, np.array(blocks_centroid[block_ind, -3:] )- np.array(pt_coor)) > 0:
                                    bar_dir = [-1*element for element in bar_dir]
                                
                                bars_data = np.append( bars_data, np.array([[int(block_ind +1) , int(pt_pos +1), float(bar_area), float(bar_fy), *bar_dir  ]] ),axis = 0   )  
            
                rs.DeleteObjects(exploded_curve)

    rs.EnableRedraw(True)






lc = 1

lc_faces_val = lc*np.ones((len(contacts_ind)), dtype = np.float64)
print(len(contacts_ind))

for ind, lc_face  in enumerate(lc_type):
    if lc_face != -1:
        print(plans_rep)
        face_ind = faces_rep.index( plans_rep[ind] )
        contact_ind = contacts_ind.index(face_ind + 1)
        print(contact_ind)
        lc_faces_val[contact_ind] = float(lc_face)



contacts_ind = np.array(contacts_ind, dtype = np.int32)
contacts_nodes = np.array(contacts_nodes, dtype = np.int32)
contacts_FE = np.array(contacts_FE, dtype = np.int32)


contacts_ind = np.array(contacts_ind, dtype = np.int32)
faces_nodes = np.array(faces_nodes, dtype = np.int32)
faces_FE = np.array(faces_FE, dtype = np.int32)

print("lccccccccccc", lc_faces_val)



output = mesh.print_hello_pyth( contacts_ind,  np.array(blocks, dtype = np.float64) , np.array(points[:, 1:4 ], dtype = np.float64),  faces_FE, faces_nodes, np.array(blocks_centroid, dtype = np.float64) , \
                       np.array(local_ref,dtype = np.float64), lc_faces_val )




supports_nodes = set()

print("lennnnnnnnn", len(output[4]))
for pos in supports_pos:
    for tri in output[4][pos]:
            supports_nodes.update({ tri[0], tri[1], tri[2] } )

supports_nodes = list(supports_nodes)

all_coh_nodes = set()
all_coh_val = []


for ind, coh in enumerate(coh_type):
    if coh != -1:
        print(faces_rep)
        faces_pos = faces_rep.index( plans_rep[ind] )
        contact_pos = contacts_ind_sorted.index( faces_pos + 1 )

        contact_nodes = set()
        for tri in output[4][contact_pos]:
            contact_nodes.update( {tri[0], tri[1], tri[2] })
        
        all_coh_nodes.update(contact_nodes) 
        all_coh_val += [coh for i in range(len(contact_nodes))] 


all_coh_nodes = list(all_coh_nodes)



all_phi_nodes = set()
all_phi_val = []


for ind, phi in enumerate(phi_type):
    if phi != -1:
        
        faces_pos = faces_rep.index( plans_rep[ind] )
        contact_pos = contacts_ind_sorted.index( faces_pos + 1 )

        contact_nodes = set()
        for tri in output[4][contact_pos]:
            contact_nodes.update( {tri[0], tri[1], tri[2] })
        
        all_phi_nodes.update(contact_nodes) 
        all_phi_val += [phi for i in range(len(contact_nodes))] 


all_phi_nodes = list(all_phi_nodes)

print(all_phi_nodes)
print(all_phi_val)

bars_data = np.unique(bars_data,axis = 0)
#old part 
print("typeeee",phi_type)

def problem_building_up(blocks_centroid,points,faces,blocks,local_ref,blocks_att):
    kc = 0
    all_mesh_points = np.empty((0,6))
    all_unique_points = np.empty((0,6))
    
    all_triangles = np.empty((0,3))
    all_triangles_coor = []

    noc_triangles = np.empty((0,3))
    noc_triangles_coor = []
    
    equilibrium_data = []
    equilibrium_line = []
    equilibrium_col = []
    equilibrium_matrix_row = []
    equilibrium_matrix_col = []

    all_supports_ind = []
   
    faces_ind_read = []
    equilibrium_matrix = []
    fig = go.Figure()
    
    equilibrium_matrix_row = output[1]
    equilibrium_matrix_col = output[2]
    equilibrium_matrix = output[0]

    from scipy.linalg import block_diag


    
    if False:
        #PORTIOLLI PART
        all_mesh_points = np.empty((0,6))
        equilibrium_data = []
        equilibrium_line = []
        equilibrium_col = []
        equilibrium_matrix_row = []
        equilibrium_matrix_col = []
    
    
        equilibrium_matrix = []
        for block_ind, block_centroid in enumerate(blocks_centroid[:,-3:]):
            
            contact_pos = np.where( (blocks[:,0] == block_ind +1) & ( (blocks[:,4] == 1 ) | (blocks[:,4] == 2 )  ) )
            
            if len(contact_pos) != 0:
                contact_ind = [ blocks[element][1] for element in contact_pos[0]]
                contact_ind_old = [ blocks[element][3] for element in contact_pos[0]]
                
                
                
                m_fx =  []
                m_fy =[]
                m_fz = []
                equilibrium_fx = []
                equilibrium_fy = []
                equilibrium_fz = []
                
                points_number = all_mesh_points.shape[0]
                
                for ind, face in  enumerate(contact_ind):
                    
                        
                    face_data = np.where(faces[:,0] == int(face))
                    
                    face_points_ind = faces[face_data,2]
                    
                    face_points_ind = [int(element)-1 for element in face_points_ind[0]]
                    face_points = points[face_points_ind]
                    face_points_coor = face_points[:,-3:]
                   
                    face_points_coor = np.array([[float(element) for element in line]  for line in face_points_coor ]  )
                    
                    
                    face_mesh(face_points_coor)
                    mesh_points = face_mesh_read()[0]
                  
                    #mesh_points = np.array([[np.round(element,30) for element in row ] for row in mesh_points])
                    mesh_points = face_points_coor # no mesh effect
                    
                    
                    local_ref_face = local_ref[ np.where(local_ref[:,0] == contact_ind_old[ind])[0] ][0]
                    
                    #global_ref_face = [  float(local_ref_face[2]) +  float(local_ref_face[5]) + float(local_ref_face[8]), \
                     #                 float(local_ref_face[3]) +  float(local_ref_face[6]) + float(local_ref_face[9]), \
                      #                  float(local_ref_face[1]) +  float(local_ref_face[4]) + float(local_ref_face[7])]
                    
                    equilibrium_fx += [  float(local_ref_face[4]) , float(local_ref_face[7]),  float(local_ref_face[1])  ]* mesh_points.shape[0] 
                    equilibrium_fy += [   float(local_ref_face[5]) , float(local_ref_face[8]), float(local_ref_face[2]) ,]* mesh_points.shape[0]
                    equilibrium_fz += [  float(local_ref_face[6]) , float(local_ref_face[9]), float(local_ref_face[3])]* mesh_points.shape[0]
                    
                    
                    
                    x_local = [ float(local_ref_face[4]) ,  float(local_ref_face[5]) , float(local_ref_face[6])]
                    y_local = [ float(local_ref_face[7]) ,  float(local_ref_face[8]) , float(local_ref_face[9])]
                    normal =   [ float(local_ref_face[1]) ,  float(local_ref_face[2]) , float(local_ref_face[3])]
                    
                    for point in mesh_points:
                        x_local_moment = calculate_moment(block_centroid, x_local, point)
                        y_local_moment = calculate_moment(block_centroid, y_local, point)
                        normal_moment = calculate_moment(block_centroid, normal, point)
                        
                        
                        
                        
                        
                        
                        #m_fx +=  [x_local_moment[0], y_local_moment[0], normal_moment[0] ] 
                        #m_fy +=  [x_local_moment[1]*0, 0*y_local_moment[1]*0, normal_moment[1]*0 ]
                        #m_fz +=  [x_local_moment[2], y_local_moment[2], normal_moment[2] ]
                        
                        m_fx +=  [x_local_moment[0], y_local_moment[0], normal_moment[0] ] 
                        m_fy +=  [x_local_moment[1], y_local_moment[1],normal_moment[1] ]
                        m_fz +=  [x_local_moment[2], y_local_moment[2], normal_moment[2] ]
                        
                    
                    
                
                        
                    #non_zero_elements = [x for x in A if x != 0]
                    #indices_of_non_zero_elements = [i for i, x in enumerate(A) if x != 0]
    
                    part_to_add = np.array([[str(block_ind +1), str(contact_ind[ind]), str(contact_ind_old[ind])] for _ in range(mesh_points.shape[0])])
                    mesh_points = np.append(mesh_points, part_to_add,axis = 1)
                    
                    all_mesh_points = np.append( all_mesh_points, mesh_points,axis = 0)
                    
        
        
        
        
            equilibrium_block = np.array([equilibrium_fx, equilibrium_fy,equilibrium_fz, m_fx, m_fy, m_fz])
            from scipy.sparse import csr_matrix
            # Convert the array to a sparse matrix
            
            sparse_matrix = csr_matrix(equilibrium_block)
            
            # Get the non-zero elements and their indices
            nonzero_elements = sparse_matrix.data
            row_indices, col_indices = sparse_matrix.nonzero()
    
            
            equilibrium_matrix = equilibrium_matrix + list(nonzero_elements)
            
            
            row_indices = [element + 6*block_ind for element in row_indices]
            
            col_indices = [element + 3*points_number for element in col_indices]
            
            
            equilibrium_matrix_row += row_indices
            equilibrium_matrix_col += col_indices


    #END PORTIOLLI
    
    
    
    
   
    
    nb = int((max(equilibrium_matrix_row) + 1)/6)
    

    
    # Specify the number of repetitions
    n = all_mesh_points.shape[0]
    
    equilibrium_matrix_row = output[1]
    equilibrium_matrix_col = output[2]
    equilibrium_matrix = output[0]
    
   
    
    equilibrium_matrix_row = [int(element) for element in equilibrium_matrix_row]
    #equilibrium_matrix_row.insert(50,3)
 
    equilibrium_matrix_col = [int(element) for element in equilibrium_matrix_col]
    
    #equilibrium_matrix_col.insert(50,3*n)
    
    print(len(equilibrium_matrix_row))

    # Print the result
    live_load = [0]*6*nb

    #B_attributes = blocks_attributes(file)
    n = int((max(equilibrium_matrix_col)+1)/3)
    

    if True:
        for block_attributes in blocks_att:
            for ind in range(1,len(block_attributes),2):
                key = block_attributes[ind]
              
                
                val = block_attributes[ind+1]

                

                #first_ind = equilibrium_matrix_row.index( (block_attributes[0])*6 - 6 + key - 1 )
                
                live_load[(block_attributes[0])*6 - 6 + key - 1] =  float(val)
                
                equilibrium_matrix.append( float(val))
                equilibrium_matrix_row.append(  (block_attributes[0])*6 - 6 + key - 1  )
                equilibrium_matrix_col.append(3*n)
    
    
    
    
    
    #equilibrium_matrix.insert(0, 1)
   
       
    
    equilibrium_matrix_row = [int(element) for element in equilibrium_matrix_row]
    #equilibrium_matrix_row.insert(0,0)
 
    equilibrium_matrix_col = [int(element) for element in equilibrium_matrix_col]
    


    c = np.zeros((3*n +1,1))
    #c =  spmatrix( [-1], [3*n], [0], (3*n +1, 1 ) )     
    c[3*n] = -1
    c = matrix(c) # dead load vector x (-1)
    
    b = np.zeros((6*nb,1 )) 
    
    density = 1/7 # blocks density
   
    for block_ind, block_volume in enumerate(blocks_volume):
        b[6*block_ind +2] = density*block_volume[1] 
        #b[6*block_ind +3] = 3*1.4



    

  


    
    all_supports_ind = np.array(all_supports_ind)
    
    import sys
    import mosek
    import math
    # Since the actual value of Infinity is ignores, we define it solely
    # for symbolic purposes:
    inf = 0.0

    aval = equilibrium_matrix
    acol = equilibrium_matrix_col
    arow = equilibrium_matrix_row

    avalk = equilibrium_matrix
    acolk = equilibrium_matrix_col
    arowk = equilibrium_matrix_row
    #static approach
    bars_cap = []
   
    if True:

        for ind, point_coor in enumerate(bpts_inter):
            blocks_data = bars_data[ np.where(bars_data[:,1] == ind +1)[0] ]
            var_col = max(acol) + 1
         
            bars_cap.append( float(blocks_data[0][2]) * float(blocks_data[0][3] ) )
            
            for block in blocks_data:
                block_ind = int(block[0]) - 1
                block_centroid = blocks_centroid[block_ind][-3:]
                bar_dir = block[-3:]
                bar_dir = [float(element) for element in bar_dir]
               
                bar_moment = calculate_moment( block_centroid, bar_dir, point_coor )
              
                bar_tor = list(bar_dir) + list(bar_moment)
           
                nzero_pos = [ind for ind in range(6) if np.round(bar_tor[ind],5) != 0]
              
                for pos in nzero_pos:
                  
                    arow.append( block_ind*6 + pos )
                    acol.append( var_col )
                    aval.append( float(bar_tor[pos]) )
                
                
            
            
            
            
        
    
    
    

    n = int((max(equilibrium_matrix_col)+1)/3)
    

    if True:
        with mosek.Task() as task:
            task = mosek.Task() 
            # Attach a printer to the task
            task.set_Stream(mosek.streamtype.log, streamprinter)
    
            nef = n 
          
            bkx = []
            cohs = []
            #bkx = [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.up]*nef   #variables lower bounded, if free just put .fr
            for i in range(nef):
                if i not in supports_nodes:
                    bkx += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.up] #variables lower bounds
                else:
                    bkx += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.fr] #variables lower bounds

            bkx +=  [mosek.boundkey.lo]*1 + [mosek.boundkey.ra]*len(bars_cap)

            blx = []
            bux = []
            for i in range(nef):
                blx += [-inf, -inf , -20 ] #variables lower bounds
                
            blx +=  [0]*1 +  [ -0*element for element in bars_cap] 
            
            for i in range(nef):
                if i not in supports_nodes:
                    bux += [-inf, -inf , 0 ] #variables lower bounds
                else:
                    bux += [-inf, -inf , -inf ] #variables lower bounds
                
            bux += [inf]*1 + list(bars_cap) #variables upper bounds
        
           
            #ne = nb*6 + max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1 #num of equalities
            
            bkc = [mosek.boundkey.fx]*nb*6 #type of constraints (equalities here)
            blc = list(b)    #b for equalities
            buc = list(b) #b for equalities
          
            
            if False:
                bkc = [mosek.boundkey.fx]*ne #type of constraints (equalities here)
    
                
                blc = list(b) + [0]*(max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1)   #b for equalities
                buc = list(b) + [0]*(max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1) #b for equalities
            
            c = np.zeros((3*nef +1,1))
            #c =  spmatrix( [-1], [3*n], [0], (3*n +1, 1 ) )     
            c[3*nef] = -1
            
            z = [element for element in list(c)] + [0]*len(bars_cap) #objective function 

            
            eq_dense =[]
            aval = list(aval)
            arow = list(arow)
            acol = list(acol)
         
            #aval[13] *= -1
            #aval[17] *= -1
            
            
            
            if False:
                aval = equilibrium_matrix.copy() + contact_eq_element_x.copy() + contact_eq_element_y.copy() + contact_eq_element_z.copy()
                    
                arow = equilibrium_matrix_row.copy()
                init_ind = max(arow) +1
                arow += [element + init_ind   for element in contact_eq_row_x]
                init_ind = max(arow) +1
                arow += [element +  init_ind for element in contact_eq_row_y]
                init_ind = max(arow) +1
                arow += [element + init_ind for element in contact_eq_row_z]
                
                
                
                acol = equilibrium_matrix_col + contact_eq_col_x + contact_eq_col_y + contact_eq_col_z
    
            # equalities for kinematic problem
            arowk = acol.copy()
            acolk = arow.copy()
            avalk = aval.copy()

            ####### a remettre
            if False:
                for i in range(3*nef):
                    index = arowk.index(i)
                    arowk.insert(index, i)
                    acolk.insert(index, 6*nb + i)
                    avalk.insert(index, 1)
                    if ((i+1) % 3 == 0):
                        index = arowk.index(i)
                        arowk.insert(index, i)
                        acolk.insert(  index, 6*nb + 3*n + int((i+1)/3) - 1 )
                        avalk.insert(index, 1)

             ###### till here        
      
           
            num_lists = 3*nef + 1 + len(bars_cap)  # Adjust the number of inner lists as needed
    
          
            asub_mos = [[] for _ in range(num_lists)]
            aval_mos  = [[] for _ in range(num_lists)]
            
            for i,j in enumerate(acol):
                asub_mos[int(j)].append(arow[i])
                aval_mos[int(j)].append(aval[i])
            
            numvar = len(bkx)
            numcon = len(bkc)
            NUMANZ = 4 
    
         # Append 'numcon' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(numcon)
    
            #Append 'numvar' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(numvar)
    
            for j in range(numvar):
                # Set the linear term c_j in the objective.
                task.putcj(j, z[j])
                # Set the bounds on variable j
                # blx[j] <= x_j <= bux[j]
                task.putvarbound(j, bkx[j], blx[j], bux[j])
            
            
            for j in range(len(aval_mos)):
                
               
                # Input column j of A
                task.putacol(j,     asub_mos[j],aval_mos[j])            # Non-zero Values of column j.
            
            
            print("supp ind",all_supports_ind)
        # Input the affine conic constraints
            # Create a matrix F such that F * x = [x(3),x(0),x(1),x(4),x(5),x(2)] 
            task.appendafes(3*nef)
            print("nefffff", nef)
            print("suppporttt",supports_nodes)
            list_con = []
            list_coef = []
            friction_ang = 33
            for ind in range(nef):
                if (ind +1) not in supports_nodes:
                    friction_ang = 33
                    c = 0 # input value
                    
                    if (ind+1) in all_phi_nodes:
                        friction_ang = float(all_phi_val[all_phi_nodes.index(ind + 1)])
                     
                        
                    

                    if (ind+1) in all_coh_nodes:
                        c = all_coh_val[all_coh_nodes.index(ind+1)]
                        
                    cohs += [-c,0,0]
                    list_con += [3*ind + 2 , 3*ind,3*ind + 1]
                   
                    list_coef += [-math.tan(math.radians(friction_ang)),1.0,1.0]
                    
            print(list_coef)    
            task.putafefentrylist(range(3*nef - len(supports_nodes) ),                      # Rows
                                list_con ,            # Columns 
                                list_coef  )          #coefficients
    
            # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
            
            
            coh_vect = []
            for ind in range(nef - len(supports_nodes) ):
                
                #if ind not in all_supports_ind:
                if True:
                   
                    coh = cohs[3*ind: 3*ind +3]
                   
                    coh_vect += coh
                    quadcone  = task.appendquadraticconedomain(3)
                    task.appendacc(quadcone,          [3*ind, 3*ind +1, 3*ind + 2],  coh)                    # None if there is no b for conic 
           
            
           
            
            for i in range(numcon):
                task.putconbound(i, bkc[i], blc[i], buc[i])
    
    
        # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)
    
            # Optimize the task
            task.optimize()
            #task.writedata("cqo1.ptf")
            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            prosta = task.getprosta(mosek.soltype.itr)
            solsta = task.getsolsta(mosek.soltype.itr)
    
            # Output a solution
            xx = task.getxx(mosek.soltype.itr)
            ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            #print("xxx",xx)
            print(task.getsuc(mosek.soltype.itr)[0:6*nb])
            print(task.getslc(mosek.soltype.itr)[0:6*nb])
            print(ub)
            #print("conic",task.getdoty(mosek.soltype.itr))
            
            #skc, y = task.getsolution(mosek.soltype.bas)
            if solsta == mosek.solsta.optimal:
                print("Objective: %s" % xx[-1:])
                #print("Optimal solution: %s" % dual(xx))
                # Get dual variables for the linear constraints
                
              
            elif solsta == mosek.solsta.dual_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif solsta == mosek.solsta.prim_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif mosek.solsta.unknown:
                print("Unknown solution status")
            else:
                print("Other solution status")
           # Clear or reset the task
            
            #del task  # Delete the task object
    
    
    
   
    
    print("allll phi nodessssss",all_phi_nodes)
    
        
    #non associative calc
    tol =1
    itr = 0
    alpha=0.3 
    beta=0.6
    #tolerance=0.001
    tolerance=0.0001
    
    non_ass = False
    if non_ass:
    
        while tol >tolerance: 
    
            # Check if the solver found an optimal solution
            if solsta == mosek.solsta.optimal:
          
                if itr>1:
                    normalvalues_old = normalvalues
                else:
                    normalvalues_old = np.zeros((int(nef) - len(all_supports_ind),1))
    
                statsolver_old = xx[-1:][0]    
                
                normalvalues = np.empty((0,1))
               
                
                for element in range(0,len(xx)-2,3):
                    if int(element /3) not in all_supports_ind:
                        normalvalues = np.append(normalvalues,xx[element+2]) 
            
            coh_vect_new = [0,0,0]*len(coh_vect)
            friction_coef = [1,1,1]*len(list_coef)
         
            for element in range(int(nef - len(all_supports_ind))):
               
                if itr >1:
                    coh_vect_new[3*element] = 0.00001*min(normalvalues) + coh_vect[3*element] + (1+alpha)*abs(list_coef[3*element])*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element])
                    friction_coef[3*element] = math.tan(math.radians(friction_ang)*alpha)
                    
                else:
        
                    coh_vect_new[3*element] = 0.00001*min(normalvalues) + coh_vect[3*element] + (1+alpha)*abs(list_coef[3*element])*(normalvalues[element] )
                    friction_coef[3*element] = math.tan(math.radians(friction_ang)*alpha)
                    
            #-list_coef[3*element]*alpha
            statsol_old = xx
    
    
            #mosek
            with mosek.Task() as task:
                # Attach a printer to the task
                task.set_Stream(mosek.streamtype.log, streamprinter)
        
                
             # Append 'numcon' empty constraints.
                # The constraints will initially have no bounds.
                task.appendcons(numcon)
        
                #Append 'numvar' variables.
                # The variables will initially be fixed at zero (x=0).
                task.appendvars(numvar)
        
                for j in range(numvar):
                    # Set the linear term c_j in the objective.
                    task.putcj(j, z[j])
                    # Set the bounds on variable j
                    # blx[j] <= x_j <= bux[j]
                    task.putvarbound(j, bkx[j], blx[j], bux[j])
                
                
                for j in range(len(aval_mos)):
                    
                   
                    # Input column j of A
                    task.putacol(j,      asub_mos[j],aval_mos[j])            # Non-zero Values of column j.
                
                
        
            # Input the affine conic constraints
                # Create a matrix F such that F * x = [x(3),x(0),x(1),x(4),x(5),x(2)] 
                task.appendafes(3*nef)
                
                list_con = []
                

                for ind in range(nef):
                    if ind not in all_supports_ind:
                    
                        list_con += [3*ind + 2 , 3*ind,3*ind + 1]
                        #list_coef += [-math.tan(math.radians(friction_ang)),1.0,1.0]
                    
                task.putafefentrylist(range(3*nef - len(all_supports_ind) ),                      # Rows
                                    list_con ,            # Columns 
                                    friction_coef  )          #coefficients
        
                # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
                
                #coh_vect = []
                for ind in range(nef - len(all_supports_ind) ):
                    
                    #if ind not in all_supports_ind:
                    if True:
                        quadcone  = task.appendquadraticconedomain(3)
                        task.appendacc(quadcone,                 [3*ind, 3*ind +1, 3*ind + 2],   [coh_vect_new[3*ind], coh_vect_new[3*ind +1] , coh_vect_new[3*ind +2] ] )                    # None if there is no b for conic 
           
                
                for i in range(numcon):
                    task.putconbound(i, bkc[i], blc[i], buc[i])
        
        
            # Input the objective sense (minimize/maximize)
                task.putobjsense(mosek.objsense.minimize)
        
                # Optimize the task
                task.optimize()
                #task.writedata("cqo1.ptf")
                # Print a summary containing information
                # about the solution for debugging purposes
                task.solutionsummary(mosek.streamtype.msg)
                prosta = task.getprosta(mosek.soltype.itr)
                solsta = task.getsolsta(mosek.soltype.itr)
        
                # Output a solution
                xx = task.getxx(mosek.soltype.itr)
                
                ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            
                
                #skc, y = task.getsolution(mosek.soltype.bas)
                if solsta == mosek.solsta.optimal:
                    print("Objective: %s" % xx[-1:])
                    #print("Optimal solution: %s" % dual(xx))
                    # Get dual variables for the linear constraints
                    
                  
                elif solsta == mosek.solsta.dual_infeas_cer:
                    print("Primal or dual infeasibility.\n")
                elif solsta == mosek.solsta.prim_infeas_cer:
                    print("Primal or dual infeasibility.\n")
                elif mosek.solsta.unknown:
                    print("Unknown solution status")
                else:
                    print("Other solution status")
            #end mosek
    
            if solsta == mosek.solsta.optimal:
    #define obj do not forget
                if itr >1:
                    tol = abs(xx[-1:][0] - statsolver_old)/abs(xx[-1:][0] )
                    print("tolerance", tol)
            else:
                statsol = statsol_old
                
            if alpha <= 0.0001:
                alpha = 0.0001
            else:
                alpha = alpha*0.7
    
            itr = itr+1
            # Clear or reset the task
            
            del task  # Delete the task object
                


    
    #end non ass

     #kinematic approach
    if False:
        with mosek.Task() as task:
            task = mosek.Task() 
            # Attach a printer to the task
            task.set_Stream(mosek.streamtype.log, streamprinter)
    
            if False:
                bkx = [mosek.boundkey.fr]*6*nb + [mosek.boundkey.fr]*(max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1) + \
                             [mosek.boundkey.fr,mosek.boundkey.fr,mosek.boundkey.up]*n + [mosek.boundkey.lo]*n  #variables lower bounded, if free just put .fr
                blx = [-inf ]*6*nb + [-inf]*(max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1) + \
                                        [-inf,-inf, -inf]*n + [0]*n #variables lower bounds
                bux = [inf]*6*nb + [inf]*(max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1) + \
                                 [inf, inf, 0]*n + [inf]*n  #variables upper bounds
        
            bkx = [mosek.boundkey.fr]*6*nb + [mosek.boundkey.fr,mosek.boundkey.fr,mosek.boundkey.up]*nef + [mosek.boundkey.lo]*2*nef   #variables lower bounded, if free just put .fr
            blx = [-inf ]*6*nb +  [-inf,-inf, -inf]*nef + [0]*2*nef #variables lower bounds
            bux = [inf]*6*nb +  [inf, inf, 0]*nef + [inf]*2*nef   #variables upper bounds
            
            ne = 3*nef
            bkc = [mosek.boundkey.fx]*(ne +1 ) #type of constraints (equalities here)
            blc =  [0]*ne + [1]   #b for equalities
            buc = [0]*ne + [1]   #b for equalities
            
            if False:
                nvcontact = (max(contact_eq_row_x) +1 + max(contact_eq_row_y) +1 + max(contact_eq_row_z) +1)
            z = [element[0] for element in list(b)] #objective function 
            
            if False:
                z += [0]*(  3*n + n + nvcontact )
            
            z += [0,0,-0*0.2/math.tan(math.radians(33))]*(  nef  )
            z += [0]*nef
            
            z += [0.4]*nef
            
            # equalities for kinematic problem
            arowk = acol.copy()

            acolk = arow.copy()
        
            avalk = aval.copy()


            
            
            #acolk += [4]
        
          
            #arowk += [ max(arowk) +1 ]
            
            #avalk += [1]
        
            #avalk += [0]
            if False:
                
                for i in range(3*n):
                    index = arowk.index(i)
                    arowk.insert(index, i)
                    acolk.insert(index, 6*nb + nvcontact + i)
                    avalk.insert(index, -1)
                    if ((i+1) % 3 == 0):
                        index = arowk.index(i)
                        arowk.insert(index, i)
                        acolk.insert(  index, 6*nb + nvcontact + 3*n + int((i+1)/3) - 1 )
                        avalk.insert(index, 1)
         
            for i in range(3*nef):
                index = arowk.index(i)
                arowk.insert(index, i)
                acolk.insert(index, 6*nb + i)
                avalk.insert(index, -1)
                if ((i+1) % 3 == 0):
                    index = arowk.index(i)
                    arowk.insert(index, i)
                    acolk.insert(  index, 6*nb + 3*nef + int((i+1)/3) - 1 )
                    avalk.insert(index, 1)

                    index = arowk.index(i)
                    arowk.insert(index, i)
                    acolk.insert(  index, 6*nb + 3*nef + nef + int((i+1)/3) - 1 )
                    avalk.insert(index, -1)
    
           
            
      
            if False:
                num_lists = 6*nb + nvcontact + 3*n + n  # Adjust the number of inner lists as needed
                
            num_lists = 6*nb  + 3*nef + nef + nef*1  # Adjust the number of inner lists as needed
    
          
            asub_mos = [[] for _ in range(num_lists)]
            aval_mos  = [[] for _ in range(num_lists)]
            
            for i,j in enumerate(acolk):
                asub_mos[int(j)].append(arowk[i])
                aval_mos[int(j)].append(avalk[i])
            
            numvar = len(bkx)
            numcon = len(bkc)
            NUMANZ = 4 
    
         # Append 'numcon' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(numcon)
    
            #Append 'numvar' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(numvar)
    
            for j in range(numvar):
                # Set the linear term c_j in the objective.
                task.putcj(j, z[j])
                # Set the bounds on variable j
                # blx[j] <= x_j <= bux[j]
                task.putvarbound(j, bkx[j], blx[j], bux[j])
            
          
            for j in range(len(aval_mos)):
                
               
                # Input column j of A
                task.putacol(j,  asub_mos[j],aval_mos[j])            # Non-zero Values of column j.
            
            
           
        # Input the affine conic constraints
            # Create a matrix F such that F * x = [x(3),x(0),x(1),x(4),x(5),x(2)] 
            if False:
                task.appendafes(3*n)
            
                list_con = []
                list_coef = []
                friction_ang = 33
                for ind in range(n):
                    list_con += [ 6*nb + nvcontact + 3*ind + 2 , 6*nb + nvcontact + 3*ind,6*nb + nvcontact + 3*ind + 1]
                    list_coef += [-1/math.tan(math.radians(30)),1.0,1.0]
                
                task.putafefentrylist(range(3*n),                      # Rows
                                      list_con ,            # Columns 
                                      list_coef  )          #coefficients
        
                # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
                
                coh_vect = []
                for ind in range(n):
                    coh = [-0,0, 0]
                    coh_vect += coh
                    quadcone  = task.appendquadraticconedomain(3)
                    task.appendacc(quadcone,            [3*ind, 3*ind +1, 3*ind + 2],  coh)                    # None if there is no b for conic 
            
            task.appendafes(3*nef)
           
            list_con = []
            list_coef = []
            
            friction_ang = 33
            for ind in range(nef):
                list_con += [ 6*nb +  3*ind + 2 , 6*nb + 3*ind,6*nb +  3*ind + 1]
                list_coef += [-1/math.tan(math.radians(friction_ang)),1.0,1.0]
            
            task.putafefentrylist(range(3*nef),                      # Rows
                                  list_con ,            # Columns 
                                  list_coef  )          #coefficients
    
            # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
            
            coh_vect = []
            for ind in range(nef):
                coh = [-0.5*0,0, 0]
                coh_vect += coh
                quadcone  = task.appendquadraticconedomain(3)
                task.appendacc(quadcone,      [3*ind, 3*ind +1, 3*ind + 2],   coh)                    # None if there is no b for conic 
           
            
           
            
            for i in range(numcon):
                task.putconbound(i, bkc[i], blc[i], buc[i])
    
    
        # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)
    
            # Optimize the task
            task.optimize()
            #task.writedata("cqo1.ptf")
            # Print a summary containing information
            # about the solution for debugging purposes
            task.solutionsummary(mosek.streamtype.msg)
            prosta = task.getprosta(mosek.soltype.itr)
            solsta = task.getsolsta(mosek.soltype.itr)
    
            # Output a solution
            kin_sol = task.getxx(mosek.soltype.itr)
            print("kinsol",kin_sol)
            #ub = task.getxx(mosek.soltype.itr)[0:6*nb]
            #ub = task.getsuc(mosek.soltype.itr)[0:6*nb] # get dual variables for equilibrium but why upper ? not slb
            #skc, y = task.getsolution(mosek.soltype.bas)
            if solsta == mosek.solsta.optimal:
                print("Objective: %s" % kin_sol[-1:])
                #print("Optimal solution: %s" % dual(xx))
                # Get dual variables for the linear constraints
                
              
            elif solsta == mosek.solsta.dual_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif solsta == mosek.solsta.prim_infeas_cer:
                print("Primal or dual infeasibility.\n")
            elif mosek.solsta.unknown:
                print("Unknown solution status")
            else:
                print("Other solution status")
           # Clear or reset the task
            
            #del task  # Delete the task object

    #end kin

    
    b = matrix(b)
   
    q = [4 for ind in range(n)]
    h = spmatrix( [0], [0], [3*n], (1,4*n ) )
    h = np.zeros((4*n,1))
    h = matrix(h)
    
    dims = {'l': 0, 'q': q, 's': []}
 
    #sol = solvers.conelp(c, MCoul_matrix_sparse , h , dims, equilibrium_matrix_sparse, b)
    
    #print(sol["x"])
    if False:
        # ipopt part
        ind = [element + 1 for element in range(len(equilibrium_matrix))]
        
        
        equilibrium_matrix_row = [element + 1 for element in equilibrium_matrix_row] 
        
        A_ipopt = np.array([ind, equilibrium_matrix_row,equilibrium_matrix_col,equilibrium_matrix]).T
        
        ind = [element + 1 for element in range(len(b))]
        
        b_ipopt = [element for element in b]
        
        b_ipopt = np.array([ind, b_ipopt]).T
        
        c_ipopt = np.array([[1,3*n,1]])


        A_df = pd.DataFrame(  A_ipopt,  columns=["Aline","ALineInd","AColInd","AVal"], ).set_index("Aline")

        # Create a pandas.DataFrame with data for n_min, n_max
        b_df = pd.DataFrame( b_ipopt , columns=["bline", "BVal"],).set_index("bline")

        c_df = pd.DataFrame(c_ipopt ,columns=["cline", "CInd", "CVal"],).set_index("cline")
        
        
        
        
        


        ind = [element + 1 for element in range(len(contact_eq_element_x))]
        contact_eq_row_x = [element + 1 for element in contact_eq_row_x] 
        Cx_ipopt = np.array([ind, contact_eq_row_x,contact_eq_col_x,contact_eq_element_x]).T

        
        ind = [element + 1 for element in range(len(contact_eq_element_y))]
        contact_eq_row_y = [element + 1 for element in contact_eq_row_y] 
        Cy_ipopt = np.array([ind, contact_eq_row_y,contact_eq_col_y,contact_eq_element_y]).T

        ind = [element + 1 for element in range(len(contact_eq_element_z))]
        contact_eq_row_z = [element + 1 for element in contact_eq_row_z] 
        Cz_ipopt = np.array([ind, contact_eq_row_z,contact_eq_col_z,contact_eq_element_z]).T
        
        Cx_df = pd.DataFrame(Cx_ipopt,columns=["Cxline"," CxLineInd","CxColInd","CxVal"],).set_index("Cxline")
        
        Cy_df = pd.DataFrame( Cy_ipopt,columns=["Cyline"," CyLineInd","CyColInd","CyVal"],).set_index("Cyline")
        
        Cz_df = pd.DataFrame(Cz_ipopt,columns=["Czline"," CzLineInd","CzColInd","CzVal"],).set_index("Czline")
    
    #ub = ipopt_data(A_df,b_df,c_df,Cx_df , Cy_df, Cz_df)
    #ub =[]
    x = []
    
    #stat_rep(all_unique_points,noc_triangles_coor,noc_triangles,all_triangles_coor,all_triangles,xx)
    #stat_rep1(all_unique_points,noc_triangles_coor,noc_triangles,all_triangles_coor,all_triangles,xx)
    all_unique_points = output[3]
    all_triangles = output[4]
    triangles_num = max(output[5]) 

    noc_triangles_coor = output[6]
    noc_triangles = output[7]
    Ntriangles_num = max(output[8]) 

    stat_parav(all_unique_points,noc_triangles_coor,noc_triangles,triangles_num ,all_triangles,Ntriangles_num,xx[0:3*n+1])

    return ub

def face_mesh(nodes):

    # Initialize gmsh:
    gmsh.initialize()
    
    # Def&ine nodes in Gmsh using the NumPy array:
    lc = 100
    node_tags = [gmsh.model.geo.add_point(x, y, z, lc) for x, y, z in nodes]
    
    # Define lines forming a square connected in a loop using a for loop:
    lines = []
    for i in range(len(node_tags)):
        line = gmsh.model.geo.add_line(node_tags[i], node_tags[(i + 1) % len(node_tags)])
        lines.append(line)
    
    # Define a loop connecting the lines to form a face:
    loop = gmsh.model.geo.add_curve_loop(lines)
    
    # Define a surface using the loop:
    face = gmsh.model.geo.add_plane_surface([loop])
    
    # Create the relevant Gmsh data structures from the Gmsh model:
    gmsh.model.geo.synchronize()
    
    # Generate mesh:
    gmsh.model.mesh.generate()
    
    # Write mesh data:
    #gmsh.write("C:\\Users\\mmoussa\\Desktop\\rhino_test\\square_face.msh")
    elements = gmsh.model.mesh.getElements(dim=2, tag=-1)[2][0]
    elements = [ int(element -1) for element in elements]
    points = gmsh.model.mesh.get_nodes()[1]  # Coordinates of the mesh nodes
    
    nodes = [  list(points[ind:ind +3])   for ind in range(0,len(points), 3)]

    triangles = [  list(elements [ind:ind +3])   for ind in range(0,len(elements), 3)]
    gmsh.finalize()
    return nodes,triangles


def face_mesh_read():
    
    
    # Read the mesh file
    mesh = meshio.read("C:\\Users\\mmoussa\\Desktop\\rhino_test\\square_face.msh")
    
    # Access mesh information
    points = mesh.points  # Coordinates of the mesh nodes
    cells = mesh.cells  # Connectivity information for mesh elements (e.g., triangles, tetrahedra)
    point_data = mesh.point_data  # Data associated with each mesh node
    cell_data = mesh.cell_data  # Data associated with each mesh element
    
    
    # Example: Print connectivity of the first 5 triangles
    triangle_cells = None
    for cell_block in cells:
        if "triangle" in cell_block.type.lower():
            triangle_cells = cell_block.data
            break
    
    if triangle_cells is not None:
        #print(f"Connectivity of the first 5 triangles:")
        #print(triangle_cells[:5])
    
        # Get coordinates of each point in the first 5 triangles
        for triangle in triangle_cells:
            triangle_coordinates = points[triangle]
           
    else:
        print("No triangle cells found in the mesh.")
    return points, triangle_cells


ub  = problem_building_up(blocks_centroid,points,faces,blocks,local_ref,blocks_att)
block_disp_plot(blocks_centroid,blocks,faces,points,ub)


print(faces_rep)
print(plans_rep)

