
import rhinoscriptsyntax as rs
import gmsh
import numpy as np
import importlib
import plotly.graph_objects as go
#import meshio
from cvxopt import matrix, solvers, spmatrix
#from scipy.sparse import coo_matrix
import sys
import time
from sklearn.cluster import KMeans
#import pandas as pd
import os
import vtk
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad, VtkVertex, VtkPolyLine, VtkPolygon
import mesh
import mosek
import math
from scipy.linalg import block_diag
from scipy.spatial import cKDTree
import os
import json
import locale
from urllib.request import urlopen
from datetime import datetime



def blocks_face_neighboring(blocks_centroid,face_centroid, k_neighbors):
    kdtree = cKDTree(blocks_centroid[:,-3:])
 
    distances, neighbors = kdtree.query(face_centroid, k=k_neighbors + 1)
  
    return distances, neighbors


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

def check_outward_normal(block_centroid,face_ind, face_centroid,face_normal,local_ref):

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


def stat_parav(data, calc_type, stat_sol):

    all_unique_points = data["all_unique_points"]
    all_triangles = data["all_triangles"]

 
    
    FILE_PATH = data["rhino_file_path"]  + "\\" + calc_type + "_stat"

    print("Running unstructured...")

    
    
    
    if calc_type == "ub":
        area = data["Coutput"][12]
  
        normal_load = [ stat_sol[ind+2]/area[int(ind/3)]   for ind in range(0,len(stat_sol[0:-1:]),3) ]
        shear_loadx = [ stat_sol[ind]/area[int(ind/3)]    for ind in range(0,len(stat_sol[0:-1:]),3) ]
        shear_loady = [ stat_sol[ind+1]/area[int(ind/3)]    for ind in range(0,len(stat_sol[0:-1:]),3) ]
     

    
    else:
        normal_load = [ stat_sol[ind+2]   for ind in range(0,len(stat_sol[0:-1:]),3) ]
        shear_loadx = [ stat_sol[ind]   for ind in range(0,len(stat_sol[0:-1:]),3) ]
        shear_loady = [ stat_sol[ind+1]   for ind in range(0,len(stat_sol[0:-1:]),3) ]
    
  
    all_unique_points = np.array(all_unique_points)

    # Define vertices
    x = [float(element) for element in all_unique_points[:, 0]]
    y = [float(element) for element in all_unique_points[:, 1]]
    z = [float(element) for element in all_unique_points[:, 2]]


    x = np.array(x, dtype = np.float64)
    y = np.array(y,dtype = np.float64)
    z = np.array(z,dtype = np.float64)
    

    # Define connectivity or vertices that belong to each element
    conn = []
    # Define offset of the last vertex of each element
    offset = []
    # Define cell types
    ctype = []
   
    node_ind = 0
 
    for face_triangles in all_triangles:
      
        for triangle in face_triangles:
            for node in triangle:
                conn.append(int(node-1))
                node_ind += 1


            offset.append(node_ind )

            if len(triangle) == 3:
                ctype.append(VtkTriangle.tid)
            else:
                
                ctype.append(VtkQuad.tid)
          

    conn =  np.array(conn)
    offset =  np.array(offset)
    ctype =  np.array(ctype)

    #cd = np.random.rand(2)
    cellData = {"pressure": []}

    
    # Define displacement components
    normal_stress = np.array(normal_load)
    shear_stress1 = np.array(shear_loadx)
    shear_stress2 = np.array(shear_loady)

    # Combine displacement components into separate arrays
    pointData = {"normal_stress":normal_stress, "shear_stress1": shear_stress1, "shear_stress2": shear_stress2  }
 
    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype, pointData=pointData)
    
def blocks_mesh_parav(data):

    noc_triangles_coor = data["noc_triangles_coor"]
    noc_triangles = data["noc_triangles"]
    
    FILE_PATH = data["rhino_file_path"] + "\\blocks_mesh"
   
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
    conn = []
    # Define offset of the last vertex of each element
    offset = []
    # Define cell types
    ctype = []

    tri_ind = 0
    node_ind = 0
 
    for face_triangles in noc_triangles:
      
        for triangle in face_triangles:
            for node in triangle:
                conn.append(int(node-1))
                node_ind += 1


            offset.append(node_ind )

            if len(triangle) == 3:
                ctype.append(VtkTriangle.tid)
            else:

                ctype.append(VtkQuad.tid)
          
  

    conn =  np.array(conn)
    offset =  np.array(offset)
    ctype =  np.array(ctype)

    pointData = {'d': np.ones(x.shape[0])}

    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype, pointData=pointData)


def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()



def kin_parav(data,calc_type,ub):

    blocks_centroid = data["blocks_centroid"]
    blocks = data["blocks"]
    faces = data["faces"]
    points = data["points"]


    ub = [ element for element in ub ]

    new_points = []
    initial_points = np.empty((0,3))
    points_disp = np.empty((0,3))

    FILE_PATH =  data["rhino_file_path"] + "\\" + calc_type +  "_disp"
  
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
             
                initial_points = np.append(initial_points, np.array([element]), axis = 0)
                

                conn.extend( [initial_points.shape[0] - 1 ])
                
                point_new_ref = [element[element1] - reference_point[element1] for element1 in range(3) ]
               
                
                new_point = transform_point(point_new_ref , displacement_vector, rotations)
                
               
                new_point = [new_point[element] + reference_point[element] for element in range(3) ]
                
                all_new.append(new_point)

                points_disp = np.append(points_disp, np.array( [new_point[:] - element[:] ] ), axis = 0)
             
                
            ctype.extend ([VtkQuad.tid])
            offset.extend( [initial_points.shape[0] ])
            
            new_points.append(all_new[0]) 
            
            all_new = np.array(all_new)   
            
            

            
    
    points = new_points
   
    point = np.array([[float(element) for element in row] for row in points])
    
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


##first part
# Reload the RhinoScriptSyntax module


def read_data_from_rhino(data):
  

    read_data_start = time.time()
    importlib.reload(rs)

    lc = data["lc"]
    faces_rep = []
    pts_coor = np.empty((0,3))

    block_ind = 1  # initialize blocks counter variable 
    face_ind = 1    # initialize faces counter variable 

    blocks = np.empty((0, 7), dtype = np.float64) # blocks array shape( block ind | face_unique_ind | block_name | face_non_unique_ind | face_type | cohesion_ind | 11 if the face is principal 22 if not  )
    faces = np.empty((0,3))   # shape(face_unique_ind | face_non_unique_ind | node ind )
    blocks_centroid = np.empty((0,4), dtype = np.float64 ) # shape(block ind |x | y |z)
    blocks_volume = np.empty((0,2), dtype = np.float64)  # shape( block ind |  volume )
    blocks_att = dict()    # variable to store blocks attributes          

    unsorted_faces_rep = [] # variable to store the faces representation without sorting the nodes indices  ( ex: a face whose the boundary passes through the nodes having indcides 1,3,7,5  have the representation 1#3#7#5)
    plans_data = np.empty((0,4)) # local variable to store the external surfaces data shape( surface ind |  coordinates for a single point that constitutes the boundary of the surface  )
    plan_ind = 1  # initialize surfaces counter variable 
    f_points = [] #local variable to store faces points
    keys_of_interest = ["fx", "fy", "fz", "mx", "my", "mz"] # keys for blocks loading applied on their centroids
    local_ref = np.empty((0,10),  dtype = np.float64) # shape(face_unique_ind | normal local axis(x,y,z) |   local_1 axis(x,y,z)  | local_2 axis(x,y,z) )
    
    blocks_brep = [] 
    bars = []
    supports_type = [] 

    # store some contacts data (A face belonging to more than one block is called a contact)
    contacts_ind = [] 
    contacts_nodes = []
    contacts_FE = [0]

    # store some faces data (A face belonging to at least one block is called face)

    faces_ind = [] 
    faces_nodes = []
    faces_FE = [0]

    #store some data for external surfaces attributes

    phi_type = [] # contains the friction agnle value attributated to an external surface. If no attribute is found -1 is stored
    coh_type = [] # contains the cohesion value attributated to an external surface. If no attribute is found -1 is stored
    lc_type = []  # contains the friction agnle value attributated to an external surface. If no attribute is found -1 is stored
    fc_type = []

    objs = rs.AllObjects()

    loading_keys_interest = ["px_l", "py_l", "pz_l", "pn_l", "px_d", "py_d", "pz_d", "pn_d"]

    linear_loading_key_interest = ["wx_l", "wy_l", "wz_l", "wn_l", "wx_d", "wy_d", "wz_d", "wn_d"]
    
    ponc_loading_key_interest = ["fx_l", "fy_l", "fz_l","mx_l", "my_l", "mz_l",               "fx_d", "fy_d", "fz_d" , "mx_d", "my_d", "mz_d"] 

    loaded_faces_values = []
    loaded_faces_adress = []
    blocks_adress = []
    loaded_blocks = []
    torseur_blocks_live = []
    torseur_blocks_dead = []
    loaded_poly_adress = []
    loaded_poly_values = []

    loaded_points_values = []
    loaded_points_adress = []

    for obj in objs:
        obj_type = rs.ObjectType(obj)

        if obj_type == 1: # chekc if the object is a point
            point_att = [element.lower() for element in rs.GetUserText(obj)]
            if set(point_att) & set(ponc_loading_key_interest ):
                loaded_points_values.append([0]*12)
                loaded_points_adress.append(obj)


                loading_att = set(point_att) & set(ponc_loading_key_interest )
                for att in loading_att:
                    
                    att_pos = ponc_loading_key_interest.index( att.lower() )
                   
                    if rs.GetUserText(obj, att) != " ":
                        loaded_points_values[-1][att_pos] += float( rs.GetUserText(obj, att) )


        
        if obj_type == 4: # Check if the object is a curve
            poly_att = [element.lower() for element in rs.GetUserText(obj)]
            if "fy" in poly_att and "as" in poly_att:
                bars.append(obj) 
            elif set(poly_att) & set(linear_loading_key_interest ) :
                loaded_poly_values.append([0]*8)
                loaded_poly_adress.append(obj)
                
                loading_att = set(poly_att) & set(linear_loading_key_interest )
                for att in loading_att:
                    
                    att_pos = linear_loading_key_interest.index( att.lower() )
                   
                    if rs.GetUserText(obj, att) != " ":
                        loaded_poly_values[-1][att_pos] += float( rs.GetUserText(obj, att) )



        if obj_type == 8: # Check if the object is a surface
        
            
            boundary_curve_id = rs.DuplicateSurfaceBorder(obj)
            # Get points along the boundary curve
            points = rs.CurvePoints(boundary_curve_id)

            rs.DeleteObjects(boundary_curve_id)
            points = points[0:-1:]
            
            for i, point in enumerate(points):
            
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

            if 'fc' in att_lcase:
                index = att_lcase.index('fc')
                fc_type.append(float(rs.GetUserText(obj, rs.GetUserText(obj)[index])))
                
            else:
                fc_type.append(-1)

            
            if 'lc' in att_lcase:
                index = att_lcase.index('lc')
                lc_type.append(float(rs.GetUserText(obj, rs.GetUserText(obj)[index])))
                
            else:
                lc_type.append(-1)

        


            loading_att = set(loading_keys_interest) & set(att_lcase) # surface loading
            
            if loading_att:
                loaded_faces_values.append([0]*8)
                loaded_faces_adress.append(obj)

                for att in loading_att:
                    
                    att_pos = loading_keys_interest.index( att.lower() )
                  
                    if rs.GetUserText(obj, att) != " ":
                        loaded_faces_values[-1][att_pos] += float( rs.GetUserText(obj, att) )
                        



                
    
        print("loadeddd face", loaded_faces_values)
        if obj_type in [16, 1073741824]:  # Check if the object is a solid
            blocks_adress.append(obj)

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
                    rs.DeleteObjects(boundary_curve_id)
                    
                    points = points[0:-1:]
                    

                    
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


                        local_ref = check_outward_normal(b_centroid,face_ind,face_centroid,face_normal,local_ref)
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

                    #  delete the exploded face
                    rs.DeleteObject(face)
                    rs.DeleteObject(boundary_curve_id)
                    
            
        
            b_centroid = rs.SurfaceVolumeCentroid(obj)[0]
            b_volume = rs.SurfaceVolume(obj)[0]

            blocks_volume = np.append(blocks_volume, np.array([[block_ind, b_volume ]]), axis = 0)

            blocks_centroid = np.append(blocks_centroid, np.array([[block_ind, b_centroid[0], b_centroid[1], b_centroid[2] ]]), axis = 0)
            
            atb_keys = rs.GetUserText(obj)
            
            if atb_keys:
                
                b_att = []

                for atb_key in atb_keys:
                    if atb_key in keys_of_interest:
                        key_ind = keys_of_interest.index(atb_key)
                        

                        b_att.extend([key_ind +1, rs.GetUserText(obj, atb_key) ])
                
                blocks_att[block_ind] = b_att

    
            
            block_ind += 1


    nb = block_ind - 1

    # Get the number of rows in the array
    num_rows = pts_coor.shape[0]

    # Create an array of 1-based indices
    indices_column = (np.arange(num_rows) + 1).reshape(-1, 1)

    # Add the indices column to the original array
    points = np.hstack((indices_column, pts_coor),dtype = np.float64)

    data["faces_number"] =  len(faces_rep)
  
    plans_rep = []
    supports_ind = []

    if plans_data.shape[0] != 0:
        for ind in range(plan_ind-1):
        
            plan_elements = plans_data[ np.where(plans_data[:,0] == ind +1 )[0] ]
    
            plan_pts = []
            for  element in plan_elements:
                pt_coor = element[-3:]
                pt_ind = np.where( np.all(points[:,-3:] == pt_coor, axis=1)  )[0]
                
                plan_pts.extend(pt_ind)
            
            sorted_plan_pts = sorted(set(plan_pts))
            sorted_plan_pts = [element + 1 for element in sorted_plan_pts]

            plan_rep = '#'.join(map(str, sorted_plan_pts ))
            
            plans_rep.append(plan_rep)
            
        
    for ind in contacts_ind:    blocks[np.where(blocks[:,1] == ind)[0], 4] = 1


    supports_rep = set(plans_rep) &  set(faces_rep)

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



    contacts_ind_sorted = sorted(contacts_ind)

    for face_pos, face_adress in enumerate(loaded_faces_adress):

        face_centroid = rs.SurfaceAreaCentroid(face_adress)[0]
        face_centroid = [face_centroid[0],face_centroid[1],  face_centroid[2] ]
        distances , neighboors = blocks_face_neighboring( blocks_centroid[:, -3:] , face_centroid, 40  )
        print(distances)

        real_distances_pos = np.where(distances != np.inf)[0] 

        print(real_distances_pos)

        for index in real_distances_pos:
            
            block_adress  = blocks_adress[neighboors[index]]
            inter = rs.IntersectBreps(face_adress, block_adress)
            
            
            if inter:
                if rs.IsCurveClosed(inter):
                    
                    print("enteredddd")
                    area = rs.Area(inter)
                    centroid = rs.CurveAreaCentroid(inter)[0]
                    centroid = [ centroid[0],  centroid[1], centroid[2] ]

                    crv_points = rs.CurvePoints(inter)  

                    print(crv_points)    
                    pt1 =  crv_points[0]
                    pt2 =  crv_points[1]
                    pt3 =  crv_points[-2]
                    print(pt1)
                    print(pt2)
                    print(pt3)
                    vect1 = [pt1[i]-pt2[i] for i in range(3)]
                    vect2 = [pt3[i]-pt2[i] for i in range(3)]

                    normal = np.cross(vect1, vect2)
                    normal /= np.linalg.norm(normal)
                    print(vect1, vect2 , normal)
                    test_vector = [blocks_centroid[index][i+1] - centroid[i]  for i in range(3)  ]

                    if np.dot(test_vector, normal ) >0:
                        normal = [-element for element in normal]

                    
                    Dead_load = [area*loaded_faces_values[face_pos][4], area*loaded_faces_values[face_pos][5], area*loaded_faces_values[face_pos][6] ]
                    print(Dead_load)

                    Dead_load = [Dead_load[i] + normal[i]*area*loaded_faces_values[face_pos][7] for i in range(3) ]
                        
                    print(Dead_load)
                    Live_load = [area*loaded_faces_values[face_pos][0], area*loaded_faces_values[face_pos][1], area*loaded_faces_values[face_pos][2] ]
                    Live_load = [Live_load[i] + normal[i]*area*loaded_faces_values[face_pos][3] for i in range(3) ]
                    
                    

                    moment_dead = calculate_moment(blocks_centroid[index][-3:], Dead_load  , centroid  )
                    print(moment_dead)

                    moment_live = calculate_moment(blocks_centroid[index][-3:], Live_load  , centroid  )

                    if neighboors[index] not in loaded_blocks:
                        loaded_blocks.append(neighboors[index])
                        torseur_blocks_dead.append( [*Dead_load, *moment_dead ])
                        torseur_blocks_live.append( [*Live_load, *moment_live ])
                    else:
                        block_pos = loaded_blocks.index(neighboors[index])#not sure

                        torseur_blocks_dead[block_pos] = [  element + torseur_blocks_dead[block_pos][ind]  for ind, element in enumerate([*Dead_load, *moment_dead ])  ]
                        torseur_blocks_live[block_pos] = [  element + torseur_blocks_live[block_pos][ind]  for ind, element in enumerate([*Live_load, *moment_live ]) ]


                    
                    
                    rs.DeleteObjects(inter)
    
       
    for point_pos, point_adress in enumerate(loaded_points_adress):
        
        loading_point_points = np.empty((0,3)) #mesh for poly loading and
        points_ind = []

        points_torseur_blocks_dead = []

        points_torseur_blocks_live = []

        points_loaded_blocks = []

    

        
        point_coor = [ rs.PointCoordinates(point_adress)[0], rs.PointCoordinates(point_adress)[1], rs.PointCoordinates(point_adress)[2] ]
        

        distances , neighboors = blocks_face_neighboring( blocks_centroid[:, -3:] , point_coor , 40  )

        real_distances_pos = np.where(distances != np.inf)[0] 

       
        

        for index in real_distances_pos:

            block_adress  = blocks_adress[neighboors[index]]

            faces_adress = rs.ExplodePolysurfaces(block_adress)
            

            for face in faces_adress :

                if rs.IsPointOnSurface(face,point_adress):
                   
                    
                    rs.DeleteObject(face)
                    search_pt1 = np.where( np.all(loading_point_points == point_coor , axis = 1) )[0]
                    if len(search_pt1) != 0:
                        points_ind.append( int(search_pt1[0]) + 1)
                    else:
                        loading_point_points = np.append(loading_point_points , np.array([point_coor]) ,axis = 0 )
                        points_ind.append( loading_point_points.shape[0] )
                
                    Dead_load = [ loaded_points_values[point_pos][6], loaded_points_values[point_pos][7], loaded_points_values[point_pos][8] ]
                        
                    Live_load = [ loaded_points_values[point_pos][0], loaded_points_values[point_pos][1], loaded_points_values[point_pos][2] ]
            
                    moment_dead = calculate_moment(blocks_centroid[neighboors[index]][-3:], Dead_load  , point_coor  )
                    moment_live = calculate_moment(blocks_centroid[neighboors[index]][-3:], Live_load  , point_coor  )

                    moment_dead = [ element1 + element2 for element1, element2 in zip(moment_dead, loaded_points_values[point_pos][9:12] ) ]   

                    moment_live = [ element1 + element2 for element1, element2 in zip(moment_live, loaded_points_values[point_pos][3:6] ) ]  


                    points_torseur_blocks_dead.append( [*Dead_load, *moment_dead ])
                    points_torseur_blocks_live.append( [*Live_load, *moment_live ])

                    points_loaded_blocks.append(neighboors[index])
                    
                    
                    
                    break
                else:
                 
                    rs.DeleteObject(face)

            rs.DeleteObjects(faces_adress)
                
      

        for i, point_loaded_block in enumerate(points_loaded_blocks):
            point_occ = points_ind.count( points_ind[i] )
     

            if point_loaded_block not in loaded_blocks:
                loaded_blocks.append(point_loaded_block)
                torseur_blocks_dead.append( [ element/point_occ  for element in points_torseur_blocks_dead[i]   ])
                torseur_blocks_live.append(  [ element/point_occ  for element in points_torseur_blocks_live[i]   ])
            
            else:

                block_pos = loaded_blocks.index(point_loaded_block)

                torseur_blocks_dead[block_pos] = [  element/point_occ + torseur_blocks_dead[block_pos][ind]  for ind, element in enumerate( points_torseur_blocks_dead[i] )  ]
                torseur_blocks_live[block_pos] = [  element/point_occ + torseur_blocks_live[block_pos][ind]  for ind, element in enumerate( points_torseur_blocks_live[i] ) ]



    for poly_pos, poly_adress in enumerate(loaded_poly_adress):
        
        loading_poly_points = np.empty((0,3)) #mesh for poly loading and
        polys_rep = []

        poly_torseur_blocks_dead = []

        poly_torseur_blocks_live = []

        poly_loaded_blocks = []

        points_poly = rs.CurvePoints(poly_adress)
             
        pt1 =  [round(points_poly[0][i],5) for i in range(3) ] 
        

        distances , neighboors = blocks_face_neighboring( blocks_centroid[:, -3:] , pt1, 40  )

        real_distances_pos = np.where(distances != np.inf)[0] 

        for index in real_distances_pos:

            
            block_adress  = blocks_adress[neighboors[index]]
            
            inter = rs.CurveBrepIntersect(poly_adress , block_adress   )
            
            
            if inter:
                for line in inter[0]:
            
                    points_poly = rs.CurvePoints(line)
             
                    pt1 =  [round(points_poly[0][i],5) for i in range(3) ] 
                    pt2 =  [round(points_poly[-1][i],5) for i in range(3) ]
                    
                    search_pt1 = np.where( np.all(loading_poly_points == pt1 , axis = 1) )[0]
                    if len(search_pt1) != 0:
                        pt1_pos = int(search_pt1[0]) + 1
                    else:
                        loading_poly_points = np.append(loading_poly_points , np.array([pt1]) ,axis = 0 )
                        pt1_pos = loading_poly_points.shape[0]

                    search_pt2 = np.where( np.all(loading_poly_points == pt2 , axis = 1) )[0]
                    if len(search_pt2) != 0:
                        pt2_pos = int(search_pt2[0]) + 1
                    else:
                        loading_poly_points = np.append(loading_poly_points , np.array([pt2]) ,axis = 0 )
                        pt2_pos = loading_poly_points.shape[0]

                    poly_rep = '#'.join(map(str, [min(pt1_pos, pt2_pos) ,  max(pt1_pos, pt2_pos) ] )) 

                    polys_rep.append(poly_rep)


                    line_centroid = [ 0.5*pt1[i] + 0.5*pt2[i] for i in range(3)]
                    line_length = rs.CurveLength(line)
            
                    Dead_load = [line_length*loaded_poly_values[poly_pos][4], line_length*loaded_poly_values[poly_pos][5], line_length*loaded_poly_values[poly_pos][6] ]
                    
                    Live_load = [line_length*loaded_poly_values[poly_pos][0], line_length*loaded_poly_values[poly_pos][1], line_length*loaded_poly_values[poly_pos][2] ]
          
                    moment_dead = calculate_moment(blocks_centroid[neighboors[index]][-3:], Dead_load  , line_centroid  )
                    moment_live = calculate_moment(blocks_centroid[neighboors[index]][-3:], Live_load  , line_centroid  )


                    poly_torseur_blocks_dead.append( [*Dead_load, *moment_dead ])
                    poly_torseur_blocks_live.append( [*Live_load, *moment_live ])

                    poly_loaded_blocks.append(neighboors[index])

                    rs.DeleteObjects(line)  
        
        for i, poly_loaded_block in enumerate(poly_loaded_blocks):
            poly_occ = polys_rep.count( polys_rep[i] )
      

            if poly_loaded_block not in loaded_blocks:
                loaded_blocks.append(poly_loaded_block)
                torseur_blocks_dead.append( [ element/poly_occ  for element in poly_torseur_blocks_dead[i]   ])
                torseur_blocks_live.append(  [ element/poly_occ  for element in poly_torseur_blocks_live[i]   ])
            
            else:

                block_pos = loaded_blocks.index(poly_loaded_block)

                torseur_blocks_dead[block_pos] = [  element/poly_occ + torseur_blocks_dead[block_pos][ind]  for ind, element in enumerate( poly_torseur_blocks_dead[i] )  ]
                torseur_blocks_live[block_pos] = [  element/poly_occ + torseur_blocks_live[block_pos][ind]  for ind, element in enumerate( poly_torseur_blocks_live[i] ) ]
        #jkjk



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

    ##Calculate steel rebars load and moment effect

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


    

    lc_faces_val = lc*np.ones((len(contacts_ind)), dtype = np.float64)


    for ind, lc_face  in enumerate(lc_type):
        if lc_face != -1:
        
            face_ind = faces_rep.index( plans_rep[ind] )
            contact_ind = contacts_ind.index(face_ind + 1)
        
            lc_faces_val[contact_ind] = float(lc_face)


    bars_data = np.unique(bars_data,axis = 0)

    contacts_ind = np.array(contacts_ind, dtype = np.int32)
    contacts_nodes = np.array(contacts_nodes, dtype = np.int32)
    contacts_FE = np.array(contacts_FE, dtype = np.int32)

    

    contacts_ind = np.array(contacts_ind, dtype = np.int32)
    faces_nodes = np.array(faces_nodes, dtype = np.int32)
    faces_FE = np.array(faces_FE, dtype = np.int32)

    data["points"] = points
    data["blocks"] = blocks
    data["faces"] = faces
    data["blocks_volume"] = blocks_volume
    data["blocks_centroid"] = blocks_centroid
    data["local_ref"] = local_ref
    data["lc_faces_val"] = lc_faces_val
    data["contacts_ind"] = contacts_ind
    data["faces_FE"] = faces_FE
    data["faces_nodes"] = faces_nodes

    data["faces_rep"] = faces_rep
    data["plans_rep"] = plans_rep

    data["supports_pos"] = supports_pos
    data["coh_type"] = coh_type
    data["phi_type"] = phi_type
    data["fc_type"] = fc_type
    data["blocks_att"] = blocks_att
    data["bpts_inter"] = bpts_inter
    data["bars_data"] = bars_data

    data["contacts_ind_sorted"] = contacts_ind_sorted

    data["nb"] = nb
    data["loaded_blocks"] = loaded_blocks
    data["torseur_blocks_dead"] = torseur_blocks_dead
    data["torseur_blocks_live"] = torseur_blocks_live

    read_data_end = time.time()

    calc_time = dict()
    data["calc_time"] = {"read_rhino_data_time" : read_data_end - read_data_start }

    
    
    

    return data
    


#### second part 
def generate_mesh_and_construct_matrices(data):  
    
    contacts_ind = data["contacts_ind"]
    blocks = data["blocks"]
    points = data["points"]
    faces_FE = data["faces_FE"]
    faces_nodes = data["faces_nodes"]
    blocks_centroid = data["blocks_centroid"]
    local_ref = data["local_ref"]
    lc_faces_val = data["lc_faces_val"]

    coh_type = data["coh_type"]
    supports_pos = data["supports_pos"]
    faces_rep = data["faces_rep"] 
    plans_rep = data["plans_rep"]
    contacts_ind_sorted = data["contacts_ind_sorted"]
    phi_type = data["phi_type"]
    
    #fc_type = data["fc_type"]

    tri_mesh = data["tri_mesh"]

    lc_ncn = data["lc_ncn"]

    # Call C++ function 

    points_coor = np.array(points[:,1:4] , dtype = np.float64)
   

    output = mesh.generate_mesh_and_construct_matrices( contacts_ind,  blocks , points_coor,  faces_FE, faces_nodes, blocks_centroid , \
                        local_ref, np.array(lc_faces_val, dtype = np.float64), lc_ncn, tri_mesh )

    
    data["calc_time"][" "] = " "
    data["calc_time"]["faces_mesh"] = output[-1:][0][0]
    data["calc_time"]["read_mesh"] = output[-1:][0][1]
    data["calc_time"]["Cpp_matrices_construction"] = output[-1:][0][2]

   
    supports_nodes = set()


    for pos in supports_pos:
        for tri in output[7][pos]:
                supports_nodes.update({ element for element in tri } )

    supports_nodes = list(supports_nodes)

    all_coh_nodes = []
    all_coh_val = []


    for ind, coh in enumerate(coh_type):
        if coh != -1:
            faces_pos = faces_rep.index( plans_rep[ind] )
            contact_pos = contacts_ind_sorted.index( faces_pos + 1 )

            contact_nodes = []
            for tri in output[7][contact_pos]:
                contact_nodes += [element for element in tri if element not in contact_nodes ]
            
            all_coh_nodes += contact_nodes 
            all_coh_val += [coh for i in range(len(contact_nodes))] 


    all_coh_nodes = list(all_coh_nodes)


    all_phi_nodes = []
    all_phi_val = []
 
    for ind, phi in enumerate(phi_type):
        if phi != -1:
            
            faces_pos = faces_rep.index( plans_rep[ind] )
            contact_pos = contacts_ind_sorted.index( faces_pos + 1 )

            contact_nodes = []
            for tri in output[7][contact_pos]:
                contact_nodes += [element for element in tri if element not in contact_nodes ]
            
            all_phi_nodes += contact_nodes
            all_phi_val += [phi for i in range(len(contact_nodes))] 

    all_fc_nodes = []
    all_fc_val = []

    fc_type = data["fc_type"]
    for ind, fc in enumerate(fc_type):
        if fc != -1:
            
            faces_pos = faces_rep.index( plans_rep[ind] )
            contact_pos = contacts_ind_sorted.index( faces_pos + 1 )

            contact_nodes = []
            for tri in output[7][contact_pos]:
                
                contact_nodes += [element for element in tri if element not in contact_nodes ]
            
            all_fc_nodes += contact_nodes
            all_fc_val += [fc for i in range(len(contact_nodes))] 

    


    all_fc_nodes = list(all_fc_nodes)
    all_phi_nodes = list(all_phi_nodes)

    data["supports_nodes"] = supports_nodes
    
    data["all_phi_nodes"] = all_phi_nodes
    data["all_fc_nodes"] = all_fc_nodes
    data["all_coh_nodes"] = all_coh_nodes

    data["all_phi_val"] = all_phi_val
    data["all_fc_val"] = all_fc_val
    data["all_coh_val"] = all_coh_val
    
    data["Coutput"] = output
 
    return data


### third part



# Define a stream printer to capture the solver output
def stream_printer(text,log_path):
    with open(log_path, "a") as f:
        f.write(text)


# New function to add a title/header to the log file
def add_solution_title(title,log_path):
    with open(log_path, "a") as f:
        f.write(f"\n{title}\n")  # Write the title with newlines for separation
        f.flush()  # Explicitly flush to ensure it's written immediately

def log_objective_value(task,log_path):
    obj_value = task.getprimalobj(mosek.soltype.itr)  # Assuming you're interested in the interior solution
    with open(log_path, "a") as f:
        f.write(f"Objective Value: {obj_value}\n")
        f.flush()

def solve_problem_with_solvers(data):
    loaded_blocks = data["loaded_blocks"]  
    torseur_blocks_dead= data["torseur_blocks_dead"] 
    torseur_blocks_live = data["torseur_blocks_live"] 

    fig = go.Figure()

    
    nb = data["nb"]

    output = data["Coutput"]  
    
    equilibrium_matrix_row = output[1]
    equilibrium_matrix_col = output[2]
    equilibrium_matrix = output[0]

    equilibrium_matrix_load_row = output[4]
    equilibrium_matrix_load_col = output[5]
    equilibrium_matrix_load = output[3]

    area = output[12]

 
    n = int((max(equilibrium_matrix_col)+1)/3)
    
    
    equilibrium_matrix_row = [int(element) for element in equilibrium_matrix_row]

    equilibrium_matrix_col = [int(element) for element in equilibrium_matrix_col]


    equilibrium_matrix_load_row = [int(element) for element in equilibrium_matrix_load_row]

    equilibrium_matrix_load_col = [int(element) for element in equilibrium_matrix_load_col]


    live_load = [0]*6*nb
    
    blocks_att = data["blocks_att"]
    bpts_inter = data["bpts_inter"] 

    #### blocks attributes
 

    for index, block_ind in enumerate(loaded_blocks):
        block_att = [ ]
        live = [0]*6
       
        

        if block_ind + 1 in blocks_att : # are both floats
                for att in range(0, len(blocks_att[block_ind + 1]), 2):
                    live[ blocks_att[block_ind + 1][att] - 1 ] += float(blocks_att[block_ind + 1][att + 1])
        
        
        live = [live[i] + torseur_blocks_live[index][i] for i in range(6)]

        for ind, element in enumerate(live):
            if abs(element) > 1e-10:
                block_att += [ind+1, element]
        
        blocks_att[block_ind + 1] = block_att

    for block_ind , block_attributes in blocks_att.items():
        for ind in range(0,len(block_attributes),2):
            key = block_attributes[ind]
            
            
            val = block_attributes[ind+1]
            
            #live_load[block_ind*6 - 6 + key - 1] =  float(val)
            
            equilibrium_matrix.append( float(val))
            equilibrium_matrix_row.append(  (block_ind)*6 - 6 + key - 1  )
            equilibrium_matrix_col.append(3*n)

            equilibrium_matrix_load.append( float(val))
            equilibrium_matrix_load_row.append(  (block_ind)*6 - 6 + key - 1  )
            equilibrium_matrix_load_col.append(3*n)
    
    

    c = np.zeros((3*n +1,1))
  
    c[3*n] = -1

    #c = matrix(c) 
    
    b = np.zeros((6*nb,1 )) 
    
    blocks_volume = data["blocks_volume"]
    
    density = 1.4e-5 # blocks density
    density = data["blocks_density"]
    

    weight = 0
    volume = 0
   
    for block_ind, block_volume in enumerate(blocks_volume):
        b[6*block_ind +2] = density*block_volume[1] 
        volume += block_volume[1] 
        weight += density*block_volume[1] 
        if data["seisme"]:
            equilibrium_matrix.append( density*block_volume[1]  )
            equilibrium_matrix_row.append(  (block_ind)*6    )
            equilibrium_matrix_col.append(3*n)

            equilibrium_matrix_load.append( density*block_volume[1]  )
            equilibrium_matrix_load_row.append(  (block_ind)*6   )
            equilibrium_matrix_load_col.append(3*n)

        if data["self_weight_failure"]:
            equilibrium_matrix.append( -density*block_volume[1]  )
            equilibrium_matrix_row.append(  (block_ind)*6 + 2    )
            equilibrium_matrix_col.append(3*n)

            equilibrium_matrix_load.append( -density*block_volume[1]  )
            equilibrium_matrix_load_row.append(  (block_ind)*6  + 2 )
            equilibrium_matrix_load_col.append(3*n)

    
    for ind ,block_ind in enumerate(loaded_blocks):
        
        b[6*block_ind:6*block_ind+6] = [ b[6*block_ind+i  ] - torseur_blocks_dead[ind][i]   for i in range(6)]

    
    print(torseur_blocks_dead)

 
    # sparse matrices for the stresses approach to be used in Mosek solver
    aval = equilibrium_matrix
    acol = equilibrium_matrix_col
    arow = equilibrium_matrix_row


    # sparse matrices for the loads approach to be used in Mosek solver

    aval_load = equilibrium_matrix_load
    acol_load = equilibrium_matrix_load_col
    arow_load = equilibrium_matrix_load_row

    avalk = equilibrium_matrix
    acolk = equilibrium_matrix_col
    arowk = equilibrium_matrix_row
    

    ####Steel rebars attributes

    bars_cap = []
    bars_data = data["bars_data"]
    blocks_centroid = data["blocks_centroid"]

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

                arow_load.append( block_ind*6 + pos )
                acol_load.append( var_col )
                aval_load.append( float(bar_tor[pos]) )



    # Since the actual value of Infinity is ignores, we define it solely
    # for symbolic purposes:
    inf = 0.0

    nef = n
    supports_nodes = data["supports_nodes"]
    all_phi_nodes = data["all_phi_nodes"]
    all_coh_nodes = data["all_coh_nodes"]

    all_coh_val = data["all_coh_val"]
    all_phi_val = data["all_phi_val"]

    # Get the directory of the current script
    rhino_file_dir = data["rhino_file_path"]

    # Construct the path to the log file in the same directory as the script
    log_path = os.path.join(rhino_file_dir, 'lower_upper_bounds_solutions.log')


    data["all_unique_points"] = output[6]
    data["all_triangles"] = output[7]
    #triangles_num = max(output[8]) 

    data["noc_triangles_coor"] = output[9]
    data["noc_triangles"] = output[10]
    #Ntriangles_num = max(output[11]) 

    blocks_mesh_parav( data )

    with open(log_path, "w") as f:
        f.write("")  # Clear the file


    if data["lower_bound_calc"]:
        
        prep_data_lower_bound_start = time.time()
        add_solution_title("Lower Bound Problem: ", log_path)
        with mosek.Task() as task:
            task = mosek.Task() 
            # Attach a printer to the task
            task.set_Stream(mosek.streamtype.log, streamprinter)

            bkx_lb = []
            cohs = []
            blx_lb = []
            bux_lb = []

            
            for i in range(nef):
                if (i+1) not in supports_nodes:
                    
                    fc = data["fc"]
                    if (i+1) in data["all_fc_nodes"]:
                        fc = float( data["all_fc_val"][ data["all_fc_nodes"].index(i + 1)])

                    if fc != False:    
                        blx_lb += [-inf, -inf , -abs(fc)  ] #variables lower bounds
                        bkx_lb += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.ra] #variables TYPE
                    else:
                        blx_lb += [-inf, -inf , -inf ] #variables lower bounds
                        bkx_lb += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.up] #variables TYPE
                else:
                    bkx_lb += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.fr] #variables TYPE
                    blx_lb += [-inf, -inf , -inf  ] #variables lower bounds
            
            bkx_lb +=  [mosek.boundkey.lo]*1 + [mosek.boundkey.ra]*len(bars_cap)

            
           
                
                
            blx_lb +=  [0]*1 +  [ -0*element for element in bars_cap] 

           
            ft = data["ft"] 
            for i in range(nef):
                if i not in supports_nodes:
                    bux_lb += [inf, inf , abs(ft) ] #variables upper bounds
                else:
                    bux_lb += [inf, inf , inf ] #variables upper bounds
                
            bux_lb += [inf]*1 + list(bars_cap) #variables upper bounds

            
            bkc_lb = [mosek.boundkey.fx]*nb*6 #type of constraints (equalities here)
            print("bbbbb", b )
            blc_lb = list(b)    #b for equalities
            buc_lb = list(b) #b for equalities
          
        
            
            c_lb = np.zeros((3*nef +1,1))
           
            c_lb[3*nef] = -1
            
            z_lb = [element for element in list(c_lb)] + [0]*len(bars_cap) #objective function 

            aval_lb = list(aval)
            arow_lb = list(arow)
            acol_lb = list(acol)


            ####### a remettre
            if False:
                # equalities for kinematic problem
                arowk = acol.copy()
                acolk = arow.copy()
                avalk = aval.copy()
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
      
           
            num_lists = 3*nef + 1 + len(bars_cap)  # Adjust the number of inner lists 
    
          
            asub_mos_lb = [[] for _ in range(num_lists)]
            aval_mos_lb  = [[] for _ in range(num_lists)]
            
            for i,j in enumerate(acol_lb):
                asub_mos_lb[int(j)].append(arow_lb[i])
                aval_mos_lb[int(j)].append(aval_lb[i])
            
            numvar_lb = len(bkx_lb)
            numcon_lb = len(bkc_lb)
            NUMANZ = 4 
    
         # Append 'numcon' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(numcon_lb)
    
            #Append 'numvar' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(numvar_lb)
    
            for j in range(numvar_lb):
                # Set the linear term c_j in the objective.
                task.putcj(j, z_lb[j])
                # Set the bounds on variable j
                # blx[j] <= x_j <= bux[j]
                task.putvarbound(j, bkx_lb[j], blx_lb[j], bux_lb[j])
            
            
            for j in range(len(aval_mos_lb)):
                # Input column j of A
                task.putacol(j,     asub_mos_lb[j],aval_mos_lb[j])            # Non-zero Values of column j.
            
         
        # Input the affine conic constraints
            task.appendafes(3*nef)
         
            list_con_lb = []
            list_coef_lb = []
         
            for ind in range(nef):
                if (ind +1) not in supports_nodes:
                    friction_ang = data["friction_ang"]
                    c = data["c"] # input value
  
                    if (ind+1) in all_phi_nodes:
                        friction_ang = float(all_phi_val[all_phi_nodes.index(ind + 1)])
      
                    if (ind+1) in all_coh_nodes:
                        c = all_coh_val[all_coh_nodes.index(ind+1)]
                    
                    cohs += [-c,0,0]
                    list_con_lb += [3*ind + 2 , 3*ind,3*ind + 1]
                   
                    list_coef_lb += [-math.tan(math.radians(friction_ang)),1.0,1.0]
                    
             
            task.putafefentrylist(range(3*nef - len(supports_nodes) ),                      # Rows
                                list_con_lb ,            # Columns 
                                list_coef_lb  )          #coefficients
    
            # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
            
            
            coh_vect_lb = []
            for ind in range(nef - len(supports_nodes) ):
                
                #if ind not in all_supports_ind:
                if True:
                   
                    coh = cohs[3*ind: 3*ind +3]
                   
                    coh_vect_lb += coh
                    quadcone  = task.appendquadraticconedomain(3)
                    task.appendacc(quadcone,          [3*ind, 3*ind +1, 3*ind + 2],  coh)                    # None if there is no b for conic 
           
            

            for i in range(numcon_lb):
                task.putconbound(i, bkc_lb[i], blc_lb[i], buc_lb[i])
    
    
           # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

            # Define a stream printer to capture the solver output
          

            # Attach the stream printer to the task for logging
            task.set_Stream(mosek.streamtype.log, lambda text: stream_printer(text, log_path))
            
            
            
            # Now, when you run the optimizer, detailed output will be captured in the specified log file
            prep_data_lower_bound_end = time.time()

            task.optimize()

            calc_lower_bound_end = time.time()

            log_objective_value(task,log_path)

            data["calc_time"]["prep_data_lower_bound"] =   prep_data_lower_bound_end - prep_data_lower_bound_start 
            data["calc_time"]["calc_data_lower_bound"] = calc_lower_bound_end - prep_data_lower_bound_end 


            
            # Print a summary containing information
            # about the solution for debugging purposes
            #task.solutionsummary(mosek.streamtype.msg)
            prosta = task.getprosta(mosek.soltype.itr)
            solsta = task.getsolsta(mosek.soltype.itr)
    
            # Output a solution
            xx_lb = task.getxx(mosek.soltype.itr)
            ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            #print("xxx",xx)
            #print(task.getsuc(mosek.soltype.itr)[0:6*nb])
            #print(task.getslc(mosek.soltype.itr)[0:6*nb])
            #print(ub)

            stat_parav(data, "lb",xx_lb[0:3*nef + 1])

            kin_parav(data, "lb",ub)

           

            if solsta == mosek.solsta.optimal:
                print("Objective: %s" % xx_lb[-1:])
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
          

    if data["upper_bound_calc"]:
        prep_data_upper_bound_start = time.time()
        add_solution_title("Upper Bound Problem: ", log_path)
        with mosek.Task() as task:
            task = mosek.Task() 

            
            bkx_ub = []
            cohs = []

            blx_ub = []
            bux_ub = []

            #bkx = [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.up]*nef   #variables lower bounded, if free just put .fr
            for i in range(nef):
                if (i+1) not in supports_nodes:
                    
                    fc = data["fc"]
                    if (i+1) in data["all_fc_nodes"]:
                        fc = float( data["all_fc_val"][ data["all_fc_nodes"].index(i + 1)])

                    if fc != False:    
                        blx_ub += [-inf, -inf , -abs(fc)*abs(area[i])  ] #variables lower bounds
                        bkx_ub += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.ra] #variables TYPE
                    else:
                        blx_ub += [-inf, -inf , -inf ] #variables lower bounds
                        bkx_ub += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.up] #variables TYPE
                else:
                    bkx_ub += [mosek.boundkey.fr, mosek.boundkey.fr, mosek.boundkey.fr] #variables TYPE
                    blx_ub += [-inf, -inf , -inf  ] #variables lower bounds

            bkx_ub +=  [mosek.boundkey.lo]*1 + [mosek.boundkey.ra]*len(bars_cap)
                              
            blx_ub +=  [0]*1 +  [ -0*element for element in bars_cap] 
            
            for i in range(nef):
                if (i+1) not in supports_nodes:
                    bux_ub += [inf, inf , 0 ] #variables upper bounds
                else:
                   
                    bux_ub += [inf, inf , inf ] #variables upper bounds
                
            bux_ub += [inf]*1 + list(bars_cap) #variables upper bounds
  
            
            bkc_ub = [mosek.boundkey.fx]*nb*6 #type of constraints (equalities here)
            blc_ub = list(b)    #b for equalities
            buc_ub = list(b) #b for equalities
          
        
            
            c = np.zeros((3*nef +1,1))
            #c =  spmatrix( [-1], [3*n], [0], (3*n +1, 1 ) )   
           
            c[3*nef] = -1
          
            z_ub = [element for element in list(c)] + [0]*len(bars_cap) #objective function 

            
       
            aval_ub = list(aval_load)
            arow_ub = list(arow_load)
            acol_ub = list(acol_load)
 
           
            num_lists = 3*nef + 1 + len(bars_cap)   # Adjust the number of inner lists as needed
    
          
            asub_mos_ub = [[] for _ in range(num_lists)]
            aval_mos_ub  = [[] for _ in range(num_lists)]
            
            for i,j in enumerate(acol_ub):
                asub_mos_ub[int(j)].append(arow_ub[i])
                aval_mos_ub[int(j)].append(aval_ub[i])
            

            numvar_ub = len(bkx_ub)
            numcon_ub = len(bkc_ub)
            NUMANZ = 4 
    
         # Append 'numcon' empty constraints.
            # The constraints will initially have no bounds.
            task.appendcons(numcon_ub)
    
            #Append 'numvar' variables.
            # The variables will initially be fixed at zero (x=0).
            task.appendvars(numvar_ub)
    
            for j in range(numvar_ub):
                # Set the linear term c_j in the objective.
                task.putcj(j, z_ub[j])
                task.putvarbound(j, bkx_ub[j], blx_ub[j], bux_ub[j])
            
            
            for j in range(len(aval_mos_ub)):
                
               
                # Input column j of A
                task.putacol(j,     asub_mos_ub[j],aval_mos_ub[j])            # Non-zero Values of column j.
            
         
        # Input the affine conic constraints
            # Create a matrix F such that F * x = [x(3),x(0),x(1),x(4),x(5),x(2)] 
            task.appendafes(3*nef)
         
            list_con_ub = []
            list_coef_ub = []
           
            
            for ind in range(nef):
                if (ind +1) not in supports_nodes:
                    friction_ang = data["friction_ang"]
                    c = data["c"]*area[ind] # input value
                    
                    if (ind+1) in all_phi_nodes:
                        friction_ang = float(all_phi_val[all_phi_nodes.index(ind + 1)])
      
                    if (ind+1) in all_coh_nodes:
                        c = all_coh_val[all_coh_nodes.index(ind+1)]*area[ind]
                        
                    cohs += [-c,0,0]
                    list_con_ub += [3*ind + 2 , 3*ind,3*ind + 1]
                   
                    list_coef_ub += [-math.tan(math.radians(friction_ang)),1.0,1.0]
                    
             
            task.putafefentrylist(range(3*nef - len(supports_nodes) ),                      # Rows
                                list_con_ub ,            # Columns 
                                list_coef_ub  )          #coefficients
    
        
            coh_vect_ub = []
            for ind in range(nef - len(supports_nodes) ):
                
                #if ind not in all_supports_ind:
                if True:
                   
                    coh = cohs[3*ind: 3*ind +3]
                   
                    coh_vect_ub += coh
                    quadcone  = task.appendquadraticconedomain(3)
                    task.appendacc(quadcone,          [3*ind, 3*ind +1, 3*ind + 2],  coh)                    # None if there is no b for conic 
           
            
           
            
            for i in range(numcon_ub):
                task.putconbound(i, bkc_ub[i], blc_ub[i], buc_ub[i])
    
    
        # Input the objective sense (minimize/maximize)
            task.putobjsense(mosek.objsense.minimize)

            task.set_Stream(mosek.streamtype.log, lambda text: stream_printer(text, log_path))

            prep_data_upper_bound_end = time.time()
            
            # Optimize the task
            task.optimize()

            calc_upper_bound_end = time.time()

            data["calc_time"]["prep_data_upper_bound"] =  prep_data_upper_bound_end - prep_data_upper_bound_start 
            data["calc_time"]["calc_data_upper_bound"] =  calc_upper_bound_end - prep_data_upper_bound_end 


            log_objective_value(task,log_path)
            #task.writedata("cqo1.ptf")
            # Print a summary containing information
            # about the solution for debugging purposes
            #task.solutionsummary(mosek.streamtype.msg)
            prosta = task.getprosta(mosek.soltype.itr)
            solsta = task.getsolsta(mosek.soltype.itr)
    
            # Output a solution
            xx_ub = task.getxx(mosek.soltype.itr)
            ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            

            stat_parav(data, "ub",xx_ub[0:3*nef + 1])

            kin_parav(data, "ub",ub)
            
            #skc, y = task.getsolution(mosek.soltype.bas)
            if solsta == mosek.solsta.optimal:
                print("Objective: %s" % xx_ub[3*nef])
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
          
        
    #non associative calc
    tol =1
    itr = 0
    
    alpha= data["alpha"]
    beta= data["beta"]
    #tolerance=0.001
    tolerance= data["tolerance"]
    

 
    if data["ub_non_ass_calc"]:

        
        # Construct the path to the log file in the same directory as the script
        log_path = os.path.join(rhino_file_dir, 'ub_non_ass_solutions.log')

        with open(log_path, "w") as f:
            f.write("")  # Clear the file

        

        while tol >tolerance: 
    
            # Check if the solver found an optimal solution
            if solsta == mosek.solsta.optimal:
          
                if itr>1:
                    normalvalues_old = normalvalues
                else:
                    normalvalues_old = np.zeros((int(nef) - len(supports_nodes),1))
                    normal = []
                    for element in range(0,len(xx_ub)-2,3):
                        if int(element /3 + 1) not in supports_nodes:
                            normal += [xx_ub[element+2]]

                    normal = np.array(normal).reshape(-1,1)
                    # Apply KMeans with 2 clusters to normal loads values
                    kmeans_normal = KMeans(n_clusters=2, random_state=0).fit(normal)

                    # Split values based on the KMeans labels
                    cluster_1 = [normal[i] for i in range(len(normal)) if kmeans_normal.labels_[i] == 0]
                    min_normal = min(cluster_1)

                statsolver_old = xx_ub[3*nef]   
                
                normalvalues = np.empty((0,1))
               
                
                for element in range(0,len(xx_ub)-2,3):
                    if int(element /3 + 1) not in supports_nodes:
                        normalvalues = np.append(normalvalues,xx_ub[element+2]) 
            
            coh_vect_new = [0,0,0]*len(coh_vect_ub)
            friction_coef_new = [1,1,1]*len(list_coef_ub)
       

            
            for element in range(int(nef - len(supports_nodes))):

                if abs(normalvalues[element]) > 1e-5:
                    coh_corr = normalvalues[element]
                else:
                    coh_corr = min_normal
                 


                if itr >1:
                    coh_vect_new[3*element] =  0.00001*coh_corr +  coh_vect_ub[3*element] + (1+alpha)*abs(list_coef_ub[3*element])*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element]) 
                    friction_coef_new[3*element] = -math.tan(math.radians(data["friction_ang"])*alpha)
                    
                else:
        
                    coh_vect_new[3*element] = 0.00001*coh_corr + coh_vect_ub[3*element] + (1+alpha)*abs(list_coef_ub[3*element])*(normalvalues[element] )
                    friction_coef_new[3*element] = -math.tan(math.radians(data["friction_ang"])*alpha)
                    
            #-list_coef[3*element]*alpha
            statsol_old = xx_ub

           
    
            #mosek
            with mosek.Task() as task:
                # Attach a printer to the task
                task.set_Stream(mosek.streamtype.log, streamprinter)
        
                
             # Append 'numcon' empty constraints.
                # The constraints will initially have no bounds.
                task.appendcons(numcon_ub)
        
                #Append 'numvar' variables.
                # The variables will initially be fixed at zero (x=0).
                task.appendvars(numvar_ub)
        
                for j in range(numvar_ub):
                    # Set the linear term c_j in the objective.
                    task.putcj(j, z_ub[j])
                    # Set the bounds on variable j
                    # blx[j] <= x_j <= bux[j]
                    task.putvarbound(j, bkx_ub[j], blx_ub[j], bux_ub[j])
                
                
                for j in range(len(aval_mos_ub)):
                    
                   
                    # Input column j of A
                    task.putacol(j,      asub_mos_ub[j],aval_mos_ub[j])            # Non-zero Values of column j.
                
                
        
            # Input the affine conic constraints

                task.appendafes(3*nef)
                
                list_con_ub = []
                

                for ind in range(nef):
                    if (ind+1) not in supports_nodes:
                    
                        list_con_ub += [3*ind + 2 , 3*ind,3*ind + 1]
                        #list_coef += [-math.tan(math.radians(friction_ang)),1.0,1.0]
                    
                task.putafefentrylist(range(3*nef - len(supports_nodes) ),                      # Rows
                                    list_con_ub ,            # Columns 
                                    friction_coef_new  )          #coefficients
        
                # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
       

           
                for ind in range(nef - len(supports_nodes) ):
                    
                    #if ind not in all_supports_ind:
                    if True:
                        quadcone  = task.appendquadraticconedomain(3)
                        task.appendacc(quadcone, [3*ind, 3*ind +1, 3*ind + 2],   [coh_vect_new[3*ind], coh_vect_new[3*ind +1] , coh_vect_new[3*ind +2] ] )                    # None if there is no b for conic 
           
                
                for i in range(numcon_ub):
                    task.putconbound(i, bkc_ub[i], blc_ub[i], buc_ub[i])
        
        
            # Input the objective sense (minimize/maximize)
                task.putobjsense(mosek.objsense.minimize)

                #task.putintparam(mosek.iparam.log, 0)

                # Attach the stream printer to the task for logging
                task.set_Stream(mosek.streamtype.log, lambda text: stream_printer(text, log_path))
                # Optimize the task
                task.optimize()
                with open(log_path, "a") as f:
                    f.write(f"min normal: {min(normalvalues)}\n")
                    f.write(f"min_normal: {min_normal}\n")

                    f.flush()
                #task.writedata("cqo1.ptf")
                # Print a summary containing information
                # about the solution for debugging purposes
                #task.solutionsummary(mosek.streamtype.msg)
                prosta = task.getprosta(mosek.soltype.itr)
                solsta = task.getsolsta(mosek.soltype.itr)
        
                # Output a solution
                xx_ub = task.getxx(mosek.soltype.itr)
                
                ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            
                
                #skc, y = task.getsolution(mosek.soltype.bas)
                if solsta == mosek.solsta.optimal:
                    print("Objective: %s" % xx_ub[3*nef])
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
     
                    if itr >1:
                    
                        tol = abs(xx_ub[3*nef] - statsolver_old)/abs(xx_ub[3*nef] )
                        print("tolerance", tol)
                else:
                    statsol = statsol_old
                    
                if alpha <= 0.0001:
                    alpha = 0.0001
                else:
                    alpha = alpha*0.5
        
                itr = itr+1
                # Clear or reset the task
                
                del task  # Delete the task object
                    
        stat_parav(data, "ub_non_ass",xx_ub[0:3*nef + 1])

        kin_parav(data, "ub_non_ass",ub)



    if data["lb_non_ass_calc"]:

        tol =1
        itr = 0
        
        alpha= data["alpha"]
        beta= data["beta"]
        #tolerance=0.001
        tolerance= data["tolerance"]

        
        # Construct the path to the log file in the same directory as the script
        log_path = os.path.join(rhino_file_dir, 'lb_non_ass_solutions.log')

        with open(log_path, "w") as f:
            f.write("")  # Clear the file

        

        while tol >tolerance: 
    
            # Check if the solver found an optimal solution
            if solsta == mosek.solsta.optimal:
          
                if itr>1:
                    normalvalues_old = normalvalues
                else:
                    normalvalues_old = np.zeros((int(nef) - len(supports_nodes),1))
                    normal = []
                    for element in range(0,len(xx_lb)-2,3):
                        if int(element /3 + 1) not in supports_nodes:
                            normal += [xx_lb[element+2]]

                    normal = np.array(normal).reshape(-1,1)
                    # Apply KMeans with 2 clusters to normal loads values
                    kmeans_normal = KMeans(n_clusters=2, random_state=0).fit(normal)

                    # Split values based on the KMeans labels
                    cluster_1 = [normal[i] for i in range(len(normal)) if kmeans_normal.labels_[i] == 0]
                    min_normal = min(cluster_1)

                statsolver_old = xx_lb[3*nef]   
                
                normalvalues = np.empty((0,1))
               
                
                for element in range(0,len(xx_lb)-2,3):
                    if int(element /3 + 1) not in supports_nodes:
                        normalvalues = np.append(normalvalues,xx_lb[element+2]) 
            
            coh_vect_new = [0,0,0]*len(coh_vect_lb)
            friction_coef_new = [1,1,1]*len(list_coef_lb)
       

            
            for element in range(int(nef - len(supports_nodes))):

                if abs(normalvalues[element]) > 1e-5:
                    coh_corr = normalvalues[element]
                else:
                    coh_corr = min(normalvalues)
                 


                if itr >1:
                    coh_vect_new[3*element] =  0.000001*coh_corr +  coh_vect_lb[3*element] + (1+alpha)*abs(list_coef_lb[3*element])*(beta*normalvalues[element] + (1-beta)*normalvalues_old[element]) 
                    friction_coef_new[3*element] = -math.tan(math.radians(data["friction_ang"])*alpha)
                    
                else:
        
                    coh_vect_new[3*element] = 0.000001*coh_corr + coh_vect_lb[3*element] + (1+alpha)*abs(list_coef_lb[3*element])*(normalvalues[element] )
                    friction_coef_new[3*element] = -math.tan(math.radians(data["friction_ang"])*alpha)
                    
            #-list_coef[3*element]*alpha
            statsol_old = xx_lb

           
    
            #mosek
            with mosek.Task() as task:
                # Attach a printer to the task
                task.set_Stream(mosek.streamtype.log, streamprinter)
    
             # Append 'numcon' empty constraints.
                # The constraints will initially have no bounds.
                task.appendcons(numcon_lb)
        
                #Append 'numvar' variables.
                # The variables will initially be fixed at zero (x=0).
                task.appendvars(numvar_lb)
        
                for j in range(numvar_lb):
                    # Set the linear term c_j in the objective.
                    task.putcj(j, z_lb[j])
                    # Set the bounds on variable j
                    # blx[j] <= x_j <= bux[j]
                    task.putvarbound(j, bkx_lb[j], blx_lb[j], bux_lb[j])
                
                
                for j in range(len(aval_mos_lb)):
                    
                   
                    # Input column j of A
                    task.putacol(j,      asub_mos_lb[j],aval_mos_lb[j])            # Non-zero Values of column j.
                
            # Input the affine conic constraints

                task.appendafes(3*nef)
                
                list_con_lb = []
                

                for ind in range(nef):
                    if (ind+1) not in supports_nodes:
                    
                        list_con_lb += [3*ind + 2 , 3*ind,3*ind + 1]
                      
                    
                task.putafefentrylist(range(3*nef - len(supports_nodes) ),                      # Rows
                                    list_con_lb ,            # Columns 
                                    friction_coef_new  )          #coefficients
        
                # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3 
       

           
                for ind in range(nef - len(supports_nodes) ):
                    
                    #if ind not in all_supports_ind:
                    if True:
                        quadcone  = task.appendquadraticconedomain(3)
                        task.appendacc(quadcone, [3*ind, 3*ind +1, 3*ind + 2],   [coh_vect_new[3*ind], coh_vect_new[3*ind +1] , coh_vect_new[3*ind +2] ] )                    # None if there is no b for conic 
           
               
                for i in range(numcon_lb):
                    task.putconbound(i, bkc_lb[i], blc_lb[i], buc_lb[i])
        
        
            # Input the objective sense (minimize/maximize)
                task.putobjsense(mosek.objsense.minimize)

                #task.putintparam(mosek.iparam.log, 0)

                # Attach the stream printer to the task for logging
                task.set_Stream(mosek.streamtype.log, lambda text: stream_printer(text, log_path))
                # Optimize the task
                task.optimize()
                with open(log_path, "a") as f:
                    f.write(f"min normal: {min(normalvalues)}\n")
                    f.write(f"min_normal: {min_normal}\n")

                    f.flush()
               
                prosta = task.getprosta(mosek.soltype.itr)
                solsta = task.getsolsta(mosek.soltype.itr)
        
                # Output a solution
                xx_lb = task.getxx(mosek.soltype.itr)
                
                ub = [task.getsuc(mosek.soltype.itr)[ind] - task.getslc(mosek.soltype.itr)[ind]  for ind in range(6*nb) ] # get dual variables for equilibrium but why upper ? not slb
            
                
                #skc, y = task.getsolution(mosek.soltype.bas)
                if solsta == mosek.solsta.optimal:
                    print("Objective: %s" % xx_lb[3*nef])
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
    
                if itr >1:
                   
                    tol = abs(xx_lb[3*nef] - statsolver_old)/abs(xx_lb[3*nef] )
                    print("tolerance", tol)
            else:
                statsol = statsol_old
                
            if alpha <= 0.0001:
                alpha = 0.0001
            else:
                alpha = alpha*0.5
    
            itr = itr+1
            # Clear or reset the task
            
            del task  # Delete the task object
                
        stat_parav(data, "lb_non_ass",xx_lb[0:3*nef + 1])

        kin_parav(data, "lb_non_ass",ub)
    

    # Construct the path to the data.json file
    json_file_path = os.path.join(data["script_dir"] , 'calc_time.json')
    

    with open(json_file_path, 'w') as file:
        json.dump(data["calc_time"], file, indent=4)


    return ub


def main():

    # Get the directory of the current script
    script_dir = os.path.dirname(__file__)


    # Get the path of the currently open Rhino file
    rhino_file_path = rs.DocumentPath()
  
    # Construct the path to the data.json file
    json_file_path = os.path.join(rhino_file_path, 'data.json')

    # Load the data from the JSON file
    with open(json_file_path, 'r') as file:
        data = json.load(file)

    data["script_dir"] = script_dir
    data["rhino_file_path"] = rhino_file_path

    #call main function
    #main(data)

    # Get the current date and time from the URL
    now = urlopen('http://just-the-time.appspot.com/')
    now_str = now.read().strip().decode('utf-8')
    try:
    # Try setting to 'en_US' for UNIX/Linux systems or 'English_United States' for Windows.
        locale.setlocale(locale.LC_ALL, 'en_US')  # Try 'en_US.UTF-8' if 'en_US' doesn't work
    except locale.Error:
        print("Locale 'en-US' could not be set.")

    exp_date = datetime( year = 2024, month = 8 , day = 31   )

    now_datetime = datetime.strptime(now_str, '%Y-%m-%d %H:%M:%S')

    # Compare the two dates
    if exp_date > now_datetime:
        #print("You can still use the script")
    
    
        ## read rhino file
        data = read_data_from_rhino(data)

        ## call C++ function to generate mesh and construct problem
        data = generate_mesh_and_construct_matrices(data)

        ## solve different approaches and generate paraview file

        ub = solve_problem_with_solvers(data)

        #
        #block_disp_plot(data["blocks_centroid"],data["blocks"],data["faces"],data["points"],ub)



#call main function
main()



