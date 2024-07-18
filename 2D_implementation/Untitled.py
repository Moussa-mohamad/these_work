
import rhinoscriptsyntax as rs
import numpy as np
import Rhino
import scriptcontext as sc
import os
# Function to round point coordinates to 3 decimal places
def round_point(pt,digits):
    return Rhino.Geometry.Point3d(round(pt.X, digits), round(pt.Y, digits), round(pt.Z, digits))

def main():
    weight13 = 5/9
    weight2 = 8/9

    edges_points_format = []
    active_edges = []
    points = []
    ref_blocks = []
    normal_vectors = []
    tang_vectors = []
    tang_moment_vectors = []
    normal_moment_vectors = []
    blocks = []
    edges = []
    false_neighboors = []
    true_neighboors = []
    blocks_centroid = []
    surfaces_guid =[]
    edges_var_pos = []
    blocks_egdes_guid = dict()
    surfaces_boundary = []
    eq_matrix_row = []
    eq_matrix_col = []
    eq_matrix_val = []
    active_edges_pts_coor = []
    joints_length = []
    var_pos = 0
    pos_to_modify = []
    keys_of_interest = ["fx",  "fz", "my"] # keys for blocks loading applied on their centroids
    # Get all surface objects in the document
    surfaces = rs.ObjectsByType(rs.filter.surface, True)
    nb = len(surfaces)

    curves_rhino = rs.ObjectsByType(rs.filter.curve, True)

    filtered_curves = []
    for curve in curves_rhino:
        if rs.GetUserText(curve):
            filtered_curves.append(curve)
    
    curves_rhino = filtered_curves
    
    

    curve_geometries = [rs.coercecurve(crv) for crv in curves_rhino]


    # Get bounding boxes
    bounding_boxes_curves = [geo.GetBoundingBox(True) for geo in curve_geometries]

    for i, srf in enumerate(surfaces):
        print(rs.ObjectName(srf))
        surfaces_guid += [srf]
        surface_boundary = rs.DuplicateSurfaceBorder(srf)
        surfaces_boundary += [rs.CopyObject(surface_boundary)  ]
        blocks_egdes_guid[i] = {"contact":[], "non_contact":surface_boundary}

    if not surfaces:
        print("No surfaces found.")
        return
    surface_breps = [rs.coercebrep(srf) for srf in surfaces]

    bounding_boxes = [srf.GetBoundingBox(True) for srf in surface_breps]
    
    #inter = rs.CurveCurveIntersection(curve_geometries[0] , surfaces_boundary[0] )
    # Create an R-tree
    rtree = Rhino.Geometry.RTree()


    # Add bounding boxes to the R-tree
    for i, bbox in enumerate(bounding_boxes):
        rtree.Insert( bbox,i)

    # Find surfaces within bounding boxes
    for i, bbox in enumerate(bounding_boxes):
        
        #block centroid loading
        obj = surfaces[i]
        atb_keys = rs.GetUserText(obj)
        b_att = dict()
        if atb_keys:
            for atb_key in atb_keys:
                if atb_key.lower() in keys_of_interest:
                    key_ind = keys_of_interest.index(atb_key.lower())
                    

                    b_att[key_ind +1] = float(rs.GetUserText(obj, atb_key)) 
                
       
        if not bbox.IsValid:
            continue

        # Define a callback function for the R-tree search
        def search_callback(sender, e):
            intersection_rep = [min(e.Id+1,i+1), max(e.Id+1,i+1)] # [first block ind, second block ind ]
            
     
            if e.Id != i and intersection_rep not in false_neighboors:  # Exclude self
         
                
                if intersection_rep not in true_neighboors:
                    
                    ref_blocks.extend([i+1])
                    intersects = rs.CurveCurveIntersection(surfaces_boundary[e.Id], surfaces_boundary[i]  )
                 
                    
                    curve_inter_num = 0
            
                    if intersects != None:
                        
                        for intersect in intersects:
                           
                            if intersect[0] != 1:
                                
                                add_to_contact = []
                                add_to_non_contact = []
                                part_to_split = []

                                curve_inter_num += 1

                                non_contacts_pos = []
                               
                                for non_contact_pos, non_contact in enumerate(blocks_egdes_guid[e.Id]["non_contact"]):
                                    
                                    sub_intersects = rs.CurveCurveIntersection(non_contact, surfaces_boundary[i]  )
                                  
                                    if sub_intersects != None:
                                        
                                        for sub_intersect in sub_intersects:
                                            
                                            if sub_intersect[0] != 1:
                                                
                                                intersect_pt1 = round_point(intersect[1],3)
                                               
                                                intersect_pt2 = round_point(intersect[2],3)

                                                sub_intersect_pt1 = round_point(sub_intersect[1],3)
                                                sub_intersect_pt2 =  round_point(sub_intersect[2],3)

                                                possibility1 = (intersect_pt1 == sub_intersect_pt1) and  (intersect_pt2 == sub_intersect_pt2)
                                                possibility2 = (intersect_pt1 == sub_intersect_pt2) and  (intersect_pt2 == sub_intersect_pt1)
                                            
                                                if possibility1 or possibility2:
                                                    sub_intersect_test = sub_intersect
                                                    #break 
                                        
                                        if sub_intersect[0] != 1:
                                            non_contacts_pos += [non_contact_pos]

                                            sub_intersect = sub_intersect_test
                                            
                                            split2 = []
                                            
                                            point1_param = sub_intersect[5]
                                            #clspoint_param = rs.CurveClosestPoint(non_contact,  sub_intersect[2])
                                            split1 = rs.SplitCurve(non_contact, point1_param,  delete_input = False)
                                            
                                            
                                        
                                            if split1  == None:
                                                
                                                split2 = rs.SplitCurve(non_contact, sub_intersect[6], delete_input = False)
                                                if split2 == None:
                                                   
                                                    add_to_contact += [rs.CopyObject(non_contact)]
                                            else:
                                                
                                                point2_coor = sub_intersect[2]
                                                
                                                clspoint1_param = rs.CurveClosestPoint(split1[0], point2_coor)
                                                clspoint2_param = rs.CurveClosestPoint(split1[1], point2_coor)
                                                
                                                clspoint1_coor = rs.EvaluateCurve(split1[0], clspoint1_param)
                                                clspoint2_coor = rs.EvaluateCurve(split1[1], clspoint2_param)
                                        
                                                dist1 = np.linalg.norm(point2_coor - clspoint1_coor)
                                                dist2 = np.linalg.norm(point2_coor - clspoint2_coor)
                                                
                                                
                                                if abs(dist1 - dist2) < 1e-5:
                                                    
                                                    
                                                    if rs.CurveCurveIntersection(split1[0], surfaces_boundary[i]   ) != None:
                                                        if rs.CurveCurveIntersection(split1[0], surfaces_boundary[i])[0][0] !=1:
                                                            
                                                            add_to_contact += [split1[0]]
                                                            add_to_non_contact += [split1[1]]
                                                        else:
                                                        
                                                            add_to_contact += [split1[1]]
                                                            add_to_non_contact += [split1[0]]
                                                    else:
                                                        add_to_contact += [split1[1]]
                                                        add_to_non_contact += [split1[0]]
                                                    
                                                else:
                                                    
                                                    
                                                    if dist1 <= dist2:
                                                        
                                                        add_to_non_contact += [split1[1]]
                                                        part_to_split = split1[0]
                                                        
                                                    else:
                                                        
                                                        add_to_non_contact += [split1[0]]
                                                        part_to_split = split1[1]
                                                    
                                                    part_points = rs.CurvePoints( part_to_split)
                                                  

                                                    cut_point_coor = rs.EvaluateCurve( part_to_split , sub_intersect[6])
                                                
                                                    cut_point_param = rs.CurveClosestPoint(part_to_split, cut_point_coor ) # or replace cut_point_coor with point2_coor
                                                    split2 = rs.SplitCurve(part_to_split, cut_point_param, delete_input = False)
                                                    
                                                    if split2 == None:
                                                        add_to_contact += [part_to_split]
                                                
                                            if split2 != None and len(split2) != 0:
                                                
                                                if rs.CurveCurveIntersection(split2[0], surfaces_boundary[i]  ) != None:
                                                    
                                                    if rs.CurveCurveIntersection(split2[0], surfaces_boundary[i])[0][0] !=1:
                                                     
                                                        add_to_contact += [split2[0]]
                                                        add_to_non_contact += [split2[1]]
                                                    else:
                                                        add_to_contact += [split2[1]]
                                                        add_to_non_contact += [split2[0]]

                                                else:
                                                    add_to_contact += [split2[1]]
                                                    add_to_non_contact += [split2[0]]
                                                
                                                rs.DeleteObjects(part_to_split)
                       

                                for non_contact_pos in non_contacts_pos:
                                    if len(blocks_egdes_guid[e.Id]["contact"]) != 0:
                                        rs.DeleteObjects(blocks_egdes_guid[e.Id]["non_contact"][non_contact_pos])

                                    blocks_egdes_guid[e.Id]["non_contact"].pop(non_contact_pos)
                                
                                
                                blocks_egdes_guid[e.Id]["non_contact"].extend(add_to_non_contact)
                                
                                blocks_egdes_guid[e.Id]["contact"].extend( add_to_contact)

                                


                                #for second block
                                add_to_contact = []
                                add_to_non_contact = []
                                part_to_split = []

                                curve_inter_num += 1

                                non_contacts_pos = []
                               
                                for non_contact_pos, non_contact in enumerate(blocks_egdes_guid[i]["non_contact"]):
                                    
                                    sub_intersects = rs.CurveCurveIntersection(non_contact, surfaces_boundary[e.Id]  )
                                  
                                    if sub_intersects != None:
                                        
                                        for sub_intersect in sub_intersects:
                                            
                                            if sub_intersect[0] != 1:
                                                
                                                intersect_pt1 = round_point(intersect[1],3)
                                               
                                                intersect_pt2 = round_point(intersect[2],3)

                                                sub_intersect_pt1 = round_point(sub_intersect[1],3)
                                                sub_intersect_pt2 =  round_point(sub_intersect[2],3)

                                                possibility1 = (intersect_pt1 == sub_intersect_pt1) and  (intersect_pt2 == sub_intersect_pt2)
                                                possibility2 = (intersect_pt1 == sub_intersect_pt2) and  (intersect_pt2 == sub_intersect_pt1)
                                               
                                                if possibility1 or possibility2:
                                                    sub_intersect_test = sub_intersect
                                                    #break 
                                        
                                        if sub_intersect[0] != 1:
                                            non_contacts_pos += [non_contact_pos]

                                            sub_intersect = sub_intersect_test
                                             
                                            split2 = []
                                            
                                            point1_param = sub_intersect[5]
                                   
                                            split1 = rs.SplitCurve(non_contact, point1_param,  delete_input = False)
                                            
                                            
                                            if split1  == None:
                                                
                                                
                                                split2 = rs.SplitCurve(non_contact, sub_intersect[6], delete_input = False)
                                                if split2 == None:
                                                    add_to_contact += [rs.CopyObject(non_contact)]
                                                    
        
                                            else:
                          
                                                part_points = rs.CurvePoints( split1[0])
                                           
                                                part_points = rs.CurvePoints( split1[1])
                                             
                                                  
                                                point2_coor = sub_intersect[2]
                                                
                                                clspoint1_param = rs.CurveClosestPoint(split1[0], point2_coor)
                                                clspoint2_param = rs.CurveClosestPoint(split1[1], point2_coor)
                                                
                                                clspoint1_coor = rs.EvaluateCurve(split1[0], clspoint1_param)
                                                clspoint2_coor = rs.EvaluateCurve(split1[1], clspoint2_param)
                                        
                                                dist1 = np.linalg.norm(point2_coor - clspoint1_coor)
                                                dist2 = np.linalg.norm(point2_coor - clspoint2_coor)
                                            
                                                
                                                if abs(dist1 - dist2) < 1e-5:
                                                    
                                                    
                                                    if rs.CurveCurveIntersection(split1[0], surfaces_boundary[e.Id]   ) != None:
                                                        if rs.CurveCurveIntersection(split1[0], surfaces_boundary[e.Id])[0][0] !=1:

                                                            add_to_contact += [split1[0]]
                                                            add_to_non_contact += [split1[1]]
                                                        else:
                                                            add_to_contact += [split1[1]]
                                                            add_to_non_contact += [split1[0]]
                                                    else:
                                                        add_to_contact += [split1[1]]
                                                        add_to_non_contact += [split1[0]]
                                                    
                                                else:
                                                    
                                                    
                                                    if dist1 <= dist2:
                                                        
                                                        add_to_non_contact += [split1[1]]
                                                        part_to_split = split1[0]
                                                        
                                               
                                                    else:
                                                    
                                                        add_to_non_contact += [split1[0]]
                                                        part_to_split = split1[1]
                                                    
                                           
                                                    part_points = rs.CurvePoints( part_to_split)
                                                    
                                                    cut_point_coor = rs.EvaluateCurve( part_to_split , sub_intersect[6])
                                                    cut_point_param = rs.CurveClosestPoint(part_to_split, cut_point_coor ) # or replace cut_point_coor with point2_coor

                                                    split2 = rs.SplitCurve(part_to_split, cut_point_param, delete_input = False)
                                          
                                                    if split2 == None:
                                                        add_to_contact += [part_to_split]
                              
                                       
                                            if split2 != None and len(split2) != 0:
                                                
                                               
                                                
                                                if rs.CurveCurveIntersection(split2[0], surfaces_boundary[e.Id]  ) != None:
                                                    
                                                    if rs.CurveCurveIntersection(split2[0], surfaces_boundary[e.Id])[0][0] !=1:
                                                       
                                                 
                                                        add_to_contact += [split2[0]]
                                                        add_to_non_contact += [split2[1]]
                                                    else:
                                                        add_to_contact += [split2[1]]
                                                        add_to_non_contact += [split2[0]]

                                                else:
                                                    add_to_contact += [split2[1]]
                                                    add_to_non_contact += [split2[0]]
                                                
                                                rs.DeleteObjects(part_to_split)
                                            
                             
                                for non_contact_pos in non_contacts_pos:
                                    if len(blocks_egdes_guid[i]["contact"]) != 0:
                                        rs.DeleteObjects(blocks_egdes_guid[i]["non_contact"][non_contact_pos])

                                    blocks_egdes_guid[i]["non_contact"].pop(non_contact_pos)
                                
                                blocks_egdes_guid[i]["non_contact"].extend(add_to_non_contact)
                                
                                blocks_egdes_guid[i]["contact"].extend( add_to_contact)

            
                           
                    else:
                        false_neighboors.append(intersection_rep)
                    
                    if curve_inter_num == 0:
                            false_neighboors.append(intersection_rep)
                    else:
                        true_neighboors.append(intersection_rep) 
        
        # Perform the R-tree search
        rtree.Search(bbox, search_callback)

        
        
        
        
        block_centroid_r = rs.SurfaceAreaCentroid(surfaces_guid[i])[0]
        
        blocks_centroid.append([coor for coor in block_centroid_r])
    

        result =  is_point_inside_curve(block_centroid_r,surfaces_boundary[i])

        if result == 1:
            is_inside = True
        elif result == 0:
            is_inside = False
        elif result == -1:
            print("Error: Not sure, recheck this case.")
            return None
        else:
            print("Error: block centroid on the curve")
            return None


        bounding_boxes_block_curves = [bounding_boxes[i]] + bounding_boxes_curves

        # Create an R-tree
        rtree1 = Rhino.Geometry.RTree()

        # Add bounding boxes to the R-tree
        for j, bbox in enumerate(bounding_boxes_block_curves):
            rtree1.Insert(bbox, j)

        #block_curve_inter = []

        blocki_bbox = surfaces_boundary[i]

        all_non_contacts = []

        block_egdes_guid = blocks_egdes_guid[i]
  

        for block_non_contact in block_egdes_guid["non_contact"]:
            all_non_contacts += [  element  for element in rs.ExplodeCurves(block_non_contact) ] 
          
            rs.DeleteObjects(block_non_contact)

        blocks_egdes_guid[i]["non_contact"] = all_non_contacts

        block_egdes_guid = blocks_egdes_guid[i]
        

        # Define a callback function for the R-tree search
        def search_callback1(sender1, e ):
            if e.Id != 0:
                add_to_contact = []
                for block_non_contact in block_egdes_guid["non_contact"]:
                
                    inter = rs.CurveCurveIntersection(   curve_geometries[e.Id - 1 ] , block_non_contact   )
              
                    if inter != None:
                        if inter[0][0] == 2:
                            add_to_contact += [block_non_contact]
                
                blocks_egdes_guid[i]["contact"] += add_to_contact
                for element in add_to_contact:

                    blocks_egdes_guid[i]["non_contact"].remove(element)
                                            
                            
        rtree1.Search(bounding_boxes_block_curves[0] , search_callback1)
        
        
        
        block_egdes_guid = blocks_egdes_guid[i]
        
     
        edges_pos = []
        div_num = 2
       

        fx_block_row = []
        fx_block_col = [] 
        fx_block_val = [ ]

        fz_block_row = []
        fz_block_col = [] 
        fz_block_val = [ ]

        my_block_row = []
        my_block_col = [] 
        my_block_val = [ ]
       
        for index, contact in enumerate(block_egdes_guid["contact"]):
           
            curves = rs.ExplodeCurves(contact)

            for curve in curves:
                factor = -1
                is_line = rs.IsLine(curve)
                if not is_line:
                    curve_points =  rs.DivideCurve(curve, div_num)
                    mid_point = curve_points[int(div_num/2)] 
                
                 
                else:
                 
                    curve_points =  rs.CurvePoints(curve)
                    mid_point = round_point( curve_points[0],3)*0.5 + round_point( curve_points[1],3)*0.5 
             
                    pt1 = curve_points[0]
                    pt2 = curve_points[1]
               
                    joint = rs.AddLine(pt1, pt2 )

                    joint_length = rs.CurveLength(joint)
                    rs.DeleteObject(joint)
            
                intersection_type, intersection_points = find_intersection_line_curve(mid_point, block_centroid_r, surfaces_boundary[i] )
             
                if intersection_type == "curve" :
         
                    print("Error: normal vector is normal to the line joining the centroids of the edge and block")
                    return None
                else:
                    A = is_inside
                    B = (len(intersection_points) % 2 == 0)
                    
                    if (A and B) or ( (not A) and ( not B)):
                        Case = "1"
                    else:
                         Case = "2"
                
            
                points_pos = []
                points_coor = []
             
                for point in curve_points:

                    point = round_point( point,3)
                    points, point_pos = check_add_point(points, [point[0], point[1], point[2] ] )
                    points_pos += [point_pos]
                    points_coor += [point]
                
                #points_pos, indices = sort_and_get_indices(points_pos)
                
                if points_pos[0] > points_pos[-1]:
                    perm = points_pos[0]
                    points_pos[0] = points_pos[-1]
                    points_pos[-1] = perm

                    perm = points_coor[0]
                    points_coor[0] = points_coor[-1]
                    points_coor[-1] = perm

                if points_pos[0] == points_pos[-1]:
                    print("Error last and first points are similar")
                    return None



                print(edges_points_format)

                if points_pos not in edges_points_format:
                    
                    normal_vectors.append([])
                    tang_vectors.append([])
                    normal_moment_vectors.append([])
                    tang_moment_vectors.append([])

                    if is_line:
                        joints_length.append(joint_length)
                        active_edges_pts_coor.append([mid_point[0], mid_point[1], mid_point[2] ])
                        normal = round_point(compute_normal_at_point( surfaces_boundary[i], mid_point ),2)

                        dot_product =  round(calculate_dot_product_rhino(normal, block_centroid_r, mid_point), 2)

                        if dot_product >0:
                            if Case =="1":
                                factor = 1
                        elif dot_product < 0:
                            if Case == "2":
                                factor = 1
                        else:
                          
                            print("Error :normal vector is normal to the line joining the centroids of the edge and block" )
                            return None
                        normal = [factor*element for element in normal ] 
                        
                        normal_vectors[-1].extend(normal)
                        tang_vector = np.cross( np.array([0,1,0]) ,np.array(normal)  )

                        tang_vectors[-1].extend(list(tang_vector))

                        
                        for ind, element in enumerate([tang_vector[0], normal[0] ]):
                            if abs(element) > 1e-3:
                                fx_block_row += [3*i  ]
                                fx_block_col += [var_pos + ind] 
                                fx_block_val += [float(element)]


                        for ind, element in enumerate([tang_vector[2], normal[2] ]):
                            if abs(element) > 1e-3 != 0:
                                fz_block_row += [3*i + 1 ]
                                fz_block_col += [var_pos + ind ] 
                                fz_block_val += [float(element)]
      
                        
                        normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), np.array([element for element in mid_point])) 
                        
                        normal_moment_vectors[-1].extend(normal_moment)
                        
                        if abs(normal_moment[1]) > 1e-4:
                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos + 1 ] 
                            my_block_val += [float(normal_moment[1])]

                      
                        tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), np.array([element for element in mid_point])) 
                        
                        tang_moment_vectors[-1].extend(tang_moment)

                        if  abs(tang_moment[1]) > 1e-4:
                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos  ] 
                            my_block_val += [float(tang_moment[1])]
                        

                        my_block_row += [3*i + 2 ]
                        my_block_col += [var_pos + 2  ] 
                        my_block_val += [1]

                        edges_var_pos.append(var_pos)
                        #eq_matrix_val.append(i)
                        var_pos += 3
                       
                     

                    else:
                        if False:
                            joints_length.append([])
                            active_edges_pts_coor.append([])

                            points_coor = apply_permutation(points_coor, indices)

                            fx_dict = dict()
                            fz_dict = dict()
                            my_dict = dict()
                            #edges_var_pos.append([])
                            for point_ind, point_coor in enumerate(points_coor[:-1]):
                                
                                pt1 = point_coor
                                pt2 = points_coor[point_ind+1]

                                joint = rs.AddLine(pt1, pt2 )

                                joint_length = rs.CurveLength(joint)
                                joints_length[-1].append(joint_length)

                                normal = round_point(compute_normal_at_point( joint, point_coor ) , 2)
                        
                                dot_product = round( calculate_dot_product_rhino(normal, block_centroid_r, point_coor) , 2)
                            
                                if dot_product >0:
                                    if Case =="1":
                                        factor = 1
                                elif dot_product <0:
                                    if Case == "2":
                                        factor = 1
                                else:
                                    print("Error :normal vector is normal to the line joining the centroids of the edge and block" )
                                    return None
                                
                                normal = [factor*element for element in normal ] 
                                
                                normal_vectors[-1].append([normal[0], normal[1], normal[2]])
                        
                                tang_vector = np.cross( np.array([0,1,0]) , np.array(normal)  )

                                tang_vectors[-1].append(list(tang_vector))

                                for ind, element in enumerate([tang_vector[0], normal[0] ]):
                                    if element != 0:
                                        if var_pos + point_ind*2 + ind in fx_dict: 
                                            #eq_matrix_row += [3*i ]
                                            fx_dict[var_pos + point_ind*2 + ind] += joint_length*0.5*element*(weight13*2 + 0.5*weight2)
                                        else:
                                            fx_dict[var_pos + point_ind*2 + ind] = joint_length*0.5*element*(weight13*2 + 0.5*weight2)


                                for ind, element in enumerate([tang_vector[1], normal[1] ]):
                                    if element != 0:
                                        if var_pos + point_ind*2 + ind in fz_dict: 
                                            #eq_matrix_row += [3*i ]
                                            fz_dict[var_pos + point_ind*2 + ind] += joint_length*0.5*element*(weight13*2 + 0.5*weight2)
                                        else:
                                            fz_dict[var_pos + point_ind*2 + ind] = joint_length*0.5*element*(weight13*2 + 0.5*weight2)

                                pt1_coor = np.array([element for element in pt1])
                                pt3_coor = np.array([element for element in pt2])
                                pt2_coor = pt1_coor*0.5 + pt3_coor*0.5
                                active_edges_pts_coor[-1].append(list(pt1_coor))

                                pt1_normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), pt1_coor)[1]
                                pt2_normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), pt2_coor)[1]
                                pt3_normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), pt3_coor)[1]

                                normal_moment = joint_length*0.5*(pt1_normal_moment*weight13 + pt3_normal_moment*weight13 + pt2_normal_moment*weight2)

                                if normal_moment !=0:
                                    if var_pos + point_ind*2 + 1 in my_dict: 
                                        #eq_matrix_row += [3*i ]
                                        my_dict[var_pos + point_ind*2 + 1 ] += normal_moment
                                    else:
                                        my_dict[var_pos + point_ind*2 + 1 ] = normal_moment


                                pt1_tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), pt1_coor)[1]
                                pt2_tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), pt2_coor)[1]
                                pt3_tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), pt3_coor)[1]

                                tang_moment = joint_length*0.5*(pt1_tang_moment*weight13 + pt3_tang_moment*weight13 + pt2_tang_moment*weight2)

                                if tang_moment !=0:
                                    if var_pos + point_ind*2  in my_dict: 
                                        #eq_matrix_row += [3*i ]
                                        my_dict[var_pos + point_ind*2  ] += tang_moment
                                    else:
                                        my_dict[var_pos + point_ind*2  ] = tang_moment

                            active_edges_pts_coor[-1].append(list(pt3_coor))

                            for element in fx_dict:
                                eq_matrix_col.append(element)
                                eq_matrix_row.append(3*i)
                                eq_matrix_val.append(fx_dict[element])
                            
                            for element in fz_dict:
                                eq_matrix_col.append(element)
                                eq_matrix_row.append(3*i +1)
                                eq_matrix_val.append(fz_dict[element])

                            for element in my_dict:
                                eq_matrix_col.append(element)
                                eq_matrix_row.append(3*i +2 )
                                eq_matrix_val.append(my_dict[element])



                            edges_var_pos.extend([var_pos])
                    
                            var_pos += 2*int(len(points_coor))
                        
                        #points_coor = apply_permutation(points_coor, indices)
                        
                        edges_var_pos.append([])
                        
                        for point_ind, point_coor in enumerate(points_coor[:-1]):
                            
                            pt1 = point_coor
                            pt2 = points_coor[point_ind+1]
                            

                            pt1_coor = np.array([element for element in pt1])
                            pt3_coor = np.array([element for element in pt2])
                            mid_point = pt1_coor*0.5 + pt3_coor*0.5
                            active_edges_pts_coor.append([mid_point[0], mid_point[1], mid_point[2] ])

                            joint = rs.AddLine(pt1, pt2 )

                            joint_length = rs.CurveLength(joint)
                            joints_length.append(joint_length)

                            normal = round_point(compute_normal_at_point( joint, point_coor ) , 2)
                    
                            dot_product = calculate_dot_product_rhino(normal, block_centroid_r, point_coor) 
                        
                            if dot_product >0:
                                if Case =="1":
                                    factor = 1
                            elif dot_product <0:
                                if Case == "2":
                                    factor = 1
                            else:
                                
                                print("Error :normal vector is normal to the line joining the centroids of the edge and block" )
                                return None
                            
                            normal = [factor*element for element in normal ] 
                            
                            normal_vectors[-1].append([normal[0], normal[1], normal[2]])
                    
                            tang_vector = np.cross( np.array([0,1,0]) , np.array(normal)  )

                            tang_vectors[-1].append(list(tang_vector))

                            for ind, element in enumerate([tang_vector[0], normal[0] ]):
                                if abs(element) > 1e-3:
                                    fx_block_row += [3*i  ]
                                    fx_block_col += [var_pos + ind] 
                                    fx_block_val += [float(element)]


                            for ind, element in enumerate([tang_vector[2], normal[2] ]):
                                if abs(element) > 1e-3 != 0:
                                    fz_block_row += [3*i + 1 ]
                                    fz_block_col += [var_pos + ind ] 
                                    fz_block_val += [float(element)]

                            normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), np.array([element for element in mid_point])) 
                        
                            normal_moment_vectors[-1].extend(normal_moment)
                            
                            if abs(normal_moment[1]) > 1e-4:
                                my_block_row += [3*i + 2 ]
                                my_block_col += [var_pos + 1 ] 
                                my_block_val += [float(normal_moment[1])]

                        
                            tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), np.array([element for element in mid_point])) 
                            
                            tang_moment_vectors[-1].extend(tang_moment)

                            if  abs(tang_moment[1]) > 1e-4:
                                my_block_row += [3*i + 2 ]
                                my_block_col += [var_pos  ] 
                                my_block_val += [float(tang_moment[1])]
                            

                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos + 2  ] 
                            my_block_val += [1]

                            edges_var_pos[-1].append(var_pos)
                            #eq_matrix_val.append(i)
                            var_pos += 3

                    if len(points_pos) <= 2:
                        edges_points_format.append(points_pos)
                        edges.append(points_pos)
                        edge_pos = len(edges) - 1
                        active_edges.append(edge_pos)
                    else:
                        edges.append([])
                        edge_pos = len(edges) - 1
                        active_edges.append(edge_pos)
                        for index, point_pos in enumerate(points_pos[:-1]):
                            edges[-1].append([ point_pos, points_pos[index + 1] ])
                            
                        edges_points_format.append(points_pos)

                else:
                    
    

                    if is_line:
                        
                        edge_pos = edges.index(points_pos)
                        edge_active_pos = active_edges.index(edge_pos) 
                        var_pos_found =  edges_var_pos[edge_active_pos]

                        tang_vector = [-element for element in tang_vectors[edge_active_pos] ]
                        normal = [ -element for element in normal_vectors[edge_active_pos] ]


                        for ind, element in enumerate([tang_vector[0], normal[0] ]):
                            if abs(element) > 1e-3:
                                fx_block_row += [3*i  ]
                                fx_block_col += [var_pos_found + ind] 
                                fx_block_val += [float(element)]


                        for ind, element in enumerate([tang_vector[2], normal[2] ]):
                            if abs(element) > 1e-3:
                                fz_block_row += [3*i + 1 ]
                                fz_block_col += [var_pos_found + ind ] 
                                fz_block_val += [float(element)]
                        
               
                        
                        normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), np.array([element for element in mid_point])) 
                        
                        if abs(normal_moment[1]) > 1e-5:
                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos_found + 1 ] 
                            my_block_val += [float(normal_moment[1])]

                      
                        tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), np.array([element for element in mid_point])) 
                        
                        if abs(tang_moment[1]) > 1e-5:
                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos_found  ] 
                            my_block_val += [float(tang_moment[1])]
                        

                        my_block_row += [3*i + 2 ]
                        my_block_col += [var_pos_found + 2  ] 
                        my_block_val += [1]

                      

                    else:
                 

                        edge_pos = edges_points_format.index(points_pos)
               
                        edge_active_pos = active_edges.index(edge_pos)

                        
                    

                        for point_ind, point_coor in enumerate(points_coor[:-1]):
                      
                           
                            #edge_active_pos = active_edges.index(edge_pos) 
                            var_pos_found =  edges_var_pos[edge_active_pos][point_ind] 

                            tang_vector = [-element for element in tang_vectors[edge_active_pos][point_ind] ]
                            normal = [ -element for element in normal_vectors[edge_active_pos][point_ind] ]
                    
                   

                            for ind, element in enumerate([tang_vector[0], normal[0] ]):
                                if abs(element) > 1e-3:
                                    fx_block_row += [3*i  ]
                                    fx_block_col += [var_pos_found + ind] 
                                    fx_block_val += [float(element)]


                            for ind, element in enumerate([tang_vector[2], normal[2] ]):
                                if abs(element) > 1e-3 != 0:
                                    fz_block_row += [3*i + 1 ]
                                    fz_block_col += [var_pos_found + ind ] 
                                    fz_block_val += [float(element)]

                            normal_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(normal), np.array([element for element in mid_point])) 
                        
                            normal_moment_vectors[-1].extend(normal_moment)
                            
                            if abs(normal_moment[1]) > 1e-4:
                                my_block_row += [3*i + 2 ]
                                my_block_col += [var_pos_found + 1 ] 
                                my_block_val += [float(normal_moment[1])]

                        
                            tang_moment = calculate_moment(np.array([element for element in block_centroid_r]), np.array(tang_vector), np.array([element for element in mid_point])) 
                            
                            tang_moment_vectors[-1].extend(tang_moment)

                            if  abs(tang_moment[1]) > 1e-4:
                                my_block_row += [3*i + 2 ]
                                my_block_col += [var_pos_found  ] 
                                my_block_val += [float(tang_moment[1])]
                            

                            my_block_row += [3*i + 2 ]
                            my_block_col += [var_pos_found + 2  ] 
                            my_block_val += [1]


                
                edges_pos += [edge_pos]
                #rs.DeleteObjects(curve)
                

            #rs.DeleteObject(contact)  
        
        if 1 in b_att:
            fx_block_row += [3*i   ]
            fx_block_col += [-1] 
            fx_block_val += [ float(b_att[1]) ]

            pos_to_modify.append(  len(eq_matrix_row) + len(fx_block_row) - 1 )
        
        if 2 in b_att:
            fz_block_row += [3*i + 1  ]
            fz_block_col += [-1] 
            fz_block_val += [ float(b_att[2]) ]

            pos_to_modify.append(  len(eq_matrix_row) + len(fx_block_row) + len(fz_block_row) - 1 )

        if 3 in b_att:
            my_block_row += [3*i + 2  ]
            my_block_col += [-1] 
            my_block_val += [ float(b_att[3]) ]
            pos_to_modify.append(  len(eq_matrix_row) + len(fx_block_row) + len(fz_block_row) + len(my_block_row) - 1 )

        eq_matrix_col += fx_block_col + fz_block_col + my_block_col
        eq_matrix_row += fx_block_row + fz_block_row + my_block_row
        eq_matrix_val += fx_block_val + fz_block_val + my_block_val

        blocks.append(edges_pos)
        

        edges_pos = []
        for non_contact in block_egdes_guid["non_contact"]:
            curves = rs.ExplodeCurves(non_contact)
            for curve in curves:
                
                points_pos = []
                
                for point in rs.CurvePoints(curve):
                    point = round_point( point,3)
                    points, point_pos = check_add_point(points, [point[0], point[1], point[2] ] )
                    points_pos += [point_pos]
                
                
                points_pos = sorted(points_pos)
                
                if points_pos not in edges:
                    edges.append(points_pos)
                    edges_points_format.append(points_pos)
                    edge_pos = len(edges) - 1
                else:
                    edge_pos = edges.index(points_pos)
                edges_pos += [edge_pos]
                #rs.DeleteObjects(curve)
                

                
                #rs.DeleteObjects(non_contact)
            

        blocks.append(edges_pos)
        
        # Get bounding boxes
       
        if False:
            bounding_boxes_block_curves = [bounding_boxes[i]] + bounding_boxes_curves

            # Create an R-tree
            rtree1 = Rhino.Geometry.RTree()

            # Add bounding boxes to the R-tree
            for j, bbox in enumerate(bounding_boxes_block_curves):
                rtree1.Insert(bbox, j)

            block_curve_inter = []

            blocki_bbox = surfaces_boundary[i]

         
            # Define a callback function for the R-tree search
            def search_callback1(sender1, e ):
                
                if e.Id != 0:
              
                    inter = rs.CurveCurveIntersection(   curve_geometries[e.Id - 1 ] , surfaces_boundary[i]   )
                    if inter != None:
                        for curve_inter in inter:
                            if curve_inter[0] == 2:
                        
                                pt1 = curve_inter[1]
                                pt1 = [ round(element,3) for element in pt1 ] 
                                
                                pt2 = curve_inter[2]
                                pt2 = [ round(element,3) for element in pt2 ] 
                        
                                pt1_found, point1_pos = check_point(points, pt1)

                                pt2_found, point2_pos = check_point(points, pt2)
                        
                                if pt1_found == False or pt2_found == False:
                                    print("Error: no intersection between a curve and the blocks")
                                    return None
                                elif point1_pos == point2_pos:
                                    print("Error: all block's edges are blocked")
                                    return None
                                else:
                                    block_curve_inter.append([min(point1_pos,point2_pos), max(point1_pos,point2_pos)])
            
            
            rtree1.Search(bounding_boxes_block_curves[0] , search_callback1)

            rs.DeleteObjects(blocki_bbox)
       
            for element in block_curve_inter:
                for non_contact_edge in  blocks[-1]:
                    
                    edge = edges[non_contact_edge]

                    if edge[0] == element[0] and edge[-1] == element[-1]:
                        blocks[-1].remove(non_contact_edge)
                        blocks[-2].append(non_contact_edge)

    # Get the current directory in Rhino
    current_directory = rs.DocumentPath()

    
    n = max(eq_matrix_col)
    n = n + 1

    for pos in pos_to_modify:
        eq_matrix_col[pos] = n

    print("blockssssss", blocks)

    print(edges)


    rs.DeleteObjects(surfaces_boundary)
    # Define the file path
    file_path = os.path.join(current_directory, 'matrices.txt')
  
    # Write the lists to the file
    with open(file_path, 'w') as file:
        file.write(' '.join(map(str, eq_matrix_val)) + '\n')
        file.write('##' + '\n' )  # Empty line as separator
        file.write(' '.join(map(str, eq_matrix_row)) + '\n')
        file.write('##' + '\n')  # Empty line as separator
        file.write(' '.join(map(str, eq_matrix_col)) + '\n')
        file.write('##' + '\n')  # Empty line as separator
        file.write(' '.join(map(str, joints_length)) + '\n')
        file.write('##' + '\n')  # Empty line as separator
        file.write(str(nb) + '\n')
        file.write('##' + '\n')  # Empty line as separator
        file.write(str(n+1) + '\n')
        file.write('##' + '\n')  # Empty line as separator

        for sublist in blocks:
            # Convert each sublist to a string with elements separated by commas
            line = ' '.join(map(str, sublist))
            # Write the line to the file followed by a newline character
            file.write(line + '\n')
          
    
        file.write('##' + '\n')  # Empty line as separator

        for sublist in edges:
            if isinstance(sublist[0], list):  # Check if the first element is a list (nested)
                flattened = [str(item) for inner_list in sublist for item in inner_list]
                line = ' '.join(flattened)
            else:
                line = ' '.join(map(str, sublist))
            file.write(line + '\n')
        file.write('##' + '\n')  # Empty line as separator

        for sublist in points:
            # Convert each sublist to a string with elements separated by commas
            line = ' '.join(map(str, sublist))
            # Write the line to the file followed by a newline character
            file.write(line + '\n')
        file.write('##' + '\n')  # Empty line as separator
        print("centroid", blocks_centroid)
        for sublist in blocks_centroid:
            # Convert each sublist to a string with elements separated by commas
            line = ' '.join(map(str, sublist))
            # Write the line to the file followed by a newline character
            file.write(line + '\n')
        file.write('##' + '\n')  # Empty line as separator


        # Use a lambda to capture i1 and pass it to the callback
  #      if i1 == 0:
 #           rtree1.Search(bbox, lambda sender1, e: search_callback(sender1, e, i1, block_curve_inter ))

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

def sort_and_get_indices(vector):
    """Sorts the vector and returns the sorted vector along with the original indices in sorted order."""
    sorted_vector_with_indices = sorted(enumerate(vector), key=lambda x: x[1])
    sorted_vector = [x[1] for x in sorted_vector_with_indices]
    indices = [x[0] for x in sorted_vector_with_indices]
    return sorted_vector, indices

def apply_permutation(vector, indices):
    """Reorganizes the vector based on the given indices."""
    return [vector[i] for i in indices]

def calculate_dot_product_rhino(normal, pt1, pt2):
    # Ensure the points are valid
    point1 = rs.coerce3dpoint(pt1)
    point2 = rs.coerce3dpoint(pt2)
    
    if not point1 or not point2:
        print("Invalid points")
        return None

    # Define the vector from pt1 to pt2
    vector = point2 - point1

    # Ensure the normal vector is valid
    normal_vector = rs.coerce3dvector(normal)
    if not normal_vector:
        print("Invalid normal vector")
        return None

    # Calculate the dot product
    dot_product = normal_vector*vector

    return dot_product


   
def compute_normal_at_point(curve_id, point):
    # Ensure the input curve is valid
    curve = rs.coercecurve(curve_id)
    if not curve:
        print("Invalid curve")
        return  None

    # Find the parameter on the curve closest to the given point
    param = rs.CurveClosestPoint(curve_id, point)
    if param is None:
        print("Failed to find the parameter on the curve closest to the point.")
        return None, None

    # Get the point on the curve at the found parameter
    point_on_curve = rs.EvaluateCurve(curve_id, param)

    # Compute the tangent vector at the given parameter
    tangent = rs.CurveTangent(curve_id, param)
    if not tangent:
        print("Failed to compute tangent vector.")
        return  None

    # Compute the normal vector (assuming the curve lies in the XY plane)
    normal = Rhino.Geometry.Vector3d(-tangent.Z, 0, tangent.X)

    # Normalize the normal vector
    normal.Unitize()
    return  normal

def is_point_inside_curve(point, curve_id):
    # Ensure the input curve is valid
    curve = rs.coercecurve(curve_id)
    if not curve:
        print("Invalid curve")
        return None

    # Ensure the curve is closed
    if not rs.IsCurveClosed(curve_id):
        print("The curve is not closed")
        return None

    # Check if the point is inside the closed curve
    result = rs.PointInPlanarClosedCurve(point, curve_id)

    return result
           
def find_intersection_line_curve(line_start, line_end, curve_id):
    # Create the line from the two points
    line_id = rs.AddLine(line_start, line_end)
    if not line_id:
        print("Failed to create the line.")
        return None

    # Find the intersection between the line and the curve
    intersections = rs.CurveCurveIntersection(line_id, curve_id)

    # Extract intersection points
    intersection_type = "points"
    intersection_points = []
    if intersections:
        for intersection in intersections:
            if intersection[0] == 1:  # Point intersection
                intersection_points.append(intersection[1])
            elif intersection[0] == 2:  # Overlap
                print(intersection)
                intersection_type = "curve"
                intersection_points.append(intersection[1])
                intersection_points.append(intersection[2])

    # Delete the temporary line
    rs.DeleteObject(line_id)

    return intersection_type,intersection_points            


def check_add_point(array, point):
    # Check if the point already exists in the array
    point_exists = [a == point for a in array]
    
    if any(point_exists):
        # If the point exists, find its index
        position = point_exists.index(True)
        return array, position
    else:
        # If the point does not exist, add it to the array
        array.append(point)
        # Return the updated array and the new point's position
        return array, len(array) - 1

def check_point(array, point):
    # Check if the point already exists in the array
    point_exists = [a == point for a in array]
    
    if any(point_exists):
        # If the point exists, find its index
        position = point_exists.index(True)
        return True, position
    else:
        # If the point does not exist, add it to the array
        return False, len(array) - 1
   
     
if __name__ == "__main__":
    main()

