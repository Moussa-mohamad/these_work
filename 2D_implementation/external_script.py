import os
from amplpy import AMPL
import numpy as np
import vtk
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad, VtkVertex, VtkPolyLine, VtkPolygon

# Helper function to read a section until an empty line is encountered
def read_until_empty_line(start_index, lines):
    section = []
    index = start_index
    while index < len(lines):
        line = lines[index].strip()
        if line == '##':
            break
        section.append(line)
        index += 1
    return section, index
    

def rotate_y(theta):
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])


def transform_point_2d(point, displacement, angle):
    
    rotation_matrix = rotate_y(angle)

    point = np.dot(rotation_matrix, point)
    # Apply displacement
    point = np.array(point) + np.array(displacement)
    return point.tolist()



script_dir = os.path.dirname(__file__)
file_path = os.path.join(script_dir, 'matrices.txt')

# Read the lists from the file
with open(file_path, 'r') as file:
    lines = file.readlines()

section, index = read_until_empty_line(0, lines)
matrices_coef = [list(map(float, line.split(' '))) for line in section][0]

# Read the second list
section, index = read_until_empty_line(index + 1, lines)
matrices_row = [list(map(int, line.split(' '))) for line in section][0]


section, index = read_until_empty_line(index + 1, lines)
matrices_col = [list(map(int, line.split(' '))) for line in section][0]

section, index = read_until_empty_line(index + 1, lines)
joints_length = [list(map(float, line.split(' '))) for line in section][0]
print(joints_length)

section, index = read_until_empty_line(index + 1, lines)

if section:
    nb = int(section[0])

section, index = read_until_empty_line(index + 1, lines)

if section:
    n = int(section[0])

section, index = read_until_empty_line(index + 1, lines)

for line in section:
    print(line)
blocks = [list(map(int, line.split(' '))) if line != '' else [] for line in section]

section, index = read_until_empty_line(index + 1, lines)
edges = [list(map(int, line.split(' '))) for line in section]

section, index = read_until_empty_line(index + 1, lines)
points = [list(map(float, line.split(' '))) for line in section]

section, index = read_until_empty_line(index + 1, lines)
blocks_centroid = [list(map(float, line.split(' '))) for line in section]
blocks_centroid = np.array(blocks_centroid)
print("List 1:", matrices_coef)
print("List 2:", matrices_row)
print("List 3:", matrices_col)
print("List 4:", joints_length)
print("List 5:", nb)
print("List 6:", n)
print("List 7:", blocks)
print("List 7:", edges)
print("list9:", points)
print("list9:", blocks_centroid)

for i in range(40,41,2):
    # Create an AMPL instance
    ampl = AMPL()
   
    # Set the solver to use
    solver = "highs" 
    ampl.set_option("solver", "ipopt")


    # Define the AMPL model with placeholders for sparse matrix
    ampl.eval(r"""
    param n;
    param nc;
    set Aline;
    param AColInd {Aline} >= 1;
    param ALineInd {Aline} >= 1;
    param AVal {Aline};
    param joints_length{1..nc};
    

    set Bline;
    param Bval{Bline};

    var x {1..n};
    param nb;

    """)

   


    # Convert to 1-based indexing for AMPL
    
    matrices_row = [index + 1 for index in matrices_row]
    #ampl_matrices_col = [index + 1 for index in matrices_col]

    # Assign data to AMPL parameters
    ampl.param['nb'] = nb
    ampl.param['n'] = n
    ampl.param['nc'] = int((n-1)/3)
    ampl.set['Bline'] = range(1, 3*nb + 1)
    print(joints_length)
    print({i+1: joints_length[i] for i in range(len(joints_length))})
    ampl.param['joints_length'] = {i+1: joints_length[i] for i in range(len(joints_length))}
    
    ampl.param['Bval'] = {1:0,2:1,3:0, 4:0,5:1,6:0 , 7:0,8:1,9:0}
    ampl.set['Aline'] = range(1, len(matrices_coef) + 1)
    ampl.param['AColInd'] = {i+1: matrices_col[i] for i in range(len(matrices_col))}
    ampl.param['ALineInd'] = {i+1: matrices_row[i] for i in range(len(matrices_row))}
    ampl.param['AVal'] = {i+1: matrices_coef[i] for i in range(len(matrices_coef))}


    ampl.eval(r"""

    param A_lnum:= max {i in Aline} ALineInd[i];
    param ind ;
    param initind;
    param allinitind{1..A_lnum};
    param allind{1..A_lnum};
    let ind :=1;
    
    var alpha1{1..nc};
    var alpha2{1..nc};
    
    param fc := 200;
    param ft := 0;
    
  

    display A_lnum;
    for {i in 1..A_lnum} {
                let initind := ind; 
                repeat while ALineInd[ind] == i   { if ind+1 > card(Aline) then break; else {let ind := ind +1;} 
                  } 
                  if ind+1 > card(Aline) then let allind[i] := ind; else let allind[i] := ind-1;
                  let allinitind[i] := initind;

                  }
    
    for { i in Aline}
        display AColInd[i];
        
    subject to equalities {i in 1..3*nb}:
    sum {j in allinitind[i]..allind[i]} AVal[j]*x[AColInd[j]+1] = Bval[i]; 
    
    subject to cap_ten {i in 1..nc}:
         x[3*i -1] >= 0;
    
    subject to sliding_top {i in 1..nc}:
         x[3*i -2] <= x[3*i -1]*0.64; 
    
    subject to sliding_bot {i in 1..nc}:
         x[3*i -2] >= -x[3*i -1]*0.64;

subject to rocking_top {i in 1..nc}:
         x[3*i] <= x[3*i -1]*joints_length[i]/2;

subject to rocking_bot {i in 1..nc}:
         x[3*i] >= -x[3*i -1]*joints_length[i]/2;        
    
      subject to crushing_bot {i in 1..nc}:
         -x[3*i]*joints_length[i]/2 - x[3*i-1]*((fc-ft)*joints_length[i]/(2*(ft+fc)) -x[3*i-1]/(2*(fc+ft) )) - 0.5*fc*ft*joints_length[i]*joints_length[i]/(fc+ft) + alpha2[i] = 0; 
     
     subject to crushing_top {i in 1..nc}:
         x[3*i]*joints_length[i]/2 - x[3*i-1]*((fc-ft)*joints_length[i]/(2*(ft+fc)) -x[3*i-1]/(2*(fc+ft) )) - 0.5*fc*ft*joints_length[i]*joints_length[i]/(fc+ft) + alpha1[i] = 0; 
     
     
     
    subject to inequalities3 {i in 1..nc}:
            alpha1[i] >= 0;

            
    subject to inequalities4 {i in 1..nc}:
    alpha2[i] >= 0;
    
    maximize obj:
        x[n];


    """)




    # Solve the problem
    ampl.solve()

    # Get the results
    x = ampl.get_variable('x').get_values().to_pandas()
    obj_value = ampl.get_objective('obj').value()

    # Print the results
    print("Solution:")
    print(x)
    print(f"Objective value: {obj_value}")
    
    
    constraint1 = ampl.get_constraint("equalities")
    constraint1_values = constraint1.getValues()
    ub = np.empty(3*nb)
    print(nb)
    ind = 0
    for i in constraint1_values:
        print(i)
        ub[ind] = i[1]
        ind = ind+1
    print(ub)







ub = [ element for element in ub ]
new_points = []
initial_points = np.empty((0,3))
points_disp = np.empty((0,3))

FILE_PATH =  "C:\\Users\\mmoussa\\Desktop\\rhino_test\\2D_implementation\\stat" +  "_disp"

# Define connectivity or vertices that belong to each element
conn = []
# Define offset of the last vertex of each element
offset = []
# Define cell types
ctype = []

for ind, block_centroid in enumerate(blocks_centroid):

  
    edges_ind = blocks[2*ind] + blocks[2*ind + 1]
     
    for edge_ind in edges_ind:
        edge_pts_ind = edges[edge_ind]
        
        points_coor = []
        #points_ind = [ int(element) -1 for element in points_ind] 
        for point_ind in edge_pts_ind:
            points_coor.append(points[point_ind])
            #points_coor.astype(float)
            #points_coor = np.array([[float(element) for element in row] for row in points_coor])
        

        all_new = []
        
        reference_point = block_centroid
   
        
        displacement_vector = [ ub[3*ind], 0, ub[3*ind +1] ]
  
        
        angle = ub[3*ind + 2]
        
        for element in points_coor:
         
            initial_points = np.append(initial_points, np.array([element]), axis = 0)
            

            conn.extend( [initial_points.shape[0] - 1 ])
            
            point_new_ref = [element[element1] - reference_point[element1] for element1 in range(3) ]
           
            
            new_point = transform_point_2d(point_new_ref , displacement_vector, angle)
            
           
            new_point = np.array([new_point[element] + reference_point[element] for element in range(3) ])
            
            all_new.append(new_point)

            points_disp = np.append(points_disp, np.array( [new_point[:] - np.array(element)[:] ] ), axis = 0)
         
            
        ctype.extend ([VtkPolyLine.tid])
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
