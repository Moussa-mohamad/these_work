import rhinoscriptsyntax as rs
import Rhino
import scriptcontext as sc
curve = []
curve = rs.ObjectsByType(rs.filter.curve, True)

curve_split =  rs.DivideCurve(curve, 100)
print(len(curve_split))

convex = curve[0]
#concave = curve[1]

points = rs.CurvePoints(convex)   

print(points)

import rhinoscriptsyntax as rs
import Rhino

def compute_normal_2d_curve(curve_id, param):
    # Get the point on the curve at the given parameter
    point = rs.EvaluateCurve(curve_id, param)
    
    # Compute the tangent vector at the given parameter
    tangent = rs.CurveTangent(curve_id, param)
    
    if not tangent:
        print("Failed to compute tangent vector.")
        return None
    
    # Compute the normal vector as a perpendicular vector to the tangent
    normal = Rhino.Geometry.Vector3d(-tangent.Y, tangent.X, 0)
    
    # Normalize the normal vector
    normal.Unitize()
    
    return point, normal

# Example usage
curve_id = rs.GetObject("Select a 2D curve", rs.filter.curve)
if curve_id:
    param = rs.GetReal("Enter the parameter along the curve", 0.5)
    if param is not None:
        point, normal = compute_normal_2d_curve(curve_id, param)
        if normal:
            print(f"Normal vector at parameter {param}: {normal}")
            rs.AddLine(point, point + normal)  # Visualize the normal vector
        else:
            print("Failed to compute normal vector.")
    else:
        print("Invalid parameter.")
else:
    print("No curve selected.")


import gmsh
import numpy as np

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("NURBS Curve")

# Example control points from Rhino (replace with your actual control points)
control_points = []

for point in points:
    control_points.append([point[0], point[1], point[2]])

gmsh_points = []
for i, cp in enumerate(control_points):
    gmsh_points.append(gmsh.model.occ.addPoint(cp[0], cp[1], cp[2]))

# Create a B-spline curve using the control points
gmsh_curve = gmsh.model.occ.addBSpline(gmsh_points)

# Synchronize the model
gmsh.model.occ.synchronize()

# Add physical groups to ensure they are meshed
curve_tag = gmsh.model.occ.getEntities(dim=1)[0][1]
gmsh.model.addPhysicalGroup(1, [curve_tag], tag=-1)
gmsh.model.setPhysicalName(1, curve_tag, "NURBS Curve")

# Set meshing options to use curvature
# Enabling curvature-based meshing with a scale factor
gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)


# Generate the mesh
gmsh.model.mesh.generate(1)

# Save the mesh to a file
gmsh.write("C:\\Users\\mmoussa\\Desktop\\rhino_test\\2D_implementation\\nurbs_curve.msh")

# Finalize Gmsh
gmsh.finalize()
