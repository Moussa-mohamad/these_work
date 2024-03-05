#include "meshKernel.hxx"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <gmsh.h_cwrap>
#include <chrono>
#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>






std::pair<std::vector<int>, std::vector<int>> extractNonZeroWithIndices(const std::vector<int>& source) {
    std::vector<int> nonZeroElements;
    std::vector<int> indices;
    for (size_t i = 0; i < source.size(); ++i) {
        if (source[i] != 0) {
            nonZeroElements.push_back(source[i]);
            indices.push_back(i);
        }
    }
    return std::make_pair(nonZeroElements, indices);
}

Eigen::Vector3d calculate_moment(const Eigen::Vector3d& point_coords,
    const Eigen::Vector3d& force_coords,
    const Eigen::Vector3d& applied_point_coords) {
    Eigen::Vector3d r = -point_coords + applied_point_coords;
    Eigen::Vector3d F = force_coords;

    // Calculate the cross product to get the moment
    Eigen::Vector3d moment = r.cross(F);

    return moment;
}

double calculate_triangle_area(const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2,
    const Eigen::Vector3d& v3) {
    // Calculate two vectors representing two sides of the triangle
    Eigen::Vector3d side1 = v2 - v1;
    Eigen::Vector3d side2 = v3 - v1;
    
    // Calculate the cross product of the two sides (normal to the triangle)
    Eigen::Vector3d normal = side1.cross(side2);
    
    // The magnitude of the cross product is twice the area of the triangle
    double area = 0.5 * normal.norm();

    return area;
}

void barycentric_coordinates(Eigen::Vector3d pt, Eigen::Matrix3d& triangle_coords,
    double* alpha, double* beta, double* gamma) {
   
    
    // Extract triangle vertices coordinates
    
    Eigen::Vector3d pt1 = triangle_coords.row(0);
    Eigen::Vector3d pt2 =  triangle_coords.row(1);
    Eigen::Vector3d pt3 = triangle_coords.row(2);
   

    double area_total = calculate_triangle_area(pt1, pt2, pt3);
    double area1 = calculate_triangle_area(pt, pt2, pt3);
    double area2 = calculate_triangle_area(pt, pt1, pt3);
    double area3 = calculate_triangle_area(pt, pt1, pt2);
    //std::cout << area_total << "," << area1 << "," << area2  << "," << area3 << std::endl;
    // Calculate barycentric coordinates
    *alpha = area1 / area_total;
    *beta = area2 / area_total;
    *gamma = area3 / area_total;
}

std::tuple<std::vector<double>, std::vector<double>> integral_gauss(Eigen::Matrix3d& vertices,
    Eigen::Matrix3d& local_ref,
    const Eigen::Vector3d& block_centroid,
    const std::vector<double>& rhos) {
    // Define Gauss points positions in normalized coordinates
    
    double s1 = 1.0 / 6.0;
    double t1 = 1.0 / 6.0;
    
    double s2 = 2.0 / 3.0;
    double t2 = 1.0 / 6.0;
    
    double s3 = 1.0 / 6.0;
    double t3 = 2.0 / 3.0;

    Eigen::Vector3d x_local = local_ref.row(0);
    Eigen::Vector3d y_local = local_ref.row(1);
    Eigen::Vector3d normal = local_ref.row(2);
    

    


    // Define first Gauss point
    Eigen::Vector3d pt1(s1 * (vertices(1, 0) - vertices(0, 0)) + t1 * (vertices(2, 0) - vertices(0, 0)) + vertices(0, 0),
        s1 * (vertices(1, 1) - vertices(0, 1)) + t1 * (vertices(2, 1) - vertices(0, 1)) + vertices(0, 1),
        s1 * (vertices(1, 2) - vertices(0, 2)) + t1 * (vertices(2, 2) - vertices(0, 2)) + vertices(0, 2) );

    Eigen::Vector3d pt2(s2 * (vertices(1, 0) - vertices(0, 0)) + t2 * (vertices(2, 0) - vertices(0, 0)) + vertices(0, 0),
        s2 * (vertices(1, 1) - vertices(0, 1)) + t2 * (vertices(2, 1) - vertices(0, 1)) + vertices(0, 1),
        s2 * (vertices(1, 2) - vertices(0, 2)) + t2 * (vertices(2, 2) - vertices(0, 2)) + vertices(0, 2) );

    Eigen::Vector3d pt3( s3 * (vertices(1, 0) - vertices(0, 0)) + t3 * (vertices(2, 0) - vertices(0, 0)) + vertices(0, 0),
        s3 * (vertices(1, 1) - vertices(0, 1)) + t3 * (vertices(2, 1) - vertices(0, 1)) + vertices(0, 1),
        s3 * (vertices(1, 2) - vertices(0, 2)) + t3 * (vertices(2, 2) - vertices(0, 2)) + vertices(0, 2) );
    
    // Calculate barycentric coordinates for each Gauss point
    //double (*coefs)[3][3] = {};
    std::array<std::array<double, 3>, 3> coefs;

    barycentric_coordinates(pt1, vertices, &coefs[0][0], &coefs[1][0], &coefs[2][0]);
    
    
    barycentric_coordinates(pt2, vertices, &coefs[0][1], &coefs[1][1], &coefs[2][1]);


    barycentric_coordinates(pt3, vertices, &coefs[0][2], &coefs[1][2], &coefs[2][2]);



    
    // Calculate triangle area
    double triangle_area = calculate_triangle_area(vertices.row(0), vertices.row(1), vertices.row(2));
    double jac = 2 * triangle_area; //jacobian of the transformation
    // Define weight of each Gauss vertex
    double w = 1.0 / 6.0;

    // Calculate load terms
    double load_term1 = 2 * triangle_area * w * (coefs[0][0] + coefs[0][1] + coefs[0][2]);
    double load_term2 = 2 * triangle_area * w * (coefs[1][0] + coefs[1][1] + coefs[1][2]);
    double load_term3 = 2 * triangle_area * w * (coefs[2][0] + coefs[2][1] + coefs[2][2]);

    std::vector<double> load_terms = { load_term1, load_term2, load_term3 };

  
 
    // Define functions to calculate moments at each Gauss point
    auto calculate_pt_moment = [&](const Eigen::Vector3d& pt) {
        return calculate_moment(block_centroid, x_local, pt),
            calculate_moment(block_centroid, y_local, pt),
            calculate_moment(block_centroid, normal, pt);
        };


    std::vector<double> moment_terms(27, 0.0); // Initialize the vector with zeros
    std::vector<int> non_zero_indices(27, 0.0); // Vector to store indices of non-zero elements

    double moments[3][3][3] = { { { calculate_moment(block_centroid, x_local, pt1)[0], calculate_moment(block_centroid, x_local, pt1)[1], calculate_moment(block_centroid, x_local, pt1)[2]},
                                  { calculate_moment(block_centroid, x_local, pt2)[0], calculate_moment(block_centroid, x_local, pt2)[1], calculate_moment(block_centroid, x_local, pt2)[2] },
                                  { calculate_moment(block_centroid, x_local, pt3)[0], calculate_moment(block_centroid, x_local, pt3)[1], calculate_moment(block_centroid, x_local, pt3)[2] } },
                                { { calculate_moment(block_centroid, y_local, pt1)[0], calculate_moment(block_centroid, y_local, pt1)[1], calculate_moment(block_centroid, y_local, pt1)[2]},
                                  { calculate_moment(block_centroid, y_local, pt2)[0], calculate_moment(block_centroid, y_local, pt2)[1], calculate_moment(block_centroid, y_local, pt2)[2] },
                                  { calculate_moment(block_centroid, y_local, pt3)[0], calculate_moment(block_centroid, y_local, pt3)[1], calculate_moment(block_centroid, y_local, pt3)[2] }},
                                { { calculate_moment(block_centroid, normal, pt1)[0], calculate_moment(block_centroid, normal, pt1)[1], calculate_moment(block_centroid, normal, pt1)[2]},
                                  { calculate_moment(block_centroid, normal, pt2)[0], calculate_moment(block_centroid, normal, pt2)[1], calculate_moment(block_centroid, normal, pt2)[2] },
                                  { calculate_moment(block_centroid, normal, pt3)[0], calculate_moment(block_centroid, normal, pt3)[1], calculate_moment(block_centroid, normal, pt3)[2] } } };



    //Eigen::Vector3d pt1_x_moment, pt1_y_moment, pt1_z_moment;
    //// Calculate moments at each Gauss point
    //pt1_x_moment = calculate_moment(block_centroid, x_local, pt1);
    //pt1_y_moment = calculate_moment(block_centroid, y_local, pt1);
    //pt1_z_moment = calculate_moment(block_centroid, normal, pt1);

    //Eigen::Vector3d pt2_x_moment, pt2_y_moment, pt2_z_moment;
    //// Calculate moments at each Gauss point
    //pt2_x_moment = calculate_moment(block_centroid, x_local, pt2);
    //pt2_y_moment = calculate_moment(block_centroid, y_local, pt2);
    //pt2_z_moment = calculate_moment(block_centroid, normal, pt2);

    //Eigen::Vector3d pt3_x_moment, pt3_y_moment, pt3_z_moment;

    //pt3_x_moment = calculate_moment(block_centroid, x_local, pt3);
    //pt3_y_moment = calculate_moment(block_centroid, y_local, pt3);
    //pt3_z_moment = calculate_moment(block_centroid, normal, pt3);
    //
    //int index = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                double value = jac * w * (coefs[k][0] * moments[i][0][j] + coefs[k][1] * moments[i][1][j] + coefs[k][2] * moments[i][2][j]);
                //if (abs(value) > 1e-5 ) {
                moment_terms[9 * k + 3 * j + i] = value;
                    //non_zero_indices[index] = 9 * k + 3 * j + i; // Store the index of the non-zero element
                    //index++; 
               // }
               
            }
        }
    }


    return std::make_tuple(load_terms, moment_terms);
}
struct Result {
    double value1;
    double value2;
};

// Template function to print vectors of any type
template<typename T>
void printVector(const std::vector<T>& vec) {
    std::cout << "Vector elements:" << std::endl;
    for (const auto& element : vec) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}


Output print_hello_c(double* blocks, int brows, int bcols, int blocks_num, double* nodes, int* active_faces, int active_faces_num, int nrows, int ncols, int* faces_FEpts, int* Faces_nodes, int faces_num, double* blocks_centroid, double* c_local_ref, double* lc_cn) {
    // Modify the value of the first element of the array

    std::unordered_map<int, int> active_faces_map;
    std::unordered_map<int, int> face_pts_map;

    std::unordered_map<int, double> lc_active_faces_map;
    

    for (int i = 0; i < active_faces_num; ++i)
        active_faces_map.insert({ active_faces[i],i + 1 });

    for (int i = 0; i < active_faces_num; ++i)
        lc_active_faces_map.insert({ active_faces[i], lc_cn[i] });




    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();



    gmsh::initialize();


    //std::vector<int> indices 
    std::vector<size_t> nodes_tags;

    // Create lines and get tags
    std::vector<int> line_tags;
    // Create points and get tags
    std::vector<size_t> node_tags;
    //double lc = 0.1; // Mesh element size

    for (int j = 0; j < faces_num; j += 1)
    {
        double lc = 1;

        if (active_faces_map.find(j + 1) != active_faces_map.end())
            lc = lc_active_faces_map[j + 1];
            


        for (size_t i = faces_FEpts[j]; i < faces_FEpts[j + 1]; i += 1) {

            int face_node = Faces_nodes[i] - 1;

            int tag = gmsh::model::geo::addPoint(nodes[3 * face_node], nodes[3 * face_node + 1], nodes[3 * face_node + 2], lc);

            node_tags.push_back(tag);

            nodes_tags.push_back(tag);
            //std::cout << tag << std::endl;
        }



        for (int i = 0; i < faces_FEpts[j + 1] - faces_FEpts[j]; i += 1) {



            int start = node_tags[i];

            int end = node_tags[(i + 1) % (faces_FEpts[j + 1] - faces_FEpts[j])]; // To close the loop

            int tag = gmsh::model::geo::addLine(start, end);
            line_tags.push_back(tag);


        }

        // Create curve loop


        int curve_loop_tag = gmsh::model::geo::addCurveLoop(line_tags);

        // Create plane surface

        int surface_tag = gmsh::model::geo::addPlaneSurface({ curve_loop_tag });

        node_tags.clear(); // Clears the content of the original nodes_tags vector
        line_tags.clear(); // Clears the content of the original nodes_tags vector

    }
    gmsh::model::geo::synchronize();



    // We finally generate and save the mesh:
    gmsh::model::mesh::generate(2);
    gmsh::write("C:\\work_dev\\build\\mesh\\Release\\contact_mesh.msh");

    // Stop measuring time
    auto end = std::chrono::high_resolution_clock::now();
    std::vector<double> Points_coords;
    std::vector<double> pcoords;
    gmsh::model::mesh::getNodes(nodes_tags, Points_coords, pcoords);  // Coordinates of the mesh nodes



    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;

    std::vector<double> Contacts_Points_Coords;
    std::vector<double> NContacts_Points_Coords;

    std::vector<int> ContactsTriNodes;
    std::vector<int> ContactsTriNum;

    int ContactTriNum = 0;
    ContactsTriNum.push_back(0);

    std::vector<int> NContactsTriNodes;
    std::vector<int> NContactsTriNum;
    int Ncontact_pt_ind = 1;
    int contact_pt_ind = 1;
    
    int NContactTriNum = 0;
    NContactsTriNum.push_back(0);
    int index = 100;
    for (int i = 0; i < faces_num; ++i)
    {
        std::cout << faces_num << std::endl;
        face_pts_map.clear();

        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2, i + 1);




        if (active_faces_map.find(i + 1) != active_faces_map.end())
        {
            active_faces_map[i + 1] = ContactsTriNum.size();

            for (const auto& tags : nodeTags) {

                for (const auto& tag : tags)
                {
                    if (face_pts_map.find(tag) == face_pts_map.end())
                    {

                        Contacts_Points_Coords.push_back(Points_coords[3 * (tag - 1)]);
                        Contacts_Points_Coords.push_back(Points_coords[3 * (tag - 1) + 1]);
                        Contacts_Points_Coords.push_back(Points_coords[3 * (tag - 1) + 2]);
                        ContactsTriNodes.push_back(contact_pt_ind);
                        face_pts_map.insert({ tag , contact_pt_ind });

                        contact_pt_ind++;
                    }
                    else
                    {
                        ContactsTriNodes.push_back(face_pts_map[tag]);
                    }

                    ContactTriNum += 1;
                }
            }

            ContactsTriNum.push_back(ContactTriNum / 3);
        }

        else
        {
            printf("entereeeed");
            for (const auto& tags : nodeTags) {

                for (const auto& tag : tags)
                {
                    if (face_pts_map.find(tag) == face_pts_map.end())
                    {
                        NContacts_Points_Coords.push_back(Points_coords[3 * (tag - 1)]);
                        NContacts_Points_Coords.push_back(Points_coords[3 * (tag - 1) + 1]);
                        NContacts_Points_Coords.push_back(Points_coords[3 * (tag - 1) + 2]);

                        NContactsTriNodes.push_back(Ncontact_pt_ind);
                        face_pts_map.insert({ tag , Ncontact_pt_ind });

                        Ncontact_pt_ind++;
                    }
                    else
                    {
                        NContactsTriNodes.push_back(face_pts_map[tag]);
                    }

                    NContactTriNum += 1;
                }
            }

            NContactsTriNum.push_back(NContactTriNum / 3);

        }

    }
        gmsh::finalize();



        // Calculate the duration
        std::chrono::duration<double> duration = end - start;

        // Print the duration in seconds
        std::cout << "Operation time: " << duration.count() << " seconds" << std::endl;



        std::vector<double> rhos = { 1, 1, 1 }; // Example rhos vector, change as needed




        Eigen::Matrix3d Tri_nodes_coor;
        Eigen::Vector3d normal; // Example normal vector
        Eigen::Vector3d x_local; // Example x_local vector
        Eigen::Vector3d y_local; // Example y_local vector
        Eigen::Vector3d block_centroid; // Example block centroid

        Eigen::Matrix3d local_ref; // Example block centroid

        // Call the integral_gauss function
        std::vector<double> loaded_terms_ref;
        std::vector<double> moment_terms_ref;


        std::tuple<std::vector<double>, std::vector<double>> terms_coefs;

        std::unordered_map<int, double> fxeq_coefs; // Key-value pairs
        std::unordered_map<int, double> fyeq_coefs; // Key-value pairs
        std::unordered_map<int, double> fzeq_coefs; // Key-value pairs
        std::unordered_map<int, double> mxeq_coefs; // Key-value pairs
        std::unordered_map<int, double> myeq_coefs; // Key-value pairs
        std::unordered_map<int, double> mzeq_coefs; // Key-value pairs

        std::vector<double> eq_coefs;
        std::vector<int> eq_cols;
        std::vector<int> eq_rows;

        for (int i = 0; i < brows; ++i)
        {
            int face_ind = blocks[i * bcols + 1];
            int block_ind = blocks[i * bcols];
            int face_rk = blocks[i * bcols + 6];

            int sign = 1;

            if (face_rk == 22)
            {
                sign = -1;
            }

            if (active_faces_map.find(face_ind) != active_faces_map.end())

            {


                local_ref << sign * c_local_ref[10 * (face_ind - 1) + 4], sign* c_local_ref[10 * (face_ind - 1) + 5], sign* c_local_ref[10 * (face_ind - 1) + 6],
                    sign* c_local_ref[10 * (face_ind - 1) + 7], sign* c_local_ref[10 * (face_ind - 1) + 8], sign* c_local_ref[10 * (face_ind - 1) + 9],
                    sign* c_local_ref[10 * (face_ind - 1) + 1], sign* c_local_ref[10 * (face_ind - 1) + 2], sign* c_local_ref[10 * (face_ind - 1) + 3];



                block_centroid << blocks_centroid[4 * block_ind - 3], blocks_centroid[4 * block_ind - 2], blocks_centroid[4 * block_ind - 1];

                int face_pos = active_faces_map[face_ind];


                for (int Tri_ind = ContactsTriNum[face_pos - 1]; Tri_ind < ContactsTriNum[face_pos]; ++Tri_ind)
                {
                    int Triangle_nodes[3] = { ContactsTriNodes[3 * Tri_ind], ContactsTriNodes[3 * Tri_ind + 1], ContactsTriNodes[3 * Tri_ind + 2] };



                    Tri_nodes_coor << Contacts_Points_Coords[3 * Triangle_nodes[0] - 3], Contacts_Points_Coords[3 * Triangle_nodes[0] - 2], Contacts_Points_Coords[3 * Triangle_nodes[0] - 1],
                        Contacts_Points_Coords[3 * Triangle_nodes[1] - 3], Contacts_Points_Coords[3 * Triangle_nodes[1] - 2], Contacts_Points_Coords[3 * Triangle_nodes[1] - 1],
                        Contacts_Points_Coords[3 * Triangle_nodes[2] - 3], Contacts_Points_Coords[3 * Triangle_nodes[2] - 2], Contacts_Points_Coords[3 * Triangle_nodes[2] - 1];




                    terms_coefs = integral_gauss(Tri_nodes_coor, local_ref, block_centroid, rhos);

                    loaded_terms_ref = std::get<0>(terms_coefs);
                    moment_terms_ref = std::get<1>(terms_coefs);
                    // fx calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 0; col < 3; ++col)
                        {
                            if (abs(loaded_terms_ref[node] * local_ref(col, 0)) > 1e-10)
                            {
                                //std::cout <<"locallll" << local_ref(col, 0) << std::endl;
                                auto it = fxeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col);
                                if (it == fxeq_coefs.end())

                                    fxeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col ,loaded_terms_ref[node] * local_ref(col, 0) });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    fxeq_coefs[(Triangle_nodes[node] - 1) * 3 + col] += loaded_terms_ref[node] * local_ref(col, 0);
                            }
                        }
                    }

                    // fy calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 0; col < 3; ++col)
                        {
                            if (abs(loaded_terms_ref[node] * local_ref(col, 1)) > 1e-10)
                            {
                                //std::cout << "locallll" << local_ref(col, 1) << std::endl;
                                auto it = fyeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col);
                                if (it == fyeq_coefs.end())

                                    fyeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col ,loaded_terms_ref[node] * local_ref(col, 1) });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    fyeq_coefs[(Triangle_nodes[node] - 1) * 3 + col] += loaded_terms_ref[node] * local_ref(col, 1);
                            }
                        }
                    }

                    // fz calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 0; col < 3; ++col)
                        {
                            if (abs(loaded_terms_ref[node] * local_ref(col, 2)) > 1e-10)
                            {
                                //std::cout << "locallll" << local_ref(col, 2) << std::endl;
                                auto it = fzeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col);
                                if (it == fzeq_coefs.end())

                                    fzeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col ,loaded_terms_ref[node] * local_ref(col, 2) });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    fzeq_coefs[(Triangle_nodes[node] - 1) * 3 + col] += loaded_terms_ref[node] * local_ref(col, 2);
                            }
                        }
                    }

                    // mx calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 0; col < 3; ++col)
                        {
                            if (abs(moment_terms_ref[node * 9 + col]) > 1e-10)
                            {

                                auto it = mxeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col);
                                if (it == mxeq_coefs.end())

                                    mxeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col ,moment_terms_ref[node * 9 + col] });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    mxeq_coefs[(Triangle_nodes[node] - 1) * 3 + col] += moment_terms_ref[node * 9 + col];
                            }
                        }
                    }

                    // my calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 3; col < 6; ++col)
                        {
                            if (abs(moment_terms_ref[node * 9 + col]) > 1e-10)
                            {

                                auto it = myeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col - 3);
                                if (it == myeq_coefs.end())

                                    myeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col - 3 ,moment_terms_ref[node * 9 + col] });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    myeq_coefs[(Triangle_nodes[node] - 1) * 3 + col - 3] += moment_terms_ref[node * 9 + col];
                            }
                        }
                    }

                    // mz calc
                    for (int node = 0; node < 3; ++node)
                    {
                        for (int col = 6; col < 9; ++col)
                        {
                            if (abs(moment_terms_ref[node * 9 + col]) > 1e-10)
                            {

                                auto it = mzeq_coefs.find((Triangle_nodes[node] - 1) * 3 + col - 6);
                                if (it == mzeq_coefs.end())

                                    mzeq_coefs.insert({ (Triangle_nodes[node] - 1) * 3 + col - 6 ,moment_terms_ref[node * 9 + col] });
                                //fxeq_rows.push_back((block_ind - 1) * 6);

                                else

                                    mzeq_coefs[(Triangle_nodes[node] - 1) * 3 + col - 6] += moment_terms_ref[node * 9 + col];
                            }
                        }
                    }

                }



                for (const auto& pair : fxeq_coefs)
                {


                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6);
                }
                for (const auto& pair : fyeq_coefs)
                {
                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6 + 1);
                }
                for (const auto& pair : fzeq_coefs)
                {
                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6 + 2);
                }
                for (const auto& pair : mxeq_coefs)
                {
                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6 + 3);
                }
                for (const auto& pair : myeq_coefs)
                {
                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6 + 4);
                }

                for (const auto& pair : mzeq_coefs)
                {
                    eq_coefs.push_back(pair.second);
                    eq_cols.push_back(pair.first);
                    eq_rows.push_back((block_ind - 1) * 6 + 5);
                }

                fxeq_coefs.clear();
                fyeq_coefs.clear();


                fzeq_coefs.clear();

                mxeq_coefs.clear();
                myeq_coefs.clear();
                mzeq_coefs.clear();

            }

        }

        
        // Allocate memory for the array
        double* eq_coefs_pt = new double[eq_coefs.size() + 1];
        int* eq_cols_pt = new int[eq_cols.size() + 1];
        int* eq_rows_pt = new int[eq_rows.size() + 1];


        // Copy elements from vector to array
        for (size_t i = 0; i < eq_coefs.size(); ++i) {
            eq_coefs_pt[i] = eq_coefs[i];
        }

        for (size_t i = 0; i < eq_cols.size(); ++i) {
            eq_cols_pt[i] = eq_cols[i];
        }
        for (size_t i = 0; i < eq_rows.size(); ++i) {
            eq_rows_pt[i] = eq_rows[i];
        }


        
        // Assign addresses of values to ptr1 and ptr2
        Output A;
        A.sparse_dim = eq_coefs.size();
        A.eq_coefs = eq_coefs_pt;
        A.eq_cols = eq_cols_pt;
        A.eq_rows = eq_rows_pt;



        std::cout << Contacts_Points_Coords.size() << std::endl;
        double* Points_coords_pt = new double[Contacts_Points_Coords.size() + 1];
        int* TriNodes_pt = new int[ContactsTriNodes.size() + 1];
        int* FacesTriNum_pt = new int[ContactsTriNum.size() + 1];

        for (size_t i = 0; i < Contacts_Points_Coords.size(); ++i) {
            Points_coords_pt[i] = Contacts_Points_Coords[i];
        }

        for (size_t i = 0; i < ContactsTriNodes.size(); ++i) {
            TriNodes_pt[i] = ContactsTriNodes[i];
        }

        for (size_t i = 0; i < ContactsTriNum.size(); ++i) {
            FacesTriNum_pt[i] = ContactsTriNum[i];
        }

        A.ContactsPointsCoords = Points_coords_pt;
        A.ContactsTriNodes = TriNodes_pt;
        A.ContactsTriNum = FacesTriNum_pt;
        A.Contacts_pts_num = int((Contacts_Points_Coords.size()) / 3);

        
        // non contact elements

        std::cout << NContacts_Points_Coords.size() << std::endl;

        double* NPoints_coords_pt = new double[NContacts_Points_Coords.size() + 1];
        int* NTriNodes_pt = new int[NContactsTriNodes.size() + 1];
        int* NFacesTriNum_pt = new int[NContactsTriNum.size() + 1];

        for (size_t i = 0; i < NContacts_Points_Coords.size(); ++i) {
            NPoints_coords_pt[i] = NContacts_Points_Coords[i];
        }

        for (size_t i = 0; i < NContactsTriNodes.size(); ++i) {
            NTriNodes_pt[i] = NContactsTriNodes[i];
        }

        for (size_t i = 0; i < NContactsTriNum.size(); ++i) {
            NFacesTriNum_pt[i] = NContactsTriNum[i];
        }

        A.NContactsPointsCoords = NPoints_coords_pt;
        A.NContactsTriNodes = NTriNodes_pt;
        A.NContactsTriNum = NFacesTriNum_pt;
        A.NContacts_pts_num = int((NContacts_Points_Coords.size()) / 3);

        

        return A;

    
}
