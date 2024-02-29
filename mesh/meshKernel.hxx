
#ifdef _MSC_VER
#define DllExport __declspec(dllexport)
#else
#define DllExport __attribute__((visibility("default")))
#endif

#include <vector>
#include <Eigen/Dense>
Eigen::Vector3d calculate_moment(const Eigen::Vector3d& point_coords,
    const Eigen::Vector3d& force_coords,
    const Eigen::Vector3d& applied_point_coords);

std::pair<std::vector<int>, std::vector<int>> extractNonZeroWithIndices(const std::vector<int>& source);

Eigen::Vector3d calculate_moment(const Eigen::Vector3d& point_coords,
    const Eigen::Vector3d& force_coords,
    const Eigen::Vector3d& applied_point_coords);

double calculate_triangle_area(const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2,
    const Eigen::Vector3d& v3) ;

void barycentric_coordinates(Eigen::Vector3d pt, Eigen::Matrix3d& triangle_coords,
    double* alpha, double* beta, double* gamma) ;

std::tuple<std::vector<double>, std::vector<double>> integral_gauss(Eigen::Matrix3d& vertices,
    Eigen::Matrix3d& local_ref,
    const Eigen::Vector3d& block_centroid,
    const std::vector<double>& rhos);

// Template function to print vectors of any type
template<typename T>
void printVector(const std::vector<T>& vec);

struct Output {
    int sparse_dim;
    double* eq_coefs;
    int* eq_cols;
    int* eq_rows;
    double* mesh_nodes;
    int* TriNodes;
    int* FacesTriNum;

    //// Destructor to release memory allocated for ptr2
    //~Pointers() {
    //    delete[] ptr2;
    // 
    //}
};



Output print_hello_c(double* blocks, int brows, int bcols, int blocks_num, double* nodes, int* active_faces, int active_faces_num, int nrows, int ncols, int* faces_FEpts, int* Faces_nodes, int faces_num, double* blocks_centroid, double* c_local_ref);
