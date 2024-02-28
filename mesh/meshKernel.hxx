
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

struct Pointers {
    int ptr1;
    double *ptr2;
};



Pointers print_hello_c(double* blocks, int brows, int bcols, int blocks_num, double* nodes, int* active_faces, int active_faces_num, int nrows, int ncols, int* faces_FEpts, int faces_num, double* blocks_centroid, double* c_local_ref);
