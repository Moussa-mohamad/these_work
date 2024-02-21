#include <gmsh.h>
#include <vector>
#include <cmath>

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);
    gmsh::model::add("t2");

    // Define points for each face
    std::vector<std::vector<std::vector<double>>> faces_nodes = {
        {{0,0,0},{1,0,0}, {1,1,0}, {0,1,0}},
        {{0,0,1},{1,0,1}, {1,1,1}, {0,1,1}},
        {{0,0,0},{0,1,0}, {0,1,1}, {0,0,1}},
        {{1,0,0},{1,1,0}, {1,1,1}, {1,0,1}}
    };

    // Initialize gmsh
    gmsh::initialize();

    // Define nodes and lines for each face
    double lc = 0.05;
    for(const auto& face_nodes : faces_nodes) {
        std::vector<std::size_t> node_tags;
        for(const auto& node : face_nodes)
            node_tags.push_back(gmsh::model::geo::addPoint(node[0], node[1], node[2], lc));

        std::vector<std::size_t> lines;
        for(std::size_t i = 0; i < node_tags.size(); ++i)
            lines.push_back(gmsh::model::geo::addLine(node_tags[i], node_tags[(i + 1) % node_tags.size()]));

        std::size_t loop = gmsh::model::geo::addCurveLoop(lines);
        gmsh::model::geo::addPlaneSurface({loop});
    }

    // Synchronize model
    gmsh::model::geo::synchronize();

    // Generate mesh
    gmsh::model::mesh::generate();

    // Write mesh data
    gmsh::write("square_face.msh");

    // Finalize gmsh
    gmsh::finalize();

    return 0;
}
