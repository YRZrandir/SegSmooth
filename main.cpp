#include <filesystem>
#include <iostream>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/version.h>
#include "SegSmooth.h"

int main(int argc, char* argv[])
{
#ifdef _DEBUG
    std::cout << "Mode: Debug" << std::endl;
#endif
    std::cout << "CGAL Version: " << CGAL_STR(CGAL_VERSION) << std::endl;

    std::cout << "working dir=" << std::filesystem::current_path() << std::endl;

    Polyhedron scanmesh;
    CGAL::IO::read_OBJ("../../test/mesh1.obj", scanmesh);
    LoadLabels(scanmesh, "../../test/mesh1.json");
    SmoothSegmentation(scanmesh);
    scanmesh.WriteOBJ("../../test/out.obj");

    return 0;
}