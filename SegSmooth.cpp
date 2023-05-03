#include "SegSmooth.h"
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/postprocess.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <nlohmann/json.hpp>

static bool gVerbose = true;

std::pair<std::vector<Point_3>, std::vector<Triangle>> PolyhedronToVF( const Polyhedron& m )
{
    std::vector<Point_3> vertices;

    std::unordered_map<hVertex, unsigned> idmap;
    int count = 0;
    for(auto hv : CGAL::vertices(m))
    {
        idmap[hv] = count;
        vertices.push_back(hv->point());
        count++; 
    }

    std::vector<Triangle> triangles;
    for(auto hf : CGAL::faces(m))
    {
        unsigned i0 = idmap[hf->halfedge()->vertex()];
        unsigned i1 = idmap[hf->halfedge()->next()->vertex()];
        unsigned i2 = idmap[hf->halfedge()->prev()->vertex()];
        triangles.emplace_back(i0, i1, i2);
    }

    return {vertices, triangles};
}

void LoadLabels( Polyhedron& mesh, std::string path )
{
    using namespace nlohmann;
    std::ifstream label_ifs( path );
    json data = json::parse( label_ifs );
    if (data.find( "labels" ) == data.end())
    {
        std::cout << "Invalid Json" << std::endl;
        return;
    }
    std::vector<int> labels = data["labels"].get<std::vector<int>>();
    if(labels.size() != mesh.size_of_vertices())
    {
        std::cout << "number of labels != number of vertices" << std::endl;
        return;
    }
    
    int count = 0;
    for(auto hv = mesh.vertices_begin(); hv != mesh.vertices_end(); hv++)
    {
        hv->_label = labels[count++];
        if(hv->_label == 100)
            hv->_label = 0;
    }

    for(auto hf : CGAL::faces(mesh))
    {
        int l0 = hf->halfedge()->vertex()->_label;
        int l1 = hf->halfedge()->next()->vertex()->_label;
        int l2 = hf->halfedge()->prev()->vertex()->_label;
        hf->_label = std::max(l0, std::max(l1, l2));
    }
}

void SmoothSegmentation(Polyhedron& mesh)
{
    std::unordered_set<hFacet> faces_to_process;
    for(auto hf : CGAL::faces(mesh))
    {
        int l0 = hf->halfedge()->vertex()->_label;
        int l1 = hf->halfedge()->next()->vertex()->_label;
        int l2 = hf->halfedge()->prev()->vertex()->_label;
        if(l0 == l1 && l1 == l2)
            continue;
        faces_to_process.insert(hf);
    }

    CGAL::Polygon_mesh_processing::angle_and_area_smoothing(mesh, CGAL::parameters::number_of_iterations(3));
}
