#include "SegSmooth.h"
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <nlohmann/json.hpp>
#include <omp.h>

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
        //hf->_label = std::max(l0, std::max(l1, l2));
        if(l0 == l1 && l1 == l2)
            hf->_label = l0;
        else if(l0 == l1 && l1 != l2)
            hf->_label = l0;
        else if(l0 == l2 && l1 != l2)
            hf->_label = l0;
        else if(l1 == l2 && l0 != l2)
            hf->_label = l1;
        else
            hf->_label = l0;
    }
}

void SmoothSegmentation(Polyhedron& mesh)
{
    using AABBPrimitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using AABBTraits = CGAL::AABB_traits<KernelEpick, AABBPrimitive>;
    using AABBTree = CGAL::AABB_tree<AABBTraits>;
    Polyhedron mesh_ori = mesh;
    AABBTree tree(CGAL::faces(mesh_ori).first, CGAL::faces(mesh_ori).second, mesh_ori);

    std::unordered_set<hFacet> faces_to_process;
    std::unordered_set<hVertex> control_vertices_set;
    for(auto hh : CGAL::halfedges(mesh))
    {
        if(hh->is_border_edge())
            continue;
        if(hh->facet()->_label != hh->opposite()->facet()->_label)
        {
            faces_to_process.insert(hh->facet());
            faces_to_process.insert(hh->opposite()->facet());
            control_vertices_set.insert(hh->vertex());
        }
    }

    std::vector<hVertex> ctrl_vertices(control_vertices_set.begin(), control_vertices_set.end());
    for(auto hv : ctrl_vertices)
    {
        hv->_is_const = true;
    }

    for(auto hf : CGAL::faces(mesh))
    {
        int l0 = hf->halfedge()->vertex()->_label;
        int l1 = hf->halfedge()->next()->vertex()->_label;
        int l2 = hf->halfedge()->prev()->vertex()->_label;
        if(l0 == l1 && l1 == l2)
            continue;
        faces_to_process.insert(hf);
    }

    for(int i = 0; i < 5; i++)
    {
        std::unordered_set<hFacet> temp = faces_to_process;
        for(auto hf : faces_to_process)
        {
            for(auto nei : CGAL::faces_around_face(hf->halfedge(), mesh))
                temp.insert(nei);
        }
        faces_to_process = temp;
    }

    std::vector<hVertex> roi_vertices;
    std::vector<hFacet> roi_faces;
    for(auto hf : faces_to_process)
    {
        roi_faces.push_back(hf);
        for(auto hv : CGAL::vertices_around_face(hf->halfedge(), mesh))
        {
            if(control_vertices_set.count(hv) == 0)
            {
                roi_vertices.push_back(hv);
            }
        }
    }
    roi_vertices.erase(std::unique(roi_vertices.begin(), roi_vertices.end()), roi_vertices.end());

    // for(auto hv : roi_vertices)
    // {
    //     hv->_label = 18;
    // }
    // for(auto hv : ctrl_vertices)
    // {
    //     hv->_label = 19;
    // }
    std::vector<hFacet> faces_to_smooth;
    for(auto hf : faces_to_process)
    {
        auto hv0 = hf->halfedge()->vertex();
        auto hv1 = hf->halfedge()->next()->vertex();
        auto hv2 = hf->halfedge()->prev()->vertex();
        if(control_vertices_set.count(hv0) == 0 && control_vertices_set.count(hv1) == 0 && control_vertices_set.count(hv2) == 0)
        {
            faces_to_smooth.push_back(hf);
        }
    }

    for(int iteration = 0; iteration < 10; iteration++)
    {
        std::vector<Point_3> newpositions(ctrl_vertices.size());
#pragma omp paralle for
        for(int i = 0; i < ctrl_vertices.size(); i++)
        {
            auto hv = ctrl_vertices[i];
            Point_3 newpos = hv->point();
            Vector_3 diffsum(0, 0, 0);
            int count = 0;
            for(auto nei : CGAL::vertices_around_target(hv, mesh))
            {
                if(control_vertices_set.count(nei) != 0)
                {
                    Vector_3 diff = nei->point() - hv->point();
                    diffsum += diff;
                    count++;
                }
            }
            newpos += diffsum / (float)count;
            newpos = tree.closest_point(newpos);
            newpositions[i] = newpos;
        }

        for(int i = 0; i < ctrl_vertices.size(); i++)
        {
            ctrl_vertices[i]->point() = newpositions[i];
        }
        // CGAL::Polygon_mesh_processing::angle_and_area_smoothing(faces_to_smooth, mesh,
        // CGAL::parameters::number_of_iterations(1).use_safety_constraints(true).use_area_smoothing(false));
        std::vector<Point_3> newpositions2(roi_vertices.size());
#pragma omp paralle for
        for(int i = 0; i < roi_vertices.size(); i++)
        {
            auto hv = roi_vertices[i];
            Point_3 newpos = hv->point();
            Vector_3 diffsum(0, 0, 0);
            int count = 0;
            for(auto nei : CGAL::vertices_around_target(hv, mesh))
            {
                diffsum += nei->point() - hv->point();
                count++;
            }
            newpos += diffsum / (float)count;
            newpos = tree.closest_point(newpos);
            newpositions2[i] = newpos;
        }
        
        for(int i = 0; i < roi_vertices.size(); i++)
        {
            roi_vertices[i]->point() = newpositions2[i];
        }
    }
}
