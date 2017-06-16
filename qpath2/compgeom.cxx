//----------------------------------------------------------------------
// compgeom.cxx : Computational geometry module for Python.
//
//                This module simply wraps a number of functions from 
//                CGAL and exposes them as a Python module. For a more
//                complete CGAL interface, check the CGAL-swig project
//                on GitHub.
// Author: Vlad Popovici
//----------------------------------------------------------------------

#include <boost/python.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <list>


namespace py = boost::python;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel>                       Polygon_set_2;


// SIMPLE_POLYGON_INTERSECTION
// Intersection of two simple polygons which may result in a list of polygons 
// without holes. This functionality is still not exposed by CGAL-swig.
//
// To simplify argument decoding, the input polygons are given as a pair
// of lists: a vertex i in P/Q has the coordinates (P/Q_x[i], P/Q_y[i]).
// The intersection may result in several polygons stored by concatenating 
// their vertices in R_x, R_y. The number of vertices in each resulting 
// component is stored in the lists R_n (whose length gives the number of
// components in the intersection).
//
int simple_polygon_intersection(py::list P_x, py::list P_y,
    py::list Q_x, py::list Q_y,
    py::list R_x, py::list R_y, py::list R_n)
{
    Polygon_2 P, Q;
    std::list<Polygon_with_holes_2> R;
    std::list<Polygon_with_holes_2>::const_iterator it;
    
    // Prepare data
    if (len(P_x) != len(P_y)) return -1;
    if (len(Q_x) != len(Q_y)) return -2;
    
    for (int k=0; k < len(P_x); ++k) {
        double x = py::extract<double>(P_x[k]);
        double y = py::extract<double>(P_y[k]);
        P.push_back(Point_2(x, y));
    }
    for (int k=0; k < len(Q_x); ++k) {
        double x = py::extract<double>(Q_x[k]);
        double y = py::extract<double>(Q_y[k]);
        Q.push_back(Point_2(x, y));
    }
    
    if (!P.is_simple() || !Q.is_simple()) return -3;

    if (!P.is_counterclockwise_oriented()) {
        P.reverse_orientation();
    }

    if (!Q.is_counterclockwise_oriented()) {
        Q.reverse_orientation();
    }

    // Intersect polygons
    if (! CGAL::do_intersect(P, Q)) {
        return 0;
    }
    
    CGAL::intersection (P, Q, std::back_inserter(R));
    
    // Prepare result
    for (it = R.begin(); it != R.end(); ++it) { // for each component...
        if ( (*it).is_unbounded() ) return -4;
        // store the outer boundary:
        const Polygon_2 b = (*it).outer_boundary();
        for (Polygon_2::Vertex_const_iterator  vit=b.vertices_begin();
            vit != b.vertices_end(); ++vit) {
                R_x.append(CGAL::to_double((*vit).x()));
                R_y.append(CGAL::to_double((*vit).y()));
        }
        R_n.append(b.size());
    }
    
    return R.size();
}


BOOST_PYTHON_MODULE(compgeom_){
    def("simple_polygon_intersection_",
        simple_polygon_intersection);
}
