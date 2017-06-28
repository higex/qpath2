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
#include <CGAL/Polygon_2_algorithms.h>
#include <list>


namespace py = boost::python;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel>                       Polygon_set_2;


// POINT_WRT_POLYGON
// Check the position of a (set of) point(s) with respect to a polygon.
//
// To simplify argument decoding, the inputs are given as separate lists of
// x- and y- coordinates. The result is a list of codes (one per tested
// point): 1 for a point inside the polygon, 0 for a point on the boundary
// and -1 for a point outside the polygon, respectively.
//
int point_wrt_polygon(py::list p_x, py::list p_y,
    py::list Q_x, py::list Q_y, py::list r)
{
    Polygon_2 Q;

    // Prepare data
    if (len(p_x) != len(p_y)) return -1;
    if (len(Q_x) != len(Q_y)) return -2;

    for (int k=0; k < len(Q_x); ++k) {
        double x = py::extract<double>(Q_x[k]);
        double y = py::extract<double>(Q_y[k]);
        Q.push_back(Point_2(x, y));
    }

    for (int k=0; k < len(p_x); ++k) {
        int rel_pos = -2;
        double x = py::extract<double>(p_x[k]);
        double y = py::extract<double>(p_y[k]);

        switch(CGAL::bounded_side_2(Q.vertices_begin(), Q.vertices_end(),
            Point_2(x, y), Kernel())) {
                case CGAL::ON_BOUNDED_SIDE:
                    rel_pos = 1;
                    break;
                case CGAL::ON_BOUNDARY:
                    rel_pos = 0;
                    break;
                case CGAL::ON_UNBOUNDED_SIDE:
                    rel_pos = -1;
                    break;
                default:
                    return -3;
            }
        r.append(rel_pos);
    }
    return 0;
}


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
    
    CGAL::intersection(P, Q, std::back_inserter(R));
    
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


// POLYGON_EQUALITY
//
// Test whether two polygons are equal (i.e. there is a permutation of the vertices of
// one polygon such that they equal the vertices of the second polygon).
//
// To simplify argument decoding, the input polygons are given as a pair
// of lists: a vertex i in P/Q has the coordinates (P/Q_x[i], P/Q_y[i]).
//
// Returns 1 for equal polygons, 0 for non-equal and negative codes for errors.
int polygon_equality(py::list P_x, py::list P_y, py::list Q_x, py::list Q_y)
{
    Polygon_2 P, Q;

    // Prepare data
    if (len(P_x) != len(P_y)) return -1;
    if (len(Q_x) != len(Q_y)) return -2;

    if (len(P_x) != len(Q_x)) return 0; // different number of vertices, polygons non-equal

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

    return static_cast<int>(P == Q);
}


BOOST_PYTHON_MODULE(compgeom_){
    def("simple_polygon_intersection_",
        simple_polygon_intersection);
    def("point_wrt_polygon_",
        point_wrt_polygon);
    def("polygon_equality_",
        polygon_equality);
}
