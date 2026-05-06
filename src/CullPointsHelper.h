#ifndef CULL_POINTS_HELPER_H
#define CULL_POINTS_HELPER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <optional>
#include <vector>

class CullPointsHelper
{
public:
    using uint = unsigned int;
    using TKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using TReal = double;
    using TVertices = CGAL::Triangulation_vertex_base_with_info_2<TReal, TKernel>;
    using TTriangulationDS = CGAL::Triangulation_data_structure_2<TVertices>;
    using TDelaunay = CGAL::Delaunay_triangulation_2<TKernel, TTriangulationDS>;
    using TPoint = TKernel::Point_2;
    using Vertex_handle = TDelaunay::Vertex_handle;
    using Face_handle = TDelaunay::Face_handle;
    using Locate_type = TDelaunay::Locate_type;
    using Vertex_circulator = TDelaunay::Vertex_circulator;
    
    struct VertIndx {
        TPoint point;
        size_t index;
        
        VertIndx(const TPoint& p, size_t i) : point(p), index(i) {}
    }; 

protected:
    double gamma;
    uint delta;
    const TDelaunay &original_tri;
    std::vector<VertIndx> orgin_points;  // Almacena puntos con sus índices
    TDelaunay reduced_tri;

public:
    CullPointsHelper(double gamma, uint delta, const TDelaunay &original_tri);
    virtual ~CullPointsHelper();

    void reduceOverThreshold();
    std::vector<TPoint> neigh(TPoint vertex);
    std::optional<size_t> getIndex(const TPoint &point) const;
    TReal avgHeight(std::vector<TPoint> neighbors);
    TReal sigmaHeight(std::vector<TPoint> neighbors);
    bool vNotEmptyOrInf(Vertex_handle vh, bool original);
    const TDelaunay& getReducedTriangulation() const;
};

#include "CullPointsHelper.hxx"

#endif