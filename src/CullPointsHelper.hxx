#ifndef CULL_POINTS_HELPER_HXX
#define CULL_POINTS_HELPER_HXX

#include "CullPointsHelper.h"
#include <iostream>
#include <vector>
#include <queue>


// Constructor: Copiar puntos y relacionarlos con índices
CullPointsHelper::CullPointsHelper(double gamma, uint delta, const TDelaunay &original_tri)
    : gamma(gamma)
    , delta(delta)
    , original_tri(original_tri)
    , reduced_tri()
{
    // Iterar sobre los vértices finitos de la triangulación original
    size_t idx = 0;
    for (auto vit = original_tri.finite_vertices_begin();
         vit != original_tri.finite_vertices_end();
         ++vit)
    {
        TPoint point = vit->point();
        
        // Guardar en orgin_points la relación punto -> índice
        orgin_points.emplace_back(point, idx);
        
        // Insertar en reduced_tri y guardar el índice en info()
        auto vertex = reduced_tri.insert(point);
        vertex->info() = vit->info();

        idx++;
    }
    
    std::cout << "CullPointsHelper inicializado con " << orgin_points.size() << " puntos" << std::endl;
}

CullPointsHelper::~CullPointsHelper() = default;

// Método para obtener el índice de un punto
std::optional<size_t> CullPointsHelper::getIndex(const TPoint &point) const
{
    // Buscar en orgin_points usando std::find_if
    auto it = std::find_if(orgin_points.begin(), orgin_points.end(),
        [&point](const VertIndx& vi) {
            return vi.point == point;
        });
    
    if (it != orgin_points.end())
    {
        return it->index;
    }
    
    return std::nullopt;
}

std::vector<CullPointsHelper::TPoint> CullPointsHelper::neigh(TPoint vertex)
{
    std::vector<TPoint> neighbors;
    std::vector<bool> m(original_tri.number_of_vertices(), false);
    std::queue<TPoint> q;
    auto first = getIndex(vertex);
    if (first == std::nullopt) return neighbors;
    m[first.value()] = true;
    q.push(vertex);
    uint depth = 0;
    while (!q.empty())
    {
        TPoint vq = q.front();
        q.pop();
        neighbors.push_back(vq);
        
        std::vector<TPoint> local_neighbors;
        auto vertex_handle =reduced_tri.nearest_vertex(vq);
        // No se encuentra el vértice o es inválido
        if (vertex_handle == Vertex_handle() 
            || reduced_tri.is_infinite(vertex_handle)) return neighbors;
        Vertex_circulator circulator = reduced_tri.incident_vertices(vertex_handle);
        Vertex_circulator done = circulator;
        
        do {
            Vertex_handle neighbor = circulator;

            if (!reduced_tri.is_infinite(neighbor)) {
                local_neighbors.push_back(neighbor->point());
            }
        } while(++circulator != done);
        depth ++;

        for (auto nIt = local_neighbors.begin(); nIt != local_neighbors.end(); ++nIt) {
            auto index = getIndex(*nIt);
            if (index == std::nullopt) return neighbors;
            if (!m[index.value()] && depth <= this->delta) {
                m[index.value()] = true;
                q.push(*nIt);
            }
        }
    }

    return neighbors;
}

CullPointsHelper::TReal 
CullPointsHelper::avgHeight(std::vector<CullPointsHelper::TPoint> neighbors) 
{
    if (neighbors.empty()) return 0.0;
    
    TReal sum = 0.0;
    for (const auto& point : neighbors) 
    {
        // Buscar el vértice en la triangulación original
        auto vertex = original_tri.nearest_vertex(point);
        if (vNotEmptyOrInf(vertex, true))
        {
            sum += vertex->info();  // info() contiene la altura
        }
    }
    
    return sum / neighbors.size();
}

CullPointsHelper::TReal 
CullPointsHelper::sigmaHeight(std::vector<CullPointsHelper::TPoint> neighbors) 
{
    if (neighbors.empty()) return 0.0;
    
    TReal sum = 0.0;
    TReal squareSum = 0.0;
    int count = 0;
    
    for (const auto& point : neighbors) 
    {
        // Buscar el vértice en la triangulación original
        auto vertex = original_tri.nearest_vertex(point);
        if (vNotEmptyOrInf(vertex, true))
        {
            TReal height = vertex->info();
            sum += height;
            squareSum += height * height;
            count++;
        }
    }
    
    if (count == 0) return 0.0;
    
    // Fórmula: sigma = sqrt( (sum(x^2) - (sum(x))^2/n) / n )
    TReal avg = sum / count;
    TReal sigma = (squareSum / count) - (avg * avg);
    
    // Asegurar que no haya errores numéricos (varianza negativa muy pequeña)
    if (sigma < 0 && sigma > -1e-10) sigma = 0.0;
    
    return std::sqrt(sigma);
}

void CullPointsHelper::reduceOverThreshold()
{
    std::vector<TPoint> to_remove;

    for (auto vIt = original_tri.finite_vertices_begin();
         vIt != original_tri.finite_vertices_end(); ++vIt)
    {
        TPoint p = vIt->point();
        auto vertex_handle = reduced_tri.nearest_vertex(p);
        if (!vNotEmptyOrInf(vertex_handle, false)) continue;

        std::vector<TPoint> v_delta = neigh(p);
        TReal avg_z   = avgHeight(v_delta);
        TReal sigma_z = sigmaHeight(v_delta);
        TReal threshold = gamma * sigma_z;

        for (const auto& v : v_delta) {
            auto vertex = reduced_tri.nearest_vertex(v);
            if (!vNotEmptyOrInf(vertex, false)) continue;

            if (std::abs(vertex->info() - avg_z) >= threshold) {
                to_remove.push_back(v); 
            }
        }
    }

    // Antes del loop de eliminación
    std::sort(to_remove.begin(), to_remove.end(), [](const TPoint& a, const TPoint& b){
        return std::tie(a.x(), a.y()) < std::tie(b.x(), b.y());
    });
    to_remove.erase(std::unique(to_remove.begin(), to_remove.end()), to_remove.end());
    // Eliminar sin duplicados
    for (const auto& p : to_remove) {
        auto vh = reduced_tri.nearest_vertex(p);
        // Verificar que siga existiendo y sea el punto correcto
        if (vNotEmptyOrInf(vh, false) && vh->point() == p) {
            reduced_tri.remove(vh);
        }
    }
}

bool CullPointsHelper::vNotEmptyOrInf(Vertex_handle vh, bool original) {
    if (original)
        return vh != Vertex_handle() && !original_tri.is_infinite(vh);
    return vh != Vertex_handle() && !reduced_tri.is_infinite(vh);
}

const CullPointsHelper::TDelaunay& 
CullPointsHelper::getReducedTriangulation() const { return reduced_tri; }
#endif
