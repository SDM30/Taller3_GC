#ifndef CULL_POINTS_HELPER_HXX
#define CULL_POINTS_HELPER_HXX

#include "CullPointsHelper.h"
#include <iostream>
#include <vector>
#include <queue>

/**
 * @brief Constructor: Copia los puntos de la triangulación original a la reducida
 *        y crea el mapeo entre puntos e índices para referencia futura
 */
CullPointsHelper::CullPointsHelper(double gamma, uint delta, const TDelaunay &original_tri)
    : gamma(gamma)
    , delta(delta)
    , original_tri(original_tri)
    , reduced_tri()
{
    size_t idx = 0;
    // Recorrer todos los vértices finitos de la triangulación original
    for (auto vit = original_tri.finite_vertices_begin();
         vit != original_tri.finite_vertices_end();
         ++vit)
    {
        TPoint point = vit->point();
        
        // Almacenar la relación punto -> índice original
        orgin_points.emplace_back(point, idx);
        
        // Insertar el punto en la triangulación reducida
        // y preservar la información de altura en info()
        auto vertex = reduced_tri.insert(point);
        vertex->info() = vit->info();

        idx++;
    }
    
    std::cout << "CullPointsHelper inicializado con " << orgin_points.size() << " puntos" << std::endl;
}

CullPointsHelper::~CullPointsHelper() = default;

/**
 * @brief Busca el índice original de un punto en el vector orgin_points
 * @return Índice si se encuentra, std::nullopt si no existe
 */
std::optional<size_t> CullPointsHelper::getIndex(const TPoint &point) const
{
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

/**
 * @brief Implementa el predicado neigh(v, delta) del algoritmo
 * 
 * Realiza una búsqueda en anchura (BFS) sobre la triangulación reducida
 * para encontrar todos los vecinos hasta profundidad delta desde vertex
 * 
 * @param vertex Punto de partida para la búsqueda
 * @return Vector con todos los puntos vecinos encontrados
 */
std::vector<CullPointsHelper::TPoint> CullPointsHelper::neigh(TPoint vertex)
{
    std::vector<TPoint> neighbors;                           // Resultado final
    std::vector<bool> m(original_tri.number_of_vertices(), false); // Marcas de visitados
    std::queue<TPoint> q;                                   // Cola para BFS
    auto first = getIndex(vertex);
    if (first == std::nullopt) return neighbors;            // Vértice no encontrado
    
    // Inicializar BFS desde el vértice dado
    m[first.value()] = true;
    q.push(vertex);
    uint depth = 0;
    
    // BFS: explorar hasta agotar la cola
    while (!q.empty())
    {
        TPoint vq = q.front();
        q.pop();
        neighbors.push_back(vq);                            // Agregar al resultado
        
        std::vector<TPoint> local_neighbors;
        auto vertex_handle = reduced_tri.nearest_vertex(vq);
        
        // Validar que el vértice existe y no es infinito
        if (vertex_handle == Vertex_handle() 
            || reduced_tri.is_infinite(vertex_handle)) return neighbors;
        
        // Obtener vecinos incidentes usando circulador de CGAL
        Vertex_circulator circulator = reduced_tri.incident_vertices(vertex_handle);
        Vertex_circulator done = circulator;
        
        do {
            Vertex_handle neighbor = circulator;
            // Solo considerar vecinos finitos
            if (!reduced_tri.is_infinite(neighbor)) {
                local_neighbors.push_back(neighbor->point());
            }
        } while(++circulator != done);
        
        depth++;

        // Explorar vecinos locales si no exceden la profundidad delta
        for (auto nIt = local_neighbors.begin(); nIt != local_neighbors.end(); ++nIt) {
            auto index = getIndex(*nIt);
            if (index == std::nullopt) return neighbors;
            if (!m[index.value()] && depth <= this->delta) {
                m[index.value()] = true;                     // Marcar como visitado
                q.push(*nIt);                                // Encolar para explorar
            }
        }
    }

    return neighbors;
}

/**
 * @brief Calcula el promedio de las alturas de los vecinos
 *        La altura se obtiene de info() de cada vértice
 */
CullPointsHelper::TReal 
CullPointsHelper::avgHeight(std::vector<CullPointsHelper::TPoint> neighbors) 
{
    if (neighbors.empty()) return 0.0;
    
    TReal sum = 0.0;
    for (const auto& point : neighbors) 
    {
        auto vertex = original_tri.nearest_vertex(point);
        if (vNotEmptyOrInf(vertex, true))
        {
            sum += vertex->info();  // info() almacena la coordenada z
        }
    }
    
    return sum / neighbors.size();
}

/**
 * @brief Calcula la desviación estándar de alturas usando la fórmula:
 *        sigma = sqrt(E[X²] - E[X]²)
 *        donde X son las alturas de los vecinos
 */
CullPointsHelper::TReal 
CullPointsHelper::sigmaHeight(std::vector<CullPointsHelper::TPoint> neighbors) 
{
    if (neighbors.empty()) return 0.0;
    
    TReal sum = 0.0;
    TReal squareSum = 0.0;
    int count = 0;
    
    for (const auto& point : neighbors) 
    {
        auto vertex = original_tri.nearest_vertex(point);
        if (vNotEmptyOrInf(vertex, true))
        {
            TReal height = vertex->info();
            sum += height;               // Acumular para la media
            squareSum += height * height; // Acumular para E[X²]
            count++;
        }
    }
    
    if (count == 0) return 0.0;
    
    // Varianza = E[X²] - (E[X])²
    TReal avg = sum / count;
    TReal sigma = (squareSum / count) - (avg * avg);
    
    // Corrección numérica para evitar valores negativos muy pequeños
    if (sigma < 0 && sigma > -1e-10) sigma = 0.0;
    
    return std::sqrt(sigma);
}

/**
 * @brief Algoritmo principal de reducción de puntos
 * 
 * Para cada vértice de la triangulación original:
 * 1. Obtiene vecinos hasta profundidad delta
 * 2. Calcula media y desviación estándar de alturas
 * 3. Calcula umbral = gamma * sigma
 * 4. Elimina vértices cuya diferencia de altura supere el umbral
 */
void CullPointsHelper::reduceOverThreshold()
{
    std::vector<TPoint> to_remove;  // Lista de puntos a eliminar

    // Iterar sobre todos los vértices del terreno original
    for (auto vIt = original_tri.finite_vertices_begin();
         vIt != original_tri.finite_vertices_end(); ++vIt)
    {
        TPoint p = vIt->point();
        auto vertex_handle = reduced_tri.nearest_vertex(p);
        if (!vNotEmptyOrInf(vertex_handle, false)) continue;

        // Paso 1: Obtener vecindario
        std::vector<TPoint> v_delta = neigh(p);
        
        // Paso 2: Calcular estadísticas del vecindario
        TReal avg_z    = avgHeight(v_delta);
        TReal sigma_z  = sigmaHeight(v_delta);
        
        // Paso 3: Determinar umbral local
        TReal threshold = gamma * sigma_z;

        // Log de depuración: mostrar información del vértice actual
        std::cout << "v=(" << p.x() << "," << p.y() << ")"
                  << " vecinos=" << v_delta.size()
                  << " avg_z=" << avg_z
                  << " sigma_z=" << sigma_z
                  << " threshold=" << threshold << std::endl;

        // Paso 4: Marcar puntos que exceden el umbral
        for (const auto& v : v_delta) {
            auto vertex = reduced_tri.nearest_vertex(v);
            if (!vNotEmptyOrInf(vertex, false)) continue;
            TReal diff = std::abs(vertex->info() - avg_z);

            std::cout << "  candidato info=" << vertex->info()
                      << " diff=" << diff
                      << (diff >= threshold ? " -> ELIMINAR" : "") << std::endl;

            if (diff >= threshold)
                to_remove.push_back(v);
        }
    }

    std::cout << "Puntos a eliminar: " << to_remove.size() << std::endl;

    // Eliminar duplicados de la lista de puntos a remover
    std::sort(to_remove.begin(), to_remove.end(), [](const TPoint& a, const TPoint& b){
        return std::tie(a.x(), a.y()) < std::tie(b.x(), b.y());
    });
    to_remove.erase(std::unique(to_remove.begin(), to_remove.end()), to_remove.end());

    // Eliminar puntos de la triangulación reducida
    for (const auto& p : to_remove) {
        auto vh = reduced_tri.nearest_vertex(p);
        if (vNotEmptyOrInf(vh, false) && vh->point() == p)
            reduced_tri.remove(vh);
    }

    std::cout << "Puntos restantes: " << reduced_tri.number_of_vertices() << std::endl;
}

/**
 * @brief Verifica que un manejador de vértice sea válido y no represente
 *        el vértice infinito en la triangulación correspondiente
 */
bool CullPointsHelper::vNotEmptyOrInf(Vertex_handle vh, bool original) {
    if (original)
        return vh != Vertex_handle() && !original_tri.is_infinite(vh);
    return vh != Vertex_handle() && !reduced_tri.is_infinite(vh);
}

/**
 * @brief Retorna la triangulación resultante después del proceso de reducción
 */
const CullPointsHelper::TDelaunay& 
CullPointsHelper::getReducedTriangulation() const { return reduced_tri; }

#endif