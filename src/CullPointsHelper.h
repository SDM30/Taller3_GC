#ifndef CULL_POINTS_HELPER_H
#define CULL_POINTS_HELPER_H

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_ds_circulators_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <optional>
#include <vector>

/**
 * @brief Clase auxiliar para reducción de puntos en triangulaciones de Delaunay
 *
 * Implementa un algoritmo que elimina vértices basándose en un umbral de
 * altura, calculado a partir de la media y desviación estándar de los vecinos
 * cercanos.
 */
class CullPointsHelper {
public:
  // ===================== DEFINICIONES DE TIPO =====================
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

  /**
   * @brief Estructura para asociar puntos con índices originales
   */
  struct VertIndx {
    TPoint point; // Coordenadas del punto
    size_t index; // Índice original en la triangulación

    VertIndx(const TPoint &p, size_t i) : point(p), index(i) {}
  };

protected:
  // ===================== ATRIBUTOS DEL ALGORITMO =====================
  double gamma;                       // Factor de escala para el umbral
  uint delta;                         // Tamaño del vecindario (saltos)
  const TDelaunay &original_tri;      // Referencia a la triangulación original
  std::vector<VertIndx> orgin_points; // Mapeo punto->índice original
  TDelaunay reduced_tri;              // Triangulación resultante (reducida)

public:
  /**
   * @brief Constructor: inicializa el helper con parámetros y triangulación
   * @param gamma Factor de escala para calcular el umbral
   * @param delta Tamaño del vecindario a considerar
   * @param original_tri Triangulación de Delaunay de entrada
   */
  CullPointsHelper(double gamma, uint delta, const TDelaunay &original_tri);
  virtual ~CullPointsHelper();

  /**
   * @brief Ejecuta el algoritmo principal de reducción de puntos
   */
  void reduceOverThreshold();

  /**
   * @brief Obtiene los vecinos de un vértice dentro de un radio delta
   * @param vertex Punto central para buscar vecinos
   * @return Vector con los puntos vecinos encontrados
   */
  std::vector<TPoint> neigh(TPoint vertex);

  /**
   * @brief Recupera el índice original de un punto
   * @param point Punto a buscar
   * @return Índice si existe, std::nullopt en caso contrario
   */
  std::optional<size_t> getIndex(const TPoint &point) const;

  /**
   * @brief Calcula la altura promedio de un conjunto de puntos
   * @param neighbors Lista de puntos vecinos
   * @return Promedio de alturas
   */
  TReal avgHeight(std::vector<TPoint> neighbors);

  /**
   * @brief Calcula la desviación estándar de alturas
   * @param neighbors Lista de puntos vecinos
   * @return Desviación estándar de las alturas
   */
  TReal sigmaHeight(std::vector<TPoint> neighbors);

  /**
   * @brief Verifica si un manejador de vértice es válido y finito
   * @param vh Manejador del vértice a verificar
   * @param original true para verificar en triangulación original
   * @return true si el vértice es válido y finito
   */
  bool vNotEmptyOrInf(Vertex_handle vh, bool original);

  /**
   * @brief Obtiene la triangulación resultante después de la reducción
   * @return Referencia constante a la triangulación reducida
   */
  const TDelaunay &getReducedTriangulation() const;
};

#endif