# Taller 3 - Reducción de Puntos en Relieves Poliedrales

## 🚀 Características

- **Algoritmo de vecindad**: Búsqueda en anchura (BFS) sobre la triangulación
- **Umbral adaptativo**: Se ajusta localmente según la variabilidad del terreno
- **Preservación topológica**: Reconstruye la triangulación de Delaunay tras eliminaciones
- **Eficiencia**: Complejidad temporal $O(n \cdot d^\delta)$ donde $d$ es el grado promedio

## 📦 Dependencias

- **CGAL** (Computational Geometry Algorithms Library) ≥ 5.0
- **C++17** o superior
- CMake ≥ 3.10 (para compilación)

## 🔧 Instalación

### Usando CMake

```bash
git clone https://github.com/usuario/CullPointsHelper.git
cd CullPointsHelper
mkdir build && cd build
cmake ..
make