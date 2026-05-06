# Taller 3 - Reducción de Puntos en Relieves Poliedrales

## Dependencias

- **CGAL** (Computational Geometry Algorithms Library) ≥ 5.0
- **C++17** o superior
- CMake ≥ 3.10 (para compilación)

## Instalación

### Usando CMake

```bash
git clone https://github.com/usuario/CullPointsHelper.git
cd CullPointsHelper
mkdir build && cd build
cmake ..
make
``` 
## Uso
```bash
./heightmap <input.png> <output> <delta> <gamma>
```