## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

# -- Basic configuration
cmake_minimum_required(VERSION 3.10)
project(Taller3)
include("./baseConf.cmake")

include_directories(
  ${CMAKE_BINARY_DIR}/lib
  ${CMAKE_SOURCE_DIR}/lib
  )

# -- Packages
find_package(CGAL REQUIRED COMPONENTS Core)
# Buscar libpng
find_package(PNG REQUIRED)

# -- Executables
add_executable(heightmap src/heightmap.cxx src/CullPointsHelper.cxx)
target_link_libraries(heightmap CGAL::CGAL)
target_link_libraries(heightmap PNG::PNG)

## eof - CMakeLists.txt
