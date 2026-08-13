# Analysis tooling for the 2-D NSCBC tests (Exec/RegTests/NSCBC-COVO).
#
#   cmake -S . -B build -DAMReX_DIR=<amrex-install>/lib/cmake/AMReX
#   cmake --build build
#   ./build/fielddump <plotfile> pressure out.dat
#   python3 metrics.py circularity <p_amb> <c> out.dat
#
# Needs only a 2-D AMReX build; nothing from PeleC or PelePhysics.
cmake_minimum_required(VERSION 3.20)
project(nscbc2d CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(AMReX REQUIRED)
add_executable(fielddump fielddump.cpp)
target_link_libraries(fielddump PRIVATE AMReX::amrex_2d)
