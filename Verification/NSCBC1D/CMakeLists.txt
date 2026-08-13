# Standalone 1-D verification driver for Source/NSCBC.H.
#
# Deliberately independent of the PeleC application build: it compiles the
# production boundary-condition header against a minimal 1-D solver so that a
# failure can be attributed to the boundary condition rather than to the AMReX
# plumbing around it.  See the header comment in nscbc1d.cpp.
#
#   cmake -S . -B build -DAMReX_DIR=<amrex-install>/lib/cmake/AMReX \
#         -DPELEC_DIR=../.. -DPELE_MECHANISM=air
#   cmake --build build && ./build/nscbc1d sweep
cmake_minimum_required(VERSION 3.20)
project(nscbc1d CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PELEC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "PeleC root")
set(PELE_MECHANISM "air" CACHE STRING "PelePhysics mechanism")
find_package(AMReX REQUIRED)

set(PP "${PELEC_DIR}/Submodules/PelePhysics")
add_executable(nscbc1d
  nscbc1d.cpp
  ${PP}/Mechanisms/${PELE_MECHANISM}/mechanism.cpp)
target_include_directories(nscbc1d PRIVATE
  ${PELEC_DIR}/Source
  ${PP}/Source
  ${PP}/Source/Eos
  ${PP}/Source/Transport
  ${PP}/Source/Utility
  ${PP}/Source/Utility/Utilities
  ${PP}/Source/Utility/Filter
  ${PP}/Source/Utility/BlackBoxFunction
  ${PP}/Source/Utility/PMF
  ${PP}/Source/Utility/TurbInflow
  ${PP}/Source/Utility/Diagnostics
  ${PP}/Mechanisms/${PELE_MECHANISM})
target_compile_definitions(nscbc1d PRIVATE
  USE_FUEGO_EOS USE_CONSTANT_TRANSPORT NUM_ADV=0 NUM_AUX=0)
target_link_libraries(nscbc1d PRIVATE AMReX::amrex_1d)
