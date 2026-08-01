cmake_minimum_required(VERSION 3.15)
project(xfemm_python LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Reuse the native targets instead of maintaining a Python-only solver build.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cfemm
                 ${CMAKE_CURRENT_BINARY_DIR}/cfemm EXCLUDE_FROM_ALL)

pybind11_add_module(_xfemm src/bindings.cpp)
target_link_libraries(_xfemm PRIVATE fsolver fpproc)
target_include_directories(_xfemm PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../cfemm/libfemm
    ${CMAKE_CURRENT_SOURCE_DIR}/../cfemm/fmesher
    ${CMAKE_CURRENT_SOURCE_DIR}/../cfemm/fsolver)
target_include_directories(_xfemm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../cfemm/fpproc)

install(TARGETS _xfemm LIBRARY DESTINATION xfemm RUNTIME DESTINATION xfemm)
