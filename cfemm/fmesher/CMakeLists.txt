# allow overriding the default behaviour for finding a suitable triangle:
set(TRIANGLE_PROVIDER auto CACHE STRING "Control which sources of the triangle library are considered.")
set_property(CACHE TRIANGLE_PROVIDER PROPERTY STRINGS auto force-builtin force-external)

if (NOT "${TRIANGLE_PROVIDER}" STREQUAL "force-builtin")
    find_package(Triangle)
endif()

# Tangle is the production mesher for new instancing work.  Upstream currently
# ships as two C++ sources rather than as a CMake project, so define a small,
# pinned library target here while preserving a local-source mode for offline
# development and release builds.
set(TANGLE_PROVIDER auto CACHE STRING
    "Select the Tangle source provider (auto, local, or fetch).")
set_property(CACHE TANGLE_PROVIDER PROPERTY STRINGS auto local fetch)
set(TANGLE_SOURCE_DIR "" CACHE PATH
    "Path to a dcm3c/tangle source checkout (required by TANGLE_PROVIDER=local).")
set(XFEMM_TANGLE_REVISION a808c624ec0584569e43593662f54890b602c6af CACHE STRING
    "Pinned dcm3c/tangle revision used by xfemm.")
set(XFEMM_TEST_MESHER_BACKEND Triangle CACHE STRING
    "Default mesher backend used by the standard test build (Triangle or Tangle).")
set_property(CACHE XFEMM_TEST_MESHER_BACKEND PROPERTY STRINGS Triangle Tangle)
string(TOLOWER "${XFEMM_TEST_MESHER_BACKEND}" _xfemm_test_mesher_backend)
if(NOT _xfemm_test_mesher_backend MATCHES "^(triangle|tangle)$")
    message(FATAL_ERROR
        "Invalid XFEMM_TEST_MESHER_BACKEND='${XFEMM_TEST_MESHER_BACKEND}'; expected Triangle or Tangle")
endif()

if(NOT TANGLE_PROVIDER MATCHES "^(auto|local|fetch)$")
    message(FATAL_ERROR
        "Invalid TANGLE_PROVIDER='${TANGLE_PROVIDER}'; expected auto, local, or fetch")
endif()

if(TANGLE_PROVIDER STREQUAL "local" AND NOT TANGLE_SOURCE_DIR)
    message(FATAL_ERROR
        "TANGLE_PROVIDER=local requires TANGLE_SOURCE_DIR to name a Tangle checkout")
endif()

if(TANGLE_PROVIDER STREQUAL "local")
    set(tangle_SOURCE_DIR "${TANGLE_SOURCE_DIR}")
elseif(TANGLE_PROVIDER STREQUAL "auto" AND TANGLE_SOURCE_DIR)
    set(tangle_SOURCE_DIR "${TANGLE_SOURCE_DIR}")
else()
    include(FetchContent)
    FetchContent_Declare(tangle_upstream
        GIT_REPOSITORY https://github.com/dcm3c/tangle.git
        GIT_TAG ${XFEMM_TANGLE_REVISION}
        GIT_SHALLOW FALSE
        GIT_PROGRESS TRUE
    )
    FetchContent_GetProperties(tangle_upstream)
    if(NOT tangle_upstream_POPULATED)
        FetchContent_Populate(tangle_upstream)
    endif()
    set(tangle_SOURCE_DIR "${tangle_upstream_SOURCE_DIR}")
endif()

foreach(_tangle_file tangle.cpp float256.cpp tangle_mesh.h LICENSE)
    if(NOT EXISTS "${tangle_SOURCE_DIR}/${_tangle_file}")
        message(FATAL_ERROR
            "Tangle source '${tangle_SOURCE_DIR}' is missing ${_tangle_file}")
    endif()
endforeach()

add_library(tangle STATIC
    "${tangle_SOURCE_DIR}/tangle.cpp"
    "${tangle_SOURCE_DIR}/float256.cpp"
)
add_library(Tangle::tangle ALIAS tangle)
target_compile_definitions(tangle PRIVATE TANGLE_AS_LIBRARY)
target_include_directories(tangle SYSTEM PUBLIC
    $<BUILD_INTERFACE:${tangle_SOURCE_DIR}>
)
message(STATUS "Using Tangle sources from ${tangle_SOURCE_DIR}")

if (Triangle_FOUND)
    message("Found triangle version ${Triangle_VERSION}...")
else()
    if (NOT "${TRIANGLE_PROVIDER}" STREQUAL "force-external")
        message("Using builtin version of triangle...")
        add_library(triangle STATIC
            triangle/triangle.c
            )
        # a bit of a dirty hack to make the link line nicer:
        add_library(Triangle::triangle-api ALIAS triangle)
        target_include_directories(triangle PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/triangle> $<INSTALL_INTERFACE:include>)
    else()
        message(FATAL_ERROR "Cache variable TRIANGLE_PROVIDER set to 'force-external', but Triangle was not found!")
    endif()
endif()

add_library(fmesher STATIC
    fmesher.cbp
    fmesher.cpp
    nosebl.cpp
    TangleMeshConverter.cpp
    TangleMesherBackend.cpp
    TriangleMesher.cpp
    TriangleMesherBackend.cpp
    writepoly.cpp
    )
target_link_libraries(fmesher PUBLIC femm PRIVATE Triangle::triangle-api Tangle::tangle)
target_include_directories(fmesher PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
if(_xfemm_test_mesher_backend STREQUAL "tangle")
    target_compile_definitions(fmesher PUBLIC XFEMM_DEFAULT_MESHER_BACKEND_TANGLE)
endif()
message(STATUS "Standard tests use the ${XFEMM_TEST_MESHER_BACKEND} mesher backend")

add_executable(fmesher-bin
    main.cpp
    )

set_target_properties(fmesher-bin
    PROPERTIES OUTPUT_NAME fmesher)

target_link_libraries(fmesher-bin fmesher Triangle::triangle-api)

add_subdirectory(test)

install(
    TARGETS fmesher-bin
    RUNTIME DESTINATION bin
    COMPONENT "cli")
install(
    FILES "${tangle_SOURCE_DIR}/LICENSE"
    DESTINATION doc/xfemm
    RENAME LICENSE-tangle.txt
    COMPONENT "cli")
# vi:expandtab:tabstop=4 shiftwidth=4:
