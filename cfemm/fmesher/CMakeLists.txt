# allow overriding the default behaviour for finding a suitable triangle:
set(TRIANGLE_PROVIDER auto CACHE STRING "Control which sources of the triangle library are considered.")
set_property(CACHE TRIANGLE_PROVIDER PROPERTY STRINGS auto force-builtin force-external)

if (NOT "${TRIANGLE_PROVIDER}" STREQUAL "force-builtin")
    find_package(Triangle)
endif()

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
    TriangleMesher.cpp
    TriangleMesherBackend.cpp
    writepoly.cpp
    )
target_link_libraries(fmesher PUBLIC femm PRIVATE Triangle::triangle-api)
target_include_directories(fmesher PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)

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
# vi:expandtab:tabstop=4 shiftwidth=4:
