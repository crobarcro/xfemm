add_library(fsolver STATIC 
    AnalysisSession.cpp
    FSolverAnalysisBackend.cpp
    fsolver.cpp
    harmonic2d.cpp
    harmonicaxi.cpp
    static2d.cpp
    staticaxi.cpp
    )
target_include_directories(fsolver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(fsolver PUBLIC femm fmesher)

add_executable(fsolver-bin
    main.cpp
    )

set_target_properties(fsolver-bin
    PROPERTIES OUTPUT_NAME fsolver)

target_link_libraries(fsolver-bin fsolver)

add_subdirectory(test)

install(
    TARGETS fsolver-bin
    RUNTIME DESTINATION bin
    COMPONENT "cli")
# vi:expandtab:tabstop=4 shiftwidth=4:
