add_library(fpproc STATIC
    fpproc.cpp
    makemask.cpp
    CMPointVals.cpp
    CPostProcMElement.cpp
    MagneticSolutionSnapshot.cpp
    InMemorySolution.cpp
    )
target_include_directories(fpproc PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(fpproc PUBLIC femm fsolver)

add_executable(fpproc-test
    main.cpp
    )

target_link_libraries(fpproc-test fpproc)

add_executable(fpproc-snapshot-test test/snapshot_parity_test.cpp)
target_link_libraries(fpproc-snapshot-test PRIVATE fpproc)
add_test(NAME fpproc_snapshot_parity
    COMMAND fpproc-snapshot-test ${CMAKE_BINARY_DIR}/fsolver/test/Temp.ans)
set_tests_properties(fpproc_snapshot_parity PROPERTIES LABELS "magnetics;postprocessor;unit")
set_tests_properties(fpproc_snapshot_parity PROPERTIES DEPENDS fsolver_Temp.solve)

add_test(NAME fpproc_snapshot_age_parity
    COMMAND fpproc-snapshot-test
            ${CMAKE_BINARY_DIR}/femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.result.ans)
set_tests_properties(fpproc_snapshot_age_parity PROPERTIES
    DEPENDS femmcli_antiperiodicBC_AGE_TorqueBenchmark.lua
    LABELS "magnetics;postprocessor;unit;airgap")
install(
    TARGETS fpproc-test
    RUNTIME DESTINATION bin
    COMPONENT "cli (debug)"
    CONFIGURATIONS "Debug")
# vi:expandtab:tabstop=4 shiftwidth=4:
