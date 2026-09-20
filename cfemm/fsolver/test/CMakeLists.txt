if(WIN32)
    set(NEWLINE_NATIVE WIN32)
else()
    set(NEWLINE_NATIVE UNIX)
endif()
## test_fsolver(<name> <premeshed>)
# Add tests to solve <name>.fem, and to compare the output file against a stored version:
# 1. Run fmesher (if <premeshed> is false)
# 2. Run fsolver <name>
# 3. Add disabled check to compare <name>.ans against <name>.ans.check
function(test_fsolver name premeshed)
    # convert check file to native line endings:
    configure_file("${CMAKE_CURRENT_LIST_DIR}/${name}.fem" "${CMAKE_CURRENT_BINARY_DIR}/${name}.fem" @ONLY NEWLINE_STYLE ${NEWLINE_NATIVE})
    configure_file("${CMAKE_CURRENT_LIST_DIR}/${name}.ans.check" "${CMAKE_CURRENT_BINARY_DIR}/${name}.ans.check" @ONLY NEWLINE_STYLE ${NEWLINE_NATIVE})
    if(premeshed)
        # copy files at test time because the solver deletes them
        # (otherwise you would need to re-run cmake after each ctest run)
        add_test(NAME fsolver_${name}.setup
            COMMAND "${CMAKE_COMMAND}"
            -E copy
            "${CMAKE_CURRENT_LIST_DIR}/${name}.edge"
            "${CMAKE_CURRENT_LIST_DIR}/${name}.ele"
            "${CMAKE_CURRENT_LIST_DIR}/${name}.node"
            "${CMAKE_CURRENT_LIST_DIR}/${name}.pbc"
            #"${CMAKE_CURRENT_LIST_DIR}/${name}.poly"
            .
            )
    else()
        add_test(NAME fsolver_${name}.mesh
            COMMAND fmesher-bin "${name}.fem"
            )
        set_tests_properties(fsolver_${name}.mesh PROPERTIES
            LABELS "magnetics;mesher"
            )
    endif()
    add_test(NAME fsolver_${name}.solve
        COMMAND fsolver-bin "${name}"
        )
    set_tests_properties(fsolver_${name}.solve PROPERTIES
        LABELS "magnetics;solver"
        )
    if(premeshed)
        set_tests_properties(fsolver_${name}.solve PROPERTIES DEPENDS fsolver_${name}.setup)
    else()
        set_tests_properties(fsolver_${name}.solve PROPERTIES DEPENDS fsolver_${name}.mesh)
    endif()

    # DISABLED property was introduced in cmake 3.9
    # skip test alltogether in earlier versions
    if(NOT "${CMAKE_VERSION}" VERSION_LESS 3.9)
        add_test(NAME fsolver_${name}.check
            COMMAND "${CMAKE_COMMAND}"
            -E compare_files "${name}.ans" "${name}.ans.check"
            )
        set_tests_properties(fsolver_${name}.check PROPERTIES
            DEPENDS fsolver_${name}.solve
            LABELS "magnetics"
            )
        # disable hair-trigger comparison by default
        set_tests_properties(fsolver_${name}.check PROPERTIES DISABLED TRUE)
    endif()
endfunction()

test_fsolver(Temp TRUE)
test_fsolver(Temp1 FALSE)

add_executable(analysis_session_test analysis_session_test.cpp)
target_link_libraries(analysis_session_test PRIVATE fsolver)
add_test(NAME analysis_session COMMAND analysis_session_test)
set_tests_properties(analysis_session PROPERTIES LABELS "magnetics;solver;unit")

add_executable(analysis_session_instanced_test analysis_session_instanced_test.cpp)
target_link_libraries(analysis_session_instanced_test PRIVATE fsolver)
add_test(NAME analysis_session_instanced COMMAND analysis_session_instanced_test)
set_tests_properties(analysis_session_instanced PROPERTIES
    LABELS "magnetics;solver;instanced;unit")

add_executable(analysis_session_machine_test analysis_session_machine_test.cpp)
target_link_libraries(analysis_session_machine_test PRIVATE fsolver)
add_test(NAME analysis_session_machine COMMAND analysis_session_machine_test)
set_tests_properties(analysis_session_machine PROPERTIES
    LABELS "magnetics;solver;instanced;unit")

add_executable(tiled_model_airgap_test tiled_model_airgap_test.cpp)
target_link_libraries(tiled_model_airgap_test PRIVATE fsolver)
add_test(NAME tiled_model_airgap COMMAND tiled_model_airgap_test)
set_tests_properties(tiled_model_airgap PROPERTIES
    LABELS "solver;Tangle;tiled;instanced;unit")

add_executable(native_static2d_test native_static2d_test.cpp)
target_link_libraries(native_static2d_test PRIVATE fsolver)
add_test(NAME native_static2d COMMAND native_static2d_test)
set_tests_properties(native_static2d PROPERTIES
    LABELS "magnetics;solver;instanced;native;unit")

add_executable(instanced_solver_benchmark instanced_solver_benchmark.cpp)
target_link_libraries(instanced_solver_benchmark PRIVATE fsolver)

add_executable(periodic_pbc_without_age_test periodic_pbc_without_age_test.cpp)
target_link_libraries(periodic_pbc_without_age_test PRIVATE fsolver)
add_test(NAME periodic_pbc_without_age
    COMMAND periodic_pbc_without_age_test "${CMAKE_CURRENT_LIST_DIR}/Temp")
set_tests_properties(periodic_pbc_without_age PROPERTIES
    LABELS "magnetics;solver;unit;periodic")
# vi:expandtab:tabstop=4 shiftwidth=4:
