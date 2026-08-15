function(test_fmesher file)
    set(test_dir "${CMAKE_CURRENT_BINARY_DIR}/${file}.work")
    file(MAKE_DIRECTORY "${test_dir}")
    configure_file("${CMAKE_CURRENT_LIST_DIR}/${file}" "${test_dir}/${file}" COPYONLY)
    add_test(NAME fmesher_${file}
        COMMAND fmesher-bin "${test_dir}/${file}"
        )
    set_tests_properties(fmesher_${file} PROPERTIES
        WORKING_DIRECTORY "${test_dir}"
        LABELS "magnetics;mesher"
        )
endfunction()

test_fmesher(Temp.fem)
test_fmesher(split_seg_err_test.fem)

add_test(NAME fmesher_selected_backend
    COMMAND fmesher-bin --version)
set_tests_properties(fmesher_selected_backend PROPERTIES
    PASS_REGULAR_EXPRESSION "default mesher backend ${XFEMM_TEST_MESHER_BACKEND}"
    LABELS "mesher;${XFEMM_TEST_MESHER_BACKEND}")
# vi:expandtab:tabstop=4 shiftwidth=4:

add_executable(fmesher_backend_test backend_test.cpp)
target_link_libraries(fmesher_backend_test PRIVATE fmesher)
foreach(backend Triangle Tangle)
    add_test(NAME fmesher_${backend}_backend_solver_mesh
        COMMAND fmesher_backend_test "${backend}"
            "${CMAKE_CURRENT_LIST_DIR}/Temp.fem"
            "${CMAKE_CURRENT_LIST_DIR}/../../femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.fem")
    set_tests_properties(fmesher_${backend}_backend_solver_mesh PROPERTIES
        LABELS "magnetics;mesher;${backend}")
endforeach()
