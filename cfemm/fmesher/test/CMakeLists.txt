# only a function test, not correctness :-|
function(test_fmesher file)
    add_test(NAME fmesher_${file}.setup
        COMMAND "${CMAKE_COMMAND}"
        -E copy "${CMAKE_CURRENT_LIST_DIR}/${file}" .
        )
    add_test(NAME fmesher_${file}
        COMMAND fmesher-bin "${file}"
        )
    set_tests_properties(fmesher_${file} PROPERTIES
        DEPENDS fmesher_${file}.setup
        LABELS "magnetics;mesher"
        )
endfunction()

test_fmesher(Temp.fem)
test_fmesher(split_seg_err_test.fem)
# vi:expandtab:tabstop=4 shiftwidth=4:

add_executable(fmesher_backend_test backend_test.cpp)
target_link_libraries(fmesher_backend_test PRIVATE fmesher)
add_test(NAME fmesher_backend_solver_mesh
    COMMAND fmesher_backend_test "${CMAKE_CURRENT_LIST_DIR}/Temp.fem")
set_tests_properties(fmesher_backend_solver_mesh PROPERTIES LABELS "magnetics;mesher")
