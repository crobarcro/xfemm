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
    PASS_REGULAR_EXPRESSION "default mesher backend ${XFEMM_MESHER_BACKEND}"
    LABELS "mesher;${XFEMM_MESHER_BACKEND}")
# vi:expandtab:tabstop=4 shiftwidth=4:

add_executable(fmesher_backend_test backend_test.cpp)
target_link_libraries(fmesher_backend_test PRIVATE fmesher)

add_executable(tangle_mesh_converter_test tangle_mesh_converter_test.cpp)
target_link_libraries(tangle_mesh_converter_test PRIVATE fmesher Tangle::tangle)
add_test(NAME fmesher_tangle_mesh_converter COMMAND tangle_mesh_converter_test)
set_tests_properties(fmesher_tangle_mesh_converter PROPERTIES
    LABELS "mesher;Tangle;unit")

add_executable(tangle_backend_execution_test tangle_backend_execution_test.cpp)
target_link_libraries(tangle_backend_execution_test PRIVATE fmesher Tangle::tangle)
add_test(NAME fmesher_tangle_backend_execution COMMAND tangle_backend_execution_test)
set_tests_properties(fmesher_tangle_backend_execution PROPERTIES
    LABELS "mesher;Tangle;unit")

foreach(backend Triangle Tangle)
    add_test(NAME fmesher_${backend}_backend_solver_mesh
        COMMAND fmesher_backend_test "${backend}"
            "${CMAKE_CURRENT_LIST_DIR}/Temp.fem"
            "${CMAKE_CURRENT_LIST_DIR}/../../femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.fem")
    set_tests_properties(fmesher_${backend}_backend_solver_mesh PROPERTIES
        LABELS "magnetics;mesher;${backend};regression"
        TIMEOUT 120)
endforeach()

# Floating-point contraction must stay off where Triangle's input and Triangle
# itself are computed (libfemm and the bundled triangle target). Each probe
# below is one a*b+c compiled with that target's own compile options; the
# control probe is compiled with contraction forced on, to show whether this
# compiler/CPU fuses at all (arm64 Clang and aarch64 GCC do, plain x86-64
# does not).
if(TARGET triangle AND CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$"
        AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    add_library(fp_probe_triangle OBJECT fp_contract_probe.c)
    target_compile_definitions(fp_probe_triangle PRIVATE XFEMM_FP_PROBE_NAME=xfemm_fp_probe_triangle)
    target_compile_options(fp_probe_triangle PRIVATE $<TARGET_PROPERTY:triangle,COMPILE_OPTIONS>)

    add_library(fp_probe_femm OBJECT fp_contract_probe.cpp)
    target_compile_definitions(fp_probe_femm PRIVATE XFEMM_FP_PROBE_NAME=xfemm_fp_probe_femm)
    target_compile_options(fp_probe_femm PRIVATE $<TARGET_PROPERTY:femm,COMPILE_OPTIONS>)

    add_library(fp_probe_control OBJECT fp_contract_probe.c)
    target_compile_definitions(fp_probe_control PRIVATE XFEMM_FP_PROBE_NAME=xfemm_fp_probe_control)
    target_compile_options(fp_probe_control PRIVATE -ffp-contract=fast)

    add_executable(fp_contract_test fp_contract_test.cpp
        $<TARGET_OBJECTS:fp_probe_triangle>
        $<TARGET_OBJECTS:fp_probe_femm>
        $<TARGET_OBJECTS:fp_probe_control>)
    add_test(NAME fmesher_fp_contract_off COMMAND fp_contract_test)
    set_tests_properties(fmesher_fp_contract_off PROPERTIES
        LABELS "mesher;Triangle;unit;regression"
        TIMEOUT 30)
endif()
