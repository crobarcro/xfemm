add_executable(cspars_multiply_test cspars_multiply_test.cpp)
target_link_libraries(cspars_multiply_test PRIVATE femm)

add_test(NAME cspars_multiply COMMAND cspars_multiply_test)
set_tests_properties(cspars_multiply PROPERTIES LABELS "libfemm;unit")

add_executable(material_curve_test material_curve_test.cpp)
target_link_libraries(material_curve_test PRIVATE femm)

add_test(NAME material_curve COMMAND material_curve_test)
set_tests_properties(material_curve PROPERTIES LABELS "libfemm;unit")

add_executable(solver_mesh_validator_test solver_mesh_validator_test.cpp)
target_link_libraries(solver_mesh_validator_test PRIVATE femm)

add_test(NAME solver_mesh_validator COMMAND solver_mesh_validator_test)
set_tests_properties(solver_mesh_validator PROPERTIES LABELS "libfemm;mesh;unit")

add_executable(instanced_mesh_test instanced_mesh_test.cpp)
target_link_libraries(instanced_mesh_test PRIVATE femm)

add_test(NAME instanced_mesh_materializer COMMAND instanced_mesh_test)
set_tests_properties(instanced_mesh_materializer PROPERTIES
    LABELS "libfemm;mesh;instanced;unit")

add_executable(instanced_mesh_hash_test instanced_mesh_hash_test.cpp)
target_link_libraries(instanced_mesh_hash_test PRIVATE femm)

add_test(NAME instanced_mesh_hashing COMMAND instanced_mesh_hash_test)
set_tests_properties(instanced_mesh_hashing PROPERTIES
    LABELS "libfemm;mesh;instanced;unit")

add_executable(instanced_mesh_airgap_test instanced_mesh_airgap_test.cpp)
target_link_libraries(instanced_mesh_airgap_test PRIVATE femm)

add_test(NAME instanced_mesh_airgap COMMAND instanced_mesh_airgap_test)
set_tests_properties(instanced_mesh_airgap PROPERTIES
    LABELS "libfemm;mesh;instanced;unit")
