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

add_executable(copy_transform_test copy_transform_test.cpp)
target_link_libraries(copy_transform_test PRIVATE femm)

add_test(NAME copy_transform COMMAND copy_transform_test)
set_tests_properties(copy_transform PROPERTIES LABELS "libfemm;geometry;unit;regression")
