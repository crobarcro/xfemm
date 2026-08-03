add_executable(cspars_multiply_test cspars_multiply_test.cpp)
target_link_libraries(cspars_multiply_test PRIVATE femm)

add_test(NAME cspars_multiply COMMAND cspars_multiply_test)
set_tests_properties(cspars_multiply PROPERTIES LABELS "libfemm;unit")

add_executable(material_curve_test material_curve_test.cpp)
target_link_libraries(material_curve_test PRIVATE femm)

add_test(NAME material_curve COMMAND material_curve_test)
set_tests_properties(material_curve PROPERTIES LABELS "libfemm;unit")
