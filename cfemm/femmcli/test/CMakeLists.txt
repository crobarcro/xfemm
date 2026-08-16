# required for cmake <= 3.4:
include(CMakeParseArguments)

if(WIN32)
    set(NEWLINE_NATIVE WIN32)
else()
    set(NEWLINE_NATIVE UNIX)
endif()

## test_lua(<name> [LABELS "a;b;c"] [WORKING_DIRECTORY "dir"])
# Add a lua test for <name>.lua.
function(test_lua testname)
    cmake_parse_arguments(test_lua
        "" # options
        "WORKING_DIRECTORY" # oneValueArgs
        "LABELS" # multiValueArgs
        "${ARGN}"
        )
    add_test(NAME ${testname}.lua
        COMMAND femmcli-bin --lua-base-dir "${CMAKE_CURRENT_LIST_DIR}/../debug" --lua-script "${CMAKE_CURRENT_LIST_DIR}/${testname}.lua"
        )
    if(test_lua_WORKING_DIRECTORY)
        set_tests_properties(${testname}.lua PROPERTIES
            WORKING_DIRECTORY "${test_lua_WORKING_DIRECTORY}"
            )
    endif()
    set_tests_properties(${testname}.lua PROPERTIES
        LABELS "lua;${test_lua_LABELS}"
        ENVIRONMENT "XFEMM_TEST_MESHER_BACKEND=${XFEMM_TEST_MESHER_BACKEND}"
        )
endfunction()

## test_lua_setup(<name> <files>...)
# Copy files required by <name>.lua.
function(test_lua_setup testname)
    foreach(file IN LISTS ARGN)
        configure_file("${CMAKE_CURRENT_LIST_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}" @ONLY NEWLINE_STYLE ${NEWLINE_NATIVE})
    endforeach()
endfunction()

## test_lua_check(<name> <suffix> <file> ...)
# Add a test <name>.check.<suffix> that compares an output file of <name>.lua with <file>.
# This kind of test is really prone to false positives.
function(test_lua_check testname suffix file)
    # convert check file to native line endings:
    configure_file("${CMAKE_CURRENT_LIST_DIR}/${file}.check" "${CMAKE_CURRENT_BINARY_DIR}/${file}.check" @ONLY NEWLINE_STYLE ${NEWLINE_NATIVE})
    if(testname STREQUAL "femmcli_femfile" AND suffix STREQUAL "ans")
        add_test(NAME ${testname}.check.${suffix}
            COMMAND "${CMAKE_COMMAND}" -DSOLUTION_FILE=${file}
                -P "${CMAKE_CURRENT_LIST_DIR}/validate_solution.cmake")
    else()
        add_test(NAME ${testname}.check.${suffix}
            COMMAND "${CMAKE_COMMAND}" -E compare_files "${file}" "${file}.check")
    endif()
    set_tests_properties(${testname}.check.${suffix} PROPERTIES DEPENDS ${testname}.lua)
endfunction()

############################################################

### base tests:
test_lua(femmcli_chdir WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")
test_lua(femmcli_compatmode)
test_lua(femmcli_complex)
test_lua(femmcli_pureLua)
test_lua(femmcli_trace)

### magnetics tests:
test_lua(femmcli_femfile LABELS "magnetics;solver")
test_lua_setup(femmcli_femfile "femmcli_femfile.fem")
test_lua_check(femmcli_femfile ans "femmcli_femfile.result.ans")
test_lua_check(femmcli_femfile fem "femmcli_femfile.result.fem")
test_lua(femmcli_fpproc LABELS "magnetics;postprocessor")
test_lua_setup(femmcli_fpproc "femmcli_fpproc.fem")
test_lua(femmcli_matlib LABELS "magnetics")
test_lua(femmcli_moverotate_magdir LABELS "magnetics")
test_lua_check(femmcli_moverotate_magdir fem "femmcli_moverotate_magdir.result.fem")
test_lua_check(femmcli_matlib fem "femmcli_matlib.result.fem")
test_lua(femmcli_TorqueBenchmark LABELS "magnetics;postprocessor;fromWiki")
test_lua_setup(femmcli_TorqueBenchmark "femmcli_TorqueBenchmark.fem")
test_lua(femmcli_antiperiodicBC_flux LABELS "magnetics;postprocessor")
test_lua_setup(femmcli_antiperiodicBC_flux "femmcli_antiperiodicBC_flux.fem")
test_lua(femmcli_antiperiodicBC_AGE_TorqueBenchmark LABELS "magnetics;postprocessor;fromWiki")
test_lua_setup(femmcli_antiperiodicBC_AGE_TorqueBenchmark "femmcli_antiperiodicBC_AGE_TorqueBenchmark.fem")

### electrostatics tests:
test_lua(femmcli_epproc LABELS "electrostatics;postprocessor")
test_lua_setup(femmcli_epproc "femmcli_epproc.fee")

### heatflow tests:
test_lua(femmcli_hpproc LABELS "heatflow;postprocessor")
test_lua_setup(femmcli_hpproc "femmcli_hpproc.feh")

# vi:expandtab:tabstop=4 shiftwidth=4:
