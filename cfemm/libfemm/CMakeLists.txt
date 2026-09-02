# lua with complex numbers
add_library(luacomplex
    liblua/femmcomplex.cpp
    liblua/llimits.h
    liblua/lopcodes.h
    liblua/lua.h
    liblua/luadebug.h
    liblua/lualib.h
    liblua/lapi.cpp
    liblua/lauxlib.cpp
    liblua/lbaselib.cpp
    liblua/lcode.cpp
    liblua/ldblib.cpp
    liblua/ldebug.cpp
    liblua/ldo.cpp
    liblua/lfunc.cpp
    liblua/lgc.cpp
    liblua/liolib.cpp
    liblua/llex.cpp
    liblua/lmathlib.cpp
    liblua/lmem.cpp
    liblua/lobject.cpp
    liblua/lparser.cpp
    liblua/lstate.cpp
    liblua/lstring.cpp
    liblua/lstrlib.cpp
    liblua/ltable.cpp
    liblua/ltests.cpp
    liblua/ltm.cpp
    liblua/lundump.cpp
    liblua/lvm.cpp
    liblua/lzio.cpp
    )
target_include_directories(luacomplex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/liblua> $<INSTALL_INTERFACE:include/liblua>)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/femmversion.h.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/femmversion.h"
    @ONLY
    )

add_library(femm
    femmconstants.cpp
    femmenums.cpp
    CArcSegment.cpp
    CBlockLabel.cpp
    CBoundaryProp.cpp
    CCircuit.cpp
    CCommonPoint.cpp
    CElement.cpp
    CAirGapElement.cpp
    CliTools.cpp
    CMaterialProp.cpp
    CMeshNode.cpp
    CNode.cpp
    CPointProp.cpp
    CSegment.cpp
    cspars.cpp
    cuthill.cpp
    feasolver.cpp
    FemmProblem.cpp
    FemmReader.cpp
    FemmStateBase.cpp
    femmversion.cpp
    fparse.cpp
    fullmatrix.cpp
    IntPoint.cpp
    locationTools.cpp
    LuaInstance.cpp
    MatlibReader.cpp
    PostProcessor.cpp
    spars.cpp
    stringTools.cpp
    linsolve/backend_factory.cpp
    mesh/SolverMeshValidator.cpp
    )
if(XFEMM_USE_PETSc)
    target_sources(femm PRIVATE linsolve/PetscLinearSystemBackend.cpp)
endif()
target_include_directories(femm
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
target_link_libraries(femm PUBLIC luacomplex)

# ---------------------------------------------------------------------------
# Optional PETSc linear solver backend
# ---------------------------------------------------------------------------
option(XFEMM_USE_PETSc "Use the PETSc linear solver backend" OFF)
if(XFEMM_USE_PETSc)
    # Locate a PETSc build with complex scalars.  The Debian/Ubuntu package
    # layout puts it in /usr/lib/petscdir/petsc<VERSION>/<arch>-complex; the
    # PETSC_DIR/PETSC_ARCH environment variables are honoured as well.
    set(_petsc_hints)
    file(GLOB _petsc_complex_dirs "/usr/lib/petscdir/petsc*/x86_64-linux-gnu-complex")
    list(APPEND _petsc_hints ${_petsc_complex_dirs})
    if(DEFINED ENV{PETSC_DIR})
        if(DEFINED ENV{PETSC_ARCH})
            list(APPEND _petsc_hints "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
        else()
            list(APPEND _petsc_hints "$ENV{PETSC_DIR}")
        endif()
    endif()
    find_path(PETSC_INCLUDE_DIR petsc.h PATHS ${_petsc_hints} PATH_SUFFIXES include NO_DEFAULT_PATH)
    find_library(PETSC_LIBRARY petsc_complex PATHS ${_petsc_hints} PATH_SUFFIXES lib NO_DEFAULT_PATH)
    if(NOT PETSC_INCLUDE_DIR OR NOT PETSC_LIBRARY)
        message(FATAL_ERROR
            "XFEMM_USE_PETSc was requested but no complex-scalar PETSc build was found. "
            "Install the PETSc complex development package (libpetsc-complex-dev) or set "
            "PETSC_DIR/PETSC_ARCH.")
    endif()
    message(STATUS "Using PETSc (complex) headers: ${PETSC_INCLUDE_DIR}")
    message(STATUS "Using PETSc (complex) library: ${PETSC_LIBRARY}")
    # PETSc headers require MPI headers.
    find_package(MPI QUIET COMPONENTS CXX)
    if(MPI_CXX_FOUND)
        set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
        set(PETSC_LINK_LIBRARIES ${PETSC_LIBRARY} ${MPI_CXX_LIBRARIES})
    else()
        set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR})
        set(PETSC_LINK_LIBRARIES ${PETSC_LIBRARY})
    endif()
    target_include_directories(femm PRIVATE ${PETSC_INCLUDE_DIRS})
    target_compile_definitions(femm PRIVATE XFEMM_USE_PETSc)
    target_link_libraries(femm PUBLIC ${PETSC_LINK_LIBRARIES})
endif()

add_subdirectory(test)
# vi:expandtab:tabstop=4 shiftwidth=4:
