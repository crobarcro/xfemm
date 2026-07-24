# XFEMM Repository Guide for AI Agents

## Project Overview

**xfemm** is a cross-platform, open-source magnetics finite element solver providing:
- Direct C++ command-line interfaces to a FEMM-based magnetics solver
- Native MATLAB/Octave bindings for problem definition and post-processing
- MEX interface for compiled C++ solver integration with MATLAB/Octave

**Primary Purpose**: Enable programmable access to high-quality 2D planar and axisymmetric finite element analysis for magnetics, particularly from MATLAB/Octave environments.

**Key Technologies**:
- C++14 with STL (core solver, mesher, post-processor)
- MATLAB/Octave with M-files and MEX interfaces
- Lua scripting (embedded in C++ for automation)
- CMake build system
- Triangle mesh generator (embedded version 1.6)

**Repository Statistics**:
- Language composition: C++ (48.3%), Edge (18.6%), C (16.2%), MATLAB (14.7%), Lua (1.2%), CMake (0.5%)
- Open issues: ~9
- Created: September 2020
- License: Aladdin Free Public License (original FEMM), Apache 2.0 (xfemm M-files)

---

## Directory Structure

### `/cfemm/` - Core C++ Finite Element System

**Purpose**: Complete C++ implementation of FEMM's magnetics solver, mesher, and post-processor. All binaries compile from here.

**Key Subdirectories**:

#### `/cfemm/libfemm/` - Shared Libraries
Core data structures and algorithms shared by all components:

- **PostProcessor.h/cpp** - Base class (`PostProcessor`, interface `PProcIface`) for reading and analyzing FEM solutions
  - Mesh element and node management
  - Field calculations at points
  - Contour integration
  - Selection and masking operations
  
- **FemmProblem.h/cpp** - Complete problem definition container
  - Node (`CNode`), line segment (`CSegment`), arc segment (`CArcSegment`) management
  - Block label management (`CBlockLabel`)
  - Material, boundary, and circuit properties
  - Geometric operations: rotation, scaling, mirroring, translation
  - Intersection detection and undo/redo support

- **libfemm/liblua/** - Embedded Lua 4.0 interface
  - Provides scripting access to solver functions
  - Custom extensions for FEMM-specific types (complex numbers)

- **Other key headers**:
  - `CMaterialProp.h` - Material properties for blocks
  - `CBoundaryProp.h` - Boundary conditions for lines/arcs
  - `CCircuit.h` - Circuit definitions
  - `CMeshNode.h`, `CElement.h` - Mesh primitives
  - `femmconstants.h` - Physical and mathematical constants
  - `femmcomplex.h` - Complex number support
  - `femmenums.h` - Enumerations for problem types, coordinate systems, etc.

#### `/cfemm/fmmesher/` - Mesh Generation
Creates triangular meshes from problem geometry:
- Input: `.fem` file problem description
- Output: Mesh data written to files
- Uses Triangle library (embedded)

#### `/cfemm/fsolver/` - Magnetostatics Solver
Solves 2D magnetostatic problems:
- Planar and axisymmetric geometries
- Linear and nonlinear materials
- Boundary conditions and loads

#### `/cfemm/esolver/` - Electrostatics Solver
Solves 2D electrostatics/current flow problems

#### `/cfemm/hsolver/` - Heat Flow Solver
Solves 2D heat transfer problems

#### `/cfemm/fpproc/`, `/cfemm/epproc/`, `/cfemm/hpproc/` - Post-Processors
Read solution files (`.ans`) and compute derived quantities (fields, forces, etc.):
- **fpproc** - Magnetics post-processor
- **epproc** - Electrostatics/current flow post-processor
- **hpproc** - Heat flow post-processor

#### `/cfemm/femmcli/` - Lua Command Line Interface
High-level interface wrapping solvers in Lua:
- Enables scriptable problem solving
- Used by MATLAB interface
- Example: `femmcli script.lua`

#### `/cfemm/CMakeLists.txt`
- Orchestrates entire build system
- Defines version from `/cfemm/VERSION` file
- Compiles all executables: `fmesher`, `fsolver`, `esolver`, `hsolver`, `fpproc`, `epproc`, `hpproc`
- Output binaries to `/cfemm/bin/`
- C++ standard: C++14
- Compiler flags: `-Wall -Wextra -Wpedantic`
- Debug flags available: `DEBUG_FEMMLUA`, `DEBUG_FEMMCLI`, `DEBUG_PARSER`

---

### `/mfemm/` - MATLAB/Octave Interface

**Purpose**: High-level, user-friendly interface to xfemm from MATLAB/Octave.

**Key Subdirectories**:

#### `/mfemm/cores/` - MEX Interfaces (C++ compiled functions)
Bridges MATLAB/Octave to C++ solvers and post-processors:
- Compiled into `.mexw64`, `.mexa64`, etc. depending on platform
- Wraps `fmesher`, `fsolver`, and post-processor libraries
- Requires `mex -setup C++` before compilation

#### `/mfemm/classes/` - MATLAB Class System
Object-oriented wrappers:
- Problem definition classes
- Solution classes
- Property classes (materials, boundaries)

#### `/mfemm/examples/` - Tutorials and Examples
Heavily commented M-files demonstrating:
- Problem creation and geometry setup
- Solver invocation
- Post-processing and visualization
- Pattern: `femm_*_tutorial_example.m`

#### `/mfemm/postproc/` - Post-Processing Functions
M-file implementations for analyzing solutions:
- Field extraction functions
- Visualization (contours, vector plots)
- Derived quantity calculations

#### `/mfemm/preproc/` - Pre-Processing Functions
M-file utilities for problem setup:
- Geometry creation helpers
- Material property management
- Boundary condition application

#### `/mfemm/README.txt`
- Package documentation and licensing (Apache 2.0 for M-files)
- Installation instructions reference: `help mfemm_setup`

#### `/mfemm/INSTALLATION.txt`
- Detailed installation steps
- Key function: `mfemm_setup('ForceMexRecompile', true)` to recompile MEX files
- Requires: MATLAB/Octave with C++ compiler configured

#### `/mfemm/mfemm_setup.m` - Setup Script (primary entry point)
Handles:
- Path configuration
- MEX compilation from source
- Validation and testing
- Parameters: `ForceMexRecompile`, `RunTests`, `DoDebug`, `Verbose`

---

### `/test/` - Test Suite

**Purpose**: Validation tests for solver correctness and regression detection.

- Tests cover mesher, solver, and post-processor components
- Executed by CI/CD pipeline (GitHub Actions)
- Platform coverage: Windows (MinGW), Linux

---

### Root-Level Configuration Files

- **CMakeLists.txt** (top-level) - Main build orchestrator
- **README.md** - User-facing introduction and compilation instructions
- **CHANGELOG.txt** - Version history and feature changes
- **release.sh**, **test_release.sh** - Release automation scripts

---

## Key Data Structures

### Problem Definition (`cfemm/libfemm/`)

```
FemmProblem (container)
├── nodelist (vector<CNode>)
├── linelist (vector<CSegment>)
├── arclist (vector<CArcSegment>)
├── labellist (vector<CBlockLabel>)
├── blockproplist (vector<CMaterialProp>)
├── lineproplist (vector<CBoundaryProp>)
├── circproplist (vector<CCircuit>)
└── nodeproplist (vector<CPointProp>)
```

### Mesh Representation

```
PostProcessor (read solution)
├── meshnodes (vector<CMeshNode>)
└── meshelems (vector<CElement>)
```

### Solution Types

- **Electromagnetics**: Magnetic field (A), permeability, force/torque
- **Electrostatics**: Electric potential (V), permittivity, capacitance
- **Current Flow**: Voltage potential, conductivity, current
- **Heat Flow**: Temperature, conductivity, heat transfer

---

## Build System

### Compilation (C++ Components)

**Linux/Mac**:
```bash
cd cfemm
cmake .
make
# Binaries in cfemm/bin/
```

**Windows** (MinGW or MSVC):
```bash
cd cfemm
cmake . -G "MinGW Makefiles"  # or -G "Visual Studio 16 2019"
cmake --build .
```

**CMake Features**:
- Optional external Triangle support (auto-detected or via `-DCMAKE_PREFIX_PATH`)
- Optional debug symbols
- CPack support for binary distribution
- Version management from `cfemm/VERSION` file

### MEX Compilation (MATLAB/Octave)

**Primary entry point**: `mfemm_setup.m`

```matlab
mfemm_setup('ForceMexRecompile', true)  % Forces recompilation of MEX files
```

**Prerequisites**:
- MATLAB R2018a+ or Octave 5+
- C++ compiler configured: `mex -setup C++`
- Compiled C++ libraries from `/cfemm/` must be available

---

## File I/O Formats

### Problem Definition Files (`.fem`)

- Human-readable text format
- Define geometry (nodes, segments, arcs, block labels)
- Material and boundary properties
- Solver settings

### Solution Files (`.ans`)

- Binary format (endian-dependent)
- Contains mesh and field solutions
- Read by post-processors (`fpproc`, `epproc`, `hpproc`)
- Accessible via MATLAB via MEX interface

### Lua Scripts

- User-definable FEMM commands via Lua syntax
- Executed by `femmcli` or embedded in other tools
- Example: `mi_moverotate(cx, cy, angle, group_id)` - geometric transformation

---

## Common Workflows

### Use Case 1: Standalone Solver (Command Line)

1. Create `.fem` file (MATLAB/Octave or manually)
2. Run: `fmesher problem.fem`
3. Run: `fsolver problem.fem`
4. Outputs: solution file (`problem.ans`)
5. Post-process: `fpproc problem.ans` or custom script

### Use Case 2: MATLAB/Octave Scripting

1. Call `mfemm_setup()` (once per session)
2. Create problem using M-file functions
3. Invoke solver via MEX interface
4. Extract and visualize results using post-processing functions

### Use Case 3: Automated Design Optimization

1. Parameterize geometry in MATLAB
2. Loop: modify geometry → solve → extract metrics
3. Optimize design parameters

---

## Common Issues & Known Limitations

### Known Issues (from GitHub)

1. **Issue #22**: Compilation failures on Windows, MEX setup issues
   - `mfemm_setup` may fail to recompile MEX files
   - Requires manual C++ compiler setup: `mex -setup C++`
   - Solution documentation is outdated

2. **Issue #19**: Lua `mi_moverotate()` in femmcli doesn't rotate magnetization angles
   - Affects geometric transformations of magnets
   - Works correctly in native FEMM GUI but not CLI

### Limitations

- 2D problems only (planar or axisymmetric)
- No direct 3D support
- Triangle mesh generator is older version (1.6 from 2005)
- Lua version is outdated (4.0)
- Windows compilation primarily tested with MinGW

---

## Code Quality & Style

- **Language**: C++ with STL (C++14 standard)
- **Naming**: CamelCase for classes (e.g., `CNode`, `CSegment`), snake_case for functions
- **Memory**: Modern C++ smart pointers (`unique_ptr`, `shared_ptr`)
- **Testing**: Regression test suite in `/test/`
- **Compiler Warnings**: Strict flags (`-Wall -Wextra -Wpedantic`)
- **Documentation**: Doxygen-compatible comments in headers

---

## Extending the Repository

### Adding a New Solver Type

1. Create solver class in `/cfemm/<new_solver>/`
2. Inherit from common base interface if available
3. Implement: meshing, matrix assembly, solving
4. Add post-processor in `/cfemm/<new_pproc>/` for result analysis
5. Add CMake target in `/cfemm/CMakeLists.txt`
6. Implement MEX wrapper in `/mfemm/cores/` for MATLAB access

### Adding MATLAB Functions

1. Create `.m` file in appropriate `/mfemm/preproc/`, `/mfemm/postproc/`, or `/mfemm/classes/`
2. Add comprehensive help documentation (viewable via `help function_name`)
3. Add example usage to `/mfemm/examples/`
4. Update `/mfemm/README.txt` if significant new feature

### Modifying Core Solver

1. Edit `/cfemm/libfemm/` headers and implementations
2. Rebuild: `cd cfemm && cmake . && make`
3. Recompile MEX if interface changed: `mfemm_setup('ForceMexRecompile', true)`
4. Run regression tests: test suite in `/test/`

---

## Key Contributors & Licenses

- **Original FEMM**: David Meeker (Aladdin Free Public License)
- **xfemm C++ Port**: Richard Crozier, Johannes Zarl-Zierl
- **M-files (mfemm)**: Apache 2.0
- **Triangle**: Jonathan Shewchuk (Carnegie Mellon University)
- **Lua**: PUC-Rio (embedded version)

---

## References & External Resources

- **FEMM Official**: http://www.femm.info
- **Citation**: Crozier, R, Mueller, M., "A New MATLAB and Octave Interface to a Popular Magnetics Finite Element Code", ICEM 2016
- **External Triangle**: https://github.com/wo80/Triangle (newer, modular version)

---

## Quick Debugging Checklist

When investigating issues:

1. **Compilation fails**: Check C++ compiler version (C++14 required), CMake version (3.0+)
2. **MEX compile fails**: Run `mex -setup C++`, verify compiler path, check `/mfemm/INSTALLATION.txt`
3. **Solver crashes**: Enable debug flags in CMake, check `.fem` file consistency via `FemmProblem::consistencyCheckOK()`
4. **Post-processor fails**: Verify `.ans` file format matches solver type (magnetics vs. electrostatics)
5. **Lua errors**: Enable `DEBUG_FEMMLUA` in CMake, check Lua syntax in `.lua` scripts
6. **Geometry issues**: Use `FemmProblem::enforcePSLG()` to detect and fix overlaps/intersections

---

**Document Version**: 1.0 (2026-07-24)  
For detailed API documentation, refer to header files in `/cfemm/libfemm/` and built-in MATLAB help via `help <function_name>`.
