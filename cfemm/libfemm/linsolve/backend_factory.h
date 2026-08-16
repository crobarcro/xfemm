/*
   Backend factory: selects the linear system backend implementation at
   runtime.

   The available implementations are:
     - LegacyLinearSystemBackend (the original linked-list solver, kept as
       the validation oracle);
     - PetscLinearSystemBackend (optional, enabled at build time with
       -DXFEMM_USE_PETSc=ON).

   The selection can be overridden with the XFEMM_SOLVER_BACKEND environment
   variable; unknown values fall back to the legacy backend.
*/

#ifndef XFEMM_LINSOLVE_BACKEND_FACTORY_H
#define XFEMM_LINSOLVE_BACKEND_FACTORY_H

#include "LegacyLinearSystemBackend.h"
#include "LinearSystemBackend.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

namespace femm {

enum class BackendKind { Legacy, Petsc };

/** Backend selected by the XFEMM_SOLVER_BACKEND environment variable
 *  (defaults to Legacy). */
inline BackendKind default_backend_kind()
{
    BackendKind kind = BackendKind::Legacy;
    const char *env = std::getenv("XFEMM_SOLVER_BACKEND");
    if (env != nullptr && *env != '\0')
    {
        std::string value(env);
        for (char &c : value)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (value == "legacy")
            kind = BackendKind::Legacy;
        else if (value == "petsc")
            kind = BackendKind::Petsc;
        else
            std::fprintf(stderr,
                         "xfemm: unknown XFEMM_SOLVER_BACKEND '%s'; using 'legacy'\n",
                         env);
    }
    return kind;
}

/** Create a backend of the given kind.
 *
 *  Returns nullptr if the requested kind is not available in this build
 *  (e.g. BackendKind::Petsc without -DXFEMM_USE_PETSc).  Defined in
 *  backend_factory.cpp with explicit instantiations for double and CComplex.
 */
template <typename Scalar>
std::unique_ptr<LinearSystemBackend<Scalar>> create_backend(BackendKind kind);

} // namespace femm

#endif // XFEMM_LINSOLVE_BACKEND_FACTORY_H
