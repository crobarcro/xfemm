/*
   Backend factory implementation: explicit instantiations of create_backend
   for the real (double) and complex (CComplex) scalar types.
*/

#include "backend_factory.h"

#ifdef XFEMM_USE_PETSc
#include "PetscLinearSystemBackend.h"
#endif

namespace femm {

template <typename Scalar>
std::unique_ptr<LinearSystemBackend<Scalar>> create_backend(BackendKind kind)
{
    switch (kind)
    {
    case BackendKind::Legacy:
        return std::make_unique<LegacyLinearSystemBackend<Scalar>>();
#ifdef XFEMM_USE_PETSc
    case BackendKind::Petsc:
        return std::make_unique<PetscLinearSystemBackend<Scalar>>();
#endif
    }
    return nullptr;
}

template std::unique_ptr<LinearSystemBackend<double>> create_backend<double>(BackendKind);
template std::unique_ptr<LinearSystemBackend<CComplex>> create_backend<CComplex>(BackendKind);

} // namespace femm
