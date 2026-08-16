#ifndef FEMM_FMESHER_TANGLEMESHCONVERTER_H
#define FEMM_FMESHER_TANGLEMESHCONVERTER_H

#include "mesh/Meshing.h"

struct Mesh;

namespace fmesher {

/**
 * Convert Tangle's value mesh into xfemm's solver transfer representation.
 *
 * Tangle coordinates use the input problem's length unit.  The caller supplies
 * the number of metres per input unit so the returned SolverMesh always obeys
 * its metre-based coordinate contract.
 */
femm::mesh::MeshResult convertTangleMesh(const ::Mesh &source,
                                         double lengthScaleToMetres);

} // namespace fmesher

#endif
