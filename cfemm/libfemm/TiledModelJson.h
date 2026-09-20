#ifndef FEMM_TILEDMODELJSON_H
#define FEMM_TILEDMODELJSON_H

#include "TiledModel.h"

#include <string>
#include <vector>

namespace femm {
namespace tiled {

/**
 * Parse a tiled magnetic problem from its JSON representation. Structural
 * errors (bad JSON, missing/unknown names, out-of-range values) are appended to
 * \p errors. Semantic tiling checks are separate: call validateTiledModel()
 * after a successful load.
 *
 * Returns true when the document parsed and every reference resolved.
 */
bool loadTiledModelJson(const std::string &text, TiledModel &out,
                        std::vector<TiledDiagnostic> &errors);

/** Serialize a tiled model to the JSON representation. */
std::string saveTiledModelJson(const TiledModel &model);

} // namespace tiled
} // namespace femm

#endif
