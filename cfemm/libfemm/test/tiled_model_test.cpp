#include "TiledModel.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

using namespace femm;
using namespace femm::tiled;

constexpr double Pi = 3.141592653589793238462643383279502884;

#define REQUIRE(condition)                                                          \
    do {                                                                            \
        if (!(condition)) {                                                         \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << #condition << '\n'; \
            return false;                                                           \
        }                                                                           \
    } while (false)

bool hasCategory(const TiledValidationResult &result, TiledDiagnosticCategory category)
{
    for (const auto &diagnostic : result.diagnostics)
        if (diagnostic.category == category)
            return true;
    return false;
}

TiledModel makeBaseModel()
{
    TiledModel model;
    model.depth = 0.0889;

    auto seam = std::make_unique<CMBoundaryProp>();
    seam->BdryName = "seam";
    seam->BdryFormat = 4; // periodic
    model.boundaryProps.push_back(std::move(seam));

    auto outer = std::make_unique<CMBoundaryProp>();
    outer->BdryName = "outer";
    outer->BdryFormat = 0; // Dirichlet
    model.boundaryProps.push_back(std::move(outer));

    auto gap = std::make_unique<CMBoundaryProp>();
    gap->BdryName = "sliding-gap";
    gap->BdryFormat = 6; // periodic AGE
    model.boundaryProps.push_back(std::move(gap));

    auto material = std::make_unique<CMMaterialProp>();
    material->BlockName = "copper";
    material->mu_x = material->mu_y = 1;
    model.materialProps.push_back(std::move(material));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    model.circuitProps.push_back(std::move(circuit));
    return model;
}

/** An annular sector tile with radial seams at 0 and sectorDeg. */
Tile annularTile(double sectorDeg, std::size_t count, Closure closure)
{
    const double sector = sectorDeg * Pi / 180.0;
    const double ri = 0.05;
    const double ro = 0.06;

    Tile tile;
    tile.name = "slot";
    tile.repeat.count = count;
    tile.repeat.closure = closure;
    tile.seamBoundary = "seam";
    tile.airGapBoundary = "sliding-gap";

    tile.geometry.nodes = {
        {ri, 0.0, 0, 0},
        {ro, 0.0, 0, 0},
        {ri * std::cos(sector), ri * std::sin(sector), 0, 0},
        {ro * std::cos(sector), ro * std::sin(sector), 0, 0},
    };
    // Radial seams (boundary 0 = "seam"), outer/inner arcs (boundary 1 = "outer").
    tile.geometry.segments = {{0, 1, -1, 0, false, 0}, {2, 3, -1, 0, false, 0}};
    tile.geometry.arcs = {
        {1, 3, 10.0, 2.0, 1, false, 0},
        {2, 0, 10.0, 2.0, 1, false, 0},
    };
    tile.geometry.labels = {{"coil", 0.055, 0.004, 0, 0, 1, 0.0, "", 0.0, 0, false}};
    return tile;
}

bool testClosedRingDerivesPitch()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 36, Closure::Closed));

    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(result.tiles.size() == 1);
    REQUIRE(std::fabs(result.tiles[0].pitchDegrees - 10.0) < 1e-9);
    REQUIRE(std::fabs(result.tiles[0].totalAngleDegrees - 360.0) < 1e-6);
    REQUIRE(!hasCategory(result, TiledDiagnosticCategory::SeamNotRadial));
    REQUIRE(!hasCategory(result, TiledDiagnosticCategory::GeometryOutsideSector));
    REQUIRE(!hasCategory(result, TiledDiagnosticCategory::ClosureCountMismatch));
    return true;
}

bool testOpenSector()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 6, Closure::Periodic));
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(result.succeeded());
    REQUIRE(std::fabs(result.tiles[0].totalAngleDegrees - 60.0) < 1e-6);
    return true;
}

bool testOverlapRejected()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 37, Closure::Closed));
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::ClosureCountMismatch));
    return true;
}

bool testNonRadialSeamRejected()
{
    TiledModel model = makeBaseModel();
    Tile tile = annularTile(10.0, 36, Closure::Closed);
    // Move the second node off the radial line through the centre.
    tile.geometry.nodes[1].y = 0.002;
    model.tiles.push_back(tile);
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::SeamNotRadial));
    return true;
}

bool testGeometryOutsideSectorRejected()
{
    TiledModel model = makeBaseModel();
    Tile tile = annularTile(10.0, 36, Closure::Closed);
    tile.geometry.nodes.push_back({0.055 * std::cos(0.4), 0.055 * std::sin(0.4), 0, 0});
    model.tiles.push_back(tile);
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::GeometryOutsideSector));
    return true;
}

bool testOverrideLengthMismatch()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 6, Closure::Periodic));
    TileLabelOverride override;
    override.tile = "slot";
    override.label = "coil";
    override.magDir = {0.0, 180.0}; // wrong length (count is 6)
    model.overrides.push_back(override);
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::OverrideLengthMismatch));
    return true;
}

bool testUnknownOverrideLabel()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 6, Closure::Periodic));
    TileLabelOverride override;
    override.tile = "slot";
    override.label = "does-not-exist";
    model.overrides.push_back(override);
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::UnknownLabel));
    return true;
}

bool testBuildTileProblem()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 6, Closure::Periodic));
    const std::unique_ptr<FemmProblem> problem = buildTileProblem(model, 0);
    REQUIRE(problem != nullptr);
    REQUIRE(problem->getTitle().empty());
    REQUIRE(problem->nodelist.size() == 4);
    REQUIRE(problem->linelist.size() == 2);
    REQUIRE(problem->arclist.size() == 2);
    REQUIRE(problem->labellist.size() == 1);
    REQUIRE(problem->lineproplist.size() == 3);
    REQUIRE(problem->blockproplist.size() == 1);
    REQUIRE(problem->circproplist.size() == 1);
    REQUIRE(problem->nodelist[0]->BoundaryMarker == 0);
    REQUIRE(problem->linelist[0]->BoundaryMarker == 0);
    REQUIRE(problem->labellist[0]->BlockType == 0);
    REQUIRE(problem->labellist[0]->InCircuit == 0);
    return true;
}

bool testCouplingConsistency()
{
    TiledModel model = makeBaseModel();
    model.tiles.push_back(annularTile(10.0, 36, Closure::Closed));
    Tile rotor = annularTile(10.0, 36, Closure::Closed);
    rotor.name = "pole";
    model.tiles.push_back(rotor);
    // Centre mismatch and non-AGE boundary.
    model.couplings.push_back({"slot", "pole", "outer", 1.0, 0.0, 0.06, 0.07});
    const TiledValidationResult result = validateTiledModel(model);
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::CouplingBoundaryMismatch));
    REQUIRE(hasCategory(result, TiledDiagnosticCategory::CouplingCentreMismatch));
    return true;
}

} // namespace

int main()
{
    const bool ok = testClosedRingDerivesPitch() && testOpenSector() &&
                    testOverlapRejected() && testNonRadialSeamRejected() &&
                    testGeometryOutsideSectorRejected() && testOverrideLengthMismatch() &&
                    testUnknownOverrideLabel() && testBuildTileProblem() &&
                    testCouplingConsistency();
    if (!ok) {
        std::cerr << "tiled model tests FAILED\n";
        return 1;
    }
    std::cout << "tiled model tests passed\n";
    return 0;
}
