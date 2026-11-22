#include "TiledModelJson.h"

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

TiledModel makeModel()
{
    TiledModel model;
    model.depth = 0.0889;
    model.problemType = PLANAR;
    model.lengthUnits = LengthMeters;

    auto seam = std::make_unique<CMBoundaryProp>();
    seam->BdryName = "slot-seam";
    seam->BdryFormat = 4;
    model.boundaryProps.push_back(std::move(seam));

    auto outer = std::make_unique<CMBoundaryProp>();
    outer->BdryName = "outer";
    outer->BdryFormat = 0;
    model.boundaryProps.push_back(std::move(outer));

    auto gap = std::make_unique<CMBoundaryProp>();
    gap->BdryName = "sliding-gap";
    gap->BdryFormat = 6;
    model.boundaryProps.push_back(std::move(gap));

    auto steel = std::make_unique<CMMaterialProp>();
    steel->BlockName = "steel";
    steel->mu_x = steel->mu_y = 1000;
    steel->Bdata = {0.0, 1.0};
    steel->Hdata = {CComplex(0, 0), CComplex(100, 0)};
    steel->BHpoints = 2;
    model.materialProps.push_back(std::move(steel));

    auto copper = std::make_unique<CMMaterialProp>();
    copper->BlockName = "copper";
    copper->mu_x = copper->mu_y = 1;
    model.materialProps.push_back(std::move(copper));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    model.circuitProps.push_back(std::move(circuit));

    const double sector = 10.0 * Pi / 180.0;
    const double ri = 0.05;
    const double ro = 0.06;
    Tile tile;
    tile.name = "slot";
    tile.repeat.count = 6;
    tile.repeat.closure = Closure::Periodic;
    tile.seamBoundary = "slot-seam";
    tile.airGapBoundary = "sliding-gap";
    tile.geometry.nodes = {
        {ri, 0.0, 0, 0},
        {ro, 0.0, 0, 0},
        {ri * std::cos(sector), ri * std::sin(sector), 0, 0},
        {ro * std::cos(sector), ro * std::sin(sector), 0, 0},
    };
    tile.geometry.segments = {{0, 1, -1, 0, false, 0}, {2, 3, -1, 0, false, 0}};
    tile.geometry.arcs = {
        {1, 3, 10.0, 2.0, 1, false, 0},
        {2, 0, 10.0, 2.0, 1, false, 0},
    };
    tile.geometry.labels = {
        {"coil-a", 0.055, 0.004, 1, 0, 200, 0.0, "", 0.0, 0, false},
        {"coil-b", 0.055, 0.006, 1, 0, -200, 0.0, "", 0.0, 0, false},
    };
    model.tiles.push_back(tile);

    TileLabelOverride override;
    override.tile = "slot";
    override.label = "coil-a";
    override.magDir = {0.0, 180.0, 0.0, 180.0, 0.0, 180.0};
    override.magDirMode = MagnetisationMode::Absolute;
    override.circuit = {0, -1, 0, -1, 0, -1};
    override.turnScale = {1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    model.overrides.push_back(override);
    return model;
}

bool testRoundTrip()
{
    const TiledModel model = makeModel();
    const std::string text = saveTiledModelJson(model);
    REQUIRE(!text.empty());

    TiledModel loaded;
    std::vector<TiledDiagnostic> errors;
    REQUIRE(loadTiledModelJson(text, loaded, errors));
    REQUIRE(errors.empty());

    REQUIRE(loaded.tiles.size() == 1);
    REQUIRE(loaded.tiles[0].name == "slot");
    REQUIRE(loaded.tiles[0].repeat.count == 6);
    REQUIRE(loaded.tiles[0].repeat.closure == Closure::Periodic);
    REQUIRE(loaded.tiles[0].seamBoundary == "slot-seam");
    REQUIRE(loaded.tiles[0].airGapBoundary == "sliding-gap");
    REQUIRE(loaded.tiles[0].geometry.nodes.size() == 4);
    REQUIRE(loaded.tiles[0].geometry.segments.size() == 2);
    REQUIRE(loaded.tiles[0].geometry.arcs.size() == 2);
    REQUIRE(loaded.tiles[0].geometry.labels.size() == 2);
    REQUIRE(loaded.tiles[0].geometry.labels[0].name == "coil-a");
    REQUIRE(loaded.tiles[0].geometry.labels[0].material == 1);
    REQUIRE(loaded.tiles[0].geometry.labels[0].circuit == 0);
    REQUIRE(loaded.tiles[0].geometry.segments[0].boundary == 0);

    REQUIRE(loaded.materialProps.size() == 2);
    REQUIRE(loaded.materialProps[0]->BlockName == "steel");
    REQUIRE(loaded.materialProps[0]->BHpoints == 2);
    REQUIRE(loaded.boundaryProps.size() == 3);
    REQUIRE(loaded.circuitProps.size() == 1);

    REQUIRE(loaded.overrides.size() == 1);
    REQUIRE(loaded.overrides[0].label == "coil-a");
    REQUIRE(loaded.overrides[0].magDirMode == MagnetisationMode::Absolute);
    REQUIRE(loaded.overrides[0].magDir.size() == 6);
    REQUIRE(loaded.overrides[0].magDir[1] == 180.0);
    REQUIRE(loaded.overrides[0].circuit[1] == -1);
    REQUIRE(loaded.overrides[0].circuit[2] == 0);
    REQUIRE(loaded.overrides[0].turnScale[2] == -1.0);

    // The reloaded model must still satisfy the tiling semantics.
    const TiledValidationResult validation = validateTiledModel(loaded);
    REQUIRE(validation.succeeded());
    REQUIRE(std::fabs(validation.tiles[0].pitchDegrees - 10.0) < 1e-9);
    REQUIRE(std::fabs(validation.tiles[0].totalAngleDegrees - 60.0) < 1e-6);
    return true;
}

bool testMalformedJson()
{
    TiledModel model;
    std::vector<TiledDiagnostic> errors;
    REQUIRE(!loadTiledModelJson("{ this is not json }", model, errors));
    REQUIRE(!errors.empty());
    return true;
}

bool testWrongVersion()
{
    TiledModel model;
    std::vector<TiledDiagnostic> errors;
    const std::string text =
        "{\"xfemm\":{\"format\":\"tiled-magnetic\",\"version\":99}}";
    REQUIRE(!loadTiledModelJson(text, model, errors));
    bool found = false;
    for (const auto &error : errors)
        found = found || error.category == TiledDiagnosticCategory::UnsupportedVersion;
    REQUIRE(found);
    return true;
}

bool testUnknownReference()
{
    const TiledModel model = makeModel();
    std::string text = saveTiledModelJson(model);
    const std::string needle = "\"material\": \"copper\"";
    const std::size_t position = text.find(needle);
    REQUIRE(position != std::string::npos);
    text.replace(position, needle.size(), "\"material\": \"unobtainium\"");
    TiledModel loaded;
    std::vector<TiledDiagnostic> errors;
    REQUIRE(!loadTiledModelJson(text, loaded, errors));
    return true;
}

} // namespace

int main()
{
    const bool ok = testRoundTrip() && testMalformedJson() && testWrongVersion() &&
                    testUnknownReference();
    if (!ok) {
        std::cerr << "tiled model JSON tests FAILED\n";
        return 1;
    }
    std::cout << "tiled model JSON tests passed\n";
    return 0;
}
