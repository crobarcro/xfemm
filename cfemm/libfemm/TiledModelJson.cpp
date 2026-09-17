#include "TiledModelJson.h"

#include "Json.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace femm {
namespace tiled {
namespace {

using json::Value;

constexpr int TiledFormatVersion = 1;
constexpr const char *TiledFormatName = "tiled-magnetic";

std::string lower(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

void addError(std::vector<TiledDiagnostic> &errors, TiledDiagnosticCategory category,
              const std::string &tile, const std::string &object, const std::string &message)
{
    errors.push_back({category, tile, object, message});
}

bool isString(const Value &value, const char *key, std::string &out)
{
    const Value *found = value.find(key);
    if (!found || !found->isString())
        return false;
    out = found->stringValue();
    return true;
}

bool isNumber(const Value &value, const char *key, double &out)
{
    const Value *found = value.find(key);
    if (!found || !found->isNumber())
        return false;
    out = found->numberValue();
    return true;
}

bool isInt(const Value &value, const char *key, int &out)
{
    double number = 0.0;
    if (!isNumber(value, key, number))
        return false;
    out = static_cast<int>(std::lround(number));
    return true;
}

bool isSize(const Value &value, const char *key, std::size_t &out)
{
    double number = 0.0;
    if (!isNumber(value, key, number) || number < 0.0)
        return false;
    out = static_cast<std::size_t>(std::lround(number));
    return true;
}

bool isBool(const Value &value, const char *key, bool &out)
{
    const Value *found = value.find(key);
    if (!found || !found->isBool())
        return false;
    out = found->boolValue();
    return true;
}

bool readComplex(const Value &value, CComplex &out)
{
    if (value.isNumber()) {
        out = CComplex(value.numberValue(), 0.0);
        return true;
    }
    if (value.isArray() && value.arrayItems().size() == 2 &&
        value.arrayItems()[0].isNumber() && value.arrayItems()[1].isNumber()) {
        out = CComplex(value.arrayItems()[0].numberValue(), value.arrayItems()[1].numberValue());
        return true;
    }
    return false;
}

Value writeComplex(const CComplex &value)
{
    Value array = Value::array();
    array.push(Value::number(value.re));
    array.push(Value::number(value.im));
    return array;
}

bool parseLengthUnit(const std::string &name, LengthUnit &out)
{
    const std::string key = lower(name);
    if (key == "inches") out = LengthInches;
    else if (key == "millimeters") out = LengthMillimeters;
    else if (key == "centimeters") out = LengthCentimeters;
    else if (key == "meters") out = LengthMeters;
    else if (key == "mils") out = LengthMils;
    else if (key == "micrometers") out = LengthMicrometers;
    else return false;
    return true;
}

const char *lengthUnitName(LengthUnit unit)
{
    switch (unit) {
    case LengthInches: return "inches";
    case LengthMillimeters: return "millimeters";
    case LengthCentimeters: return "centimeters";
    case LengthMeters: return "meters";
    case LengthMils: return "mils";
    case LengthMicrometers: return "micrometers";
    }
    return "meters";
}

bool parseClosure(const std::string &name, Closure &out)
{
    const std::string key = lower(name);
    if (key == "closed") out = Closure::Closed;
    else if (key == "periodic") out = Closure::Periodic;
    else if (key == "antiperiodic") out = Closure::Antiperiodic;
    else return false;
    return true;
}

const char *closureName(Closure closure)
{
    switch (closure) {
    case Closure::Closed: return "closed";
    case Closure::Periodic: return "periodic";
    case Closure::Antiperiodic: return "antiperiodic";
    }
    return "closed";
}

bool parseBoundaryType(const std::string &name, int &format)
{
    const std::string key = lower(name);
    if (key == "dirichlet" || key == "prescribeda") format = 0;
    else if (key == "smallskin" || key == "smallskindepth") format = 1;
    else if (key == "mixed") format = 2;
    else if (key == "periodic") format = 4;
    else if (key == "antiperiodic") format = 5;
    else if (key == "age") format = 6;
    else if (key == "age-antiperiodic" || key == "ageantiperiodic") format = 7;
    else return false;
    return true;
}

const char *boundaryTypeName(int format)
{
    switch (format) {
    case 0: return "dirichlet";
    case 1: return "smallskin";
    case 2: return "mixed";
    case 4: return "periodic";
    case 5: return "antiperiodic";
    case 6: return "age";
    case 7: return "age-antiperiodic";
    default: return "unknown";
    }
}

std::string boundaryName(const TiledModel &model, int index)
{
    if (index < 0 || static_cast<std::size_t>(index) >= model.boundaryProps.size())
        return {};
    return model.boundaryProps[static_cast<std::size_t>(index)]
        ? model.boundaryProps[static_cast<std::size_t>(index)]->BdryName
        : std::string();
}

std::string materialName(const TiledModel &model, int index)
{
    if (index < 0 || static_cast<std::size_t>(index) >= model.materialProps.size())
        return {};
    return model.materialProps[static_cast<std::size_t>(index)]
        ? model.materialProps[static_cast<std::size_t>(index)]->BlockName
        : std::string();
}

std::string circuitName(const TiledModel &model, int index)
{
    if (index < 0 || static_cast<std::size_t>(index) >= model.circuitProps.size())
        return {};
    return model.circuitProps[static_cast<std::size_t>(index)]
        ? model.circuitProps[static_cast<std::size_t>(index)]->CircName
        : std::string();
}

int indexOf(const std::vector<std::string> &names, const std::string &name)
{
    for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i] == name)
            return static_cast<int>(i);
    return -1;
}

class Loader
{
public:
    Loader(TiledModel &model, std::vector<TiledDiagnostic> &errors)
        : m_model(model), m_errors(errors) {}

    bool load(const Value &root)
    {
        if (!root.isObject()) {
            addError(m_errors, TiledDiagnosticCategory::FileFormatError, {}, {},
                     "top-level JSON value must be an object");
            return false;
        }
        if (!loadHeader(root))
            return false;
        loadProblem(root.find("problem"));
        loadPointProps(root.find("pointProps"));
        loadBoundaries(root.find("boundaries"));
        loadMaterials(root.find("materials"));
        loadCircuits(root.find("circuits"));
        loadTiles(root.find("tiles"));
        loadCouplings(root.find("couplings"));
        loadOverrides(root.find("overrides"));
        return m_ok;
    }

private:
    TiledModel &m_model;
    std::vector<TiledDiagnostic> &m_errors;
    bool m_ok = true;

    void fail(TiledDiagnosticCategory category, const std::string &tile,
              const std::string &object, const std::string &message)
    {
        m_ok = false;
        addError(m_errors, category, tile, object, message);
    }

    bool loadHeader(const Value &root)
    {
        const Value *header = root.find("xfemm");
        if (!header || !header->isObject()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, {},
                 "missing xfemm format header");
            return false;
        }
        std::string format;
        if (!isString(*header, "format", format) || format != TiledFormatName) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, {},
                 "unexpected format name; expected tiled-magnetic");
            return false;
        }
        int version = 0;
        if (!isInt(*header, "version", version) || version != TiledFormatVersion) {
            fail(TiledDiagnosticCategory::UnsupportedVersion, {}, {},
                 "unsupported tiled format version");
            return false;
        }
        return true;
    }

    void loadProblem(const Value *problem)
    {
        if (!problem || !problem->isObject()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, {}, "missing problem section");
            return;
        }
        std::string text;
        if (isString(*problem, "type", text)) {
            if (lower(text) == "planar") m_model.problemType = PLANAR;
            else if (lower(text) == "axisymmetric") m_model.problemType = AXISYMMETRIC;
            else fail(TiledDiagnosticCategory::FileFormatError, {}, "problem.type",
                      "unknown problem type: " + text);
        }
        if (isString(*problem, "coords", text)) {
            if (lower(text) == "cartesian") m_model.coords = CART;
            else if (lower(text) == "polar") m_model.coords = POLAR;
            else fail(TiledDiagnosticCategory::FileFormatError, {}, "problem.coords",
                      "unknown coordinate system: " + text);
        }
        if (isString(*problem, "units", text) &&
            !parseLengthUnit(text, m_model.lengthUnits))
            fail(TiledDiagnosticCategory::FileFormatError, {}, "problem.units",
                 "unknown length unit: " + text);
        isNumber(*problem, "depth", m_model.depth);
        isNumber(*problem, "frequency", m_model.frequency);
        isNumber(*problem, "precision", m_model.precision);
        isNumber(*problem, "minAngle", m_model.minAngle);
    }

    void loadPointProps(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "pointProps",
                 "pointProps must be an array");
            return;
        }
        for (const Value &entry : array->arrayItems()) {
            auto prop = std::make_unique<CMPointProp>();
            std::string name;
            isString(entry, "name", name);
            prop->PointName = name;
            const Value *a = entry.find("A");
            if (a && !readComplex(*a, prop->A))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "point property A must be a number or [re,im]");
            const Value *j = entry.find("J");
            if (j && !readComplex(*j, prop->J))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "point property J must be a number or [re,im]");
            m_model.pointProps.push_back(std::move(prop));
        }
    }

    void loadBoundaries(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "boundaries",
                 "boundaries must be an array");
            return;
        }
        for (const Value &entry : array->arrayItems()) {
            auto prop = std::make_unique<CMBoundaryProp>();
            std::string name;
            if (!isString(entry, "name", name)) {
                fail(TiledDiagnosticCategory::FileFormatError, {}, "boundaries",
                     "boundary property is missing a name");
                continue;
            }
            prop->BdryName = name;
            std::string type;
            if (!isString(entry, "type", type) ||
                !parseBoundaryType(type, prop->BdryFormat)) {
                int format = 0;
                if (isInt(entry, "format", format))
                    prop->BdryFormat = format;
                else
                    fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                         "unknown boundary type: " + type);
            }
            isNumber(entry, "A0", prop->A0);
            isNumber(entry, "A1", prop->A1);
            isNumber(entry, "A2", prop->A2);
            isNumber(entry, "phi", prop->phi);
            isNumber(entry, "mu", prop->Mu);
            isNumber(entry, "sigma", prop->Sig);
            const Value *c0 = entry.find("c0");
            if (c0 && !readComplex(*c0, prop->c0))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "boundary c0 must be a number or [re,im]");
            const Value *c1 = entry.find("c1");
            if (c1 && !readComplex(*c1, prop->c1))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "boundary c1 must be a number or [re,im]");
            isNumber(entry, "innerAngle", prop->InnerAngle);
            isNumber(entry, "outerAngle", prop->OuterAngle);
            m_model.boundaryProps.push_back(std::move(prop));
        }
    }

    void loadMaterials(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "materials",
                 "materials must be an array");
            return;
        }
        for (const Value &entry : array->arrayItems()) {
            auto prop = std::make_unique<CMMaterialProp>();
            std::string name;
            if (!isString(entry, "name", name)) {
                fail(TiledDiagnosticCategory::FileFormatError, {}, "materials",
                     "material is missing a name");
                continue;
            }
            prop->BlockName = name;
            isNumber(entry, "muX", prop->mu_x);
            isNumber(entry, "muY", prop->mu_y);
            isNumber(entry, "Hc", prop->H_c);
            isNumber(entry, "sigma", prop->Cduct);
            isNumber(entry, "lamD", prop->Lam_d);
            isNumber(entry, "phiH", prop->Theta_hn);
            isNumber(entry, "phiHx", prop->Theta_hx);
            isNumber(entry, "phiHy", prop->Theta_hy);
            isInt(entry, "lamType", prop->LamType);
            isNumber(entry, "lamFill", prop->LamFill);
            isInt(entry, "strands", prop->NStrands);
            isNumber(entry, "wireD", prop->WireD);
            const Value *j = entry.find("J");
            if (j && !readComplex(*j, prop->J))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "material J must be a number or [re,im]");
            const Value *bh = entry.find("bh");
            if (bh && bh->isArray()) {
                for (const Value &point : bh->arrayItems()) {
                    if (point.isArray() && point.arrayItems().size() == 2 &&
                        point.arrayItems()[0].isNumber() && point.arrayItems()[1].isNumber()) {
                        prop->Bdata.push_back(point.arrayItems()[0].numberValue());
                        prop->Hdata.push_back(CComplex(point.arrayItems()[1].numberValue(), 0.0));
                    }
                }
                prop->BHpoints = static_cast<int>(prop->Bdata.size());
            }
            m_model.materialProps.push_back(std::move(prop));
        }
    }

    void loadCircuits(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "circuits",
                 "circuits must be an array");
            return;
        }
        for (const Value &entry : array->arrayItems()) {
            auto circuit = std::make_unique<CMCircuit>();
            std::string name;
            if (!isString(entry, "name", name)) {
                fail(TiledDiagnosticCategory::FileFormatError, {}, "circuits",
                     "circuit is missing a name");
                continue;
            }
            circuit->CircName = name;
            isInt(entry, "type", circuit->CircType);
            const Value *amps = entry.find("amps");
            if (amps && !readComplex(*amps, circuit->Amps))
                fail(TiledDiagnosticCategory::FileFormatError, {}, name,
                     "circuit amps must be a number or [re,im]");
            m_model.circuitProps.push_back(std::move(circuit));
        }
    }

    int resolveName(const Value &entry, const char *key,
                    const std::vector<std::string> &names,
                    TiledDiagnosticCategory category, const std::string &context,
                    const std::string &kind)
    {
        std::string name;
        if (!isString(entry, key, name) || name.empty())
            return -1;
        const int index = indexOf(names, name);
        if (index < 0)
            fail(category, context, name, "unknown " + kind + ": " + name);
        return index;
    }

    std::vector<std::string> boundaryNames() const
    {
        std::vector<std::string> names;
        for (const auto &prop : m_model.boundaryProps)
            names.push_back(prop ? prop->BdryName : std::string());
        return names;
    }
    std::vector<std::string> materialNames() const
    {
        std::vector<std::string> names;
        for (const auto &prop : m_model.materialProps)
            names.push_back(prop ? prop->BlockName : std::string());
        return names;
    }
    std::vector<std::string> circuitNames() const
    {
        std::vector<std::string> names;
        for (const auto &prop : m_model.circuitProps)
            names.push_back(prop ? prop->CircName : std::string());
        return names;
    }

    void loadTiles(const Value *array)
    {
        if (!array || !array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "tiles",
                 "tiles must be an array");
            return;
        }
        const std::vector<std::string> boundaries = boundaryNames();
        const std::vector<std::string> materials = materialNames();
        const std::vector<std::string> circuits = circuitNames();

        for (const Value &entry : array->arrayItems()) {
            Tile tile;
            isString(entry, "name", tile.name);
            isString(entry, "seamBoundary", tile.seamBoundary);
            isString(entry, "airGapBoundary", tile.airGapBoundary);

            const Value *repeat = entry.find("repeat");
            if (repeat && repeat->isObject()) {
                std::string kind;
                if (isString(*repeat, "kind", kind) && lower(kind) != "rotation")
                    fail(TiledDiagnosticCategory::InvalidRepeatKind, tile.name, {},
                         "only rotation repeat is supported");
                const Value *center = repeat->find("center");
                if (center && center->isArray() && center->arrayItems().size() == 2) {
                    if (center->arrayItems()[0].isNumber())
                        tile.repeat.centerXMetres = center->arrayItems()[0].numberValue();
                    if (center->arrayItems()[1].isNumber())
                        tile.repeat.centerYMetres = center->arrayItems()[1].numberValue();
                }
                isSize(*repeat, "count", tile.repeat.count);
                std::string closure;
                if (isString(*repeat, "closure", closure) &&
                    !parseClosure(closure, tile.repeat.closure))
                    fail(TiledDiagnosticCategory::FileFormatError, tile.name, "repeat.closure",
                         "unknown closure: " + closure);
            }

            const Value *geometry = entry.find("geometry");
            if (!geometry || !geometry->isObject()) {
                fail(TiledDiagnosticCategory::FileFormatError, tile.name, {},
                     "tile is missing a geometry object");
                m_model.tiles.push_back(std::move(tile));
                continue;
            }

            const Value *nodes = geometry->find("nodes");
            if (nodes && nodes->isArray())
                for (const Value &node : nodes->arrayItems()) {
                    TileNode out;
                    isNumber(node, "x", out.x);
                    isNumber(node, "y", out.y);
                    out.boundaryMarker = resolveName(node, "boundary", boundaries,
                                                     TiledDiagnosticCategory::UnknownBoundary,
                                                     tile.name, "boundary");
                    isInt(node, "group", out.group);
                    tile.geometry.nodes.push_back(out);
                }

            const Value *segments = geometry->find("segments");
            if (segments && segments->isArray())
                for (const Value &segment : segments->arrayItems()) {
                    TileSegment out;
                    isInt(segment, "n0", out.n0);
                    isInt(segment, "n1", out.n1);
                    isNumber(segment, "maxSideLength", out.maxSideLength);
                    out.boundary = resolveName(segment, "boundary", boundaries,
                                               TiledDiagnosticCategory::UnknownBoundary,
                                               tile.name, "boundary");
                    isBool(segment, "hidden", out.hidden);
                    isInt(segment, "group", out.group);
                    tile.geometry.segments.push_back(out);
                }

            const Value *arcs = geometry->find("arcs");
            if (arcs && arcs->isArray())
                for (const Value &arc : arcs->arrayItems()) {
                    TileArc out;
                    isInt(arc, "n0", out.n0);
                    isInt(arc, "n1", out.n1);
                    isNumber(arc, "arcLength", out.arcLength);
                    isNumber(arc, "maxSegDegrees", out.maxSegDegrees);
                    out.boundary = resolveName(arc, "boundary", boundaries,
                                               TiledDiagnosticCategory::UnknownBoundary,
                                               tile.name, "boundary");
                    isBool(arc, "hidden", out.hidden);
                    isInt(arc, "group", out.group);
                    tile.geometry.arcs.push_back(out);
                }

            const Value *labels = geometry->find("labels");
            if (labels && labels->isArray())
                for (const Value &label : labels->arrayItems()) {
                    TileLabel out;
                    isString(label, "name", out.name);
                    isNumber(label, "x", out.x);
                    isNumber(label, "y", out.y);
                    out.material = resolveName(label, "material", materials,
                                               TiledDiagnosticCategory::UnknownMaterial,
                                               tile.name, "material");
                    out.circuit = resolveName(label, "circuit", circuits,
                                              TiledDiagnosticCategory::UnknownCircuit,
                                              tile.name, "circuit");
                    isInt(label, "turns", out.turns);
                    isNumber(label, "magDir", out.magDir);
                    isString(label, "magDirFctn", out.magDirFctn);
                    isNumber(label, "maxArea", out.maxArea);
                    isInt(label, "group", out.group);
                    isBool(label, "hole", out.hole);
                    tile.geometry.labels.push_back(std::move(out));
                }

            m_model.tiles.push_back(std::move(tile));
        }
    }

    void loadCouplings(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "couplings",
                 "couplings must be an array");
            return;
        }
        for (const Value &entry : array->arrayItems()) {
            TileCoupling coupling;
            isString(entry, "inner", coupling.innerTile);
            isString(entry, "outer", coupling.outerTile);
            isString(entry, "boundary", coupling.boundary);
            const Value *center = entry.find("center");
            if (center && center->isArray() && center->arrayItems().size() == 2) {
                if (center->arrayItems()[0].isNumber())
                    coupling.centerXMetres = center->arrayItems()[0].numberValue();
                if (center->arrayItems()[1].isNumber())
                    coupling.centerYMetres = center->arrayItems()[1].numberValue();
            }
            isNumber(entry, "innerRadius", coupling.innerRadiusMetres);
            isNumber(entry, "outerRadius", coupling.outerRadiusMetres);
            m_model.couplings.push_back(std::move(coupling));
        }
    }

    void loadOverrides(const Value *array)
    {
        if (!array || array->isNull())
            return;
        if (!array->isArray()) {
            fail(TiledDiagnosticCategory::FileFormatError, {}, "overrides",
                 "overrides must be an array");
            return;
        }
        const std::vector<std::string> circuits = circuitNames();
        for (const Value &entry : array->arrayItems()) {
            TileLabelOverride override;
            isString(entry, "tile", override.tile);
            isString(entry, "label", override.label);

            const Value *magDir = entry.find("magDir");
            if (magDir) {
                if (magDir->isArray()) {
                    for (const Value &value : magDir->arrayItems())
                        if (value.isNumber())
                            override.magDir.push_back(value.numberValue());
                } else if (magDir->isObject()) {
                    std::string mode;
                    if (isString(*magDir, "mode", mode)) {
                        if (lower(mode) == "absolute")
                            override.magDirMode = MagnetisationMode::Absolute;
                        else if (lower(mode) != "delta")
                            fail(TiledDiagnosticCategory::FileFormatError, override.tile,
                                 override.label, "unknown magDir mode: " + mode);
                    }
                    const Value *values = magDir->find("values");
                    if (values && values->isArray())
                        for (const Value &value : values->arrayItems())
                            if (value.isNumber())
                                override.magDir.push_back(value.numberValue());
                }
            }

            const Value *circuit = entry.find("circuit");
            if (circuit && circuit->isArray())
                for (const Value &value : circuit->arrayItems()) {
                    if (!value.isString()) {
                        override.circuit.push_back(-1);
                        continue;
                    }
                    const std::string &name = value.stringValue();
                    override.circuit.push_back(name.empty() ? -1 : indexOf(circuits, name));
                }

            const Value *turnScale = entry.find("turnScale");
            if (turnScale && turnScale->isArray())
                for (const Value &value : turnScale->arrayItems())
                    if (value.isNumber())
                        override.turnScale.push_back(value.numberValue());

            m_model.overrides.push_back(std::move(override));
        }
    }
};

} // namespace

bool loadTiledModelJson(const std::string &text, TiledModel &out,
                        std::vector<TiledDiagnostic> &errors)
{
    Value root;
    std::string error;
    if (!json::Value::parse(text, root, error)) {
        addError(errors, TiledDiagnosticCategory::FileFormatError, {}, {}, error);
        return false;
    }
    out = TiledModel();
    Loader loader(out, errors);
    return loader.load(root);
}

namespace {

Value tileToJson(const TiledModel &model, const Tile &tile)
{
    Value entry = Value::object();
    entry.set("name", Value::string(tile.name));

    Value geometry = Value::object();
    Value nodes = Value::array();
    for (const TileNode &node : tile.geometry.nodes) {
        Value value = Value::object();
        value.set("x", Value::number(node.x));
        value.set("y", Value::number(node.y));
        if (node.boundaryMarker >= 0)
            value.set("boundary", Value::string(boundaryName(model, node.boundaryMarker)));
        if (node.group != 0)
            value.set("group", Value::number(node.group));
        nodes.push(std::move(value));
    }
    geometry.set("nodes", std::move(nodes));

    Value segments = Value::array();
    for (const TileSegment &segment : tile.geometry.segments) {
        Value value = Value::object();
        value.set("n0", Value::number(segment.n0));
        value.set("n1", Value::number(segment.n1));
        value.set("maxSideLength", Value::number(segment.maxSideLength));
        if (segment.boundary >= 0)
            value.set("boundary", Value::string(boundaryName(model, segment.boundary)));
        if (segment.hidden)
            value.set("hidden", Value::boolean(true));
        if (segment.group != 0)
            value.set("group", Value::number(segment.group));
        segments.push(std::move(value));
    }
    geometry.set("segments", std::move(segments));

    Value arcs = Value::array();
    for (const TileArc &arc : tile.geometry.arcs) {
        Value value = Value::object();
        value.set("n0", Value::number(arc.n0));
        value.set("n1", Value::number(arc.n1));
        value.set("arcLength", Value::number(arc.arcLength));
        value.set("maxSegDegrees", Value::number(arc.maxSegDegrees));
        if (arc.boundary >= 0)
            value.set("boundary", Value::string(boundaryName(model, arc.boundary)));
        if (arc.hidden)
            value.set("hidden", Value::boolean(true));
        if (arc.group != 0)
            value.set("group", Value::number(arc.group));
        arcs.push(std::move(value));
    }
    geometry.set("arcs", std::move(arcs));

    Value labels = Value::array();
    for (const TileLabel &label : tile.geometry.labels) {
        Value value = Value::object();
        value.set("name", Value::string(label.name));
        value.set("x", Value::number(label.x));
        value.set("y", Value::number(label.y));
        if (label.material >= 0)
            value.set("material", Value::string(materialName(model, label.material)));
        if (label.circuit >= 0)
            value.set("circuit", Value::string(circuitName(model, label.circuit)));
        value.set("turns", Value::number(label.turns));
        value.set("magDir", Value::number(label.magDir));
        if (!label.magDirFctn.empty())
            value.set("magDirFctn", Value::string(label.magDirFctn));
        value.set("maxArea", Value::number(label.maxArea));
        if (label.group != 0)
            value.set("group", Value::number(label.group));
        if (label.hole)
            value.set("hole", Value::boolean(true));
        labels.push(std::move(value));
    }
    geometry.set("labels", std::move(labels));
    entry.set("geometry", std::move(geometry));

    Value repeat = Value::object();
    repeat.set("kind", Value::string("rotation"));
    Value center = Value::array();
    center.push(Value::number(tile.repeat.centerXMetres));
    center.push(Value::number(tile.repeat.centerYMetres));
    repeat.set("center", std::move(center));
    repeat.set("count", Value::number(static_cast<double>(tile.repeat.count)));
    repeat.set("closure", Value::string(closureName(tile.repeat.closure)));
    entry.set("repeat", std::move(repeat));

    entry.set("seamBoundary", Value::string(tile.seamBoundary));
    if (!tile.airGapBoundary.empty())
        entry.set("airGapBoundary", Value::string(tile.airGapBoundary));
    return entry;
}

} // namespace

std::string saveTiledModelJson(const TiledModel &model)
{
    Value root = Value::object();
    Value header = Value::object();
    header.set("format", Value::string(TiledFormatName));
    header.set("version", Value::number(TiledFormatVersion));
    root.set("xfemm", std::move(header));

    Value problem = Value::object();
    problem.set("type", Value::string(model.problemType == AXISYMMETRIC ? "axisymmetric"
                                                                        : "planar"));
    problem.set("coords", Value::string(model.coords == POLAR ? "polar" : "cartesian"));
    problem.set("units", Value::string(lengthUnitName(model.lengthUnits)));
    problem.set("depth", Value::number(model.depth));
    problem.set("frequency", Value::number(model.frequency));
    problem.set("precision", Value::number(model.precision));
    problem.set("minAngle", Value::number(model.minAngle));
    root.set("problem", std::move(problem));

    Value pointProps = Value::array();
    for (const auto &prop : model.pointProps) {
        if (!prop)
            continue;
        Value entry = Value::object();
        entry.set("name", Value::string(prop->PointName));
        entry.set("A", writeComplex(prop->A));
        entry.set("J", writeComplex(prop->J));
        pointProps.push(std::move(entry));
    }
    root.set("pointProps", std::move(pointProps));

    Value boundaries = Value::array();
    for (const auto &prop : model.boundaryProps) {
        if (!prop)
            continue;
        Value entry = Value::object();
        entry.set("name", Value::string(prop->BdryName));
        entry.set("type", Value::string(boundaryTypeName(prop->BdryFormat)));
        entry.set("format", Value::number(prop->BdryFormat));
        entry.set("A0", Value::number(prop->A0));
        entry.set("A1", Value::number(prop->A1));
        entry.set("A2", Value::number(prop->A2));
        entry.set("phi", Value::number(prop->phi));
        entry.set("mu", Value::number(prop->Mu));
        entry.set("sigma", Value::number(prop->Sig));
        entry.set("c0", writeComplex(prop->c0));
        entry.set("c1", writeComplex(prop->c1));
        entry.set("innerAngle", Value::number(prop->InnerAngle));
        entry.set("outerAngle", Value::number(prop->OuterAngle));
        boundaries.push(std::move(entry));
    }
    root.set("boundaries", std::move(boundaries));

    Value materials = Value::array();
    for (const auto &prop : model.materialProps) {
        if (!prop)
            continue;
        Value entry = Value::object();
        entry.set("name", Value::string(prop->BlockName));
        entry.set("muX", Value::number(prop->mu_x));
        entry.set("muY", Value::number(prop->mu_y));
        entry.set("Hc", Value::number(prop->H_c));
        entry.set("J", writeComplex(prop->J));
        entry.set("sigma", Value::number(prop->Cduct));
        entry.set("lamD", Value::number(prop->Lam_d));
        entry.set("phiH", Value::number(prop->Theta_hn));
        entry.set("phiHx", Value::number(prop->Theta_hx));
        entry.set("phiHy", Value::number(prop->Theta_hy));
        entry.set("lamType", Value::number(prop->LamType));
        entry.set("lamFill", Value::number(prop->LamFill));
        entry.set("strands", Value::number(prop->NStrands));
        entry.set("wireD", Value::number(prop->WireD));
        Value bh = Value::array();
        for (std::size_t i = 0; i < prop->Bdata.size() && i < prop->Hdata.size(); ++i) {
            Value point = Value::array();
            point.push(Value::number(prop->Bdata[i]));
            point.push(Value::number(prop->Hdata[i].re));
            bh.push(std::move(point));
        }
        entry.set("bh", std::move(bh));
        materials.push(std::move(entry));
    }
    root.set("materials", std::move(materials));

    Value circuits = Value::array();
    for (const auto &prop : model.circuitProps) {
        if (!prop)
            continue;
        Value entry = Value::object();
        entry.set("name", Value::string(prop->CircName));
        entry.set("amps", writeComplex(prop->Amps));
        entry.set("type", Value::number(prop->CircType));
        circuits.push(std::move(entry));
    }
    root.set("circuits", std::move(circuits));

    Value tiles = Value::array();
    for (const Tile &tile : model.tiles)
        tiles.push(tileToJson(model, tile));
    root.set("tiles", std::move(tiles));

    Value couplings = Value::array();
    for (const TileCoupling &coupling : model.couplings) {
        Value entry = Value::object();
        entry.set("inner", Value::string(coupling.innerTile));
        entry.set("outer", Value::string(coupling.outerTile));
        entry.set("boundary", Value::string(coupling.boundary));
        Value center = Value::array();
        center.push(Value::number(coupling.centerXMetres));
        center.push(Value::number(coupling.centerYMetres));
        entry.set("center", std::move(center));
        entry.set("innerRadius", Value::number(coupling.innerRadiusMetres));
        entry.set("outerRadius", Value::number(coupling.outerRadiusMetres));
        couplings.push(std::move(entry));
    }
    root.set("couplings", std::move(couplings));

    Value overrides = Value::array();
    for (const TileLabelOverride &override : model.overrides) {
        Value entry = Value::object();
        entry.set("tile", Value::string(override.tile));
        entry.set("label", Value::string(override.label));
        if (!override.magDir.empty()) {
            Value magDir = Value::object();
            magDir.set("mode", Value::string(
                                   override.magDirMode == MagnetisationMode::Absolute
                                       ? "absolute"
                                       : "delta"));
            Value values = Value::array();
            for (double value : override.magDir)
                values.push(Value::number(value));
            magDir.set("values", std::move(values));
            entry.set("magDir", std::move(magDir));
        }
        if (!override.circuit.empty()) {
            Value circuit = Value::array();
            for (int index : override.circuit)
                circuit.push(Value::string(circuitName(model, index)));
            entry.set("circuit", std::move(circuit));
        }
        if (!override.turnScale.empty()) {
            Value turnScale = Value::array();
            for (double value : override.turnScale)
                turnScale.push(Value::number(value));
            entry.set("turnScale", std::move(turnScale));
        }
        overrides.push(std::move(entry));
    }
    root.set("overrides", std::move(overrides));

    return root.dump(2) + "\n";
}

} // namespace tiled
} // namespace femm
