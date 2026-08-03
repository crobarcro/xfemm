#include "CMaterialProp.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void expectInvalid(femm::CMMaterialProp &material, const std::string &messagePart)
{
    try {
        material.GetSlopes();
        assert(false && "invalid material curve was accepted");
    } catch (const std::invalid_argument &error) {
        assert(std::string(error.what()).find(messagePart) != std::string::npos);
    }
}

femm::CMMaterialProp curve(std::initializer_list<double> b,
                           std::initializer_list<CComplex> h)
{
    femm::CMMaterialProp result;
    result.Bdata.assign(b);
    result.Hdata.assign(h);
    result.BHpoints = static_cast<int>(result.Bdata.size());
    return result;
}

} // namespace

int main()
{
    auto onePoint = curve({0.0}, {{0.0, 0.0}});
    expectInvalid(onePoint, "zero points or at least two");

    auto duplicateB = curve({0.0, 0.0}, {{0.0, 0.0}, {100.0, 0.0}});
    expectInvalid(duplicateB, "strictly increasing");

    auto descendingB = curve({1.0, 0.5}, {{0.0, 0.0}, {100.0, 0.0}});
    expectInvalid(descendingB, "strictly increasing");

    auto nanB = curve({0.0, std::numeric_limits<double>::quiet_NaN()},
                      {{0.0, 0.0}, {100.0, 0.0}});
    expectInvalid(nanB, "non-finite");

    auto infiniteH = curve({0.0, 1.0},
                           {{0.0, 0.0}, {std::numeric_limits<double>::infinity(), 0.0}});
    expectInvalid(infiniteH, "non-finite");

    auto dc = curve({0.0, 1.0}, {{0.0, 0.0}, {100.0, 0.0}});
    dc.GetSlopes();
    assert(dc.slope.size() == 2);

    auto ac = curve({0.0, 1.0}, {{0.0, 0.0}, {100.0, 25.0}});
    ac.GetSlopes();
    assert(ac.slope.size() == 2);
    assert(std::isfinite(ac.slope[0].re) && std::isfinite(ac.slope[0].im));

    std::istringstream malformed(
        "<BeginBlock>\n<BlockName> = \"bad\"\n<BHPoints> = 1\n0 0\n<EndBlock>\n");
    std::ostringstream errors;
    femm::CMSolverMaterialProp::fromStream(malformed, errors);
    assert(malformed.fail());
    assert(errors.str().find("zero points or at least two") != std::string::npos);
}
