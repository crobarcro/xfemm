#include "mesh/InstancedMesh.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

using namespace femm::mesh;

#define REQUIRE(condition)                                                            \
    do {                                                                              \
        if (!(condition)) {                                                           \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << #condition << '\n';   \
            return false;                                                             \
        }                                                                             \
    } while (false)

MeshTemplate makeTriangleTemplate()
{
    MeshTemplate meshTemplate;
    meshTemplate.localMesh.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    meshTemplate.localMesh.elements.push_back({{{0, 1, 2}}, 1});
    meshTemplate.localMesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    meshTemplate.seams = {{"edge", {0, 1}}};
    return meshTemplate;
}

InstancedMesh makeInstanced()
{
    InstancedMesh instanced;
    instanced.templates.push_back(makeTriangleTemplate());
    MeshInstance first;
    MeshInstance second;
    second.transform.translationXMetres = 1.0;
    instanced.instances.push_back(first);
    instanced.instances.push_back(second);
    return instanced;
}

bool testDeterminismAndSeparation()
{
    InstancedMesh first = makeInstanced();
    InstancedMesh second = makeInstanced();
    REQUIRE(templateTopologyIdentity(first) == templateTopologyIdentity(second));
    REQUIRE(instanceLayoutIdentity(first) == instanceLayoutIdentity(second));
    REQUIRE(materializedTopologyIdentity(first.materialize().mesh) ==
            materializedTopologyIdentity(second.materialize().mesh));

    // A transform change alters layout and materialised topology only.
    InstancedMesh moved = first;
    moved.instances[1].transform.translationXMetres = 5.0;
    REQUIRE(templateTopologyIdentity(moved) == templateTopologyIdentity(first));
    REQUIRE(instanceLayoutIdentity(moved) != instanceLayoutIdentity(first));
    REQUIRE(materializedTopologyIdentity(moved.materialize().mesh) !=
            materializedTopologyIdentity(first.materialize().mesh));

    // A seam change alters template and materialised topology only.
    InstancedMesh reseamed = first;
    reseamed.templates[0].seams[0].orderedNodes = {0, 2};
    REQUIRE(templateTopologyIdentity(reseamed) != templateTopologyIdentity(first));
    REQUIRE(instanceLayoutIdentity(reseamed) == instanceLayoutIdentity(first));

    // A region attribute change alters template and materialised topology.
    InstancedMesh reattributed = first;
    reattributed.templates[0].localMesh.elements[0].regionAttribute = 7;
    REQUIRE(templateTopologyIdentity(reattributed) != templateTopologyIdentity(first));
    REQUIRE(instanceLayoutIdentity(reattributed) == instanceLayoutIdentity(first));
    REQUIRE(materializedTopologyIdentity(reattributed.materialize().mesh) !=
            materializedTopologyIdentity(first.materialize().mesh));

    // Repeating a materialisation is byte-stable.
    REQUIRE(materializedTopologyIdentity(first.materialize().mesh) ==
            materializedTopologyIdentity(first.materialize().mesh));
    return true;
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, bool (*)()>> tests = {
        {"determinism-and-separation", testDeterminismAndSeparation},
    };
    bool success = true;
    for (const auto &test : tests) {
        if (!test.second()) {
            std::cerr << "FAILED: " << test.first << '\n';
            success = false;
        }
    }
    if (success)
        std::cout << "instanced mesh hash tests passed\n";
    return success ? 0 : 1;
}
