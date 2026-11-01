#include "FemmProblem.h"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

enum class CopyTransform
{
    Mirror,
    Rotate,
    Translate
};

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

const char *transformName(CopyTransform transform)
{
    switch (transform)
    {
    case CopyTransform::Mirror: return "mirror";
    case CopyTransform::Rotate: return "rotate";
    case CopyTransform::Translate: return "translate";
    }
    return "unknown";
}

void addSelectedObject(femm::FemmProblem &problem, femm::EditMode mode)
{
    if (mode == femm::EditMode::EditNodes || mode == femm::EditMode::EditGroup)
    {
        auto node = std::make_unique<femm::CNode>(2., 1.);
        node->IsSelected = true;
        problem.nodelist.push_back(std::move(node));
        return;
    }

    if (mode == femm::EditMode::EditLines)
    {
        problem.nodelist.push_back(std::make_unique<femm::CNode>(2., 0.));
        problem.nodelist.push_back(std::make_unique<femm::CNode>(3., 0.));
        auto line = std::make_unique<femm::CSegment>();
        line->n0 = 0;
        line->n1 = 1;
        line->IsSelected = true;
        problem.linelist.push_back(std::move(line));
        return;
    }

    if (mode == femm::EditMode::EditArcs)
    {
        problem.nodelist.push_back(std::make_unique<femm::CNode>(2., 0.));
        problem.nodelist.push_back(std::make_unique<femm::CNode>(3., 0.));
        auto arc = std::make_unique<femm::CArcSegment>();
        arc->n0 = 0;
        arc->n1 = 1;
        arc->ArcLength = 30.;
        arc->MaxSideLength = 5.;
        arc->IsSelected = true;
        problem.arclist.push_back(std::move(arc));
        return;
    }

    if (mode == femm::EditMode::EditLabels)
    {
        auto label = std::make_unique<femm::CMBlockLabel>();
        label->x = 2.;
        label->y = 1.;
        label->IsSelected = true;
        problem.labellist.push_back(std::move(label));
        return;
    }

    throw std::runtime_error("invalid edit mode in copy-transform test");
}

void applyCopy(femm::FemmProblem &problem, CopyTransform transform, femm::EditMode mode, int copies = 1)
{
    switch (transform)
    {
    case CopyTransform::Mirror:
        require(copies == 1, "mirror copy only supports one copy");
        problem.mirrorCopy(10., -1., 10., 1., mode);
        break;
    case CopyTransform::Rotate:
        problem.rotateCopy(CComplex(100., 100.), 30., copies, mode);
        break;
    case CopyTransform::Translate:
        problem.translateCopy(10., 3., copies, mode);
        break;
    }
}

std::size_t selectedContainerSize(const femm::FemmProblem &problem, femm::EditMode mode)
{
    switch (mode)
    {
    case femm::EditMode::EditNodes:
    case femm::EditMode::EditGroup:
        return problem.nodelist.size();
    case femm::EditMode::EditLines:
        return problem.linelist.size();
    case femm::EditMode::EditArcs:
        return problem.arclist.size();
    case femm::EditMode::EditLabels:
        return problem.labellist.size();
    case femm::EditMode::Invalid:
        break;
    }
    throw std::runtime_error("invalid edit mode in copy-transform test");
}

void testSingleCopy(CopyTransform transform, femm::EditMode mode)
{
    femm::FemmProblem problem(femm::FileType::MagneticsFile);
    addSelectedObject(problem, mode);
    applyCopy(problem, transform, mode);

    require(selectedContainerSize(problem, mode) == 2,
            std::string(transformName(transform)) + " copy must create exactly one object copy");
}

void testMultipleNodeCopies(CopyTransform transform)
{
    femm::FemmProblem problem(femm::FileType::MagneticsFile);
    addSelectedObject(problem, femm::EditMode::EditNodes);
    applyCopy(problem, transform, femm::EditMode::EditNodes, 2);

    require(problem.nodelist.size() == 3,
            std::string(transformName(transform)) + " copy must not recursively copy newly appended nodes");
}

} // namespace

int main()
{
    const CopyTransform transforms[] = {
        CopyTransform::Mirror,
        CopyTransform::Rotate,
        CopyTransform::Translate
    };
    const femm::EditMode modes[] = {
        femm::EditMode::EditNodes,
        femm::EditMode::EditLines,
        femm::EditMode::EditArcs,
        femm::EditMode::EditLabels,
        femm::EditMode::EditGroup
    };

    for (CopyTransform transform : transforms)
        for (femm::EditMode mode : modes)
            testSingleCopy(transform, mode);

    testMultipleNodeCopies(CopyTransform::Rotate);
    testMultipleNodeCopies(CopyTransform::Translate);

    return 0;
}
