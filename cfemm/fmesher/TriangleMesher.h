#ifndef TRIANGLEMESHER_H
#define TRIANGLEMESHER_H

#include <memory>
#include <string>
#include <vector>

namespace femm { class CNode; class CSegment; class FemmProblem; }

namespace fmesher {

enum class PointMarkerInfo { None, FromProblem };
enum class SegmentMarkerInfo { FromCnt, FromProblem };

struct RawTriangulation {
    struct Point { double x, y; int marker; };
    struct Triangle { std::vector<int> corners; std::vector<double> attributes; };
    struct Edge { int first, second, marker; };

    std::vector<Point> points;
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
    int cornersPerTriangle = 0;
    int attributesPerTriangle = 0;
};

/** C++ boundary around Triangle. Triangle's C types are intentionally hidden. */
class TriangleMesher {
    using nodelist_t = std::vector<std::unique_ptr<femm::CNode>>;
    using linelist_t = std::vector<std::unique_ptr<femm::CSegment>>;
public:
    TriangleMesher();
    ~TriangleMesher();
    TriangleMesher(const TriangleMesher &) = delete;
    TriangleMesher &operator=(const TriangleMesher &) = delete;

    bool initPointsWithMarkers(const nodelist_t &, const femm::FemmProblem &, PointMarkerInfo);
    bool initSegmentsWithMarkers(const linelist_t &, const femm::FemmProblem &, SegmentMarkerInfo);
    bool initHolesAndRegions(const femm::FemmProblem &, bool forceMaxMeshArea, double defaultMeshSize);
    int triangulate(bool verbose);
    std::string triangulateParams(bool verbose = false) const;
    bool writePolyFile(std::string filename, std::string comment) const;
    RawTriangulation rawTriangulation() const;

    int (*WarnMessage)(const char *, ...);
    int (*TriMessage)(const char *, ...);
    void setMinAngle(double);
    void suppressExteriorSteinerPoints();
    void suppressUnusedVertices();

private:
    class Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace fmesher
#endif
