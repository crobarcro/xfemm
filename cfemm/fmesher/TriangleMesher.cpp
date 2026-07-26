/* Triangle backend implementation. */
#include "TriangleMesher.h"
#include "fparse.h"
#include "FemmProblem.h"
#include "CNode.h"
#include "CSegment.h"
#include "triangle.h"
#ifndef XFEMM_BUILTIN_TRIANGLE
#include "triangle_api.h"
#endif
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#ifndef REAL
#define REAL double
#endif
using namespace std;
using namespace femm;
namespace fmesher {
namespace {
#ifdef XFEMM_BUILTIN_TRIANGLE
void initialize(triangulateio &io)
#else
void initialize(triangleio &io)
#endif
{
    io = {};
}
}
class TriangleMesher::Impl {
public:
#ifdef XFEMM_BUILTIN_TRIANGLE
    triangulateio in;
    triangulateio out;
#else
    triangleio in;
    context *ctx;
#endif
    double m_minAngle = 0.;
    bool m_suppressExteriorSteinerPoints = false;
    bool m_suppressUnusedVertices = false;
    Impl() {
        initialize(in);
#ifdef XFEMM_BUILTIN_TRIANGLE
        initialize(out);
#else
        ctx = triangle_context_create();
#endif
    }
    ~Impl() {
        free(in.pointlist); free(in.pointattributelist); free(in.pointmarkerlist);
        free(in.regionlist); free(in.segmentlist); free(in.segmentmarkerlist); free(in.holelist);
#ifdef XFEMM_BUILTIN_TRIANGLE
        free(out.pointlist); free(out.pointattributelist); free(out.pointmarkerlist);
        free(out.trianglelist); free(out.triangleattributelist); free(out.trianglearealist);
        free(out.neighborlist); free(out.segmentlist); free(out.segmentmarkerlist);
        free(out.edgelist); free(out.edgemarkerlist);
#else
        triangle_context_destroy(ctx);
#endif
    }
};
TriangleMesher::TriangleMesher() : WarnMessage(&PrintWarningMsg), TriMessage(nullptr), impl(new Impl) {}
TriangleMesher::~TriangleMesher() = default;
bool TriangleMesher::initPointsWithMarkers(const TriangleMesher::nodelist_t &nodelst, const FemmProblem &problem, PointMarkerInfo info)
{
    // calling this method on an already initialized object would leak memory
    if (impl->in.numberofpoints!=0)
    {
        WarnMessage("initPointsWithMarkers called twice!\n");
        return false;
    }

    impl->in.numberofpoints = nodelst.size();

    impl->in.pointlist = (REAL *) malloc(impl->in.numberofpoints * 2 * sizeof(REAL));
    if (!impl->in.pointlist) {
        WarnMessage("Point list for triangulation is null!\n");
        return false;
    }

    for(int i=0; i < impl->in.numberofpoints; i++)
    {
        impl->in.pointlist[2*i] = nodelst[i]->x;
        impl->in.pointlist[2*i+1] = nodelst[i]->y;
    }

    // Initialise the pointmarkerlist
    impl->in.pointmarkerlist = (int *) malloc(impl->in.numberofpoints * sizeof(int));
    if (!impl->in.pointmarkerlist) {
        WarnMessage("Point marker list for triangulation is null!\n");
        return false;
    }

    // write impl->out node marker list
    for(int i=0; i<impl->in.numberofpoints; i++)
    {
        int t=0;
        if (info==PointMarkerInfo::FromProblem)
        {
            for(int j=0; j<(int)problem.nodeproplist.size(); j++)
                if(problem.nodeproplist[j]->PointName==nodelst[i]->BoundaryMarkerName)
                    t = j + 2;

            if (problem.filetype != femm::FileType::MagneticsFile)
            {
                // include conductor number;
                for(int j = 0; j < (int)problem.circproplist.size(); j++)
                {
                    // add the conductor number using a mask
                    if(problem.circproplist[j]->CircName == nodelst[i]->InConductorName)
                        t += ((j+1) * 0x10000);
                }
            }
        }

        impl->in.pointmarkerlist[i] = t;
    }
    return true;
}

bool TriangleMesher::initSegmentsWithMarkers(const TriangleMesher::linelist_t &linelst, const FemmProblem &problem, SegmentMarkerInfo info)
{
    // calling this method on an already initialized object would leak memory
    if (impl->in.numberofsegments!=0)
    {
        WarnMessage("initSegmentsWithMarkers called twice!\n");
        return false;
    }

    impl->in.numberofsegments = linelst.size();

    // Initialise the segmentlist
    impl->in.segmentlist = (int *) malloc(2 * impl->in.numberofsegments * sizeof(int));
    if (!impl->in.segmentlist) {
        WarnMessage("Segment list for triangulation is null!\n");
        return false;
    }
    // Initialise the segmentmarkerlist
    impl->in.segmentmarkerlist = (int *) malloc(impl->in.numberofsegments * sizeof(int));
    if (!impl->in.segmentmarkerlist) {
        WarnMessage("Segment marker list for triangulation is null!\n");
        return false;
    }

    // build the segmentlist
    for(int i=0; i<impl->in.numberofsegments; i++)
    {
        impl->in.segmentlist[2*i] = linelst[i]->n0;
        impl->in.segmentlist[2*i+1] = linelst[i]->n1;
    }

    // now build the segment marker list
    // construct the segment list
    for(int i=0; i<impl->in.numberofsegments; i++)
    {
        int t=0;
        if (info==SegmentMarkerInfo::FromProblem)
        {
            for(int j=0; j <(int)problem.lineproplist.size(); j++)
            {
                if (problem.lineproplist[j]->BdryName == linelst[i]->BoundaryMarkerName)
                {
                    t = -(j+2);
                }
            }

            if (problem.filetype != femm::FileType::MagneticsFile)
            {
                // include conductor number;
                for (int j=0; j <(int)problem.circproplist.size(); j++)
                {
                    if (problem.circproplist[j]->CircName == linelst[i]->InConductorName)
                    {
                        t -= ((j+1) * 0x10000);
                    }
                }
            }
        } else {
            t = -(linelst[i]->cnt+2);
        }
        impl->in.segmentmarkerlist[i] = t;
    }
    return true;
}

bool TriangleMesher::initHolesAndRegions(const FemmProblem &problem, bool forceMaxMeshArea, double defaultMeshSize)
{
    // calling this method on an already initialized object would leak memory
    if (impl->in.numberofholes!=0)
    {
        WarnMessage("initHolesAndRegions called twice!\n");
        return false;
    }

    impl->in.numberofholes = problem.countHoles();
    if(impl->in.numberofholes > 0)
    {
        impl->in.holelist = (REAL *) malloc(impl->in.numberofholes * 2 * sizeof(REAL));
        if (!impl->in.holelist) {
            WarnMessage("Hole list for triangulation is null!\n");
            return false;
        }

        // Construct the holes array
        int k=0;
        for(const auto &label: problem.labellist)
        {
            // we search through the block list looking for blocks that have
            // the tag <No Mesh>
            if(label->isHole())
            {
#ifdef DEBUG
                {
                    char buf[1028];
                    SNPRINTF (buf, sizeof(buf), "Adding hole (at (%g,%g)) to triangle input hole list\n",
                              label->x,label->y);
                    WarnMessage(buf);
                }
#endif // DEBUG
                impl->in.holelist[k++] = label->x;
                impl->in.holelist[k++] = label->y;
            }
        }
    }

    impl->in.numberofregions = problem.labellist.size() - impl->in.numberofholes;
    impl->in.regionlist = (REAL *) malloc(impl->in.numberofregions * 4 * sizeof(REAL));
    if (!impl->in.regionlist) {
        WarnMessage("Region list for triangulation is null!\n");
        return false;
    }

    int j=0;
    int k=0;
    for(const auto & label: problem.labellist)
    {
        if(!label->isHole())
        {
            impl->in.regionlist[j] = label->x;
            impl->in.regionlist[j+1] = label->y;
            impl->in.regionlist[j+2] = k + 1; // Regional attribute (for whole mesh).
#ifdef DEBUG
            {
                char buf[1028];
                SNPRINTF (buf, sizeof(buf), "Adding region (at (%g,%g)) with attribute value %g to triangle input region list\n",
                          label->x, label->y, impl->in.regionlist[j+2]);
                WarnMessage(buf);
            }
#endif // DEBUG
            // Note(ZaJ): this is the code that was used impl->in the periodic bc triangulation:
            //  if (label->MaxArea>0 && (label->MaxArea<defaultMeshSize))
            //      impl->in.regionlist[j+3] = label->MaxArea;  // Area constraint
            //  else
            //      impl->in.regionlist[j+3] = defaultMeshSize;
            // ... which is equivalent to the code below (if forceMaxMeshArea is true).
            // ... the code below is a copy of the nonperiodic case (if forceMaxMeshArea is set to problem->DoForceMaxMeshArea)

            // Area constraint
            if (label->MaxArea <= 0)
            {
                // if no mesh size has been specified use the default
                impl->in.regionlist[j+3] = defaultMeshSize;
            }
            else if ((label->MaxArea > defaultMeshSize) && (forceMaxMeshArea))
            {
                // if the user has specied that FEMM should choose an
                // upper mesh size limit, regardles of their choice,
                // and their choice is less than that limit, change it
                // to that limit
                impl->in.regionlist[j+3] = defaultMeshSize;
            }
            else
            {
                // Use the user's choice of mesh size
                impl->in.regionlist[j+3] = label->MaxArea;
            }

            j += 4;
            k++;
        }
    }
    return true;
}

int TriangleMesher::triangulate(bool verbose)
{
    std::string triArgs = triangulateParams(verbose);
    // this is a mess, but building the string with std::string is more flexible than sprintf
    // (and the triangulate api is ancient)
    char cmdline[512];
    sprintf(cmdline, "%s",triArgs.c_str());

#ifdef XFEMM_BUILTIN_TRIANGLE
    int tristatus = ::triangulate(cmdline, &impl->in, &impl->out, (struct triangulateio *) nullptr, this->TriMessage);
    if (tristatus!=0)
    {
        std::string msg = "Call to triangulate failed with status code: " + to_string(tristatus) +"\n";
        WarnMessage(msg.c_str());
        return tristatus;
    }
#else
    // parse options
    int tristatus = triangle_context_options(impl->ctx, cmdline);
    if (tristatus != TRI_OK)
    {
        WarnMessage("Invalid option string for triangle!\n");
        return tristatus;
    }
    // Triangulate the polygon.
    tristatus = triangle_mesh_create(impl->ctx, &impl->in);
    if (tristatus != TRI_OK)
    {
        std::string msg = "Call to triangulate failed with status code: " + to_string(tristatus) +"\n";
        WarnMessage(msg.c_str());
        return tristatus;
    }
#endif
    return 0;
}

string TriangleMesher::triangulateParams(bool verbose) const
{
    // An explaination of the input parameters used for Triangle
    //
    // -p Triangulates a Planar Straight Line Graph, i.e. list of segments.
    // -P Suppresses the output .poly file.
    // -q Quality mesh generation with no angles smaller than specified impl->in the following number
    // -e Outputs a list of edges of the triangulation.
    // -A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
    // -a Imposes a maximum triangle area constraint.
    // -z Numbers all items starting from zero (rather than one)
    // -I Suppresses mesh iteration numbers
    // -j prevents duplicated input vertices, or vertices `eaten' by holes,
    //    from appearing impl->in the output .node file.  Thus, if two input vertices
    //    have exactly the same coordinates, only the first appears impl->in the
    //    output.
    // -Y Suppresses the creation of Steiner points on the exterior boundary.
    //
    // See http://www.cs.cmu.edu/~quake/triangle.switch.html for more info
    std::string triArgs = "-pPq" + to_string(impl->m_minAngle) + "eAaz" + (verbose?"":"Q") + "I";
    if (impl->m_suppressUnusedVertices)
        triArgs += "j";
    if (impl->m_suppressExteriorSteinerPoints)
        triArgs += "Y";

    return triArgs;
}

bool TriangleMesher::writePolyFile(string filename, std::string comment) const
{
    std::ofstream polyFile (filename);
    // set floating point precision once for the whole stream
    polyFile << std::setprecision(17);
    // when filling to a width, adjust to the left
    polyFile.setf(std::ios::left);

    polyFile << impl->in.numberofpoints << "\t2\t0\t1\n";
    for (int i=0; i < impl->in.numberofpoints; i++)
    {
        polyFile << i << "\t" << impl->in.pointlist[2*i] << "\t" << impl->in.pointlist[2*i+1] << "\t" << impl->in.pointmarkerlist[i] << "\n";
    }

    polyFile << impl->in.numberofsegments << "\t1\n";
    for (int i=0; i < impl->in.numberofsegments; i++)
    {
        polyFile << i << "\t" << impl->in.segmentlist[2*i] << "\t" << impl->in.segmentlist[2*i+1] << "\t" << impl->in.segmentmarkerlist[i] <<"\n";
    }

    polyFile << impl->in.numberofholes << "\n";
    for (int i=0; i < impl->in.numberofholes; i++)
    {
        polyFile << i << "\t" << impl->in.holelist[2*i] << "\t" << impl->in.holelist[2*i+1] << "\n";
    }

    polyFile << impl->in.numberofregions << "\n";
    for (int i=0; i < impl->in.numberofregions; i++)
    {
        int j=4*i;
        polyFile << i << "\t"
                 << impl->in.regionlist[j] << "\t"
                 << impl->in.regionlist[j+1] << "\t"
                 << impl->in.regionlist[j+2] << "\t"
                 << impl->in.regionlist[j+3] << "\n";
    }

    polyFile << "# " << comment << "\n";
    return true;
}

void TriangleMesher::setMinAngle(double value)
{
    impl->m_minAngle = value;
}

void TriangleMesher::suppressExteriorSteinerPoints()
{
    impl->m_suppressExteriorSteinerPoints = true;
}

void TriangleMesher::suppressUnusedVertices()
{
    impl->m_suppressUnusedVertices = true;
}




femm::mesh::RawMesh TriangleMesher::rawTriangulation() const
{
    femm::mesh::RawMesh result;
#ifdef XFEMM_BUILTIN_TRIANGLE
    if (impl->out.numberofcorners != 3)
        throw runtime_error("Triangle produced non-linear elements");
    for (int i=0; i<impl->out.numberofpoints; ++i)
        result.points.push_back({impl->out.pointlist[2*i], impl->out.pointlist[2*i+1], impl->out.pointmarkerlist[i]});
    for (int i=0; i<impl->out.numberoftriangles; ++i) {
        femm::mesh::RawMesh::Triangle triangle;
        for (int j=0; j<3; ++j)
            triangle.nodes[j] = static_cast<femm::mesh::MeshIndex>(impl->out.trianglelist[i*3+j]);
        if (impl->out.numberoftriangleattributes > 0)
            triangle.regionAttribute = impl->out.triangleattributelist[i*impl->out.numberoftriangleattributes];
        result.triangles.push_back(triangle);
    }
    for (int i=0; i<impl->out.numberofedges; ++i)
        result.edges.push_back({static_cast<femm::mesh::MeshIndex>(impl->out.edgelist[2*i]),
                                static_cast<femm::mesh::MeshIndex>(impl->out.edgelist[2*i+1]),
                                impl->out.edgemarkerlist[i]});
#else
    // The external API deliberately keeps its mesh opaque.  Its stable writer
    // API is used as an interchange boundary, without exposing it publicly.
    if (triangle_check_mesh(impl->ctx) != 0) {
        WarnMessage("Mesh has topological inconsistencies!\n");
        throw runtime_error("Triangle returned an inconsistent mesh");
    }
    auto read = [](FILE *file, auto parser) { rewind(file); parser(file); fclose(file); };
    FILE *file = tmpfile();
    if (!file || triangle_write_nodes(impl->ctx, file) != TRI_OK) throw runtime_error("Could not extract Triangle nodes");
    read(file, [&result](FILE *f) { int count, dim, attrs, markers; fscanf(f,"%d%d%d%d",&count,&dim,&attrs,&markers); for(int i=0;i<count;++i){ int id,marker=0; double x,y; fscanf(f,"%d%lf%lf",&id,&x,&y); for(int j=0;j<attrs;++j){double ignored;fscanf(f,"%lf",&ignored);} if(markers)fscanf(f,"%d",&marker); result.points.push_back({x,y,marker}); }});
    file = tmpfile();
    if (!file || triangle_write_edges(impl->ctx, file) != TRI_OK) throw runtime_error("Could not extract Triangle edges");
    read(file, [&result](FILE *f) { int count,markers; fscanf(f,"%d%d",&count,&markers); for(int i=0;i<count;++i){int id,a,b,marker=0;fscanf(f,"%d%d%d",&id,&a,&b);if(markers)fscanf(f,"%d",&marker);result.edges.push_back({static_cast<femm::mesh::MeshIndex>(a),static_cast<femm::mesh::MeshIndex>(b),marker});}});
    file = tmpfile();
    if (!file || triangle_write_elements(impl->ctx, file) != TRI_OK) throw runtime_error("Could not extract Triangle elements");
    read(file, [&result](FILE *f) { int count,corners,attributes; fscanf(f,"%d%d%d",&count,&corners,&attributes); if(corners!=3)throw runtime_error("Triangle produced non-linear elements"); for(int i=0;i<count;++i){int id;fscanf(f,"%d",&id);femm::mesh::RawMesh::Triangle t;for(int j=0;j<3;++j){int n;fscanf(f,"%d",&n);t.nodes[j]=static_cast<femm::mesh::MeshIndex>(n);}for(int j=0;j<attributes;++j){double a;fscanf(f,"%lf",&a);if(j==0)t.regionAttribute=a;}result.triangles.push_back(t);}});
#endif
    return result;
}

} // namespace fmesher
