/*
   This code is a modified version of an algorithm
   forming part of the software program Finite
   Element Method Magnetics (FEMM), authored by
   David Meeker. The original software code is
   subject to the Aladdin Free Public Licence
   version 8, November 18, 1999. For more information
   on FEMM see www.femm.info. This modified version
   is not endorsed in any way by the original
   authors of FEMM.

   This software has been modified to use the C++
   standard template libraries and remove all Microsoft (TM)
   MFC dependent code to allow easier reuse across
   multiple operating system platforms.

   Date Modified: 2011 - 11 - 10
   By: Richard Crozier
   Contact: richard.crozier@yahoo.co.uk

   Additional changes:
   Copyright 2016-2017 Johannes Zarl-Zierl <johannes.zarl-zierl@jku.at>
   Contributions by Johannes Zarl-Zierl were funded by Linz Center of
   Mechatronics GmbH (LCM)
*/


// implementation of various incarnations of calls
// to triangle from the FMesher class
#include "fmesher.h"
#include "fparse.h"
#include "IntPoint.h"
#include "femmconstants.h"
#include "CCommonPoint.h"
#include "CAirGapElement.h"
//extern "C" {
#include "triangle.h"
#ifndef XFEMM_BUILTIN_TRIANGLE
#include "triangle_api.h"
#endif
//}

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <malloc.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#ifndef REAL
#define REAL double
#endif

#ifndef LineFraction
#define LineFraction 500.0
#endif

// Default mesh size is the diagonal of the geometry's
// bounding box divided by BoundingBoxFraction
#ifndef BoundingBoxFraction
#define BoundingBoxFraction 100.0
#endif

#define toDegrees(x) ((Im(x)>=0) ? arg(x) : (arg(x) + 2.*PI))*(180./PI)

using namespace std;
using namespace femm;
using namespace fmesher;
using namespace femmsolver;

namespace {

enum class PointMarkerInfo {
    None ///< Use zero / Don't store information in marker list
    , FromProblem ///< Generate marker info using the problem descripton
};
enum class SegmentMarkerInfo {
    FromCnt ///< Generate marker infor from cnt field of segments.
    , FromProblem ///< Generate marker info using the problem descripton
};

/**
 * @brief The TriangulateHelper class encapsulates the interface to triangle,
 * so that the rest of the code doesn't have to deal with changes in its api.
 *
 * All memory that is allocated by its member functions is freed in the destructor.
 * Don't call initialization functions more than once.
 */
class TriangulateHelper {
    using nodelist_t = std::vector<std::unique_ptr<CNode> >;
    using linelist_t = std::vector<std::unique_ptr<CSegment> >;
public:
    TriangulateHelper();
    ~TriangulateHelper();

    /**
     * @brief Build a point list and point marker list as input for triangle.
     * The point marker list is later used to make the connection of meshed nodes back to original nodes.
     * @param nodelst
     * @param problem
     * @return \c true on success, \c false on (allocation) error
     */
    bool initPointsWithMarkers(const nodelist_t &nodelst, const FemmProblem &problem, PointMarkerInfo info);
    /**
     * @brief Build a segment list and segment marker list as input for triangle.
     * The segment marker list is later used to make the connection of meshed arcs/lines back to original arcs/lines.
     * @param linelst
     * @param problem
     * @return \c true on success, \c false on (allocation) error
     */
    bool initSegmentsWithMarkers(const linelist_t &linelst, const FemmProblem &problem, SegmentMarkerInfo info);

    /**
     * @brief Build a list of holes and regions as input for triangle.
     * This translates the CBlockLabel info for triangle.
     * @param problem
     * @param forceMaxMeshArea if \c true, this enforces an upper bound (defaultMeshSize) for the size of regional attributes
     * @param defaultMeshSize size of regional attributes that are not valid (i.e. <=0 or over the upper bound (if enforced))
     * @return \c true on success, \c false on (allocation) error
     */
    bool initHolesAndRegions(const FemmProblem &problem, bool forceMaxMeshArea, double defaultMeshSize);

    /**
     * @brief triangulate
     * The values of minAngle and suppressExteriourSteinerPoints are applied.
     * @param verbose Verbosity of triangle
     * @return
     */
    int triangulate(bool verbose);

    /**
     * @brief Create a parameter list for consumption of triangle
     * @param verbose Verbosity of triangle
     * @return
     */
    std::string triangulateParams(bool verbose=false) const;

    /**
     * @brief Write an input file for triangle.
     * Normally, the \c.poly file is not used because triangle is called directly as a library.
     * This function is intended as a debugging aid.
     * @param filename the complete file name ending in ".poly"
     * @return \c true, if writing succeeded, \c false otherwise.
     */
    bool writePolyFile(std::string filename, std::string comment) const;
    bool writeTriangulationFiles(std::string Pathname) const;

    // pointer to function to call when issuing warning messages
    int (*WarnMessage)(const char*, ...);

    // pointer to function to use for triangle to issue warning messages
    int (*TriMessage)(const char * format, ...);

    void setMinAngle(double value);
    /**
     * @brief Suppress insertion of Steiner points on the mesh boundary.
     */
    void suppressExteriorSteinerPoints();
    /**
     * @brief Purge vertices that are not part of the final triangulation from the triangulation output.
     * Should be safe to add all the time, but I didn't test with periodic bc triangulation.
     * Therefore only used with nonperiodic triangulation.
     */
    void suppressUnusedVertices();

private:
#ifdef XFEMM_BUILTIN_TRIANGLE
    struct triangulateio in;
    struct triangulateio out;
#else
    triangleio in;
    context *ctx;
#endif
    double m_minAngle = 0.;
    bool m_suppressExteriorSteinerPoints = false;
    bool m_suppressUnusedVertices = false;
};

/**
 * @brief Initialize a triangulateio to all zero.
 * @param io
 */
#ifdef XFEMM_BUILTIN_TRIANGLE
void initialize(struct triangulateio &io)
#else
void initialize(triangleio &io)
#endif
{
    io.pointlist = nullptr;
    io.pointattributelist = nullptr;
    io.pointmarkerlist = nullptr;
    io.numberofpoints = 0;
    io.numberofpointattributes = 0;

    io.trianglelist = nullptr;
    io.triangleattributelist = nullptr;
    io.trianglearealist = nullptr;
    io.neighborlist = nullptr;
    io.numberoftriangles = 0;
    io.numberofcorners = 0;
    io.numberoftriangleattributes = 0;

    io.segmentlist = nullptr;
    io.segmentmarkerlist = nullptr;
    io.numberofsegments = 0;

    io.holelist = nullptr;
    io.numberofholes = 0;

    io.regionlist = nullptr;
    io.numberofregions = 0;

    io.edgelist = nullptr;
    io.edgemarkerlist = nullptr;
#ifdef XFEMM_BUILTIN_TRIANGLE
    io.normlist = nullptr; // only used by voronoi diagram
#endif
    io.numberofedges = 0;
}

}

double FMesher::averageLineLength() const
{
    double z=0;
    const double numLines = problem->linelist.size();
    for (const auto &line : problem->linelist)
    {
        z += problem->lengthOfLine(*line) / numLines;
    }
    return z;
}

double fmesher::defaultMeshSizeHeuristics(const std::vector<std::unique_ptr<CNode> > &nodelst, bool doSmartMesh)
{
    if (nodelst.empty())
        return -1;

    // compute minimum and maximum x/y values
    CComplex min=nodelst[0]->CC();
    CComplex max=min;
    for(const auto &node: nodelst)
    {
        if (node->x < min.re) min.re = node->x;
        if (node->y < min.im) min.im = node->y;
        if (node->x > max.re) max.re = node->x;
        if (node->y > max.im) max.im = node->y;
    }

    if (doSmartMesh)
    {
        double absdist = abs(max-min)/BoundingBoxFraction;
        return absdist * absdist;
    } else {
        return abs(max-min);
    }
}

void fmesher::discretizeInputSegments(const FemmProblem &problem, std::vector<std::unique_ptr<CNode> > &nodelst, std::vector<std::unique_ptr<CSegment> > &linelst, double dL, SegmentFilter filter)
{
    for(int i=0; i<(int)problem.linelist.size(); i++)
    {
        const CSegment &line = *problem.linelist[i];

        if (filter == SegmentFilter::OnlyUnselected && line.IsSelected )
            continue;

        const CNode &n0 = *problem.nodelist[line.n0];
        const CNode &n1 = *problem.nodelist[line.n1];
        const CComplex a0 = n0.CC();
        const CComplex a1 = n1.CC();
        // create working copy:
        CSegment segm = line;
        // use the cnt flag to carry a notation
        // of which line or arc in the input geometry a
        // particular segment is associated with
        // (this info is only used in the periodic BC triangulation, and is ignored in the nonperiodic one)
        segm.cnt = i;

        double lineLength = problem.lengthOfLine(line);

        int numParts;
        if (line.MaxSideLength == -1) {
            numParts = 1;
        }
        else{
            numParts = (unsigned int) std::ceil(lineLength/line.MaxSideLength);
        }

        if (numParts == 1) // default condition where discretization on line is not specified
        {
            if (lineLength < (3. * dL) || problem.DoSmartMesh == false)
            {
                // line is too short to add extra points
                linelst.push_back(segm.clone());
            }
            else{
//                // add extra points at a distance of dL from the ends of the line.
//                // this forces Triangle to finely mesh near corners
//                int l = (int) nodelst.size();
//
//                // first part
//                CComplex a2 = a0 + dL * (a1-a0) / abs(a1-a0);
//                CNode node1 (a2.re, a2.im);
//                nodelst.push_back(node1.clone());
//                segm.n0 = line.n0;
//                segm.n1 = l;
//                linelst.push_back(segm.clone());
//
//                // middle part
//                a2 = a1 + dL * (a0-a1) / abs(a1-a0);
//                CNode node2 (a2.re, a2.im);
//                nodelst.push_back(node2.clone());
//                segm.n0 = l;
//                segm.n1 = l + 1;
//                linelst.push_back(segm.clone());
//
//                // end part
//                segm.n0 = l + 2;
//                segm.n1 = line.n1;
//                linelst.push_back(segm.clone());

// add extra points at a distance of dL from the ends of the line.
                CComplex a2;
                CNode node;
                int l = 0;

				// this forces Triangle to finely mesh near corners
				for(int j=0;j<3;j++)
				{
					if(j==0)
					{
						a2=a0+dL*(a1-a0)/abs(a1-a0);
						node.x=a2.re; node.y=a2.im;
						l=(int) nodelst.size();
						nodelst.push_back (node.clone());
						segm.n0=line.n0;
						segm.n1=l;
						linelst.push_back (segm.clone());
					}

					if(j==1)
					{
						a2=a1+dL*(a0-a1)/abs(a1-a0);
						node.x=a2.re; node.y=a2.im;
						l=(int) nodelst.size ();
						nodelst.push_back(node.clone());
						segm.n0=l-1;
						segm.n1=l;
						linelst.push_back(segm.clone());
					}

					if(j==2)
					{
						l=(int) nodelst.size()-1;
						segm.n0=l;
						segm.n1=line.n1;
						linelst.push_back(segm.clone());
					}

				}
            }
        }
        else{
            for(int j=0; j<numParts; j++)
            {
                CComplex a2 = a0 + (a1-a0)*((double) (j+1)) / ((double) numParts);
                CNode node (a2.re, a2.im);
                if(j == 0){
                    // first part -> n0 == line.n0
                    int l=nodelst.size();
                    nodelst.push_back(node.clone());
                    segm.n0=line.n0;
                    segm.n1=l;
                    linelst.push_back(segm.clone());
                }
                else if(j == (numParts-1))
                {
                    // last part -> n1 == line.n1 ; endpoint already exists
                    int l=nodelst.size()-1;
                    segm.n0=l;
                    segm.n1=line.n1;
                    linelst.push_back(segm.clone());
                }
                else{
                    int l=nodelst.size();
                    nodelst.push_back(node.clone());
                    segm.n0=l-1;
                    segm.n1=l;
                    linelst.push_back(segm.clone());
                }
            }
        }
    }
}

void fmesher::discretizeInputArcSegments(const FemmProblem &problem, std::vector<std::unique_ptr<CNode> > &nodelst, std::vector<std::unique_ptr<CSegment> > &linelst, SegmentFilter filter)
{
    for(int i=0;i<(int)problem.arclist.size();i++)
    {

        problem.arclist[i]->mySideLength = problem.arclist[i]->MaxSideLength;

        const CArcSegment &arc = *problem.arclist[i];

        if (filter == SegmentFilter::OnlyUnselected && arc.IsSelected )
            continue;

        // smart meshing does not apply to arc segments
        assert(arc.MaxSideLength != -1);

        // create working copy:
        CSegment segm = arc;
        // use the cnt flag to carry a notation
        // of which line or arc in the input geometry a
        // particular segment is associated with
        // (this info is only used in the periodic BC triangulation, and is ignored in the nonperiodic one)
        segm.cnt=i+problem.linelist.size();

        segm.BoundaryMarkerName=arc.BoundaryMarkerName;
        if (problem.filetype != FileType::MagneticsFile)
            segm.InConductorName=arc.InConductorName; // not relevant/compatible to magnetics problems

        int numParts=(int) ceil(arc.ArcLength/arc.MaxSideLength);

        CComplex center;
        double R=0;
        problem.getCircle(arc,center,R);

        CComplex a1=exp(I*arc.ArcLength*PI/(((double) numParts)*180.));
        CComplex a2=problem.nodelist[arc.n0]->CC();

        if(numParts==1){
            linelst.push_back(segm.clone());
        }
        else for(int j=0;j<numParts;j++)
        {
            // move point along arc
            a2=(a2-center)*a1+center;
            CNode node(a2.re,a2.im);
            int l = (int)nodelst.size();
            if(j==0){
                // first part -> n0 == arc.n0
                nodelst.push_back(node.clone());
                segm.n0=arc.n0;
                segm.n1=l;
            }
            else if(j==(numParts-1))
            {
                // last part -> n1 == arc.n1 ; endpoint already exists
                segm.n0=l-1;
                segm.n1=arc.n1;
            }
            else{
                nodelst.push_back(node.clone());
                segm.n0=l-1;
                segm.n1=l;
            }
            linelst.push_back(segm.clone());
        }
    }
}

/**
 * @brief FMesher::HasPeriodicBC
 * @return \c true, if there are periodic or antiperiodic boundary conditions, \c false otherwise.
 *
 * \internal
 *  * \femm42{femm/writepoly.cpp}
 *  * \femm42{femm/bd_writepoly.cpp}
 *  * \femm42{femm/hd_writepoly.cpp}
 */
bool FMesher::HasPeriodicBC()
{
    bool hasPeriodicBC = false;
    for (const auto &bdryProp: problem->lineproplist)
    {
        if (bdryProp->isPeriodic())
            hasPeriodicBC = true;
    }
    // if flag is false, there can't be any lines
    // with periodic BC's, because no periodic boundary
    // conditions have been defined.
    if (hasPeriodicBC==false) return false;

    //now, if there are some periodic boundary conditions,
    //we have to check to see if any have actually been
    //applied to the model

    // first, test the segments
    for(int i=0;i<(int)problem->linelist.size();i++)
    {
        int k=-1;
        for(int j=0;j<(int)problem->lineproplist.size();j++)
        {
            if(problem->lineproplist[j]->BdryName==
               problem->linelist[i]->BoundaryMarkerName)
            {
                k=j;
                break;
            }
        }
        if(k>=0){
            if (problem->lineproplist[k]->isPeriodic())
            {
                return true;
            }
        }
    }

    // If we've gotten this far, we still need to check the
    // arc segments.
    for(int i=0;i<(int)problem->arclist.size();i++)
    {
        int k=-1;
        for(int j=0;j<(int)problem->lineproplist.size();j++)
        {
            if(problem->lineproplist[j]->BdryName==
               problem->arclist[i]->BoundaryMarkerName)
            {
                k=j;
                break;
            }
        }
        if(k>=0){
            if (problem->lineproplist[k]->isPeriodic())
            {
                return true;
            }
        }
    }

    // Finally, we're done.
    // we could not find periodic and/or antiperiodic boundaries.
    return false;
}


bool TriangulateHelper::writeTriangulationFiles(string PathName) const
{
    FILE *fp;
    std::string msg;
    std::string plyname;

#ifndef XFEMM_BUILTIN_TRIANGLE
    if (triangle_check_mesh(ctx)!=0)
    {
        WarnMessage("Mesh has topological inconsistencies!\n");
        return false;
    }
#endif

    // write the .node file
    plyname = PathName.substr(0, PathName.find_last_of('.')) + ".node";

    // check to see if we are ready to write a .node datafile containing
    // the nodes

    if ((fp = fopen(plyname.c_str(),"wt"))==NULL){
        WarnMessage("Couldn't write to specified .node file");
        return false;
    }

#ifdef XFEMM_BUILTIN_TRIANGLE
    if (out.numberofpoints > 0)
    {
        // <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
        fprintf(fp, "%i\t%i\t%i\t%i\n", out.numberofpoints, 2, 0, 1);
        //fprintf(fp, "%i\t%i\t%i\n", out.numberofpoints, 2, out.numberofpoints, 1);

        // <vertex #> <x> <y> [attributes] [boundary marker]
        for(int i = 0; i < (2 * out.numberofpoints) - 1; i = i + 2)
        {
            fprintf(fp, "%i\t%.17g\t%.17g\t%i\n", i/2, out.pointlist[i], out.pointlist[i+1], out.pointmarkerlist[i/2]);
        }

        fclose(fp);
    }
#else
    int status = triangle_write_nodes(ctx, fp);
    fclose(fp);
    if (status != TRI_OK)
    {
        msg = "Failed to write to specified .node file\n";
        WarnMessage(msg.c_str());
        return false;
    }
#endif

    // write the .edge file
    plyname = PathName.substr(0, PathName.find_last_of('.')) + ".edge";

    // check to see if we are ready to write an edge datafile;

    if ((fp = fopen(plyname.c_str(),"wt"))==NULL){
        msg = "Couldn't write to specified .edge file\n";
        WarnMessage(msg.c_str());
        return false;
    }
#ifdef XFEMM_BUILTIN_TRIANGLE

    if (out.numberofedges > 0)
    {
        // write number of edges, number of boundary markers, 0 or 1
        fprintf(fp, "%i\t%i\n", out.numberofedges, 1);

        // write the edges in the format
        // <edge #> <endpoint> <endpoint> [boundary marker]
        // Endpoints are indices into the corresponding .edge file.
        for(int i=0; i < 2 * (out.numberofedges) - 1; i = i + 2)
        {
            fprintf(fp, "%i\t%i\t%i\t%i\n", i/2, out.edgelist[i], out.edgelist[i+1], out.edgemarkerlist[i/2]);
        }

        fclose(fp);
    } else {
        WarnMessage("No edges to write!\n");
    }
#else
    // Note: triangle_write_edges also numbers the edges, which is required for writing the .ele file
    status = triangle_write_edges(ctx, fp);
    fclose(fp);
    if (status != TRI_OK)
    {
        msg = "Failed to write to specified .edge file\n";
        WarnMessage(msg.c_str());
        return false;
    }
#endif

    // write the .ele file
    plyname = PathName.substr(0, PathName.find_last_of('.')) + ".ele";


    // check to see if we are ready to write a .ele datafile containing
    // thr triangle elements

    if ((fp = fopen(plyname.c_str(),"wt"))==NULL){
        WarnMessage("Couldn't write to specified .ele file");
        return false;
    }

#ifdef XFEMM_BUILTIN_TRIANGLE
    if (out.numberoftriangles > 0)
    {
        // write number of triangle elements, number of corners per triangle and
        // the number of attributes per triangle
        fprintf(fp, "%i\t%i\t%i\n", out.numberoftriangles, out.numberofcorners, out.numberoftriangleattributes);

        // write the triangle info to the file with the format
        // <triangle #> <node> <node> <node> ... [attributes]
        // Endpoints are indices into the corresponding .node file.
        for(int i=0, nexttriattrib=0; i < (out.numberofcorners) * (out.numberoftriangles) - (out.numberofcorners - 1); i = i + (out.numberofcorners))
        {
            // print the triangle number
            fprintf(fp, "%i\t", i / (out.numberofcorners));

            // print the corner nodes
            for (int j = 0; j < (out.numberofcorners); j++)
            {
                fprintf(fp, "%i\t", out.trianglelist[i+j]);
            }

            // print the triangle attributes, if there are any
            if (out.numberoftriangleattributes > 0)
            {
                for(int j = 0; j < (out.numberoftriangleattributes); j++)
                {
                    fprintf(fp, "%.17g\t", out.triangleattributelist[nexttriattrib+j]);
                }

                // set the position of the next set of triangle attributes
                nexttriattrib = nexttriattrib + (out.numberoftriangleattributes);
            }

            // go to the next line
            fprintf(fp, "\n");
        }

        fclose(fp);

    }
#else
    status = triangle_write_elements(ctx, fp);
    fclose(fp);
    if (status != TRI_OK)
    {
        msg = "Failed to write to specified .ele file\n";
        WarnMessage(msg.c_str());
        return false;
    }
#endif
    return true;
}

/**
 * @brief FMesher::DoNonPeriodicBCTriangulation
 * What we do in the normal case is DoNonPeriodicBCTriangulation
 * @param PathName
 * @return
 *
 * Original function name:
 *  * \femm42{femm/writepoly.cpp,CFemmeDoc::OnWritePoly()}
 *  * \femm42{femm/bd_writepoly.cpp,CbeladrawDoc::OnWritePoly()}
 *  * \femm42{femm/hd_writepoly.cpp,ChdrawDoc::OnWritePoly()}
 */
int FMesher::DoNonPeriodicBCTriangulation(string PathName)
{
    // // if incremental permeability solution, we crib mesh from the previous problem.
    // // we can just bail out in that case.
    // if (!problem->previousSolutionFile.empty() && problem->Frequency>0)
    //     return true;

    FILE *fp;
    double dL;
    //CStdString s;
    string plyname;
    std::vector < std::unique_ptr<CNode> >       nodelst;
    std::vector < std::unique_ptr<CSegment> >    linelst;

#ifdef DEBUG
    WarnMessage("writepoly: beginning NON periodic boundary triangulation\n");
#endif // DEBUG

    nodelst.clear();
    linelst.clear();
    // calculate length used to kludge fine meshing near input node points
    dL = averageLineLength() / LineFraction;

    // copy node list as it is;
    for (const auto &node : problem->nodelist)
        nodelst.push_back(node->clone());

    problem->clearNotationTags();
    // discretize input segments
    discretizeInputSegments(*problem, nodelst, linelst, dL);

    // discretize input arc segments
    discretizeInputArcSegments(*problem, nodelst, linelst);

    // create correct output filename;
    string pn = PathName;

    // figure out a good default mesh size for block labels where
    // mesh size isn't explicitly specified
    double DefaultMeshSize = defaultMeshSizeHeuristics(nodelst, problem->DoSmartMesh);

//    for(i=0,k=0;i<blocklist.size();i++)
//        if(blocklist[i]->BlockTypeName=="<No Mesh>")
//        {
//            fprintf(fp,"%i    %.17g    %.17g\n",k,blocklist[i]->x,blocklist[i]->y);
//            k++;
//        }
//
//    // write out regional attributes
//    fprintf(fp,"%i\n",blocklist.size()-j);
//
//    for(i=0,k=0;i<blocklist.size();i++)
//        if(blocklist[i]->BlockTypeName!="<No Mesh>" && (blocklist[i].BlockType!="<Inf>"))
//        {
//            fprintf(fp,"%i    %.17g    %.17g    ",k,blocklist[i]->x,blocklist[i]->y);
//            fprintf(fp,"%i    ",k+1);
//            if (blocklist[i]->MaxArea>0)
//                fprintf(fp,"%.17g\n",blocklist[i]->MaxArea);
//            else fprintf(fp,"-1\n");
//            k++;
//        }
//    fclose(fp);

    // write out a trivial pbc file
    plyname = pn.substr(0,pn.find_last_of('.')) + ".pbc";
    if ((fp=fopen(plyname.c_str(),"wt"))==NULL){
        WarnMessage("Couldn't write to specified .pbc file");
        return -1;
    }
    fprintf(fp,"0\n0\n");
    fclose(fp);

    // **********         call triangle       ***********

    {
        TriangulateHelper triHelper;
        triHelper.WarnMessage = WarnMessage;
        triHelper.TriMessage = this->TriMessage;
        if (!triHelper.initPointsWithMarkers(nodelst,*problem, PointMarkerInfo::FromProblem))
            return -1;
        if (!triHelper.initSegmentsWithMarkers(linelst,*problem,SegmentMarkerInfo::FromProblem))
            return -1;
        if (!triHelper.initHolesAndRegions(*problem, problem->DoForceMaxMeshArea, DefaultMeshSize))
            return -1;
        triHelper.setMinAngle(std::min(problem->MinAngle+MINANGLE_BUMP,MINANGLE_MAX));
        triHelper.suppressUnusedVertices();
        if (writePolyFiles)
        {
            string plyname = PathName.substr(0, PathName.find_last_of('.')) + ".poly";
            triHelper.writePolyFile(plyname, triHelper.triangulateParams());
        }
        int tristatus = triHelper.triangulate(Verbose);
        if (tristatus != 0)
            return tristatus;

        triHelper.writeTriangulationFiles(PathName);
    }
    problem->clearNotationTags();

    return 0;
}


/**
 * \brief Call triangle twice to order segments on the boundary properly
 * for periodic or antiperiodic boundary conditions
 *
 * Original function name:
 *  * \femm42{femm/writepoly.cpp,CFemmeDoc::FunnyOnWritePoly()}
 *  * \femm42{femm/bd_writepoly.cpp,CbeladrawDoc::FunnyOnWritePoly()}
 *  * \femm42{femm/hd_writepoly.cpp,ChdrawDoc::FunnyOnWritePoly()}
 */
int FMesher::DoPeriodicBCTriangulation(string PathName)
{
    // // if incremental permeability solution, we crib mesh from the previous problem.
    // // we can just bail out in that case.
    // if (!problem->previousSolutionFile.empty() && problem->Frequency>0)
    //     return true;
    FILE *fp;
    int i, j, k, n;
    int l,n0,n1,n2;
    double z,R,dL;
    CComplex a0,a1,a2,c;
    CComplex b0,b1,b2;
    char instring[1024];
    //string s;
    string plyname;
    std::vector < std::unique_ptr<CNode> >              nodelst;
    std::vector < std::unique_ptr<CSegment> >           linelst;
    //std::vector < std::unique_ptr<CCBlockLabel> >       blocklst;
    std::vector < std::unique_ptr<CPeriodicBoundary> >  pbclst;
    std::vector < std::unique_ptr<CAirGapElement> >     agelst;
    std::vector < std::unique_ptr<CCommonPoint> >       ptlst;
    CNode node;
    CSegment segm;
    CCommonPoint pt;
    CPeriodicBoundary pbc;
    CAirGapElement age;

#ifdef DEBUG
    WarnMessage("writepoly: beginning periodic boundary triangulation\n");
#endif // DEBUG

    problem->updateUndo();

    // calculate length used to kludge fine meshing near input node points
    dL = averageLineLength() / LineFraction;

    // copy node list as it is;
    for (const auto &node : problem->nodelist)
        nodelst.push_back(node->clone());

    problem->clearNotationTags();
    // discretize input segments
    discretizeInputSegments(*problem, nodelst, linelst, dL);

    // discretize input arc segments
    discretizeInputArcSegments(*problem, nodelst, linelst);


    // create correct output filename;
    string pn = PathName;

    // figure out a good default mesh size for block labels where
    // mesh size isn't explicitly specified
    double DefaultMeshSize = defaultMeshSizeHeuristics(nodelst, problem->DoSmartMesh);

#ifdef DEBUG
    WarnMessage("writepoly: about to call triangle\n");
#endif // DEBUG

    // **********         call triangle       ***********

    {
        TriangulateHelper triHelper;
        triHelper.WarnMessage = WarnMessage;
        triHelper.TriMessage = this->TriMessage;

        if (!triHelper.initPointsWithMarkers(nodelst,*problem, PointMarkerInfo::None))
            return -1;
        if (!triHelper.initSegmentsWithMarkers(linelst,*problem,SegmentMarkerInfo::FromCnt))
            return -1;
        if (!triHelper.initHolesAndRegions(*problem, true, DefaultMeshSize))
            return -1;

        triHelper.setMinAngle(problem->MinAngle);
        if (writePolyFiles)
        {
            string plyname = PathName.substr(0, PathName.find_last_of('.')) + ".raw.poly";
            triHelper.writePolyFile(plyname, triHelper.triangulateParams());
        }
        int tristatus = triHelper.triangulate(Verbose);
        if (tristatus != 0)
            return tristatus;

        triHelper.writeTriangulationFiles(PathName);
    }

#ifdef DEBUG
    WarnMessage("writepoly: finished calling triangle\n");
#endif // DEBUG

    // So far, so good.  Now, read back in the .edge file
    // to make sure the points in the segments and arc
    // segments are ordered in a consistent way so that
    // the (anti)periodic boundary conditions can be applied.

#ifdef DEBUG
    WarnMessage("writepoly: 876\n");
#endif // DEBUG

    // read meshlines;
    plyname = pn.substr(0,pn.find_last_of('.')) + ".edge";
    if((fp=fopen(plyname.c_str(),"rt"))==NULL){
        WarnMessage("Call to triangle was unsuccessful\n");
        problem->undo();  problem->unselectAll();
        return -1;
    }
    fgets(instring,1024,fp);
    sscanf(instring,"%i",&k);
    problem->clearNotationTags();
    // use cnt again to keep a
    // tally of how many subsegments each
    // entity is sliced into.
    for(auto &arc: problem->arclist) arc->cnt=0;

#ifdef DEBUG
    WarnMessage("writepoly: 898\n");
#endif // DEBUG

    // resize initializes the new elements using the default ctor:
    ptlst.clear();
    ptlst.shrink_to_fit();
    int npt = problem->linelist.size()+problem->arclist.size();
    ptlst.reserve(npt);
    for(i=0; i<npt; i++)
        ptlst.push_back(std::unique_ptr <CCommonPoint> (new CCommonPoint()));

#ifdef DEBUG
    WarnMessage("writepoly: 910\n");
#endif // DEBUG

    for(i=0;i<k;i++)
    {
        // get the next edge from the file
        fgets(instring,1024,fp);
        // get the edge number, start and end points (n0 and n1) and the
        // segment/arc marker j
        sscanf(instring,"%i    %i    %i    %i",&l,&n0,&n1,&j);
        // if j != 0, this edge is part of a segment/arc
        if(j!=0)
        {
            // convert back to the `right' numbering
            j=-(j+2);
            assert(j>=0);

            // store a reference line that we can use to
            // determine whether or not this is a
            // boundary segment w/out re-running triangle.
            if (ptlst[j]->t==0)
            {
                ptlst[j]->t=1;
                ptlst[j]->setSortedValues(n0,n1);
            }

            if(j<(int)problem->linelist.size())
            {
                // deal with segments

                // increment cnt for the segment which the edge we are
                // examining is a part of to get a tally of how many edges
                // are a part of the segment/boundary
                problem->linelist[j]->cnt++;
                // check if the end n0 of the segment is the same node as the
                // end n1 of the edge, or if the end n1 of the segment is the
                // same node and end n0 of the edge. If so, flip the direction
                // of the segment
                if((problem->linelist[j]->n0 == n1) || (problem->linelist[j]->n1 == n0))
                {
                    // flip the end points of the segment
                    std::swap(problem->linelist[j]->n0,problem->linelist[j]->n1);
                }
            }
            else{
                // deal with arc segments;
                // Can't just flip the point order with
                // impunity in the arc segments, so we flip
                // a marker which denotes which side the
                // normal is on.

                j=j-(int)problem->linelist.size();
                problem->arclist[j]->cnt++;
                // Note(ZaJ): I assume that the second if statement should just be the else branch of the first if statement...
                if((problem->arclist[j]->n0==n1) || (problem->arclist[j]->n1==n0))
                    problem->arclist[j]->NormalDirection=false;
                if((problem->arclist[j]->n0==n0) || (problem->arclist[j]->n1==n1))
                    problem->arclist[j]->NormalDirection=true;
            }
        }
    }
    fclose(fp);

#ifdef DEBUG
    WarnMessage("writepoly: 974\n");
#endif // DEBUG

    // figure out which segments / arcsegments are on the
    // boundary and force an appropriate mesh density on
    // these based on how many divisions are in the first
    // trial meshing of the domain.

    // paw through the element list to find out how many
    // elements each reference segment appears in.  If a
    // segment is on the boundary, it ought to appear in just
    // one element.  Otherwise, it appears in two.
    plyname = pn.substr(0,pn.find_last_of('.')) + ".ele";
    if((fp=fopen(plyname.c_str(),"rt"))==NULL){
        WarnMessage("Call to triangle was unsuccessful");
        problem->undo();  problem->unselectAll();
        return -1;
    }
    fgets(instring,1024,fp);
    sscanf(instring,"%i",&k);

#ifdef DEBUG
    WarnMessage("writepoly: 996\n");
#endif // DEBUG

    for(i=0;i<k;i++)
    {
        fgets(instring,1024,fp);
        sscanf(instring,"%i    %i    %i    %i",&j,&n0,&n1,&n2);

        // Sort out the three nodes...
        if (n0>n1) { n=n0; n0=n1; n1=n; }
        if (n1>n2) { n=n1; n1=n2; n2=n; }
        if (n0>n1) { n=n0; n0=n1; n1=n; }

        // now, check to see if any of the test segments
        // are sides of this node...
        for(j=0;j<(int)ptlst.size();j++)
        {
            if ((n0==ptlst[j]->x) && (n1==ptlst[j]->y)) ptlst[j]->t--;
            if ((n0==ptlst[j]->x) && (n2==ptlst[j]->y)) ptlst[j]->t--;
            if ((n1==ptlst[j]->x) && (n2==ptlst[j]->y)) ptlst[j]->t--;
        }
    }
    fclose(fp);

#ifdef DEBUG
    WarnMessage("writepoly: 1021\n");
#endif // DEBUG

    // impose "new" mesh constraints on bdry arcs and segments....
    for(i=0; i < (int)problem->linelist.size(); i++)
    {
        if (ptlst[i]->t == 0)
        {
            // simply make the max side length equal to the
            // length of the boundary divided by the number
            // of elements that were created in the first
            // attempt at meshing
            problem->linelist[i]->MaxSideLength = problem->lengthOfLine(i) / ((double) problem->linelist[i]->cnt);
        }
    }

    for(i=0; i < (int)problem->arclist.size(); i++)
    {
        if (ptlst[i+problem->linelist.size()]->t == 0)
        {
            // alter maxsidelength, but do it in such
            // a way that it carries only 4 significant
            // digits.  There's no use in carrying double
            // precision here, because it looks crappy
            // when you open up the arc segment to see
            // its properties.
            char kludge[32];
            double newMaxSideLength;
            newMaxSideLength = problem->arclist[i]->ArcLength/((double) problem->arclist[i]->cnt);
            sprintf(kludge,"%.1e",newMaxSideLength);
            sscanf(kludge,"%lf",&newMaxSideLength);

            problem->arclist[i]->MaxSideLength = newMaxSideLength;
        }
    }

    ptlst.clear();
    ptlst.shrink_to_fit();

#ifdef DEBUG
    {
        char buf[1048];
        SNPRINTF(buf, sizeof(buf), "writepoly: lineproplist.size() is %li\n", problem->lineproplist.size());
        WarnMessage(buf);
    }
#endif // DEBUG

	// Search through defined bc's for pbcs and ages;
	// Allocate space to store their properties if they are detected
	for(i=0;i<(int)problem->lineproplist.size();i++)
	{
#ifdef DEBUG
        {
            char buf[1048];
            SNPRINTF(buf, sizeof(buf), "writepoly: searching for periodic boundaries, checking boundary %i called %s\n", i, problem->lineproplist[i]->BdryName.c_str ());
            WarnMessage(buf);
        }
#endif // DEBUG

		// pbc
		if ( (problem->lineproplist[i]->BdryFormat==4)
             || (problem->lineproplist[i]->BdryFormat==5))
        {
#ifdef DEBUG
        {
            char buf[1048];
            SNPRINTF(buf, sizeof(buf), "writepoly: %s is periodic boundary, updating pbc\n", problem->lineproplist[i]->BdryName.c_str ());
            WarnMessage(buf);
        }
#endif // DEBUG
			pbc.BdryName=problem->lineproplist[i]->BdryName;
			pbc.BdryFormat = problem->lineproplist[i]->BdryFormat - 4; // 0 for pbc, 1 for apbc

			pbc.antiPeriodic = problem->lineproplist[i]->isPeriodic(CBoundaryProp::PeriodicityType::AntiPeriodic);

			pbclst.push_back(pbc.clone());
		}

		// age
		if ( (problem->lineproplist[i]->BdryFormat==6)
             || (problem->lineproplist[i]->BdryFormat==7) )
		{
#ifdef DEBUG
            {
                char buf[1048];
                SNPRINTF(buf, sizeof(buf), "writepoly: found air gap boundary in boundary %i\n", i);
                WarnMessage(buf);
            }
#endif // DEBUG
			// only add an AGE to the list if it's actually being used
			for(j=0,k=0;j<(int)problem->arclist.size();j++)
				if (problem->arclist[j]->BoundaryMarkerName==problem->lineproplist[i]->BdryName) k++;
			if (k>1)
			{
				age.BdryName=problem->lineproplist[i]->BdryName;
				age.BdryFormat=problem->lineproplist[i]->BdryFormat-6; // 0 for pbc, 1 for apbc
				age.InnerAngle=problem->lineproplist[i]->InnerAngle;
				age.OuterAngle=problem->lineproplist[i]->OuterAngle;
				agelst.push_back(age.clone());
#ifdef DEBUG
                {
                    char buf[1048];
                    SNPRINTF(buf, sizeof(buf), "writepoly: added air gap boundary to agelist\n");
                    WarnMessage(buf);
                }
#endif // DEBUG
			}
		}
	}
#ifdef DEBUG
    {
        char buf[1048];
        SNPRINTF(buf, sizeof(buf), "writepoly: found %i air gap boundaries\n", agelst.size ());
        WarnMessage(buf);
    }
#endif // DEBUG

	// make sure all Air Gap Element arcs have the same discretization
	// for each arc associated with a particular Air Gap Element...

	// find out the total arc length and arc elements
	// corresponding ot each lineproplist entry
	for(i=0;i<(int)problem->arclist.size();i++)
	{
		if (problem->arclist[i]->BoundaryMarkerName!="<None>")
		{
			for(j=0;j<(int)agelst.size();j++)
			{

				if (problem->arclist[i]->BoundaryMarkerName==agelst[j]->BdryName)
				{
					agelst[j]->totalArcLength += problem->arclist[i]->ArcLength;
					agelst[j]->totalArcElements += problem->arclist[i]->IsSelected;

					problem->GetCircle(*(problem->arclist[i]),agelst[j]->agc,R);
					if (agelst[j]->ro==0)
					{
						agelst[j]->ri=R;
						agelst[j]->ro=R;
					}
					if (R>agelst[j]->ro) agelst[j]->ro=R;
					if (R<agelst[j]->ri) agelst[j]->ri=R;

					break;
				}
			}
		}
	}

	// cycle through AGEs and fix constituent arcs so that all arcs have the same discretization
	for (i=0;i<(int)agelst.size();i++)
	{
		if (agelst[i]->totalArcLength>0) // if the AGE is actually in play
		{
			char kludge[32];
			double myMaxSideLength,altMaxSideLength;

			myMaxSideLength=agelst[i]->totalArcLength/agelst[i]->totalArcElements;
			agelst[i]->totalArcLength/=2;	// this is now the angle spanned by the AGE

			// however, don't want long, skinny air gap elements.  Impose some limits
			// based on the inner and outer radii;
			altMaxSideLength=(360./PI)*(agelst[i]->ro-agelst[i]->ri)/(agelst[i]->ro+agelst[i]->ri);
			if (altMaxSideLength<myMaxSideLength) myMaxSideLength=altMaxSideLength;
			sprintf(kludge,"%.1e",myMaxSideLength);
			sscanf(kludge,"%lf",&myMaxSideLength);

			// apply new side length to all arcs in this AGE
			for(j=0;j<(int)problem->arclist.size();j++)
				if (problem->arclist[j]->BoundaryMarkerName==agelst[i]->BdryName)
					problem->arclist[j]->MaxSideLength=myMaxSideLength;
		}
	}

	// and perform a quick error check; AGE BCs can't be applied to segments (at least yet)
	for (i=0;i<(int)problem->linelist.size();i++)
	{
		if (problem->linelist[i]->BoundaryMarkerName!="<None>")
		{
			for(j=0;j<(int)agelst.size();j++)
			{

				if (problem->linelist[i]->BoundaryMarkerName==agelst[j]->BdryName)
				{
					WarnMessage("Can't apply Air Gap Element BCs to line segments");
					problem->undo();
					//UnselectAll();
					return -2;
				}
			}
		}
	}

        // want to impose explicit discretization only on
        // the boundary arcs and segments.  After the meshing
        // is done, spacing on boundary segments should be
        // restored to the value that was there before meshing
        // was called, but the arc segments should keep the
        // "new" MaxSideLength--this is used in other places
        // and must always be consistent with the the mesh.


    // Now, do a shitload of checking to make sure that
    // the PBCs haven't been defined by the user
    // in a messed up way.

    // First, search through defined bc's for periodic ones;
//    for(i=0;i<(int)problem->lineproplist.size();i++)
//    {
//        if (problem->lineproplist[i]->isPeriodic())
//        {
//            CPeriodicBoundary pbc;
//            pbc.BdryName=problem->lineproplist[i]->BdryName;
//            pbc.antiPeriodic = problem->lineproplist[i]->isPeriodic(CBoundaryProp::PeriodicityType::AntiPeriodic);
//            pbclst.push_back(pbc.clone());
//        }
//    }

    for(i=0;i<(int)problem->linelist.size();i++)
    {
        for(j=0;j<(int)pbclst.size();j++)
        {
            if (pbclst[j]->BdryName==problem->linelist[i]->BoundaryMarkerName)
            {
                // A pbc or apbc can only be applied to 2 segs
                // at a time.  If it is applied to multiple segs
                // at the same time, flag it and kick it out.
                if (pbclst[j]->nseg==2)
                {
                    char buf[2048];
                    SNPRINTF(buf, sizeof(buf), "An (anti)periodic BC (named \"%s\") is assigned to more than two segments", pbclst[j]->BdryName.c_str ());
                    WarnMessage(buf);
                    problem->undo();  problem->unselectAll();
                    return -1;
                }
                pbclst[j]->seg[pbclst[j]->nseg]=i;
                pbclst[j]->nseg++;
            }
        }
    }

    for(i=0;i<(int)problem->arclist.size();i++)
    {
        for(j=0;j<(int)pbclst.size();j++)
        {
            if (pbclst[j]->BdryName==problem->arclist[i]->BoundaryMarkerName)
            {
                // A pbc or apbc can only be applied to 2 arcs
                // at a time.  If it is applied to multiple arcs
                // at the same time, flag it and kick it out.
                if (pbclst[j]->narc==2)
                {
                    char buf[2048];
                    SNPRINTF(buf, sizeof(buf), "An (anti)periodic BC (named \"%s\") is assigned to more than two arcs", pbclst[j]->BdryName.c_str ());
                    WarnMessage(buf);
                    problem->undo();  problem->unselectAll();
                    return -1;
                }
                pbclst[j]->seg[pbclst[j]->narc]=i;
                pbclst[j]->narc++;
            }
        }
    }

    j=0;
    while(j<(int)pbclst.size())
    {
        // check for a bc that is a mix of arcs and segments.
        // this is an error, and it should get flagged.
        if ((pbclst[j]->nseg>0) && (pbclst[j]->narc>0))
        {
            WarnMessage("Can't mix arcs and segments for (anti)periodic BCs");
            problem->undo();  problem->unselectAll();
            return -1;
        }


        // remove any periodic BC's that aren't actually in play
        if((pbclst[j]->nseg<2) && (pbclst[j]->narc<2)) pbclst.erase(pbclst.begin()+j);
        else j++;
    }

    for(j=0;j<(int)pbclst.size();j++)
    {
        // check to see if adjoining entries are applied
        // to objects of compatible size/shape, and
        // reconcile meshing on the objects.

        // for segments:
        if(pbclst[j]->nseg>0){

            // make sure that lines are pretty much the same length
            if(fabs(problem->lengthOfLine(pbclst[j]->seg[0])
                   -problem->lengthOfLine(pbclst[j]->seg[1]))>1.e-06)
            {
                WarnMessage("(anti)periodic BCs applied to dissimilar segments");
                problem->undo();  problem->unselectAll();
                return -1;
            }

            // make sure that both lines have the same spacing
            double len1,len2,len;
            len1=problem->linelist[pbclst[j]->seg[0]]->MaxSideLength;
            len2=problem->linelist[pbclst[j]->seg[1]]->MaxSideLength;

            if(len1<=0) len1=len2;
            if(len2<=0) len2=len1;
            len=(std::min)(len1,len2);

            problem->linelist[pbclst[j]->seg[0]]->MaxSideLength=len;
            problem->linelist[pbclst[j]->seg[1]]->MaxSideLength=len;
        }

        // for arc segments:
        if(pbclst[j]->narc>0){

            // make sure that arcs are pretty much the
            // same arc length
            if(fabs(problem->arclist[pbclst[j]->seg[0]]->ArcLength
                   -problem->arclist[pbclst[j]->seg[1]]->ArcLength)>1.e-06)
            {
                WarnMessage("(anti)periodic BCs applied to dissimilar arc segments");
                problem->undo();  problem->unselectAll();
                return -1;
            }

            // make sure that both lines have the same spacing
            double len1,len2,len;
            len1=problem->arclist[pbclst[j]->seg[0]]->MaxSideLength;
            len2=problem->arclist[pbclst[j]->seg[1]]->MaxSideLength;

            len=(std::min)(len1,len2);

            problem->arclist[pbclst[j]->seg[0]]->MaxSideLength=len;
            problem->arclist[pbclst[j]->seg[1]]->MaxSideLength=len;
        }
    }

    // write out new poly and write out adjacent
    // boundary nodes in a separate .pbc file.

    // kludge things a bit and use IsSelected to denote
    // whether or not a line or arc has already been processed.
    problem->clearNotationTags();
    problem->unselectAll();
    nodelst.clear();
    nodelst.shrink_to_fit();
    linelst.clear();
    linelst.shrink_to_fit();

    // first, add in existing nodes
    for(const auto &node: problem->nodelist)
        nodelst.push_back(node->clone());

    for(n=0; n<(int)pbclst.size(); n++)
    {
        if (pbclst[n]->nseg != 0) // if this pbc is a line segment...
        {
            int s0,s1;
            CNode node0,node1;

            s0=pbclst[n]->seg[0];
            s1=pbclst[n]->seg[1];
            problem->linelist[s0]->IsSelected=true;
            problem->linelist[s1]->IsSelected=true;

            // make it so that first point on first line
            // maps to first point on second line...
            std::swap(problem->linelist[s1]->n0,problem->linelist[s1]->n1);

            // store number of sub-segments in k
            if (problem->linelist[s0]->MaxSideLength == -1)
            {
                k = 1;
            }
            else{
                a0 = problem->nodelist[problem->linelist[s0]->n0]->CC();
                a1 = problem->nodelist[problem->linelist[s0]->n1]->CC();
                b0 = problem->nodelist[problem->linelist[s1]->n0]->CC();
                b1 = problem->nodelist[problem->linelist[s1]->n1]->CC();
                z = abs(a1-a0);
                k = (int) std::ceil(z/problem->linelist[s0]->MaxSideLength);
            }

            // add segment end points to the list;
            pt.x = problem->linelist[s0]->n0;
            pt.y = problem->linelist[s1]->n0;
            pt.t = pbclst[n]->antiPeriodic;
            ptlst.push_back(pt.clone ());
            pt.x = problem->linelist[s0]->n1;
            pt.y = problem->linelist[s1]->n1;
            pt.t = pbclst[n]->antiPeriodic;
            ptlst.push_back(pt.clone ());

            if (k == 1){
                // catch the case in which the line
                // doesn't get subdivided.
                linelst.push_back(problem->linelist[s0]->clone());
                linelst.push_back(problem->linelist[s1]->clone());
            }
            else{
                segm = *problem->linelist[s0];
                for(j=0; j<k; j++)
                {
                    a2 = a0+(a1-a0)*((double) (j+1))/((double) k);
                    b2 = b0+(b1-b0)*((double) (j+1))/((double) k);
                    node0.x = a2.re; node0.y = a2.im;
                    node1.x = b2.re; node1.y = b2.im;
                    if(j==0){
                        l = nodelst.size();
                        nodelst.push_back(node0.clone());
                        segm.n0 = problem->linelist[s0]->n0;
                        segm.n1 = l;
                        linelst.push_back(segm.clone());
                        pt.x = l;

                        l = nodelst.size();
                        nodelst.push_back(node1.clone());
                        segm.n0 = problem->linelist[s1]->n0;
                        segm.n1 = l;
                        linelst.push_back(segm.clone());
                        pt.y = l;

                        pt.t = pbclst[n]->antiPeriodic;
                        ptlst.push_back(pt.clone());
                    }
                    else if(j==(k-1))
                    {
                        // last subdivision--no ptlst
                        // entry associated with this one.
                        l = nodelst.size()-2;
                        segm.n0 = l;
                        segm.n1 = problem->linelist[s0]->n1;
                        linelst.push_back(segm.clone());

                        l = nodelst.size()-1;
                        segm.n0 = l;
                        segm.n1 = problem->linelist[s1]->n1;
                        linelst.push_back(segm.clone());
                    }
                    else{
                        l = nodelst.size();

                        nodelst.push_back(node0.clone());
                        nodelst.push_back(node1.clone());

                        segm.n0 = l-2;
                        segm.n1 = l;
                        linelst.push_back(segm.clone());

                        segm.n0 = l-1;
                        segm.n1 = l+1;
                        linelst.push_back(segm.clone());

                        pt.x = l;
                        pt.y = l+1;
                        pt.t = pbclst[n]->antiPeriodic;
                        ptlst.push_back(pt.clone());
                    }
                }
            }
        }
        else{  // if this pbc is an arc segment...

            int s0,s1;
            int p0[2],p1[2];
            CNode node0,node1;
            CComplex bgn0,bgn1,c0,c1,d0,d1;
            double r0,r1;

            s0 = pbclst[n]->seg[0];
            s1 = pbclst[n]->seg[1];
            problem->arclist[s0]->IsSelected = true;
            problem->arclist[s1]->IsSelected = true;

            k = (int) ceil(problem->arclist[s0]->ArcLength/problem->arclist[s0]->MaxSideLength);
            segm.BoundaryMarkerName = problem->arclist[s0]->BoundaryMarkerName;
            if (problem->filetype != FileType::MagneticsFile)
                segm.InConductorName=problem->arclist[i]->InConductorName; // not relevant for magnetics
            problem->getCircle(*problem->arclist[s0],c0,r0);
            problem->getCircle(*problem->arclist[s1],c1,r1);

            if (problem->arclist[s0]->NormalDirection ==0){
                bgn0 = problem->nodelist[problem->arclist[s0]->n0]->CC();
                d0 = exp(I*problem->arclist[s0]->ArcLength*PI/(((double) k)*180.));
                p0[0] = problem->arclist[s0]->n0;
                p0[1] = problem->arclist[s0]->n1;
            }
            else{
                bgn0 = problem->nodelist[problem->arclist[s0]->n1]->CC();
                d0 = exp(-I*problem->arclist[s0]->ArcLength*PI/(((double) k)*180.));
                p0[0] = problem->arclist[s0]->n1;
                p0[1] = problem->arclist[s0]->n0;
            }

            if (problem->arclist[s1]->NormalDirection!=0){
                bgn1 = problem->nodelist[problem->arclist[s1]->n0]->CC();
                d1 = exp(I*problem->arclist[s1]->ArcLength*PI/(((double) k)*180.));
                p1[0] = problem->arclist[s1]->n0;
                p1[1] = problem->arclist[s1]->n1;
            }
            else{
                bgn1 = problem->nodelist[problem->arclist[s1]->n1]->CC();
                d1 = exp(-I*problem->arclist[s1]->ArcLength*PI/(((double) k)*180.));
                p1[0] = problem->arclist[s1]->n1;
                p1[1] = problem->arclist[s1]->n0;
            }

            // add arc segment end points to the list;
            pt.x=p0[0]; pt.y=p1[0]; pt.t=pbclst[n]->antiPeriodic;
            ptlst.push_back(pt.clone());
            pt.x=p0[1]; pt.y=p1[1]; pt.t=pbclst[n]->antiPeriodic;
            ptlst.push_back(pt.clone());

            if (k==1){

                // catch the case in which the line
                // doesn't get subdivided.
                segm.n0=p0[0]; segm.n1=p0[1];
                linelst.push_back(segm.clone());
                segm.n0=p1[0]; segm.n1=p1[1];
                linelst.push_back(segm.clone());
            }
            else{
                for(j=0;j<k;j++)
                {
                    bgn0=(bgn0-c0)*d0+c0;
                    node0.x=bgn0.re; node0.y=bgn0.im;

                    bgn1=(bgn1-c1)*d1+c1;
                    node1.x=bgn1.re; node1.y=bgn1.im;

                    if(j==0){
                        l=nodelst.size();
                        nodelst.push_back(node0.clone());
                        segm.n0=p0[0];
                        segm.n1=l;
                        linelst.push_back(segm.clone());
                        pt.x=l;

                        l=nodelst.size();
                        nodelst.push_back(node1.clone());
                        segm.n0=p1[0];
                        segm.n1=l;
                        linelst.push_back(segm.clone());
                        pt.y=l;

                        pt.t=pbclst[n]->antiPeriodic;
                        ptlst.push_back(pt.clone());
                    }
                    else if(j==(k-1))
                    {
                        // last subdivision--no ptlst
                        // entry associated with this one.
                        l=nodelst.size()-2;
                        segm.n0=l;
                        segm.n1=p0[1];
                        linelst.push_back(segm.clone());

                        l=nodelst.size()-1;
                        segm.n0=l;
                        segm.n1=p1[1];
                        linelst.push_back(segm.clone());
                    }
                    else{
                        l=nodelst.size();

                        nodelst.push_back(node0.clone());
                        nodelst.push_back(node1.clone());

                        segm.n0=l-2;
                        segm.n1=l;
                        linelst.push_back(segm.clone());

                        segm.n0=l-1;
                        segm.n1=l+1;
                        linelst.push_back(segm.clone());

                        pt.x=l;
                        pt.y=l+1;
                        pt.t=pbclst[n]->antiPeriodic;
                        ptlst.push_back(pt.clone());
                    }
                }

            }
        }
    }


	// Now, discretize arcs that are part of an AGE
	for(n=0;n<(int)agelst.size();n++)
	{
		std::vector <int> myVector;

		z = (agelst[n]->ro + agelst[n]->ri)/2.;

		for(i=0;i<(int)problem->arclist.size();i++)
		if((problem->arclist[i]->IsSelected==false) && (problem->arclist[i]->BoundaryMarkerName==agelst[n]->BdryName)){
			problem->arclist[i]->IsSelected=true;
			a2.Set(problem->nodelist[problem->arclist[i]->n0]->x,problem->nodelist[problem->arclist[i]->n0]->y);
			k=(int) ceil(problem->arclist[i]->ArcLength/problem->arclist[i]->MaxSideLength);
			segm.BoundaryMarker=problem->arclist[i]->BoundaryMarker;
			problem->GetCircle(*(problem->arclist[i]),c,R);
			a1=exp(I*problem->arclist[i]->ArcLength*PI/(((double) k)*180.));

			// insert the starting node
			if (R>z) // on outer radius
				myVector.push_back(problem->arclist[i]->n0);
			else	// on inner radius
				myVector.insert(myVector.begin(),problem->arclist[i]->n0);

			if(k==1){
				segm.n0=problem->arclist[i]->n0;
				segm.n1=problem->arclist[i]->n1;
				linelst.push_back(segm.clone());
			}
			else for(j=0;j<k;j++)
			{
				a2=(a2-c)*a1+c;
				node.x=a2.re; node.y=a2.im;
				if(j==0){
					l=(int) nodelst.size();
					nodelst.push_back(node.clone());
					segm.n0=problem->arclist[i]->n0;
					segm.n1=l;
					linelst.push_back(segm.clone());

					// insert newly created node
					if (R>z) // on outer radius
						myVector.push_back(l);
					else	// on inner radius
						myVector.insert(myVector.begin(),l);
				}
				else if(j==(k-1))
				{
					l=(int) nodelst.size()-1;
					segm.n0=l;
					segm.n1=problem->arclist[i]->n1;
					linelst.push_back(segm.clone());
				}
				else{
					l=(int) nodelst.size();
					nodelst.push_back(node.clone());
					segm.n0=l-1;
					segm.n1=l;
					linelst.push_back(segm.clone());

					// insert newly created node
					if (R>z) // on outer radius
						myVector.push_back(l);
					else	// on inner radius
						myVector.insert(myVector.begin(),l);

				}
			}
		}

		agelst[n]->nodeNums.resize (myVector.size() + 1);
		agelst[n]->nodeNums.shrink_to_fit ();
		agelst[n]->nodeNums[0] = static_cast<int>(myVector.size());
		std::copy(myVector.begin(), myVector.end(), agelst[n]->nodeNums.begin() + 1);
	}


    // Then, do the rest of the lines and arcs in the
    // "normal" way and write .poly file.

    // discretize input segments
    discretizeInputSegments(*problem, nodelst, linelst, dL, SegmentFilter::OnlyUnselected);

    // discretize input arc segments
    discretizeInputArcSegments(*problem, nodelst, linelst, SegmentFilter::OnlyUnselected);

    // create correct output filename;
    pn = PathName;

//    plyname=pn.Left(pn.ReverseFind('.')) + ".poly";
//
//    // check to see if we are ready to write a datafile;
//
//    if ((fp=fopen(plyname,"wt"))==NULL){
//        WarnMessage("Couldn't write to specified .poly file");
//        Undo();  UnselectAll();
//        return false;
//    }
//
//    // write out node list
//    fprintf(fp,"%i    2    0    1\n",nodelst.size());
//    for(i=0;i<nodelst.size();i++)
//    {
//        for(j=0,t=0;j<nodeproplist.size();j++)
//                if(nodeproplist[j]->PointName==nodelst[i]->BoundaryMarkerName) t=j+2;
//        fprintf(fp,"%i    %.17g    %.17g    %i\n",i,nodelst[i]->x,nodelst[i]->y,t);
//    }
//
//    // write out segment list
//    fprintf(fp,"%i    1\n",linelst.size());
//    for(i=0;i<linelst.size();i++)
//    {
//        for(j=0,t=0;j<lineproplist.size();j++)
//                if(lineproplist[j]->BdryName==linelst[i]->BoundaryMarkerName) t=-(j+2);
//        fprintf(fp,"%i    %i    %i    %i\n",i,linelst[i]->n0,linelst[i]->n1,t);
//    }

//    for(i=0,k=0;i<blocklist.size();i++)
//        if(blocklist[i]->BlockTypeName=="<No Mesh>")
//        {
//            fprintf(fp,"%i    %.17g    %.17g\n",k,blocklist[i]->x,blocklist[i]->y);
//            k++;
//        }
//
//    // write out regional attributes
//    fprintf(fp,"%i\n",blocklist.size()-j);
//
//    for(i=0,k=0;i<blocklist.size();i++)
//        if(blocklist[i]->BlockTypeName!="<No Mesh>")
//        {
//            fprintf(fp,"%i    %.17g    %.17g    ",k,blocklist[i]->x,blocklist[i]->y);
//            fprintf(fp,"%i    ",k+1);
//            if (blocklist[i]->MaxArea>0)
//                fprintf(fp,"%.17g\n",blocklist[i]->MaxArea);
//            else fprintf(fp,"-1\n");
//            k++;
//        }
//    fclose(fp);

    // Make sure to prune out any duplications in the ptlst
    for(k=0;k<(int)ptlst.size();k++) ptlst[k]->sortXY();
    k=0;
    while((k+1) < (int)ptlst.size())
    {
        j=k+1;
        while(j < (int)ptlst.size())
        {
            if((ptlst[k]->x==ptlst[j]->x) && (ptlst[k]->y==ptlst[j]->y))
                ptlst.erase(ptlst.begin()+j);
            else j++;
        }
        k++;
    }

    // used to have a check to eliminate the case where a point
    // and its companion are the same point--actually, this shouldn't
    // be a problem just to let the algorithm deal with this
    // as usual.

    // One last error check--each point must have only one companion point.
    // however, it would be possible to screw up in the definition of the BCs
    // so that this isn't the case.  Look through the points to try and catch
    // this one.
/*
    // let's let this check go away for a minute...

    for(k=0,n=false;(k+1)<ptlst.size();k++)
    {
        for(j=k+1;j<ptlst.size();j++)
        {
            if(ptlst[k]->x==ptlst[j]->x) n=true;
            if(ptlst[k]->y==ptlst[j]->y) n=true;
            if(ptlst[k]->x==ptlst[j]->y) n=true;
            if(ptlst[k]->y==ptlst[j]->x) n=true;
        }
    }
    if (n==true){
        WarnMessage("Nonphysical (anti)periodic boundary assignments");
        Undo();  UnselectAll();
        return false;
    }
*/
    // write out a pbc file containing a list of linked nodes
    plyname = pn.substr(0,pn.find_last_of('.')) + ".pbc";
    if ((fp=fopen(plyname.c_str(),"wt"))==NULL){
        WarnMessage("Couldn't write to specified .pbc file");
        problem->undo();  problem->unselectAll();
        return -1;
    }
    fprintf(fp,"%i\n", (int) ptlst.size());

    for(k=0;k<(int)ptlst.size();k++)
    {
        fprintf(fp,"%i    %i    %i    %i\n",k,ptlst[k]->x,ptlst[k]->y,ptlst[k]->t);
    }

#ifdef DEBUG
    {
        char buf[1048];
        SNPRINTF(buf, sizeof(buf), "writepoly: writing out %i air gap boundaries\n", agelst.size ());
        WarnMessage(buf);
    }
#endif // DEBUG
	fprintf(fp,"%i\n",(int) agelst.size());
	for(k=0;k<(int)agelst.size();k++)
	{
		n = agelst[k]->nodeNums[0] / 2;
		const double dtta = agelst[k]->totalArcLength / n;
		n0 = static_cast<int>(round(360.0 / dtta)); // total elements in a 360deg annular ring;
		n1 = static_cast<int>(round(360.0 / agelst[k]->totalArcLength)); // number of copied segments

		std::vector<CQuadPoint> InnerRing(n0);
		std::vector<CQuadPoint> OuterRing(n0);

		// Should do some consistency checking here;
		//   totalArcLength*n1 should equal 360
		//   no*dtta should equal 360
		//   if antiperiodic, n1 should be an even number
		//   otherwise, throw error message, clean up, and return

		// map each bdry point onto points on the ring;
		int kk;
		for(j=0,kk=0;j<n1;j++)  // do each slice
		{
			if ((agelst[k]->BdryFormat==1) && (j % 2 != 0)) dL=-1; // antiperiodic
			else dL=1;

			a1=exp(I*(j*agelst[k]->totalArcLength+agelst[k]->InnerAngle)*DEGREE);
			a2=exp(I*(j*agelst[k]->totalArcLength+agelst[k]->OuterAngle)*DEGREE);
			for(i=1;i<=n;i++)
			{
				a0=a1*(nodelst[agelst[k]->nodeNums[i]]->CC()-agelst[k]->agc); // position of the shifted mesh node
				z=toDegrees(a0)/dtta;

				InnerRing[kk].n0=agelst[k]->nodeNums[i];
				InnerRing[kk].w0=z;
				InnerRing[kk].w1=dL;

				a0=a2*(nodelst[agelst[k]->nodeNums[i+n]]->CC()-agelst[k]->agc); // position of the shifted mesh node
				z=toDegrees(a0)/dtta;

				OuterRing[kk].n0=agelst[k]->nodeNums[i+n];
				OuterRing[kk].w0=z;
				OuterRing[kk].w1=dL;

				kk++;
			}
		}

		// InnerRing and OuterRing contain a list of boundary nodes, but the aren't yet properly sorted.
		// Sort out the rings based on the angle of the points in the ring
		for(int ii=0;ii<n0;ii++)
		{
			int bDone=1;
			CQuadPoint qq;

			for(int jj=0;jj<(n0-1);jj++)
			{
				if (InnerRing[jj].w0 > InnerRing[jj+1].w0)
				{
					qq=InnerRing[jj];
					InnerRing[jj]=InnerRing[jj+1];
					InnerRing[jj+1]=qq;
					bDone=0;
				}
			}
			if (bDone) break;
		}

		for(int ii=0;ii<n0;ii++)
		{
			int bDone=1;
			CQuadPoint qq;

			for(int jj=0;jj<(n0-1);jj++)
			{
				if (OuterRing[jj].w0 > OuterRing[jj+1].w0)
				{
					qq=OuterRing[jj];
					OuterRing[jj]=OuterRing[jj+1];
					OuterRing[jj+1]=qq;
					bDone=0;
				}
			}
			if (bDone) break;
		}

		// print out AGE definition
		fprintf(fp,"\"%s\"\n",agelst[k]->BdryName.c_str ());
		fprintf(fp,"%i %.17g %.17g %.17g %.17g %.17g %.17g %.17g %i %.17g %.17g\n",
			agelst[k]->BdryFormat,agelst[k]->InnerAngle,agelst[k]->OuterAngle,
			agelst[k]->ri,agelst[k]->ro,agelst[k]->totalArcLength,
			Re(agelst[k]->agc),Im(agelst[k]->agc),n,
			InnerRing[0].w0,OuterRing[0].w0);

		for(i=0;i<=n;i++)
		{
			int p0,p1;

			p1=i; if(p1==n0) p1=0;
			p0=p1-1; if(p0<0) p0=n0+p0;

			// ring points that bracket points in the annulus mesh
			// and their sign, for the purposes of periodicity/antiperiodicity
			fprintf(fp,"%i %g %i %g %i %g %i %g\n",
				InnerRing[p0].n0, InnerRing[p0].w1,
				InnerRing[p1].n0, InnerRing[p1].w1,
				OuterRing[p0].n0, OuterRing[p0].w1,
				OuterRing[p1].n0, OuterRing[p1].w1);
		}

/*
		fprintf(fp,"%s\n",agelst[k]->BdryName);
		fprintf(fp,"%i %.17g %.17g %.17g %.17g %.17g %.17g %.17g %i\n",
			agelst[k]->BdryFormat,agelst[k]->InnerAngle,agelst[k]->OuterAngle,
			agelst[k]->ri,agelst[k]->ro,agelst[k]->totalArcLength,
			Re(agelst[k]->agc),Im(agelst[k]->agc),n);
		for(i=1;i<=n;i++)
			fprintf(fp,"%i %i\n",agelst[k]->quadNode[i],agelst[k]->quadNode[n+i]); */




	}


    fclose(fp);

    // call triangle with -Y flag.
    {
        TriangulateHelper triHelper;
        triHelper.WarnMessage = WarnMessage;
        triHelper.TriMessage = this->TriMessage;

        if (!triHelper.initPointsWithMarkers(nodelst,*problem, PointMarkerInfo::FromProblem))
            return -1;
        if (!triHelper.initSegmentsWithMarkers(linelst,*problem,SegmentMarkerInfo::FromProblem))
            return -1;
        if (!triHelper.initHolesAndRegions(*problem, true, DefaultMeshSize))
            return -1;

        triHelper.setMinAngle(std::min(problem->MinAngle+MINANGLE_BUMP,MINANGLE_MAX));
        triHelper.suppressExteriorSteinerPoints();
        if (writePolyFiles)
        {
            string plyname = PathName.substr(0, PathName.find_last_of('.')) + ".poly";
            triHelper.writePolyFile(plyname, triHelper.triangulateParams());
        }
        int tristatus = triHelper.triangulate(Verbose);
        if (tristatus != 0)
            return tristatus;

        triHelper.writeTriangulationFiles(PathName);
    }

    problem->unselectAll();

    // Now restore boundary segment discretizations that have
    // been mucked up in the process...
    problem->undoLines();
    problem->undoArcs();

    // and save the latest version of the document to make sure
    // any changes to arc discretization get propagated into
    // the solution description....
    //SaveFEMFile(pn);
    problem->saveFEMFile(pn);

    return 0;
}

TriangulateHelper::TriangulateHelper()
    : WarnMessage(&PrintWarningMsg)
    , TriMessage(nullptr)
{
    initialize(in);
#ifdef XFEMM_BUILTIN_TRIANGLE
    initialize(out); // triangle-api writes to input struct
#else
    ctx = triangle_context_create();
#endif
}

TriangulateHelper::~TriangulateHelper()
{
    if (in.pointlist) { free(in.pointlist); }
    if (in.pointattributelist) { free(in.pointattributelist); }
    if (in.pointmarkerlist) { free(in.pointmarkerlist); }
    if (in.regionlist) { free(in.regionlist); }
    if (in.segmentlist) { free(in.segmentlist); }
    if (in.segmentmarkerlist) { free(in.segmentmarkerlist); }
    if (in.holelist) { free(in.holelist); }

#ifdef XFEMM_BUILTIN_TRIANGLE
    if (out.pointlist) { free(out.pointlist); }
    if (out.pointattributelist) { free(out.pointattributelist); }
    if (out.pointmarkerlist) { free(out.pointmarkerlist); }
    if (out.trianglelist) { free(out.trianglelist); }
    if (out.triangleattributelist) { free(out.triangleattributelist); }
    if (out.trianglearealist) { free(out.trianglearealist); }
    if (out.neighborlist) { free(out.neighborlist); }
    if (out.segmentlist) { free(out.segmentlist); }
    if (out.segmentmarkerlist) { free(out.segmentmarkerlist); }
    if (out.edgelist) { free(out.edgelist); }
    if (out.edgemarkerlist) { free(out.edgemarkerlist); }
#else
    triangle_context_destroy(ctx);
#endif
}

bool TriangulateHelper::initPointsWithMarkers(const TriangulateHelper::nodelist_t &nodelst, const FemmProblem &problem, PointMarkerInfo info)
{
    // calling this method on an already initialized object would leak memory
    if (in.numberofpoints!=0)
    {
        WarnMessage("initPointsWithMarkers called twice!\n");
        return false;
    }

    in.numberofpoints = nodelst.size();

    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    if (!in.pointlist) {
        WarnMessage("Point list for triangulation is null!\n");
        return false;
    }

    for(int i=0; i < in.numberofpoints; i++)
    {
        in.pointlist[2*i] = nodelst[i]->x;
        in.pointlist[2*i+1] = nodelst[i]->y;
    }

    // Initialise the pointmarkerlist
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    if (!in.pointmarkerlist) {
        WarnMessage("Point marker list for triangulation is null!\n");
        return false;
    }

    // write out node marker list
    for(int i=0; i<in.numberofpoints; i++)
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

        in.pointmarkerlist[i] = t;
    }
    return true;
}

bool TriangulateHelper::initSegmentsWithMarkers(const TriangulateHelper::linelist_t &linelst, const FemmProblem &problem, SegmentMarkerInfo info)
{
    // calling this method on an already initialized object would leak memory
    if (in.numberofsegments!=0)
    {
        WarnMessage("initSegmentsWithMarkers called twice!\n");
        return false;
    }

    in.numberofsegments = linelst.size();

    // Initialise the segmentlist
    in.segmentlist = (int *) malloc(2 * in.numberofsegments * sizeof(int));
    if (!in.segmentlist) {
        WarnMessage("Segment list for triangulation is null!\n");
        return false;
    }
    // Initialise the segmentmarkerlist
    in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
    if (!in.segmentmarkerlist) {
        WarnMessage("Segment marker list for triangulation is null!\n");
        return false;
    }

    // build the segmentlist
    for(int i=0; i<in.numberofsegments; i++)
    {
        in.segmentlist[2*i] = linelst[i]->n0;
        in.segmentlist[2*i+1] = linelst[i]->n1;
    }

    // now build the segment marker list
    // construct the segment list
    for(int i=0; i<in.numberofsegments; i++)
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
        in.segmentmarkerlist[i] = t;
    }
    return true;
}

bool TriangulateHelper::initHolesAndRegions(const FemmProblem &problem, bool forceMaxMeshArea, double defaultMeshSize)
{
    // calling this method on an already initialized object would leak memory
    if (in.numberofholes!=0)
    {
        WarnMessage("initHolesAndRegions called twice!\n");
        return false;
    }

    in.numberofholes = problem.countHoles();
    if(in.numberofholes > 0)
    {
        in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));
        if (!in.holelist) {
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
                in.holelist[k++] = label->x;
                in.holelist[k++] = label->y;
            }
        }
    }

    in.numberofregions = problem.labellist.size() - in.numberofholes;
    in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
    if (!in.regionlist) {
        WarnMessage("Region list for triangulation is null!\n");
        return false;
    }

    int j=0;
    int k=0;
    for(const auto & label: problem.labellist)
    {
        if(!label->isHole())
        {
            in.regionlist[j] = label->x;
            in.regionlist[j+1] = label->y;
            in.regionlist[j+2] = k + 1; // Regional attribute (for whole mesh).
#ifdef DEBUG
            {
                char buf[1028];
                SNPRINTF (buf, sizeof(buf), "Adding region (at (%g,%g)) with attribute value %g to triangle input region list\n",
                          label->x, label->y, in.regionlist[j+2]);
                WarnMessage(buf);
            }
#endif // DEBUG
            // Note(ZaJ): this is the code that was used in the periodic bc triangulation:
            //  if (label->MaxArea>0 && (label->MaxArea<defaultMeshSize))
            //      in.regionlist[j+3] = label->MaxArea;  // Area constraint
            //  else
            //      in.regionlist[j+3] = defaultMeshSize;
            // ... which is equivalent to the code below (if forceMaxMeshArea is true).
            // ... the code below is a copy of the nonperiodic case (if forceMaxMeshArea is set to problem->DoForceMaxMeshArea)

            // Area constraint
            if (label->MaxArea <= 0)
            {
                // if no mesh size has been specified use the default
                in.regionlist[j+3] = defaultMeshSize;
            }
            else if ((label->MaxArea > defaultMeshSize) && (forceMaxMeshArea))
            {
                // if the user has specied that FEMM should choose an
                // upper mesh size limit, regardles of their choice,
                // and their choice is less than that limit, change it
                // to that limit
                in.regionlist[j+3] = defaultMeshSize;
            }
            else
            {
                // Use the user's choice of mesh size
                in.regionlist[j+3] = label->MaxArea;
            }

            j += 4;
            k++;
        }
    }
    return true;
}

int TriangulateHelper::triangulate(bool verbose)
{
    std::string triArgs = triangulateParams(verbose);
    // this is a mess, but building the string with std::string is more flexible than sprintf
    // (and the triangulate api is ancient)
    char cmdline[512];
    sprintf(cmdline, "%s",triArgs.c_str());

#ifdef XFEMM_BUILTIN_TRIANGLE
    int tristatus = ::triangulate(cmdline, &in, &out, (struct triangulateio *) nullptr, this->TriMessage);
    if (tristatus!=0)
    {
        std::string msg = "Call to triangulate failed with status code: " + to_string(tristatus) +"\n";
        WarnMessage(msg.c_str());
        return tristatus;
    }
#else
    // parse options
    int tristatus = triangle_context_options(ctx, cmdline);
    if (tristatus != TRI_OK)
    {
        WarnMessage("Invalid option string for triangle!\n");
        return tristatus;
    }
    // Triangulate the polygon.
    tristatus = triangle_mesh_create(ctx, &in);
    if (tristatus != TRI_OK)
    {
        std::string msg = "Call to triangulate failed with status code: " + to_string(tristatus) +"\n";
        WarnMessage(msg.c_str());
        return tristatus;
    }
#endif
    return 0;
}

string TriangulateHelper::triangulateParams(bool verbose) const
{
    // An explaination of the input parameters used for Triangle
    //
    // -p Triangulates a Planar Straight Line Graph, i.e. list of segments.
    // -P Suppresses the output .poly file.
    // -q Quality mesh generation with no angles smaller than specified in the following number
    // -e Outputs a list of edges of the triangulation.
    // -A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
    // -a Imposes a maximum triangle area constraint.
    // -z Numbers all items starting from zero (rather than one)
    // -I Suppresses mesh iteration numbers
    // -j prevents duplicated input vertices, or vertices `eaten' by holes,
    //    from appearing in the output .node file.  Thus, if two input vertices
    //    have exactly the same coordinates, only the first appears in the
    //    output.
    // -Y Suppresses the creation of Steiner points on the exterior boundary.
    //
    // See http://www.cs.cmu.edu/~quake/triangle.switch.html for more info
    std::string triArgs = "-pPq" + to_string(m_minAngle) + "eAaz" + (verbose?"":"Q") + "I";
    if (m_suppressUnusedVertices)
        triArgs += "j";
    if (m_suppressExteriorSteinerPoints)
        triArgs += "Y";

    return triArgs;
}

bool TriangulateHelper::writePolyFile(string filename, std::string comment) const
{
    std::ofstream polyFile (filename);
    // set floating point precision once for the whole stream
    polyFile << std::setprecision(17);
    // when filling to a width, adjust to the left
    polyFile.setf(std::ios::left);

    polyFile << in.numberofpoints << "\t2\t0\t1\n";
    for (int i=0; i < in.numberofpoints; i++)
    {
        polyFile << i << "\t" << in.pointlist[2*i] << "\t" << in.pointlist[2*i+1] << "\t" << in.pointmarkerlist[i] << "\n";
    }

    polyFile << in.numberofsegments << "\t1\n";
    for (int i=0; i < in.numberofsegments; i++)
    {
        polyFile << i << "\t" << in.segmentlist[2*i] << "\t" << in.segmentlist[2*i+1] << "\t" << in.segmentmarkerlist[i] <<"\n";
    }

    polyFile << in.numberofholes << "\n";
    for (int i=0; i < in.numberofholes; i++)
    {
        polyFile << i << "\t" << in.holelist[2*i] << "\t" << in.holelist[2*i+1] << "\n";
    }

    polyFile << in.numberofregions << "\n";
    for (int i=0; i < in.numberofregions; i++)
    {
        int j=4*i;
        polyFile << i << "\t"
                 << in.regionlist[j] << "\t"
                 << in.regionlist[j+1] << "\t"
                 << in.regionlist[j+2] << "\t"
                 << in.regionlist[j+3] << "\n";
    }

    polyFile << "# " << comment << "\n";
    return true;
}

void TriangulateHelper::setMinAngle(double value)
{
    m_minAngle = value;
}

void TriangulateHelper::suppressExteriorSteinerPoints()
{
    m_suppressExteriorSteinerPoints = true;
}

void TriangulateHelper::suppressUnusedVertices()
{
    m_suppressUnusedVertices = true;
}


