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

#ifndef FMESHER_H
#define FMESHER_H

#include "nosebl.h"
#include "femmcomplex.h"
#include "IntPoint.h"

#include "CArcSegment.h"
#include "CBlockLabel.h"
#include "CBoundaryProp.h"
#include "CCircuit.h"
#include "CNode.h"
#include "CPointProp.h"
#include "CSegment.h"
#include "femmenums.h"
#include "FemmProblem.h"
#include "mesh/SolverMesh.h"

#include <memory>
#include <vector>
#include <string>

#ifndef LineFraction
#define LineFraction 500.0
#endif

// Default mesh size is the diagonal of the geometry's
// bounding box divided by BoundingBoxFraction
#ifndef BoundingBoxFraction
#define BoundingBoxFraction 100.0
#endif

#ifndef _WIN32
#define _strnicmp strncasecmp
#ifndef SNPRINTF
  #define SNPRINTF std::snprintf
#endif
#else
#ifndef SNPRINTF
  #define SNPRINTF _snprintf
#endif
#endif

namespace fmesher
{

// FMesher Class

class FMesher
{

protected:


// Attributes
public:

    enum loaderrors { F_FILE_OK, F_FILE_UNKNOWN_TYPE, F_FILE_NOT_OPENED, F_FILE_MALFORMED};
    enum class Backend { Triangle, Tangle };

    explicit FMesher();
    explicit FMesher(std::string);
    explicit FMesher(std::shared_ptr<femm::FemmProblem> p);
    virtual ~FMesher() {}

    std::shared_ptr<femm::FemmProblem> problem;
    bool Verbose = true;
    bool writePolyFiles = false; ///< write .poly files when calling triangle
    Backend backend = defaultBackend();

	std::string BinDir;

	// vectors containing the mesh information
    std::vector< std::unique_ptr<femm::IntPoint> > meshline;
    std::vector< std::unique_ptr<femm::IntPoint> > greymeshline;
    std::vector< std::unique_ptr<femm::CNode> >	meshnode;

    // used to echo start of input file to output
    std::vector< std::string > probdescstrings;

// Operations
public:

    static femm::FileType GetFileType(std::string PathName);
    static constexpr Backend defaultBackend()
    {
#ifdef XFEMM_MESHER_BACKEND_TANGLE
        return Backend::Tangle;
#else
        return Backend::Triangle;
#endif
    }
    static const char *backendName(Backend selected);

    /**
     * @brief Calculate length used to kludge fine meshing near input node points
     * @return the average line length
     * \internal
     * \note This method does not exist in FEMM42. It contains code that was duplicated in the two triangulation methods.
     * \endinternal
     */
    double averageLineLength() const;
    bool SaveFEMFile(std::string PathName); ///< \deprecated

    //void downstr(char *s);

    // Core functions
	/**
	 * @brief LoadMesh
	 * @param PathName
	 * @return
	 * \internal
	 * ### FEMM sources (identical implementation):
	 * - \femm42{beladrawDoc.cpp,CbeladrawDoc::LoadMesh()}
	 * - \femm42{CDRAWDOC.cpp,CcdrawDoc::LoadMesh()}
	 * - \femm42{FemmeDoc.cpp,CFemmeDoc::LoadMesh()}
	 * - \femm42{HDRAWDOC.cpp,ChdrawDoc::LoadMesh()}
	 * \endinternal
	 */
	bool LoadMesh(std::string PathName);
	int DoNonPeriodicBCTriangulation(std::string PathName);
	int DoPeriodicBCTriangulation(std::string PathName);
	bool HasPeriodicBC();

    // pointer to function to call when issuing warning messages
    int (*WarnMessage)(const char*, ...);

    // pointer to function to use for triangle to issue warning messages
    int (*TriMessage)(const char * format, ...);

private:
    friend class TriangleMesherBackend;
    int doNonPeriodicTriangleWorkflow(std::string PathName, femm::mesh::SolverMesh &mesh);
    int doPeriodicTriangleWorkflow(std::string PathName, femm::mesh::SolverMesh &mesh);

    virtual bool Initialize(femm::FileType t);
	void addFileStr (char * q);
};

/**
 * @brief The SegmentFilter enum is used by the discretization functions.
 * @see discretizeInputSegments
 * @see discretizeInputArcSegments
 */
enum class SegmentFilter {
    AllSegments ///< Use all segments
    , OnlyUnselected ///< Only use segments that are not selected
};

/**
 * @brief Figure out a good default mesh size for block labels where mesh size isn't explicitly specified.
 * @param nodelst
 * @param doSmartMesh
 * @return a suitable mesh size, or -1 if nodelst is empty.
 *
 * \internal
 * This function contains code originally duplicated in both DoPeriodicBCTriangulation and DoNonPeriodicBCTriangulation.
 * \endinternal
 */
double defaultMeshSizeHeuristics(const std::vector<std::unique_ptr<femm::CNode>> &nodelst, bool doSmartMesh);

/**
 * @brief Create a copy of the problem's segment list where the segment length is bounded by their MaxSideLength.
 * All segments in the problem's linelist are copied into \p linelst, and additional segments are added as needed.
 *
 * Where MaxSideLength is not specified (i.e. -1), the segment is either copied verbatim (if smart meshing is disabled),
 * or further nodes are inserted at distance \p dL from the segment's existing nodes and the segment is split.
 *
 * @param problem
 * @param nodelst
 * @param linelst
 * @param dL distance to corner for smart meshing
 *
 * \internal
 * This function contains code originally duplicated in both DoPeriodicBCTriangulation and DoNonPeriodicBCTriangulation.
 * \endinternal
 */
void discretizeInputSegments(
        const femm::FemmProblem &problem,
        std::vector <std::unique_ptr<femm::CNode>> &nodelst,
        std::vector <std::unique_ptr<femm::CSegment>> &linelst,
        double dL,
        SegmentFilter filter = SegmentFilter::AllSegments
        );

/**
 * @brief Create a copy of the problem's segment list where the segment length is bounded by their MaxSideLength.
 * All segments in the problem's linelist are copied into \p linelst, and additional segments are added as needed.
 *
 * @param problem
 * @param nodelst
 * @param linelst
 *
 * \internal
 * This function contains code originally duplicated in both DoPeriodicBCTriangulation and DoNonPeriodicBCTriangulation.
 * \endinternal
 */
void discretizeInputArcSegments(
                const femm::FemmProblem &problem,
                std::vector <std::unique_ptr<femm::CNode>> &nodelst,
                std::vector <std::unique_ptr<femm::CSegment>> &linelst,
                SegmentFilter filter = SegmentFilter::AllSegments
                );

/**
 * @brief triangleVersionString exposes information about the triangle version that is used.
 * @return a string describing the version of triangle in use.
 */
std::string triangleVersionString();

} // namespace femm

#endif
