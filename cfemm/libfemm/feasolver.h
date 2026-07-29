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

   Date Modified: 2017
   By: Richard Crozier
       Johannes Zarl-Zierl
   Contact:
        richard.crozier@yahoo.co.uk
       johannes@zarl-zierl.at

   Contributions by Johannes Zarl-Zierl were funded by
    Linz Center of Mechatronics GmbH (LCM)
*/

// feasolver.h : interface of the generic feasolver class
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FEASOLVER_H
#define FEASOLVER_H

#include "femmenums.h"
#include "fparse.h"
#include "spars.h"
#include "CAirGapElement.h"
#include "CBoundaryProp.h"
#include "CCommonPoint.h"
#include "CNode.h"

#include <string>
#include <vector>

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

enum LoadMeshErr
{
    NOERROR,
    BADFEMFILE,
    BADNODEFILE,
    BADPBCFILE,
    BADELEMENTFILE,
    BADEDGEFILE,
    MISSINGMATPROPS,
    ELMLABELTOOBIG
};

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
class FEASolver
{

// Attributes
public:
    using FEASolver_type = FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>;
    using PointProp_type = PointPropT;
    using BoundaryProp_type = BoundaryPropT;
    using BlockProp_type = BlockPropT;
    using CircuitProp_type = CircuitPropT;
    using BlockLabel_type = BlockLabelT;
    using MeshElement_type = MeshElementT;

    FEASolver();
    virtual ~FEASolver();

    // General problem attributes
    double FileFormat; ///< \brief format version of the file
    double Precision;  ///< \brief Computing precision within FEMM
    double MinAngle;   ///< \brief angle restriction for triangulation [deg]
    double Depth;      ///< \brief typical length in z-direction [lfac]
    femm::LengthUnit  LengthUnits;  ///< \brief Unit for lengths. Also referred to as \em lfac.
    femm::CoordsType  Coords;  ///< \brief definition of the coordinate system
    femm::ProblemType ProblemType; ///< \brief The 2D problem is either planar or axisymmetric
    // axisymmetric external region parameters
    double extZo;  ///< \brief center of exterior [lfac], only valid for axisymmetric problems
    double extRo;  ///< \brief radius of exterior [lfac], only valid for axisymmetric problems
    double extRi;  ///< \brief radius of interior [lfac], only valid for axisymmetric problems
    std::string comment; ///< \brief Problem description

    int		ACSolver;
    bool    DoForceMaxMeshArea;
    bool    DoSmartMesh;
    bool    bMultiplyDefinedLabels;


    // CArrays containing the mesh information
    int	BandWidth;
    std::vector<MeshElementT> meshele;

    int NumNodes;
    int NumEls;

    // lists of properties
    int NumBlockProps;
    int NumPBCs;
    int NumLineProps;
    int NumPointProps;
    int NumCircProps;
    int NumBlockLabels;
    int NumAirGapElems;

    std::vector<femm::CCommonPoint>	pbclist;

    // string to hold the location of the files
    std::string PathName;

    int PrevType; ///< \brief flag indicating type of previous solution, 0 for None, 1 for Incremental or 2 for Frozen \verbatim[prevtype]\endverbatim
    std::string previousSolutionFile; ///< \brief name of a previous solution file for hsolver and fsolver incremental permeability \verbatim[prevsoln]\endverbatim

    std::vector< PointPropT > nodeproplist;
    std::vector< BoundaryPropT > lineproplist;
    std::vector< BlockPropT > blockproplist;
    std::vector< CircuitPropT > circproplist;
    std::vector< BlockLabelT > labellist;
    std::vector< femmsolver::CAirGapElement > agelist;
    std::vector< femm::CNode > nodes;
// Operations
public:

    virtual LoadMeshErr LoadMesh(bool deleteFiles=true) = 0;
    /**
     * @brief Solve the problem.
     * This is basically a convenience function that encompasses all steps as would be normally taken when running the solver with default settings.
     * Note: You need to set the PathName beforehand.
     * @return \c true on success, \c false on error.
     */
    virtual bool runSolver(bool verbose=false) = 0;

    /**
     * @brief Return an error string for the mesh error enum.
     * @param err
     * @return an error string, or an empty string if \c err is \c NOERROR
     */
    static std::string getErrorString(LoadMeshErr err);

    int Cuthill(bool deleteFiles=true);
    /** Apply Cuthill--McKee ordering from connectivity already held in memory. */
    int Cuthill(const std::vector<std::pair<std::size_t, std::size_t>> &edges);
    int SortElements();

    // pointer to function to call when issuing warning messages
    int (*WarnMessage)(const char*, ...);
    int (*PrintMessage)(const char*, ...);

    virtual void CleanUp();

protected:

    bool meshLoadedFromPrevSolution;

protected:
    /**
     * @brief LoadProblemFile
     * @param file
     * @return
     * \internal
     * ### FEMM source:
     * - \femm42{fkn/femmedoccore.cpp:CFemmeDocCore::OnOpenDocument()}
     * - \femm42{hsolv/hsolvdoc.cpp,Chsolvdoc::OnOpenDocument()}
     * - \femm42{belasolv/femmedoccore.cpp,CFemmeDocCore::OnOpenDocument()}
     * \endinternal
     */
    bool LoadProblemFile(std::string &file);
    /**
     * @brief handleToken is called by LoadProblemFile() when a token is encountered that it can not handle.
     *
     * Classes subclassing FeaSolver can override this method to handle problem-specific problem file entries.
     * If a token is understood, this method should read its text from the input stream and then return \c true.
     * If a token is not understood, the method must not change the position in the input stream, and must return \c false.
     *
     * \note The input stream contains the remaining contents of the line that \c token was found on.
     * Therefore, this method can never read beyond the current line.
     *
     * @param token the token in question
     * @param input input stream for the current line
     * @param err output stream for error messages
     * @return \c false, if the token is not handled. \c true, if it is handled
     */
    virtual bool handleToken(const std::string &token, std::istream &input, std::ostream &err);

private:

    virtual void SortNodes (std::vector<int> newnum) = 0;

};

/////////////////////////////////////////////////////////////////////////////

double Power(double x, int n);

#endif
