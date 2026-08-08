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
   By:  Emoke Szelitzky
        Tibor Szelitzky
        Richard Crozier
        Johannes Zarl-Zierl
   Contact:
        szelitzkye@gmail.com
        sztibi82@yahoo.com
        richard.crozier@yahoo.co.uk
        johannes.zarl-zierl@jku.at

   Contributions by Johannes Zarl-Zierl were funded by
   Linz Center of Mechatronics GmbH (LCM)
*/

// hsolver.h : interface of the HSolver class

#ifndef HSOLVER_H
#define HSOLVER_H

#include "feasolver.h"
#include "linsolve/LinearSystemBackend.h"
#include "CBlockLabel.h"
#include "CBoundaryProp.h"
#include "CCircuit.h"
#include "CElement.h"
#include "CNode.h"
#include "CMaterialProp.h"
#include "CPointProp.h"

#include <string>

class HSolver : public FEASolver<
        femm::CHPointProp
        , femm::CHBoundaryProp
        , femm::CHMaterialProp
        , femm::CHConductor
        , femm::CHBlockLabel
        , femmsolver::CElement
        >
{

// Attributes
public:

	HSolver();
	~HSolver();


    // General problem attributes
    double	dT; ///< \brief delta T used by hsolver \verbatim[dT]\endverbatim

    // mesh information
    femm::CNode *meshnode;

    // Vector containing previous solution for time-transient analysis, only valid when dT!=0
	double *Tprev;


// Operations
public:

    LoadMeshErr LoadMesh(bool deleteFiles=true) override;
    int LoadPrev();
    bool LoadProblemFile();
    double ChargeOnConductor(int OnConductor, femm::LinearSystemBackend<double> &L);
	int WriteResults(femm::LinearSystemBackend<double> &L);
    int AnalyzeProblem(femm::LinearSystemBackend<double> &L);
    int (*WarnMessage)(const char*, ...);

    virtual bool runSolver(bool verbose=false) override;
private:

    void MsgBox(const char* message);
    void CleanUp() override;

    // override parent class virtual method
    void SortNodes (std::vector<int> newnum) override;

    virtual bool handleToken(const std::string &token, std::istream &input, std::ostream &err) override;

};

#endif
