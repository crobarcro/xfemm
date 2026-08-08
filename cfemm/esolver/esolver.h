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

   Date Modified: 2014 - 03 - 21
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

// esolver.h : interface of the ESolver class

#ifndef ESOLVER_H
#define ESOLVER_H

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

class ESolver : public FEASolver<
        femm::CSPointProp
        , femm::CSBoundaryProp
        , femm::CSMaterialProp
        , femm::CSCircuit
        , femm::CSBlockLabel
        , femmsolver::CElement
        >
{

// Attributes
public:

	ESolver();
	~ESolver();


    // mesh information
    femm::CNode *meshnode;

// Operations
public:

    LoadMeshErr LoadMesh(bool deleteFiles=true) override;
    bool LoadProblemFile();
    double ChargeOnConductor(int conductor, femm::LinearSystemBackend<double> &L);
    int WriteResults(femm::LinearSystemBackend<double> &L);
    int AnalyzeProblem(femm::LinearSystemBackend<double> &L);
    int (*WarnMessage)(const char*, ...);

    virtual bool runSolver(bool verbose=false) override;
private:

    void MsgBox(const char* message);
    void CleanUp() override;

    // override parent class virtual method
    void SortNodes (std::vector<int> newnum) override;

    virtual bool handleToken(const std::string &, std::istream &, std::ostream &) override;

};

#endif
