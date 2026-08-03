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

// fsolver.cpp : implementation of the FSolver class
//

#include "femmcomplex.h"
#include "spars.h"
#include "fparse.h"
#include "feasolver.h"
#include "stringTools.h"

#include <assert.h>
#include <ctype.h>
#include <fstream>
#include <ios>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#ifndef _MSC_VER
#define _strnicmp strncasecmp
#endif

#define STOP_ON_UNKNOWN_TOKEN

using namespace std;
using namespace femm;

/////////////////////////////////////////////////////////////////////////////
// FEASolver construction/destruction

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::FEASolver()
    : FileFormat(-1)
    , Precision(1.e-08)
    , MinAngle(0.)
    , Depth(-1)
    , LengthUnits(LengthInches)
    , Coords(CART)
    , ProblemType(PLANAR)
    , extZo(0.0)
    , extRo(0.0)
    , extRi(0.0)
    , comment()
    , ACSolver(0)
    , DoForceMaxMeshArea(false)
    , DoSmartMesh(true)
    , bMultiplyDefinedLabels(false)
    , BandWidth(0)
    , meshele()
    , NumNodes(0)
    , NumEls(0)
    , NumBlockProps(0)
    , NumPBCs(0)
    , NumLineProps(0)
    , NumPointProps(0)
    , NumCircProps(0)
    , NumBlockLabels(0)
    , NumAirGapElems(0)
    , pbclist()
    , PathName()
    , PrevType(0)
    , nodeproplist()
    , lineproplist()
    , blockproplist()
    , circproplist()
    , labellist()
    , agelist()
    , nodes()
{

    // initialise the warning message box function pointer to
    // point to the PrintWarningMsg function
    WarnMessage = &PrintWarningMsg;
    PrintMessage = &PrintWarningMsg;

    meshLoadedFromPrevSolution = false;
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::~FEASolver()
{
}
template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
void FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::CleanUp()
{
    FileFormat = -1;
    Precision = 1.e-08;
    MinAngle = 0.;
    Depth = -1;
    LengthUnits = LengthInches;
    Coords = CART;
    ProblemType = PLANAR;
    extZo = 0.0;
    extRo = 0.0;
    extRi = 0.0;
    comment.clear();
    ACSolver = 0;
    DoForceMaxMeshArea = false;
    DoSmartMesh = true;
    bMultiplyDefinedLabels = false;
    BandWidth = 0;
    meshele.clear();
    NumNodes = 0;
    NumEls = 0;
    NumBlockProps = 0;
    NumPBCs = 0;
    NumLineProps = 0;
    NumPointProps = 0;
    NumCircProps = 0;
    NumBlockLabels = 0;
    pbclist.clear();
    agelist.clear();
    // *do not* remove the PathName
    //PathName.clear();
    PrevType = 0;
    nodeproplist.clear();
    lineproplist.clear();
    blockproplist.clear();
    circproplist.clear();
    labellist.clear();
    nodes.clear();
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
bool FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::LoadProblemFile(std::string &file)
{
    std::ifstream input;
    std::stringstream err;
    err >> noskipws; // don't discard whitespace from message stream

    WarnMessage ("FEASolver::LoadProblemFile\n");

    input.open(file.c_str(), std::ifstream::in);
    if (!input.is_open())
    {
        err << "Couldn't read from specified .fem file: "
            << file.c_str()
            << "\n";
        WarnMessage(err.str().c_str());
        return false;
    }

    // define some defaults
    CleanUp();

    // parse the file
#ifdef DEBUG_PARSER
    std::cout << "feasolver starting parsing file " << file << std::endl;
#endif // DEBUG_PARSER

    string line;
    string token;
    bool success = true;
    while (input.good() && success)
    {
        if (input.eof())
            break;

        std::getline(input, line);

        trim(line);
        if (line.empty())
        {
#ifdef DEBUG_PARSER
            std::cout << "current line is empty, moving on" << std::endl;
#endif // DEBUG_PARSER
            continue;
        }

#ifdef DEBUG_PARSER
        std::cout << "current line reads:" << std::endl
                  << line << std::endl;
#endif // DEBUG_PARSER

        std::istringstream lineStream(line);
        nextToken(lineStream, &token);

#ifdef DEBUG_PARSER
        std::cout << "feasolver token is: " << token << std::endl;
#endif // DEBUG_PARSER

        if( token == "[format]")
        {
            success &= expectChar(lineStream, '=',err);
            success &= parseValue(lineStream, FileFormat, err);
            continue;
        }

        // Precision
        if( token == "[precision]")
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, Precision, err);
            continue;
        }

        if( token == "[minangle]")
        {
            success &= expectChar(lineStream, '=',err);
            success &= parseValue(lineStream, MinAngle, err);
            continue;
        }

        // Depth for 2D planar problems;
        if( token == "[depth]")
        {
            success &= expectChar(lineStream, '=',err);
            success &= parseValue(lineStream, Depth, err);
            continue;
        }

        // Units of length used by the problem
        if( token == "[lengthunits]")
        {
            success &= expectChar(lineStream, '=', err);
            nextToken(lineStream, &token);

            if( token == "inches" ) LengthUnits=LengthInches;
            else if( token == "millimeters" ) LengthUnits=LengthMillimeters;
            else if( token == "centimeters" ) LengthUnits=LengthCentimeters;
            else if( token == "mils" ) LengthUnits=LengthMils;
            else if( token == "microns" ) LengthUnits=LengthMicrometers;
            else if( token == "meters" ) LengthUnits=LengthMeters;
            continue;
        }

        // Coordinates (cartesian or polar)
        if( token == "[coordinates]" )
        {
            success &= expectChar(lineStream, '=', err);
            nextToken(lineStream, &token);

            if ( token == "cartesian" ) Coords=CART;
            if ( token == "polar" ) Coords=POLAR;
            continue;
        }

        // Problem Type (planar or axisymmetric)
        if( token == "[problemtype]" )
        {
            success &= expectChar(lineStream, '=', err);
            nextToken(lineStream, &token);

            if( token == "planar" ) ProblemType=PLANAR;
            if( token == "axisymmetric" ) ProblemType=AXISYMMETRIC;
            continue;
        }

        // properties for axisymmetric external region
        if( token == "[extzo]" )
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, extZo, err);
            continue;
        }

        if( token == "[extro]" )
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, extRo, err);
            continue;
        }

        if( token == "[extri]" )
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, extRi, err);
            continue;
        }


        if( token == "[comment]" )
        {
            success &= expectChar(lineStream, '=', err);
            parseString(lineStream, &comment, err);
            continue;
        }

        // AC Solver Type
        if( token == "[acsolver]")
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, ACSolver, err);
            continue;
        }

		// Previous solution type
		if( token == "[prevtype]" )
        {
			success &= expectChar(lineStream, '=', err);
			success &= parseValue(lineStream, PrevType, err);
			continue;
		}

        if( token == "[prevsoln]" )
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseString(lineStream, &previousSolutionFile, err);
            continue;
        }

        // Option to force use of default max mesh, overriding
        // user choice
        if( token == "[forcemaxmesh]")
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, DoForceMaxMeshArea, err);
            continue;
        }

        // Option to use smart meshing
        //std::cout << "checking for dosmartmesh" << std::endl;
        if( token == "[dosmartmesh]" )
        {
            //std::cout << "found dosmartmesh" << std::endl;
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, DoSmartMesh, err);
            continue;
        }

        // Point Properties
        if( token == "[pointprops]" )
        {
            int k;
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, k, err);
            if (k>0) nodeproplist.reserve(k);
            while (input.good() && NumPointProps < k)
            {
                PointPropT next = PointPropT::fromStream(input, err);
                nodeproplist.push_back(next);
                NumPointProps++;
            }
            // message will be printed after parsing is done
            if (NumPointProps != k)
            {
                err << "Expected "<<k<<" PointProps, but got " << NumPointProps << "\n";
                break; // stop parsing
            }
            continue;
        }


        // Boundary Properties;
        if( token == "[bdryprops]" )
        {
            success &= expectChar(lineStream, '=', err);
            int k;
            success &= parseValue(lineStream, k, err);
            if (k>0) lineproplist.reserve(k);

            while (input.good() && NumLineProps < k)
            {
                BoundaryPropT next = BoundaryPropT::fromStream(input, err);
                lineproplist.push_back(next);
                NumLineProps++;
            }
            // message will be printed after parsing is done
            if (NumLineProps != k)
            {
                err << "Expected "<<k<<" BoundaryProps, but got " << NumLineProps << "\n";
                break; // stop parsing
            }
            continue;
        }


        // Block Properties;
        if( token == "[blockprops]" )
        {
            success &= expectChar(lineStream, '=', err);
            int k;
            success &= parseValue(lineStream, k, err);
            if (k>0) blockproplist.reserve(k);

            while (input.good() && NumBlockProps < k)
            {
                BlockPropT next = BlockPropT::fromStream(input, err);
                if (input.good()) {
                    blockproplist.push_back(next);
                    NumBlockProps++;
                }
            }
            // message will be printed after parsing is done
            if (NumBlockProps != k)
            {
                err << "Expected "<<k<<" BlockProps, but got " << NumBlockProps << "\n";
                break; // stop parsing
            }
            continue;
        }

        // Circuit Properties
        if( token == "[circuitprops]" || token == "[conductorprops]")
        {
            success &= expectChar(lineStream, '=', err);
            int k;
            success &= parseValue(lineStream, k, err);
            if(k>0) circproplist.reserve(k);

            while (input.good() && NumCircProps < k)
            {
                CircuitPropT next = CircuitPropT::fromStream(input, err);
                circproplist.push_back(next);
                NumCircProps++;
            }
            // message will be printed after parsing is done
            if (NumCircProps != k)
            {
                err << "Expected "<<k<<" CircuitProps, but got " << NumCircProps << "\n";
                break; // stop parsing
            }
            continue;
        }

        // read in regional attributes
        if(token == "[numblocklabels]" )
        {
            success &= expectChar(lineStream, '=', err);
            success &= parseValue(lineStream, NumBlockLabels, err);
            if (NumBlockLabels>0) labellist.reserve(NumBlockLabels);

#ifdef DEBUG_PARSER
            {
                char buf[1048]; SNPRINTF(buf, sizeof(buf), "feasolver reserved space for %i block labels, labellist length is %i\n", NumBlockLabels, (int) labellist.size ());
                WarnMessage(buf);
            }
#endif // DEBUG_PARSER

            for(int i=0; i<NumBlockLabels; i++)
            {

                BlockLabelT blk = BlockLabelT::fromStream(input, err);

                labellist.push_back(blk);
            }

            // message will be printed after parsing is done
//            if (NumBlockLabels != i)
//            {
//                err << "Expected "<<NumBlockLabels<<" BlockLabels, but got " << i << "\n";
//                break; // stop parsing
//            }
#ifdef DEBUG
            {
                char buf[1048]; SNPRINTF(buf, sizeof(buf), "feasolver: Expected %i block labels, labellist length is %i\n", NumBlockLabels, (int) labellist.size ());
                WarnMessage(buf);
            }
#endif // DEBUG
#ifdef DEBUG_PARSER
            {
                char buf[1048]; SNPRINTF(buf, sizeof(buf), "feasolver: finished reading block labels\n");
                WarnMessage(buf);
            }
#endif // DEBUG_PARSER
            continue;
        }

        if(token == "[numpoints]"
                || token == "[numsegments]"
                || token == "[numarcsegments]"
                || token == "[numholes]"
                )
        {
            success &= expectChar(lineStream, '=', err);
            // the nodes in the .fem file seem to be unused, and
            // the original femm code does not even read them here
            // -> just skip the lines here:
            int i;
            success &= parseValue(lineStream, i, err);
            while (0 < i--)
                std::getline(input,line);
            continue;
        }

        // fall-through; token was not used
        if (!handleToken(token, lineStream, err))
        {
            err << "Unknown token: " << token << "\n"
                << "Context line:\n" << line << "\n";
            success = false;
#ifdef STOP_ON_UNKNOWN_TOKEN
            // stop parsing:
            break;
#endif
        }
    }
    if (!success)
    {
        string msg = "Parse error while reading input file " + file + "!\n";
        msg += "Last token was: " + token + "\n";
        msg += "Last input line was: " + line + "\n";
        msg += err.str();
        WarnMessage(msg.c_str());
        return false;
    }

    if (!err.str().empty())
    {
        WarnMessage(err.str().c_str());
    }

#ifdef DEBUG_PARSER
    std::cout << "feasolver finished parsing file " << file << std::endl;
#endif // DEBUG_PARSER

    return true;
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
bool FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::handleToken(const string &, istream &, ostream &)
{
    // token not handled
    return false;
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
std::string FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::getErrorString(LoadMeshErr err)
{
    switch (err)
    {
    case (NOERROR):
        return std::string();
    case ( BADEDGEFILE ):
        return "problem loading mesh:\nCould not open .edge file.\n";
    case ( BADELEMENTFILE ):
        return "problem loading mesh:\nCould not open .ele file.\n";
    case( BADFEMFILE ):
        return "problem loading mesh:\nCould not open .fem file.\n";
    case( BADNODEFILE ):
        return "problem loading mesh:\nCould not open .node file.\n";
    case( BADPBCFILE ):
        return "problem loading mesh:\nCould not open .pbc file.\n";
    case( MISSINGMATPROPS ):
        return "problem loading mesh:\nMaterial properties have not been defined for all regions.\n";
    case( ELMLABELTOOBIG ):
        return "problem loading mesh:\nElemnet label number was greater than the number of labels in the problem.\n";
    }

    assert(false);
    return std::string();
}

