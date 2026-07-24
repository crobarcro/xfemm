#include "FemmProblem.h"

#include "femmconstants.h"
#include "make_unique.h"

#include <cassert>
#include <ctgmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>

#ifdef DEBUG_MEX
#include "mex.h"
#endif // DEBUG_MEX

femm::FemmProblem::~FemmProblem()
{
}

bool femm::FemmProblem::saveFEMFile(const std::string &filename) const
{
    if ( filetype == FileType::Unknown )
    {
        // can't save if file type is not set
        std::cerr << "Cannot save file because file type is unknown!\n";
        return false;
    }
    std::ofstream fem (filename);
    if ( fem.fail() )
    {
        // open failed
        std::cerr << "Opening file " << filename << " for writing failed!\n";
        return false;
    }

    writeProblemDescription(fem);
    return true;
}

void femm::FemmProblem::writeProblemDescription(std::ostream &output) const
{
    // set floating point precision once for the whole stream
    output << std::setprecision(17);
    // when filling to a width, adjust to the left
    output.setf(std::ios::left);
    // width and fill-direction are only used to emulate the output of femm42 as close as possible.

    // ***** write properties ****
    if (filetype == FileType::MagneticsFile)
        output << "[Format]      =  4.0\n";
    else
        output << "[Format]      =  1\n";
    if (filetype != FileType::ElectrostaticsFile)
    {
        output.width(12);
        output << "[Frequency]" << "  =  " << Frequency << "\n";
    }
    output.width(12);
    output << "[Precision]" << "  =  " << Precision << "\n";
    output.width(12);
    output << "[MinAngle]" << "  =  " << MinAngle << "\n";
    output.width(12);
    output << "[Depth]" << "  =  " << Depth << "\n";
    output << "[LengthUnits] =  ";
    switch (LengthUnits) {
    case femm::LengthMillimeters:
        output << "millimeters\n";
        break;
    case femm::LengthCentimeters:
        output << "centimeters\n";
        break;
    case femm::LengthMeters:
        output << "meters\n";
        break;
    case femm::LengthMils:
        output << "mils\n";
        break;
    case femm::LengthMicrometers:
        output << "microns\n";
        break;
    case femm::LengthInches:
        output << "inches\n";
        break;
    }

    if (problemType == femm::PLANAR)
    {
        output << "[ProblemType] =  planar\n";
    } else {
        output << "[ProblemType] =  axisymmetric\n";
        if ( extRo != 0 && extRi != 0)
        {
            output.width(12);
            output << "[extZo]" << "  =  " << extZo << "\n";
            output.width(12);
            output << "[extRo]" << "  =  " << extRo << "\n";
            output.width(12);
            output << "[extRi]" << "  =  " << extRi << "\n";
        }
    }
    if (Coords == femm::CART)
    {
        output << "[Coordinates] =  cartesian\n";
    } else {
        output << "[Coordinates] =  polar\n";
    }

    if (filetype == FileType::MagneticsFile)
    {
        output.width(12);
        output << "[ACSolver]" << "  =  " << ACSolver <<"\n";
    }


    output.width(12);
    output << "[PrevSoln]" << "  = \"" << previousSolutionFile << "\"\n";

    output.width(12);
    output << "[PrevType]" << "  =  " << PrevType << "\n";

    std::string commentString (comment);
    // escape line-breaks
    size_t pos = commentString.find('\n');
    while (pos != std::string::npos)
        commentString.replace(pos,1,"\\n");
    output.width(12);
    output << "[Comment]" << "  =  \"" << commentString <<"\"\n";

    // ***** write data ****
    output.width(12);
    output << "[PointProps]" << "  =  " << nodeproplist.size() <<"\n";
    for( const auto &prop: nodeproplist)
    {
        prop->toStream(output);
    }

    output.width(12);
    output << "[BdryProps]" << "  = " << lineproplist.size() <<"\n";
    for( const auto &prop: lineproplist)
    {
        prop->toStream(output);
    }

    output.width(12);
    output << "[BlockProps]" << "  = " << blockproplist.size() <<"\n";
    for( const auto &prop: blockproplist)
    {
        prop->toStream(output);
    }

    std::string circuitHeader;
    if (filetype == FileType::MagneticsFile)
    {
        circuitHeader =  "[CircuitProps]";
    } else {
        circuitHeader =  "[ConductorProps]";
    }

    output.width(12);
    output << circuitHeader << "  = " << circproplist.size() <<"\n";
    for( const auto &prop: circproplist)
    {
        prop->toStream(output);
    }

    // the solver does not even read this, but write it anyways...
    output << "[NumPoints] = " << nodelist.size() <<"\n";
    for (const auto &node: nodelist)
    {
        output << node->x
               << "\t" << node->y
               << "\t" << node->BoundaryMarker+1
               << "\t" << node->InGroup;

        if ( filetype == FileType::HeatFlowFile ||
             filetype == FileType::ElectrostaticsFile )
        {
            output << "\t" << node->InConductor+1;
        }
        output << "\n";
    }
    // the solver does not even read this, but write it anyways...
    output << "[NumSegments] = " << linelist.size() <<"\n";
    for (const auto &line: linelist)
    {
        // n0 n1 MaxSideLength BoundaryMarker Hidden Group [InConductor]
        output << line->n0 << "\t" << line->n1;
        if (line->MaxSideLength < 0)
        {
            output << "\t-1";
        } else {
            output << "\t" << line->MaxSideLength;
        }

        output << "\t" << line->BoundaryMarker+1
            << "\t" << (int)line->Hidden
            << "\t" << line->InGroup;

        if ( filetype == FileType::HeatFlowFile ||
             filetype == FileType::ElectrostaticsFile )
        {
            output << "\t" << line->InConductor+1;
        }
        output << "\n";
    }
    // the solver does not even read this, but write it anyways...
    output << "[NumArcSegments] = " << arclist.size() <<"\n";
    for (const auto &arc: arclist)
    {
        // n0 n1 ArcLength MaxSideLength BoundaryMarker Hidden Group [InConductor|MeshedSideLength]
        output << arc->n0
               << "\t" << arc->n1
               << "\t" << arc->ArcLength
               << "\t" << arc->MaxSideLength
               << "\t" << arc->BoundaryMarker+1
               << "\t" << (int)arc->Hidden
               << "\t" << arc->InGroup;


        if ( filetype == FileType::HeatFlowFile ||
             filetype == FileType::ElectrostaticsFile )
        {
            output << "\t" << arc->InConductor+1;
        }
        else if (filetype == FileType::MagneticsFile)
        {
            output << "\t" << arc->mySideLength;
        }
        output << "\n";
    }


    // write out list of holes;
    int numHoles=countHoles();

    output << "[NumHoles] = " << numHoles << "\n";
    for (const auto &label: labellist)
    {
        if(label->isHole())
            output << label->x << " " << label->y << " " << label->InGroup << "\n";
    }

    output << "[NumBlockLabels] = " << labellist.size()-numHoles << "\n";
    for (const auto &label: labellist)
    {
        if(!label->isHole())
        {
            label->toStream(output);
        }
    }
}

void femm::FemmProblem::updateLabelsFromIndex()
{
    // block labels
    for (auto &label: labellist)
    {
        if (label->hasBlockType())
            label->BlockTypeName = blockproplist[label->BlockType]->BlockName;
        if (label->isInCircuit())
            label->InCircuitName = circproplist[label->InCircuit]->CircName;
    }
    // points
    for (auto &node: nodelist)
    {
        if ( node->hasBoundaryMarker())
            node->BoundaryMarkerName = nodeproplist[node->BoundaryMarker]->PointName;
        if ( node->isInConductor())
            node->InConductorName = circproplist[node->InConductor]->CircName;
    }
    // segments
    for (auto &segm: linelist)
    {
        if (segm->hasBoundaryMarker())
            segm->BoundaryMarkerName = lineproplist[segm->BoundaryMarker]->BdryName;
        if (segm->isInConductor())
            segm->InConductorName = circproplist[segm->InConductor]->CircName;
    }
    // arc segments
    for (auto &asegm: arclist)
    {
        if (asegm->hasBoundaryMarker())
            asegm->BoundaryMarkerName = lineproplist[asegm->BoundaryMarker]->BdryName;
        if (asegm->isInConductor())
            asegm->InConductorName = circproplist[asegm->InConductor]->CircName;
    }
}

void femm::FemmProblem::updateBlockMap()
{
    int idx=0;
    blockMap.clear();
    for (const auto &bprop: blockproplist)
    {
        blockMap[bprop->BlockName] = idx++;
    }
}

void femm::FemmProblem::updateCircuitMap()
{
    int idx=0;
    circuitMap.clear();
    for (const auto &circuit: circproplist)
    {
        circuitMap[circuit->CircName] = idx++;
    }
}

void femm::FemmProblem::updateLineMap()
{
    int idx=0;
    lineMap.clear();
    for (const auto &bprop: lineproplist)
    {
        lineMap[bprop->BdryName] = idx++;
    }
}

void femm::FemmProblem::updateNodeMap()
{
    int idx=0;
    nodeMap.clear();
    for (const auto &pprop: nodeproplist)
    {
        nodeMap[pprop->PointName] = idx++;
    }
}

bool femm::FemmProblem::addArcSegment(femm::CArcSegment &asegm, double tol)
{
    // don't add if line is degenerate
    if (asegm.n0==asegm.n1)
        return false;

    // don't add if the arc is already in the list;
    for(int i=0; i<(int)arclist.size(); i++){
        if ((arclist[i]->n0==asegm.n0) && (arclist[i]->n1==asegm.n1) &&
                (fabs(arclist[i]->ArcLength-asegm.ArcLength)<1.e-02)) return false;
        // arcs are ``the same'' if start and end points are the same, and if
        // the arc lengths are relatively close (but a lot farther than
        // machine precision...
    }

    // add proposed arc to the linelist
    asegm.IsSelected = false;

    CComplex p[2];
    std::vector < CComplex > newnodes;
    // check to see if there are intersections
    for(int i=0; i<(int)linelist.size(); i++)
    {
        int j = getLineArcIntersection(*linelist[i],asegm,p);
        if (j>0)
            for(int k=0; k<j; k++)
                newnodes.push_back(p[k]);
    }
    for (int i=0; i<(int)arclist.size(); i++)
    {
        int j = getArcArcIntersection(asegm,*arclist[i],p);
        if (j>0)
            for(int k=0; k<j; k++)
                newnodes.push_back(p[k]);
    }

    // add nodes at intersections
    double t;
    if (tol==0)
    {
        if (nodelist.size()<2) t=1.e-08;
        else{
            CComplex p0,p1;
            p0 = nodelist[0]->CC();
            p1 = p0;
            for (int i=1; i<(int)nodelist.size(); i++)
            {
                if(nodelist[i]->x<p0.re) p0.re = nodelist[i]->x;
                if(nodelist[i]->x>p1.re) p1.re = nodelist[i]->x;
                if(nodelist[i]->y<p0.im) p0.im = nodelist[i]->y;
                if(nodelist[i]->y>p1.im) p1.im = nodelist[i]->y;
            }
            t = abs(p1-p0)*CLOSE_ENOUGH;
        }
    }
    else t = tol;

    for (int i=0; i<(int)newnodes.size(); i++)
        addNode(newnodes[i].re,newnodes[i].im,t);

    // add proposed arc segment;
    arclist.push_back(MAKE_UNIQUE<CArcSegment>(asegm));

    // check to see if proposed arc passes through other points;
    // if so, delete arc and create arcs that link intermediate points;
    // does this by recursive use of AddArcSegment;

    unselectAll();
    CComplex c;
    double R;
    getCircle(asegm,c,R);

    double dmin = tol;
    if (tol==0)
        dmin = fabs(R*PI*asegm.ArcLength/180.)*1.e-05;

    int k = (int)arclist.size()-1;
    for(int i=0; i<(int)nodelist.size(); i++)
    {
        if( (i!=asegm.n0) && (i!=asegm.n1) )
        {
            double d=shortestDistanceFromArc(CComplex(nodelist[i]->x,nodelist[i]->y),*arclist[k]);

            //	MsgBox("d=%g dmin=%g",d,dmin);
            // what is the purpose of this test?
            //	if (abs(nodelist[i]->CC()-nodelist[asegm.n0]->CC())<2.*dmin) d=2.*dmin;
            //	if (abs(nodelist[i]->CC()-nodelist[asegm.n1]->CC())<2.*dmin) d=2.*dmin;


            if (d<dmin){

                CComplex a0,a1,a2;
                a0.Set(nodelist[asegm.n0]->x,nodelist[asegm.n0]->y);
                a1.Set(nodelist[asegm.n1]->x,nodelist[asegm.n1]->y);
                a2.Set(nodelist[i]->x,nodelist[i]->y);
                arclist[k]->ToggleSelect();
                deleteSelectedArcSegments();

                CArcSegment newarc = asegm;
                newarc.n1 = i;
                newarc.ArcLength = arg((a2-c)/(a0-c))*180./PI;
                addArcSegment(newarc,dmin);

                newarc = asegm;
                newarc.n0 = i;
                newarc.ArcLength = arg((a1-c)/(a2-c))*180./PI;
                addArcSegment(newarc,dmin);

                i = nodelist.size();
            }
        }
    }

    return true;
}

bool femm::FemmProblem::addBlockLabel(double x, double y, double d)
{
    std::unique_ptr<CBlockLabel> pt;
    switch (filetype) {
    case FileType::MagneticsFile:
        pt = MAKE_UNIQUE<CMBlockLabel>();
        break;
    case FileType::HeatFlowFile:
        pt = MAKE_UNIQUE<CHBlockLabel>();
        break;
    case FileType::ElectrostaticsFile:
        pt = MAKE_UNIQUE<CSBlockLabel>();
        break;
    default:
        assert(false && "Unhandled file type");
        break;
    }
    pt->x = x;
    pt->y = y;

    return addBlockLabel(std::move(pt), d);
}

bool femm::FemmProblem::addBlockLabel(std::unique_ptr<femm::CBlockLabel> &&label, double d)
{
    double x = label->x;
    double y = label->y;

    // can't put a block label on top of an existing node...
    for (int i=0; i<(int)nodelist.size(); i++)
        if(nodelist[i]->GetDistance(x,y)<d) return false;

    // can't put a block label on a line, either...
    for (int i=0; i<(int)linelist.size(); i++)
        if(shortestDistanceFromSegment(x,y,i)<d) return false;

    // test to see if ``too close'' to existing node...
    bool exists=false;
    for (int i=0; i<(int)labellist.size(); i++)
        if(labellist[i]->GetDistance(x,y)<d) {
            exists=true;
            break;
        }

    // if all is OK, add point in to the node list...
    if(!exists){
        labellist.push_back(std::move(label));
    }

    return true;
}



bool femm::FemmProblem::addNode(double x, double y, double d)
{
    // create an appropriate node and call addNode on it
    std::unique_ptr<CNode> node = std::unique_ptr<CNode>(new CNode(x,y));
    return addNode(std::move(node), d);
}

bool femm::FemmProblem::addNode(std::unique_ptr<femm::CNode> &&node, double d)
{
    CComplex c,a0,a1,a2;
    double R;
    double x = node->x;
    double y = node->y;

    // test to see if ``too close'' to existing node...
    for (int i=0; i<(int)nodelist.size(); i++)
        if(nodelist[i]->GetDistance(x,y)<d) return false;

    // can't put a node on top of a block label; do same sort of test.
    for (int i=0;i<(int)labellist.size();i++)
        if(labellist[i]->GetDistance(x,y)<d) return false;

    // if all is OK, add point in to the node list...
    nodelist.push_back(std::move(node));

    // test to see if node is on an existing line; if so,
    // break into two lines;

    for(int i=0, k=(int)linelist.size(); i<k; i++)
    {
        if (fabs(shortestDistanceFromSegment(x,y,i))<d)
        {
            std::unique_ptr<CSegment> segm;
            segm = linelist[i]->clone();
            linelist[i]->n1=nodelist.size()-1;
            segm->n0=nodelist.size()-1;
            linelist.push_back(std::move(segm));
        }
    }

    // test to see if node is on an existing arc; if so,
    // break into two arcs;
    for(int i=0, k=(int)arclist.size(); i<k; i++)
    {
        if (shortestDistanceFromArc(CComplex(x,y),*arclist[i])<d)
        {
            a0.Set(nodelist[arclist[i]->n0]->x,nodelist[arclist[i]->n0]->y);
            a1.Set(nodelist[arclist[i]->n1]->x,nodelist[arclist[i]->n1]->y);
            a2.Set(x,y);
            getCircle(*arclist[i],c,R);

            std::unique_ptr<CArcSegment> asegm;
            asegm = MAKE_UNIQUE<CArcSegment>(*arclist[i]);
            arclist[i]->n1 = nodelist.size()-1;
            arclist[i]->ArcLength = arg((a2-c)/(a0-c))*180./PI;
            asegm->n0 = nodelist.size()-1;
            asegm->ArcLength = arg((a1-c)/(a2-c))*180./PI;
            arclist.push_back(std::move(asegm));
        }
    }
    return true;
}

bool femm::FemmProblem::addSegment(int n0, int n1, double tol)
{
    return addSegment(n0,n1,nullptr,tol);
}

bool femm::FemmProblem::addSegment(int n0, int n1, const femm::CSegment *parsegm, double tol)
{
    double xi,yi,t;
    CComplex p[2];
    CSegment segm;
    std::vector < CComplex > newnodes;

    // don't add if line is degenerate
    if (n0==n1) return false;

    // don't add if the line is already in the list;
    for (int i=0; i<(int)linelist.size(); i++){
        if ((linelist[i]->n0==n0) && (linelist[i]->n1==n1)) return false;
        if ((linelist[i]->n0==n1) && (linelist[i]->n1==n0)) return false;
    }

    // add proposed line to the linelist
    segm.BoundaryMarkerName="<None>";
    if (parsegm!=NULL) segm=*parsegm;
    segm.IsSelected=false;
    segm.n0=n0; segm.n1=n1;

    // check to see if there are intersections with segments
    for (int i=0; i<(int)linelist.size(); i++)
        if(getIntersection(n0,n1,i,&xi,&yi)) newnodes.push_back(CComplex(xi,yi));

    // check to see if there are intersections with arcs
    for (int i=0; i<(int)arclist.size(); i++){
        int j = getLineArcIntersection(segm,*arclist[i],p);
        if (j>0)
            for(int k=0;k<j;k++)
                newnodes.push_back(p[k]);
    }

    // add nodes at intersections
    if (tol==0)
    {
        if (nodelist.size()<2)
            t = 1.e-08;
        else{
            CComplex p0,p1;
            p0 = nodelist[0]->CC();
            p1 = p0;
            for (int i=1; i<(int)nodelist.size(); i++)
            {
                if(nodelist[i]->x<p0.re) p0.re=nodelist[i]->x;
                if(nodelist[i]->x>p1.re) p1.re=nodelist[i]->x;
                if(nodelist[i]->y<p0.im) p0.im=nodelist[i]->y;
                if(nodelist[i]->y>p1.im) p1.im=nodelist[i]->y;
            }
            t=abs(p1-p0)*CLOSE_ENOUGH;
        }
    }
    else t=tol;

    for (int i=0; i<(int)newnodes.size(); i++)
        addNode(newnodes[i].re,newnodes[i].im,t);

    // Add proposed line segment
    linelist.push_back(segm.clone());

    // check to see if proposed line passes through other points;
    // if so, delete line and create lines that link intermediate points;
    // does this by recursive use of AddSegment;
    double d,dmin;
    unselectAll();
    if (tol==0)
        dmin = abs(nodelist[n1]->CC()-nodelist[n0]->CC())*1.e-05;
    else dmin = tol;

    for (int i=0, k=linelist.size()-1; i<(int)nodelist.size(); i++)
    {
        if( (i!=n0) && (i!=n1) )
        {
            d=shortestDistanceFromSegment(nodelist[i]->x,nodelist[i]->y,k);
            if (abs(nodelist[i]->CC()-nodelist[n0]->CC())<dmin) d=2.*dmin;
            if (abs(nodelist[i]->CC()-nodelist[n1]->CC())<dmin) d=2.*dmin;
            if (d<dmin){
                linelist[k]->ToggleSelect();
                deleteSelectedSegments();
                if(parsegm==NULL)
                {
                    addSegment(n0,i,dmin);
                    addSegment(i,n1,dmin);
                }
                else{
                    addSegment(n0,i,&segm,dmin);
                    addSegment(i,n1,&segm,dmin);
                }
                i=nodelist.size();
            }
        }
    }

    return true;
}

void femm::FemmProblem::clearNotationTags()
{
    for (auto &line : linelist)
    {
        line->cnt = 0;
    }
    for (auto &arc : arclist)
    {
        arc->cnt = 0;
    }
}




// identical in fmesher, FPProc and HPProc
int femm::FemmProblem::closestArcSegment(double x, double y) const
{
    if(arclist.size()==0) return -1;

    double d0=shortestDistanceFromArc(CComplex(x,y),*arclist[0]);
    int idx=0;
    for(int i=0; i<(int)arclist.size(); i++)
    {
        double d1=shortestDistanceFromArc(CComplex(x,y),*arclist[i]);
        if(d1<d0)
        {
            d0=d1;
            idx=i;
        }
    }

    return idx;
}

int femm::FemmProblem::closestBlockLabel(double x, double y) const
{
    if(labellist.size()==0) return -1;

    int idx=0;
    double d0=labellist[0]->GetDistance(x,y);
    for(int i=0; i<(int)labellist.size(); i++)
    {
        double d1=labellist[i]->GetDistance(x,y);
        if(d1<d0)
        {
            d0=d1;
            idx=i;
        }
    }

    return idx;
}

// identical in fmesher, FPProc, and HPProc
int femm::FemmProblem::closestNode(double x, double y) const
{
    if(nodelist.size()==0) return -1;

    int idx=0;
    double d0=nodelist[0]->GetDistance(x,y);
    for(int i=0; i<(int)nodelist.size(); i++)
    {
        double d1=nodelist[i]->GetDistance(x,y);
        if(d1<d0)
        {
            d0=d1;
            idx=i;
        }
    }

    return idx;
}

// identical in fmesher, hpproc
int femm::FemmProblem::closestSegment(double x, double y) const
{
    if(linelist.size()==0) return -1;

    int idx=0;
    double d0=shortestDistanceFromSegment(x,y,0);
    for(int i=0; i<(int)linelist.size(); i++)
    {
        double d1=shortestDistanceFromSegment(x,y,i);
        if(d1<d0)
        {
            d0=d1;
            idx=i;
        }
    }

    return idx;
}

bool femm::FemmProblem::consistencyCheckOK() const
{
    using std::to_string;
    bool ok = true;
    // block labels
    for (const auto &label: labellist)
    {
        if (label->hasBlockType())
        {
            int idx = label->BlockType;
            if ( idx < 0 || idx >= (int)blockproplist.size() )
            {
                std::cerr << "Label has invalid block type index: " << to_string(idx)
                          << " (0 <= index < " << to_string(blockproplist.size()) << ")\n";
                ok = false;
            } else if (label->BlockTypeName != blockproplist[idx]->BlockName)
            {
                std::cerr << "Label has inconsistent block type name: " << label->BlockTypeName
                          << ", but block type " << to_string(idx)
                          << " has name " << blockproplist[idx]->BlockName << "\n";
                ok = false;
            }
        }
        if (label->isInCircuit())
        {
            int idx = label->InCircuit;
            if ( idx < 0 || idx >= (int)circproplist.size() )
            {
                std::cerr << "Label has invalid circuit index: " << to_string(idx)
                          << " (0 <= index < " << to_string(circproplist.size()) << ")\n";
                ok = false;
            } else if (label->InCircuitName != circproplist[idx]->CircName)
            {
                std::cerr << "Label has inconsistent circuit name: " << label->InCircuitName
                          << ", but circuit " << to_string(idx)
                          << " has name " << circproplist[idx]->CircName << "\n";
                ok = false;
            }
        }
    }
    // points
    for (const auto &node: nodelist)
    {
        if (node->hasBoundaryMarker())
        {
            int idx = node->BoundaryMarker;
            if ( idx < 0 || idx >= (int)nodeproplist.size() )
            {
                std::cerr << "Point has invalid boundary marker index: " << to_string(idx)
                          << " (0 <= index < " << to_string(nodeproplist.size()) << ")\n";
                ok = false;
            } else if (node->BoundaryMarkerName != nodeproplist[idx]->PointName)
            {
                std::cerr << "Point has inconsistent boundary marker name: " << node->BoundaryMarkerName
                          << ", but boundary marker " << to_string(idx)
                          << " has name " << nodeproplist[idx]->PointName << "\n";
                ok = false;
            }
        }
        if (node->isInConductor())
        {
            int idx = node->InConductor;
            if ( idx < 0 || idx >= (int)circproplist.size() )
            {
                std::cerr << "Point has invalid conductor index: " << to_string(idx)
                          << " (0 <= index < " << to_string(circproplist.size()) << ")\n";
                ok = false;
            } else if (node->InConductorName != circproplist[idx]->CircName)
            {
                std::cerr << "Point has inconsistent conductor name: " << node->InConductorName
                          << ", but conductor " << to_string(idx)
                          << " has name " << circproplist[idx]->CircName << "\n";
                ok = false;
            }
        }
    }
    // segments
    for (const auto &segm: linelist)
    {
        if (segm->hasBoundaryMarker())
        {
            int idx = segm->BoundaryMarker;
            if ( idx < 0 || idx >= (int)lineproplist.size() )
            {
                std::cerr << "Line has invalid boundary marker index: " << to_string(idx)
                          << " (0 <= index < " << to_string(lineproplist.size()) << ")\n";
                ok = false;
            } else if (segm->BoundaryMarkerName != lineproplist[idx]->BdryName)
            {
                std::cerr << "Line has inconsistent boundary marker name: " << segm->BoundaryMarkerName
                          << ", but boundary marker " << to_string(idx)
                          << " has name " << lineproplist[idx]->BdryName << "\n";
                ok = false;
            }
        }
        if (segm->isInConductor())
        {
            int idx = segm->InConductor;
            if ( idx < 0 || idx >= (int)circproplist.size() )
            {
                std::cerr << "Point has invalid circuit index: " << to_string(idx)
                          << " (0 <= index < " << to_string(circproplist.size()) << ")\n";
                ok = false;
            } else if (segm->InConductorName != circproplist[idx]->CircName)
            {
                std::cerr << "Point has inconsistent circuit name: " << segm->InConductorName
                          << ", but circuit " << to_string(idx)
                          << " has name " << circproplist[idx]->CircName << "\n";
                ok = false;
            }
        }
    }
    // arc segments
    for (const auto &asegm: arclist)
    {
        if (asegm->hasBoundaryMarker())
        {
            int idx = asegm->BoundaryMarker;
            if ( idx < 0 || idx >= (int)lineproplist.size() )
            {
                std::cerr << "Line has invalid boundary marker index: " << to_string(idx)
                          << " (0 <= index < " << to_string(lineproplist.size()) << ")\n";
                ok = false;
            } else if (asegm->BoundaryMarkerName != lineproplist[idx]->BdryName)
            {
                std::cerr << "Line has inconsistent boundary marker name: " << asegm->BoundaryMarkerName
                          << ", but boundary marker " << to_string(idx)
                          << " has name " << lineproplist[idx]->BdryName << "\n";
                ok = false;
            }
        }
        if (asegm->isInConductor())
        {
            int idx = asegm->InConductor;
            if ( idx < 0 || idx >= (int)circproplist.size() )
            {
                std::cerr << "Point has invalid circuit index: " << to_string(idx)
                          << " (0 <= index < " << to_string(circproplist.size()) << ")\n";
                ok = false;
            } else if (asegm->InConductorName != circproplist[idx]->CircName)
            {
                std::cerr << "Point has inconsistent circuit name: " << asegm->InConductorName
                          << ", but circuit " << to_string(idx)
                          << " has name " << circproplist[idx]->CircName << "\n";
                ok = false;
            }
        }
    }
    return ok;
}

int femm::FemmProblem::countHoles() const
{
    int numHoles = 0;
    for (const auto &label: labellist)
    {
        if(label->isHole())
            numHoles++;
    }
    return numHoles;
}

bool femm::FemmProblem::canCreateRadius(int n) const
{
    // check to see if a selected point, specified by n, can be
    // converted to a radius.  To be able to be converted to a radius,
    // the point must be an element of either 2 lines, 2 arcs, or
    // 1 line and 1 arc.

    int j=0;
    for (const auto &line: linelist)
    {
        if (n==line->n0 || n==line->n1)
            j++;
    }
    for (const auto &arc: arclist)
    {
        if (n==arc->n0 || n==arc->n1)
            j++;
    }

    return (j==2);
}

bool femm::FemmProblem::createRadius(int n, double r)
{
    // replace the node indexed by n with a radius of r

    if(r<=0)
        return false;

    int nseg=0;
    int narc=0;
    int arc[2],seg[2];

    for(int k=0; k<(int)linelist.size(); k++)
        if ((linelist[k]->n0==n) || (linelist[k]->n1==n))
            seg[nseg++]=k;
    for(int k=0; k<(int)arclist.size(); k++)
        if ((arclist[k]->n0==n) || (arclist[k]->n1==n))
            arc[narc++]=k;

    if ((narc+nseg)!=2)
        return false;

    // there are three valid cases:
    switch (nseg-narc)
    {
    case 0:  // One arc and one line
    {
        CComplex c,u,p0,p1,q,p[4],v[8],i1[8],i2[8];
        double rc,b,R[4],phi;
        CArcSegment ar;

        // inherit the boundary condition from the arc so that
        // we can apply it to the newly created arc later;
        ar.InGroup       =arclist[arc[0]]->InGroup;
        ar.BoundaryMarkerName=arclist[arc[0]]->BoundaryMarkerName;

        // get the center and radius of the circle associated with the arc;
        getCircle(*arclist[arc[0]],c,rc);

        // get the locations of the endpoints of the segment;
        p0=nodelist[n]->CC();
        if(linelist[seg[0]]->n0==n)
            p1=nodelist[linelist[seg[0]]->n1]->CC();
        else
            p1=nodelist[linelist[seg[0]]->n0]->CC();

        u=(p1-p0)/abs(p1-p0);  // unit vector along the line
        q=p0 + u*Re((c-p0)/u); // closest point on line to center of circle
        u=(q-c)/abs(q-c); // unit vector from center to closest point on line;

        p[0]=q+r*u; R[0]=rc+r;
        p[1]=q-r*u; R[1]=rc+r;
        p[2]=q+r*u; R[2]=rc-r;
        p[3]=q-r*u; R[3]=rc-r;

        int j=0;
        for(int k=0;k<4;k++)
        {
            b=R[k]*R[k]-abs(p[k]-c)*abs(p[k]-c);
            if (b>=0){
                b=sqrt(b);
                v[j++]=p[k]+I*b*(p[k]-c)/abs(p[k]-c);
                v[j++]=p[k]-I*b*(p[k]-c)/abs(p[k]-c);
            }
        }

        // locations of the centerpoints that could be for the radius that
        // we are looking for are stored in v.  We now need to paw through
        // them to find the one solution that we are after.
        u=(p1-p0)/abs(p1-p0);  // unit vector along the line
        int m=0;
        for(int k=0;k<j;k++)
        {
            i1[m]=p0 +u*Re((v[k]-p0)/u); // intersection with the line
            i2[m]=c + rc*(v[k]-c)/abs(v[k]-c); // intersection with the arc;
            v[m]=v[k];

            // add this one to the list of possibly valid solutions if
            // both of the intersection points actually lie on the arc
            if ( shortestDistanceFromArc(i2[m],*arclist[arc[0]])<(r/10000.) &&
                 shortestDistanceFromSegment(Re(i1[m]),Im(i1[m]),seg[0])<(r/10000.)
                 && abs(i1[m]-i2[m])>(r/10000.))
            {
                m++;
                if (m==2) break;
            }

        }

        if (m==0) return false;

        // But there are also special cases where there could be two answers.
        // We then pick the solution that has the center point closest to the point to be removed.
        if(m>1)
        {
            if (abs(v[0]-p0)<abs(v[1]-p0)) j=0;
            else j=1;
        }
        else j=0;	// The index of the winning case is in j....

        updateUndo();
        addNode(Re(i1[j]),Im(i1[j]),r/10000.);
        addNode(Re(i2[j]),Im(i2[j]),r/10000.);
        unselectAll();

        // delete the node that is to be replace by a radius;
        n=closestNode(Re(p0),Im(p0));
        nodelist[n]->IsSelected=true;
        deleteSelectedNodes();

        // compute the angle spanned by the new arc;
        phi=arg((i2[j]-v[j])/(i1[j]-v[j]));
        if (phi<0)
        {
            c=i2[j]; i2[j]=i1[j]; i1[j]=c;
            phi=fabs(phi);
        }

        // add in the new radius;
        ar.n0=closestNode(Re(i1[j]),Im(i1[j]));
        ar.n1=closestNode(Re(i2[j]),Im(i2[j]));
        ar.ArcLength=phi/DEG;
        addArcSegment(ar);

        return true;
    }
    case 2:  // Two lines
    {
        CComplex p0,p1,p2;
        double phi,len;
        CArcSegment ar;

        if (linelist[seg[0]]->n0==n) p1=nodelist[linelist[seg[0]]->n1]->CC();
        else p1=nodelist[linelist[seg[0]]->n0]->CC();

        if (linelist[seg[1]]->n0==n) p2=nodelist[linelist[seg[1]]->n1]->CC();
        else p2=nodelist[linelist[seg[1]]->n0]->CC();

        p0=nodelist[n]->CC();

        // get the angle between the lines
        phi=arg((p2-p0)/(p1-p0));

        // check to see if this case is degenerate
        if (fabs(phi)>(179.*DEG)) return false;

        // check to see if the points are in the wrong order
        // and fix it if they are.
        if (phi<0){
            p0=p1; p1=p2; p2=p0; p0=nodelist[n]->CC();
            std::swap(seg[0],seg[1]);
            phi=fabs(phi);
        }

        len = r/tan(phi/2.); // distance from p0 to the tangency point;

        // catch the case where the desired radius is too big to fit;
        if ((abs(p1-p0)<len) || (abs(p2-p0)<len)) return false;

        // compute the locations of the tangency points;
        p1=len*(p1-p0)/abs(p1-p0)+p0;
        p2=len*(p2-p0)/abs(p2-p0)+p0;

        // inherit the boundary condition from one of the segments
        // so that we can apply it to the newly created arc later;
        ar.BoundaryMarkerName=linelist[seg[0]]->BoundaryMarkerName;
        ar.InGroup       =linelist[seg[0]]->InGroup;

        // add new nodes at ends of radius
        updateUndo();
        addNode(Re(p1),Im(p1),len/10000.);
        addNode(Re(p2),Im(p2),len/10000.);
        unselectAll();

        // delete the node that is to be replace by a radius;
        n=closestNode(Re(p0),Im(p0));
        nodelist[n]->IsSelected=true;
        deleteSelectedNodes();

        // add in the new radius;
        ar.n0=closestNode(Re(p2),Im(p2));
        ar.n1=closestNode(Re(p1),Im(p1));
        ar.ArcLength=180.-phi/DEG;
        addArcSegment(ar);

        return true;
    }
    case -2: // Two arcs
    {
        CComplex c0,c1,c2,p[8],i1[8],i2[8];
        double a[8],b[8],c,d[8],x[8],r0,r1,r2,phi;
        CArcSegment ar;

        r0=r;
        getCircle(*arclist[arc[0]],c1,r1);
        getCircle(*arclist[arc[1]],c2,r2);
        c=abs(c2-c1);

        // solve for all of the different possible cases;
        a[0]=r1+r0; b[0]=r2+r0;
        a[1]=r1+r0; b[1]=r2+r0;
        a[2]=r1-r0; b[2]=r2-r0;
        a[3]=r1-r0; b[3]=r2-r0;
        a[4]=r1-r0; b[4]=r2+r0;
        a[5]=r1-r0; b[5]=r2+r0;
        a[6]=r1+r0; b[6]=r2-r0;
        a[7]=r1+r0; b[7]=r2-r0;
        for(int k=0;k<8;k++){
            x[k]=(b[k]*b[k]+c*c-a[k]*a[k])/(2.*c*c);
            d[k]=sqrt(b[k]*b[k]-x[k]*x[k]*c*c);
        }
        for(int k=0;k<8;k+=2)
        {
            // solve for the center point of the radius for each solution.
            p[k]  =((1-x[k])*c+I*d[k])*(c2-c1)/abs(c2-c1) + c1;
            p[k+1]=((1-x[k])*c-I*d[k])*(c2-c1)/abs(c2-c1) + c1;
        }

        c0=nodelist[n]->CC();

        int j=0;
        for(int k=0;k<8;k++)
        {
            i1[j]=c1 + r1*(p[k]-c1)/abs(p[k]-c1); // compute possible intersection points
            i2[j]=c2 + r2*(p[k]-c2)/abs(p[k]-c2); // with the arcs;
            p[j] =p[k];

            // add this one to the list of possibly valid solutions if
            // both of the intersection points actually lie on the arc
            if ( shortestDistanceFromArc(i1[j],*arclist[arc[0]])<(r0/10000.) &&
                 shortestDistanceFromArc(i2[j],*arclist[arc[1]])<(r0/10000.) &&
                 abs(i1[j]-i2[j])>(r0/10000.))
            {
                j++;
                if (j==2) break;
            }

        }

        // There could be no valid solutions...
        if (j==0) return false;

        // But there are also special cases where there could be two answers.
        // We then pick the solution that has the center point closest to the point to be removed.
        if(j>1)
        {
            if (abs(p[0]-c0)<abs(p[1]-c0)) j=0;
            else j=1;
        }
        else j=0;	// The index of the winning case is in j....

        // inherit the boundary condition from one of the segments
        // so that we can apply it to the newly created arc later;
        ar.BoundaryMarkerName=arclist[arc[0]]->BoundaryMarkerName;
        ar.InGroup=arclist[arc[0]]->InGroup;

        // add new nodes at ends of radius
        updateUndo();
        addNode(Re(i1[j]),Im(i1[j]),c/10000.);
        addNode(Re(i2[j]),Im(i2[j]),c/10000.);
        unselectAll();

        // delete the node that is to be replace by a radius;
        n=closestNode(Re(c0),Im(c0));
        nodelist[n]->IsSelected=true;
        deleteSelectedNodes();

        // compute the angle spanned by the new arc;
        phi=arg((i2[j]-p[j])/(i1[j]-p[j]));
        if (phi<0)
        {
            c0=i2[j]; i2[j]=i1[j]; i1[j]=c0;
            phi=fabs(phi);
        }

        // add in the new radius;
        ar.n0=closestNode(Re(i1[j]),Im(i1[j]));
        ar.n1=closestNode(Re(i2[j]),Im(i2[j]));
        ar.ArcLength=phi/DEG;
        addArcSegment(ar);

        return true;
    }
    }

    return false;
}

femm::EditMode femm::FemmProblem::defaultEditMode() const
{
    return d_EditMode;
}

void femm::FemmProblem::setDefaultEditMode(femm::EditMode mode)
{
    d_EditMode = mode;
}

bool femm::FemmProblem::deleteSelectedArcSegments()
{
    size_t oldsize = arclist.size();

    if (!arclist.empty())
    {
        // remove selected elements
        arclist.erase(
                    std::remove_if(arclist.begin(),arclist.end(),
                                   [](const std::unique_ptr<femm::CArcSegment>& arc){ return arc->IsSelected;} ),
                    arclist.end()
                    );
    }
    arclist.shrink_to_fit();

    return arclist.size() != oldsize;
}

bool femm::FemmProblem::deleteSelectedBlockLabels()
{
    size_t oldsize = labellist.size();

    if (!labellist.empty())
    {
        // remove selected elements
        labellist.erase(
                    std::remove_if(labellist.begin(),labellist.end(),
                                   [](const std::unique_ptr<femm::CBlockLabel>& label){ return label->IsSelected;} ),
                    labellist.end()
                    );
    }
    labellist.shrink_to_fit();

    return labellist.size() != oldsize;
}

bool femm::FemmProblem::deleteSelectedNodes()
{
    bool changed = false;

    if (nodelist.size() > 0)
    {
        int i=0;
        do
        {
            if(nodelist[i]->IsSelected!=false)
            {
                changed=true;
                // first remove all lines that contain the point;
                for (int j=0; j<(int)linelist.size(); j++)
                    if((linelist[j]->n0==i) || (linelist[j]->n1==i))
                        linelist[j]->ToggleSelect();
                deleteSelectedSegments();

                // remove all arcs that contain the point;
                for (int j=0; j<(int)arclist.size(); j++)
                    if((arclist[j]->n0==i) || (arclist[j]->n1==i))
                        arclist[j]->ToggleSelect();
                deleteSelectedArcSegments();

                // remove node from the nodelist...
                nodelist.erase(nodelist.begin()+i);

                // update lines to point to the new node numbering
                for (int j=0; j<(int)linelist.size(); j++)
                {
                    if (linelist[j]->n0>i) linelist[j]->n0--;
                    if (linelist[j]->n1>i) linelist[j]->n1--;
                }

                // update arcs to point to the new node numbering
                for (int j=0; j<(int)arclist.size(); j++)
                {
                    if (arclist[j]->n0>i) arclist[j]->n0--;
                    if (arclist[j]->n1>i) arclist[j]->n1--;
                }
            } else
                i++;
        } while (i<(int)nodelist.size());
    }

    nodelist.shrink_to_fit();
    return changed;
}

bool femm::FemmProblem::deleteSelectedSegments()
{
    size_t oldsize = linelist.size();

    if (!linelist.empty())
    {
        // remove selected elements
        linelist.erase(
                    std::remove_if(linelist.begin(),linelist.end(),
                                   [](const std::unique_ptr<femm::CSegment>& segm){ return segm->IsSelected;} ),
                    linelist.end()
                    );
    }
    linelist.shrink_to_fit();

    return linelist.size() != oldsize;
}

void femm::FemmProblem::enforcePSLG(double tol)
{
    std::vector< std::unique_ptr<CNode>> newnodelist;
    std::vector< std::unique_ptr<CSegment>> newlinelist;
    std::vector< std::unique_ptr<CArcSegment>> newarclist;
    std::vector< std::unique_ptr<CBlockLabel>> newlabellist;

    // save existing objects into new*list:
    newnodelist.swap(nodelist);
    newlinelist.swap(linelist);
    newarclist.swap(arclist);
    newlabellist.swap(labellist);

    // find out what tolerance is so that there are not nodes right on
    // top of each other;
    double d = tol;
    if (tol==0)
    {
        if (newnodelist.size()<=1)
            d = 1.e-08;
        else
        {
            CComplex p0 = newnodelist[0]->CC();
            CComplex p1 = p0;
            for (int i=1; i<(int)newnodelist.size(); i++)
            {
                if(newnodelist[i]->x<p0.re) p0.re = newnodelist[i]->x;
                if(newnodelist[i]->x>p1.re) p1.re = newnodelist[i]->x;
                if(newnodelist[i]->y<p0.im) p0.im = newnodelist[i]->y;
                if(newnodelist[i]->y>p1.im) p1.im = newnodelist[i]->y;
            }
            d = abs(p1-p0)*CLOSE_ENOUGH;
        }
    }

    // put in all of the nodes;
    for (const auto &node: newnodelist)
    {
        addNode(node->clone(), d);
    }

    // put in all of the lines;
    for (const auto &line: newlinelist)
    {
        // Note(ZaJ): the original code uses a variant of AddSegment that uses coordinates instead of node indices.
        // As far as I can tell, we can use our "default" implementation that uses node indices,
        // as long as we search the correct endpoints using ClosestNode().
        // This is necessary, because calling addNode to copy the nodes from the newnodelist does
        // not necessarily preserve all nodes and therefore inices into the newnodelist are not always
        // valid at this point.

        CComplex p0 (newnodelist[line->n0]->x, newnodelist[line->n0]->y);
        CComplex p1 (newnodelist[line->n1]->x, newnodelist[line->n1]->y);
        // using the raw pointer is ok here, because AddSegment creates a copy anyways
        addSegment(ClosestNode(p0.re,p0.im), ClosestNode(p1.re,p1.im), line.get(), d);
    }

    // put in all of the arcs;
    for (const auto &arc: newarclist)
    {
        // Note(ZaJ): the original code uses a variant of AddArcSegment that uses coordinates instead of the node indices.
        // As far as I can tell, we can use our "default" implementation that uses node indices,
        // as long as we search the correct endpoints using ClosestNode().
        // This is necessary, because calling addNode to copy the nodes from the newnodelist does
        // not necessarily preserve all nodes and therefore inices into the newnodelist are not always
        // valid at this point.

        CComplex p0 (newnodelist[arc->n0]->x, newnodelist[arc->n0]->y);
        CComplex p1 (newnodelist[arc->n1]->x, newnodelist[arc->n1]->y);
        arc->n0 = closestNode(p0.re,p0.im);
        arc->n1 = closestNode(p1.re,p1.im);
        // using the raw pointer is ok here, because AddArcSegment creates a copy anyways
        addArcSegment(*arc.get(), d);
    }

    // put in all of the block labels;
    for (auto &label: newlabellist)
    {
        addBlockLabel(std::move(label), d);
    }

    unselectAll();
}



int femm::FemmProblem::getArcArcIntersection(const femm::CArcSegment &arc0, const femm::CArcSegment &arc1, CComplex *p) const
{
    CComplex a0,a1,t,c0,c1;
    double d,l,R0,R1,z0,z1,c,tta0,tta1;
    int i=0;

    a0.Set(nodelist[arc0.n0]->x,nodelist[arc0.n0]->y);
    a1.Set(nodelist[arc1.n0]->x,nodelist[arc1.n0]->y);

    getCircle(arc1,c1,R1);
    getCircle(arc0,c0,R0);

    d=abs(c1-c0);			// distance between centers

    if ((d>R0+R1) || (d<1.e-08)) return 0;
    // directly eliminate case where there can't
    // be any crossings....

    l=sqrt((R0+R1-d)*(d+R0-R1)*(d-R0+R1)*(d+R0+R1))/(2.*d);
    c=1.+(R0/d)*(R0/d)-(R1/d)*(R1/d);
    t=(c1-c0)/d;
    tta0=arc0.ArcLength*PI/180;
    tta1=arc1.ArcLength*PI/180;

    p[i]=c0 + (c*d/2.+ I*l)*t;		// first possible intersection;
    z0=arg((p[i]-c0)/(a0-c0));
    z1=arg((p[i]-c1)/(a1-c1));
    if ((z0>0.) && (z0<tta0) && (z1>0.) && (z1<tta1)) i++;

    if(fabs(d-R0+R1)/(R0+R1)< 1.e-05)
    {
        p[i]=c0+ c*d*t/2.;
        return i;
    }

    p[i]=c0 + (c*d/2.-I*l)*t;		// second possible intersection
    z0=arg((p[i]-c0)/(a0-c0));
    z1=arg((p[i]-c1)/(a1-c1));
    if ((z0>0.) && (z0<tta0) && (z1>0.) && (z1<tta1)) i++;

    // returns the number of valid intersections found;
    // intersections are returned in the array p[];
    return i;
}

bool femm::FemmProblem::getBoundingBox(double (&x)[2], double (&y)[2]) const
{
    if (nodelist.size()<2)
        return false;

    // initial estimate
    x[0]=nodelist[0]->x;  x[1]=nodelist[0]->x;
    y[0]=nodelist[0]->y;  y[1]=nodelist[0]->y;
    // expand to fit
    for (const auto &node: nodelist)
    {
        if(node->x<x[0]) x[0]=node->x;
        if(node->x>x[1]) x[1]=node->x;
        if(node->y<y[0]) y[0]=node->y;
        if(node->y>y[1]) y[1]=node->y;
    }

    for (const auto &label: labellist)
    {
        if(label->x<x[0]) x[0]=label->x;
        if(label->x>x[1]) x[1]=label->x;
        if(label->y<y[0]) y[0]=label->y;
        if(label->y>y[1]) y[1]=label->y;
    }

    // arcs can "curve" outside of the bounding box defined by their endpoints
    for (const auto &arc: arclist)
    {
        int k=(int) ceil(arc->ArcLength/arc->MaxSideLength);
        double dt = arc->ArcLength*PI/(((double) k)*180.);
        CComplex c;
        double R;
        getCircle(*arc,c,R);
        CComplex p(nodelist[arc->n0]->x, nodelist[arc->n0]->y);
        CComplex s=exp(I*dt);
        for(int j=0; j<k; j++)
        {
            p=(p-c)*s+c;
            if(p.re<x[0]) x[0]=p.re;
            if(p.re>x[1]) x[1]=p.re;
            if(p.im<y[0]) y[0]=p.im;
            if(p.im>y[1]) y[1]=p.im;
        }
    }

    return true;
}

// identical in fmesher, FPProc and HPProc
void femm::FemmProblem::getCircle(const femm::CArcSegment &arc, CComplex &c, double &R) const
{
    // construct the coordinates of the two points on the circle
    CComplex a0(nodelist[arc.n0]->x,nodelist[arc.n0]->y);
    CComplex a1(nodelist[arc.n1]->x,nodelist[arc.n1]->y);

    // calculate distance between arc endpoints
    double d = abs(a1 - a0);

    // figure out what the radius of the circle is...

    // get unit vector pointing from a0 to a1
    CComplex t = (a1 - a0) / d;

    // convert swept angle from degrees to radians
    double tta = arc.ArcLength * PI / 180.;

    // the radius is half the chord length divided by sin of
    // half the swept angle (finding the side length of a
    // triangle formed by the two points and the centre)
    R = d / (2.*sin(tta/2.));

    // center of the arc segment's circle
    c = a0 + (d/2. + I * sqrt(R*R - d*d / 4.)) * t;
}

std::string femm::FemmProblem::getTitle() const
{
    return pathName;
}

bool femm::FemmProblem::getIntersection(int n0, int n1, int segm, double *xi, double *yi) const
{
    CComplex p0,p1,q0,q1;
    double ee,x,z;

    // Check to see if the two lines have a common endpoint
    // If they do, there can be no other intersection...
    if (n0==linelist[segm]->n0) return false;
    if (n0==linelist[segm]->n1) return false;
    if (n1==linelist[segm]->n0) return false;
    if (n1==linelist[segm]->n1) return false;

    // Get a definition of "real small" based on the lengths
    // of the lines of interest;
    p0=nodelist[linelist[segm]->n0]->CC();
    p1=nodelist[linelist[segm]->n1]->CC();
    q0=nodelist[n0]->CC();
    q1=nodelist[n1]->CC();
    ee = (std::min)(abs(p1-p0),abs(q1-q0))*1.0e-8;

    // Rotate and scale the prospective line
    q0=(q0-p0)/(p1-p0);
    q1=(q1-p0)/(p1-p0);

    // Check for cases where there is obviously no intersection
    if ((Re(q0)<=0.) && (Re(q1)<=0.)) return false;
    if ((Re(q0)>=1.) && (Re(q1)>=1.)) return false;
    if ((Im(q0)<=0.) && (Im(q1)<=0.)) return false;
    if ((Im(q0)>=0.) && (Im(q1)>=0.)) return false;

    // compute intersection
    z=Im(q0)/Im(q0-q1);

    // check to see if the line segments intersect at a point sufficiently
    // far from the segment endpoints....
    x=Re((1.0 - z)*q0 + z*q1);
    if((x < ee) || (x > (1.0 - ee))) return false;

    // return resulting intersection point
    p0 = (1.0 - z)*nodelist[n0]->CC() + z*nodelist[n1]->CC();
    *xi=Re(p0);
    *yi=Im(p0);

    return true;
}

int femm::FemmProblem::getLineArcIntersection(const femm::CSegment &seg, const femm::CArcSegment &arc, CComplex *p) const
{
    CComplex p0,p1,a0,a1,t,v,c;
    double d,l,R,z,tta;
    int i=0;

    p0.Set(nodelist[seg.n0]->x,nodelist[seg.n0]->y);
    p1.Set(nodelist[seg.n1]->x,nodelist[seg.n1]->y);
    a0.Set(nodelist[arc.n0]->x,nodelist[arc.n0]->y);
    a1.Set(nodelist[arc.n1]->x,nodelist[arc.n1]->y);
    d=abs(a1-a0);			// distance between arc endpoints

    // figure out what the radius of the circle is...
    t = (a1-a0) / d;
    tta = arc.ArcLength * PI / 180.;
    R = d / ( 2. * std::sin(tta / 2.));
    c = a0 + (d/2. + I*sqrt(R*R-d*d/4.))*t; // center of the arc segment's circle...

    // figure out the distance between line and circle's center;
    d=abs(p1-p0);
    t=(p1-p0)/d;
    v=(c-p0)/t;
    if (fabs(Im(v))>R) return 0;
    l=sqrt( R*R - Im(v)*Im(v));	// Im(v) is distance between line and center...

    if ((l/R) < 1.e-05) 		// case where line is very close to a tangent;
    {
        p[i]=p0 + Re(v)*t;		// make it be a tangent.
        R=Re((p[i]-p0)/t);
        z=arg((p[i]-c)/(a0-c));
        if ((R>0) && (R<d) && (z>0.) && (z<tta)) i++;
        return i;
    }

    p[i]=p0 + (Re(v)+l)*t;		// first possible intersection;
    R=Re((p[i]-p0)/t);
    z=arg((p[i]-c)/(a0-c));
    if ((R>0) && (R<d) && (z>0.) && (z<tta)) i++;

    p[i]=p0 + (Re(v)-l)*t;		// second possible intersection
    R=Re((p[i]-p0)/t);
    z=arg((p[i]-c)/(a0-c));
    if ((R>0) && (R<d) && (z>0.) && (z<tta)) i++;

    // returns the number of valid intersections found;
    // intersections are returned in the array p[];
    return i;
}



double femm::FemmProblem::lengthOfLine(int i) const
{
    return lengthOfLine(*linelist[i]);
}

double femm::FemmProblem::lengthOfLine(const femm::CSegment &seg) const
{
    return abs(nodelist[seg.n0]->CC()-
            nodelist[seg.n1]->CC());
}

void femm::FemmProblem::mirrorCopy(double x0, double y0, double x1, double y1, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    CComplex x=x0 + I*y0;
    CComplex p=(x1-x0) + I*(y1-y0);
    if(abs(p)==0)
        return;
    p/=abs(p);

    if (selector==EditMode::EditNodes || selector == EditMode::EditGroup)
    {
        for (const auto &node: nodelist)
        {
            if (node->IsSelected)
            {
                CComplex y (node->x,node->y);
                y = (y-x) / p;
                y = p*y.Conj()+x;
                // create copy
                std::unique_ptr<CNode> newnode = node->clone();
                // overwrite coordinates in copy
                newnode->x = y.re;
                newnode->y = y.im;
                newnode->IsSelected = false;
                nodelist.push_back(std::move(newnode));
            }
        }
    }
    if (selector == EditMode::EditLines || selector == EditMode::EditGroup)
    {
        for (const auto &line: linelist)
        {
            if (line->IsSelected)
            {
                // copy endpoints
                std::unique_ptr<CNode> n0 = nodelist[line->n0]->clone();
                CComplex y0 (n0->x,n0->y);
                y0 = (y0-x) / p;
                y0 = p*y0.Conj()+x;
                n0->x = y0.re;
                n0->y = y0.im;
                n0->IsSelected = false;

                std::unique_ptr<CNode> n1 = nodelist[line->n1]->clone();
                CComplex y1 (n1->x,n1->y);
                y1 = (y1-x) / p;
                y1 = p*y1.Conj()+x;
                n1->x = y1.re;
                n1->y = y1.im;
                n1->IsSelected = false;

                // copy line (with identical endpoints)
                std::unique_ptr<CSegment> newline = line->clone();
                newline->IsSelected = false;
                // set endpoints
                newline->n0 = (int)nodelist.size();
                nodelist.push_back(std::move(n0));
                newline->n1 = (int)nodelist.size();
                nodelist.push_back(std::move(n1));
                linelist.push_back(std::move(newline));
            }
        }
    }

    if (selector == EditMode::EditLabels || selector == EditMode::EditGroup)
    {
        for (const auto &label: labellist)
        {
            if (label->IsSelected)
            {
                std::unique_ptr<CBlockLabel> newlabel = label->clone();
                CComplex y (label->x,label->y);
                y = (y-x) / p;
                y = p*y.Conj()+x;
                newlabel->x = y.re;
                newlabel->y = y.im;
                newlabel->IsSelected = false;
                // set specific problem parameters:
                if (CMBlockLabel *ptr=dynamic_cast<CMBlockLabel*>(newlabel.get()))
                    ptr->MagDir = (180./PI)*arg(p*conj(exp(I*ptr->MagDir*PI/180.)/p));

                labellist.push_back(std::move(newlabel));
            }
        }
    }
    if (selector == EditMode::EditArcs || selector == EditMode::EditGroup)
    {
        for (const auto &arc: arclist)
        {
            if (arc->IsSelected)
            {
                // copy endpoints
                std::unique_ptr<CNode> n0 = nodelist[arc->n0]->clone();
                CComplex y0 (n0->x,n0->y);
                y0 = (y0-x) / p;
                y0 = p*y0.Conj()+x;
                n0->x = y0.re;
                n0->y = y0.im;
                n0->IsSelected = false;

                std::unique_ptr<CNode> n1 = nodelist[arc->n1]->clone();
                CComplex y1 (n1->x,n1->y);
                y1 = (y1-x) / p;
                y1 = p*y1.Conj()+x;
                n1->x = y1.re;
                n1->y = y1.im;
                n1->IsSelected = false;

                // copy arc (with identical endpoints)
                std::unique_ptr<CArcSegment> newarc = MAKE_UNIQUE<CArcSegment>(*arc);
                newarc->IsSelected = false;
                // set endpoints
                newarc->n0 = (int)nodelist.size();
                nodelist.push_back(std::move(n0));
                newarc->n1 = (int)nodelist.size();
                nodelist.push_back(std::move(n1));
                arclist.push_back(std::move(newarc));
            }
        }
    }
    enforcePSLG();
}

void femm::FemmProblem::rotateCopy(CComplex c, double dt, int ncopies, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    for(int nc=0; nc<ncopies; nc++)
    {
        // accumulated angle
        double t = ((double) (nc+1))*dt;

        CComplex z = exp(I*t*PI/180);

        if (selector==EditMode::EditNodes || selector == EditMode::EditGroup)
        {
            for (const auto &node: nodelist)
            {
                if (node->IsSelected)
                {
                    CComplex x (node->x, node->y);
                    x=(x-c)*z+c;
                    // create copy
                    std::unique_ptr<CNode> newnode = node->clone();
                    // overwrite coordinates in copy
                    newnode->x = x.re;
                    newnode->y = x.im;
                    newnode->IsSelected = false;
                    nodelist.push_back(std::move(newnode));
                }
            }
        }

        if (selector == EditMode::EditLines || selector == EditMode::EditGroup)
        {
            for (const auto &line: linelist)
            {
                if (line->IsSelected)
                {
                    // copy endpoints
                    std::unique_ptr<CNode> n0 = nodelist[line->n0]->clone();
                    CComplex x0 (n0->x,n0->y);
                    x0 = (x0-c)*z+c;
                    n0->x = x0.re;
                    n0->y = x0.im;
                    n0->IsSelected = false;

                    std::unique_ptr<CNode> n1 = nodelist[line->n1]->clone();
                    CComplex x1 (n1->x,n1->y);
                    x1 = (x1-c)*z+c;
                    n1->x = x1.re;
                    n1->y = x1.im;
                    n1->IsSelected = false;

                    // copy line (with identical endpoints)
                    std::unique_ptr<CSegment> newline = line->clone();
                    newline->IsSelected = false;
                    // set endpoints
                    newline->n0 = (int)nodelist.size();
                    nodelist.push_back(std::move(n0));
                    newline->n1 = (int)nodelist.size();
                    nodelist.push_back(std::move(n1));
                    linelist.push_back(std::move(newline));
                }
            }
        }

        if (selector == EditMode::EditArcs || selector == EditMode::EditGroup)
        {
            for (const auto &arc: arclist)
            {
                if (arc->IsSelected)
                {
                    // copy endpoints
                    std::unique_ptr<CNode> n0 = nodelist[arc->n0]->clone();
                    CComplex x0 (n0->x,n0->y);
                    x0 = (x0-c)*z+c;
                    n0->x = x0.re;
                    n0->y = x0.im;
                    n0->IsSelected = false;

                    std::unique_ptr<CNode> n1 = nodelist[arc->n1]->clone();
                    CComplex x1 (n1->x,n1->y);
                    x1 = (x1-c)*z+c;
                    n1->x = x1.re;
                    n1->y = x1.im;
                    n1->IsSelected = false;

                    // copy arc (with identical endpoints)
                    std::unique_ptr<CArcSegment> newarc = MAKE_UNIQUE<CArcSegment>(*arc);
                    newarc->IsSelected = false;
                    // set endpoints
                    newarc->n0 = (int)nodelist.size();
                    nodelist.push_back(std::move(n0));
                    newarc->n1 = (int)nodelist.size();
                    nodelist.push_back(std::move(n1));
                    arclist.push_back(std::move(newarc));
                }
            }
        }

        if (selector == EditMode::EditLabels || selector == EditMode::EditGroup)
        {
            for (const auto &label: labellist)
            {
                if (label->IsSelected)
                {
                    std::unique_ptr<CBlockLabel> newlabel = label->clone();
                    CComplex x(label->x,label->y);
                    x = (x-c)*z+c;
                    newlabel->x = x.re;
                    newlabel->y = x.im;
                    newlabel->IsSelected = false;

                    if (CMBlockLabel *ptr=dynamic_cast<CMBlockLabel*>(newlabel.get()))
                        ptr->MagDir += t;

                    labellist.push_back(std::move(newlabel));
                }
            }
        }
    }

    enforcePSLG();
}

void femm::FemmProblem::rotateMove(CComplex c, double t, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    bool processNodes = (selector == EditMode::EditNodes);

    const CComplex z = exp(I*t*PI/180);

    if(selector==EditMode::EditLines || selector==EditMode::EditGroup)
    {
        for (const auto &line: linelist)
        {
            if (line->IsSelected)
            {
                nodelist[line->n0]->IsSelected = true;
                nodelist[line->n1]->IsSelected = true;
            }
        }
        processNodes = true;
    }

    if(selector==EditMode::EditArcs || selector==EditMode::EditGroup)
    {
        for (const auto &arc: arclist)
        {
            if (arc->IsSelected)
            {
                nodelist[arc->n0]->IsSelected = true;
                nodelist[arc->n1]->IsSelected = true;
            }
        }
        processNodes = true;
    }

    if(selector==EditMode::EditLabels || selector==EditMode::EditGroup)
    {
        for (auto &label: labellist)
        {
            if (label->IsSelected)
            {
                CComplex x (label->x, label->y);
                x = (x-c)*z+c;
                label->x = x.re;
                label->y = x.im;

                // only relevant to magnetics problems:
                if (CMBlockLabel *ptr=dynamic_cast<CMBlockLabel*>(label.get()))
                    ptr->MagDir += t;
            }
        }
    }

    if(processNodes)
    {
        for (auto &node : nodelist)
        {
            if (node->IsSelected)
            {
                CComplex x(node->x,node->y);
                x = (x-c)*z+c;
                node->x = x.re;
                node->y = x.im;
            }
        }
    }
    enforcePSLG();
}

void femm::FemmProblem::scaleMove(double bx, double by, double sf, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    bool processNodes = (selector == EditMode::EditNodes);

    if (selector==EditMode::EditLines || selector==EditMode::EditGroup)
    {
        for (const auto& line: linelist)
        {
            if (line->IsSelected)
            {
                nodelist[line->n0]->IsSelected = true;
                nodelist[line->n1]->IsSelected = true;
            }
        }
        processNodes = true;
    }

    if (selector==EditMode::EditArcs || selector==EditMode::EditGroup)
    {
        for (const auto &arc: arclist)
        {
            if (arc->IsSelected)
            {
                nodelist[arc->n0]->IsSelected = true;
                nodelist[arc->n1]->IsSelected = true;
            }
        }
        processNodes = true;
    }

    if (selector==EditMode::EditLabels || selector==EditMode::EditGroup)
    {
        for (auto &label: labellist)
        {
            if (label->IsSelected)
            {
                label->x = bx+sf*(label->x - bx);
                label->y = by+sf*(label->y - by);
                label->MaxArea *= (sf*sf);
            }
        }
    }

    if (processNodes)
    {
        for (auto &node : nodelist)
        {
            if (node->IsSelected)
            {
                node->x = bx+sf*(node->x - bx);
                node->y = by+sf*(node->y - by);
            }
        }
    }
    enforcePSLG();
}









// identical in fmesher, hpproc, fpproc
double femm::FemmProblem::shortestDistanceFromSegment(double p, double q, int segm) const
{
    double x0=nodelist[linelist[segm]->n0]->x;
    double y0=nodelist[linelist[segm]->n0]->y;
    double x1=nodelist[linelist[segm]->n1]->x;
    double y1=nodelist[linelist[segm]->n1]->y;

    double t=((p-x0)*(x1-x0) + (q-y0)*(y1-y0))/
            ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

    if (t>1.) t=1.;
    if (t<0.) t=0.;

    double x2=x0+t*(x1-x0);
    double y2=y0+t*(y1-y0);

    return sqrt((p-x2)*(p-x2) + (q-y2)*(q-y2));
}



// identical in fmesher, FPProc and HPProc
double femm::FemmProblem::shortestDistanceFromArc(CComplex p, const femm::CArcSegment &arc) const
{
    double R;
    CComplex c;
    getCircle(arc,c,R);

    double d=abs(p-c);
    if(d==0) return R;

    CComplex a0(nodelist[arc.n0]->x,nodelist[arc.n0]->y);
    CComplex a1(nodelist[arc.n1]->x,nodelist[arc.n1]->y);
    CComplex t=(p-c)/d;
    double l=abs(p-c-R*t);
    double z=arg(t/(a0-c))*180/PI;
    if ((z>0) && (z<arc.ArcLength)) return l;

    z=abs(p-a0);
    l=abs(p-a1);
    if(z<l) return z;
    return l;
}

void femm::FemmProblem::translateCopy(double incx, double incy, int ncopies, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    for(int nc=0; nc<ncopies; nc++)
    {
        // accumulated offsets
        double dx=((double)(nc+1))*incx;
        double dy=((double)(nc+1))*incy;

        if (selector==EditMode::EditNodes || selector == EditMode::EditGroup)
        {
            for (const auto &node: nodelist)
            {
                if (node->IsSelected)
                {
                    // create copy
                    std::unique_ptr<CNode> newnode = node->clone();
                    // overwrite coordinates in copy
                    newnode->x += dx;
                    newnode->y += dy;
                    newnode->IsSelected = false;
                    nodelist.push_back(std::move(newnode));
                }
            }
        }

        if (selector == EditMode::EditLines || selector == EditMode::EditGroup)
        {
            for (const auto &line: linelist)
            {
                if (line->IsSelected)
                {
                    // copy endpoints
                    std::unique_ptr<CNode> n0 = nodelist[line->n0]->clone();
                    n0->x += dx;
                    n0->y += dy;
                    n0->IsSelected = false;

                    std::unique_ptr<CNode> n1 = nodelist[line->n1]->clone();
                    n1->x += dx;
                    n1->y += dy;
                    n1->IsSelected = false;

                    // copy line (with identical endpoints)
                    std::unique_ptr<CSegment> newline = line->clone();
                    newline->IsSelected = false;
                    // set endpoints
                    newline->n0 = (int)nodelist.size();
                    nodelist.push_back(std::move(n0));
                    newline->n1 = (int)nodelist.size();
                    nodelist.push_back(std::move(n1));
                    linelist.push_back(std::move(newline));
                }
            }
        }

        if (selector == EditMode::EditLabels || selector == EditMode::EditGroup)
        {
            for (const auto &label: labellist)
            {
                if (label->IsSelected)
                {
                    std::unique_ptr<CBlockLabel> newlabel = label->clone();
                    newlabel->x += dx;
                    newlabel->y += dy;
                    newlabel->IsSelected = false;

                    labellist.push_back(std::move(newlabel));
                }
            }
        }

        if (selector == EditMode::EditArcs || selector == EditMode::EditGroup)
        {
            for (const auto &arc: arclist)
            {
                if (arc->IsSelected)
                {
                    // copy endpoints
                    std::unique_ptr<CNode> n0 = nodelist[arc->n0]->clone();
                    n0->x += dx;
                    n0->y += dy;
                    n0->IsSelected = false;

                    std::unique_ptr<CNode> n1 = nodelist[arc->n1]->clone();
                    n1->x += dx;
                    n1->y += dy;
                    n1->IsSelected = false;

                    // copy arc (with identical endpoints)
                    std::unique_ptr<CArcSegment> newarc = MAKE_UNIQUE<CArcSegment>(*arc);
                    newarc->IsSelected = false;
                    // set endpoints
                    newarc->n0 = (int)nodelist.size();
                    nodelist.push_back(std::move(n0));
                    newarc->n1 = (int)nodelist.size();
                    nodelist.push_back(std::move(n1));
                    arclist.push_back(std::move(newarc));
                }
            }
        }
    }

    enforcePSLG();
}

void femm::FemmProblem::translateMove(double dx, double dy, femm::EditMode selector)
{
    assert(selector != EditMode::Invalid);
    bool processNodes = (selector == EditMode::EditNodes);

    if (selector == EditMode::EditLines || selector == EditMode::EditGroup)
    {
        // select end points of selected lines:
        for (auto &line: linelist)
        {
            if (line->IsSelected)
            {
                nodelist[line->n0]->IsSelected = true;
                nodelist[line->n1]->IsSelected = true;
            }
        }
        // make sure to translate endpoints
        processNodes = true;
    }
    if (selector == EditMode::EditArcs || selector == EditMode::EditGroup)
    {
        // select end points of selected arcs:
        for (auto &arc: arclist)
        {
            if (arc->IsSelected)
            {
                nodelist[arc->n0]->IsSelected = true;
                nodelist[arc->n1]->IsSelected = true;
            }
        }
        // make sure to translate endpoints
        processNodes = true;
    }

    if (selector == EditMode::EditLabels || selector == EditMode::EditGroup)
    {
        for (auto &lbl: labellist)
        {
            if (lbl->IsSelected)
            {
                lbl->x += dx;
                lbl->y += dy;
            }
        }
    }
    if (processNodes)
    {
        for (auto &node: nodelist)
        {
            if (node->IsSelected)
            {
                node->x += dx;
                node->y += dy;
            }
        }
    }
    enforcePSLG();
}

int femm::FemmProblem::ClosestNode(const double x, const double y) const
{
    int i,j;
    double d0,d1;

    if(nodelist.size()==0) return -1;

    j=0;
    d0=nodelist[0]->GetDistance(x,y);
    for(i=0; i<(int)nodelist.size(); i++)
    {
        d1=nodelist[i]->GetDistance(x,y);
        if(d1<d0)
        {
            d0=d1;
            j=i;
        }
    }

    return j;
}


int femm::FemmProblem::ClosestArcSegment(double x, double y) const
{
    double d0,d1;
    int i,j;

    if(arclist.size()==0) return -1;

    j=0;
    d0=ShortestDistanceFromArc(CComplex(x,y),*(arclist[0]));
    for(i=0; i<(int)arclist.size(); i++)
    {
        d1=ShortestDistanceFromArc(CComplex(x,y),*(arclist[i]));
        if(d1<d0)
        {
            d0=d1;
            j=i;
        }
    }

    return j;
}

void femm::FemmProblem::GetCircle(const CArcSegment &arc, CComplex &c, double &R) const
{
    CComplex a0,a1,t;
    double d,tta;

    a0.Set(nodelist[arc.n0]->x, nodelist[arc.n0]->y);
    a1.Set(nodelist[arc.n1]->x, nodelist[arc.n1]->y);
    d=abs(a1-a0);            // distance between arc endpoints

    // figure out what the radius of the circle is...
    t=(a1-a0)/d;
    tta=arc.ArcLength*PI/180.;
    R=d/(2.*sin(tta/2.));
    c=a0 + (d/2. + I*sqrt(R*R-d*d/4.))*t; // center of the arc segment's circle...
}

double femm::FemmProblem::ShortestDistanceFromArc(const CComplex p, const CArcSegment &arc) const
{
    double R,d,l,z;
    CComplex a0,a1,c,t;

    a0.Set(nodelist[arc.n0]->x,nodelist[arc.n0]->y);
    a1.Set(nodelist[arc.n1]->x,nodelist[arc.n1]->y);
    GetCircle(arc,c,R);
    d=abs(p-c);

    if(d==0) return R;

    t=(p-c)/d;
    l=abs(p-c-R*t);
    z=arg(t/(a0-c))*180/PI;
    if ((z>0) && (z<arc.ArcLength)) return l;

    z=abs(p-a0);
    l=abs(p-a1);
    if(z<l) return z;
    return l;
}

double femm::FemmProblem::ShortestDistanceFromSegment(double p, double q, int segm) const
{
    double t,x[3],y[3];

    x[0]=nodelist[linelist[segm]->n0]->x;
    y[0]=nodelist[linelist[segm]->n0]->y;
    x[1]=nodelist[linelist[segm]->n1]->x;
    y[1]=nodelist[linelist[segm]->n1]->y;

    t=((p-x[0])*(x[1]-x[0]) + (q-y[0])*(y[1]-y[0]))/
      ((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]));

    if (t>1.) t=1.;
    if (t<0.) t=0.;

    x[2]=x[0]+t*(x[1]-x[0]);
    y[2]=y[0]+t*(y[1]-y[0]);

    return sqrt((p-x[2])*(p-x[2]) + (q-y[2])*(q-y[2]));
}


void femm::FemmProblem::unselectAll()
{
    for(auto &node: nodelist) node->IsSelected = false;
    for(auto &line: linelist) line->IsSelected = false;
    for(auto &label: labellist) label->IsSelected = false;
    for(auto &arc: arclist) arc->IsSelected = false;
}

void femm::FemmProblem::undo()
{
    for(int i=0; i<(int)undolinelist.size(); i++)
        linelist[i].swap(undolinelist[i]);
    for(int i=0; i<(int)undoarclist.size(); i++)
        arclist[i].swap(undoarclist[i]);
    for(int i=0; i<(int)undolabellist.size(); i++)
        labellist[i].swap(undolabellist[i]);
    for(int i=0; i<(int)undonodelist.size(); i++)
        nodelist[i].swap(undonodelist[i]);
}

void femm::FemmProblem::undoLines()
{
    for(int i=0; i<(int)undolinelist.size(); i++)
        linelist[i].swap(undolinelist[i]);
}

void femm::FemmProblem::undoArcs()
{
	for(int i=0;i<(int)arclist.size();i++)
	{
		arclist[i]->mySideLength=arclist[i]->MaxSideLength;
		arclist[i]->MaxSideLength=undoarclist[i]->MaxSideLength;
	}
}

void femm::FemmProblem::updateUndo()
{
    undonodelist.clear();
    undolinelist.clear();
    undoarclist.clear();
    undolabellist.clear();

    // copy each entry
    for(const auto& node: nodelist)
        undonodelist.push_back(node->clone());
    for(const auto& line: linelist)
        undolinelist.push_back(line->clone());
    for(const auto& arc: arclist)
        undoarclist.push_back(MAKE_UNIQUE<CArcSegment>(*arc));
    for(const auto& label: labellist)
        undolabellist.push_back(label->clone());
}

femm::FemmProblem::FemmProblem(FileType ftype)
    : FileFormat(-1)
    , Frequency(0.0)
    , Precision(1.e-08)
    , MinAngle(DEFAULT_MINIMUM_ANGLE)
    , Depth(1)
    , LengthUnits(LengthInches)
    , Coords(CART)
    , problemType(PLANAR)
    , extZo(0)
    , extRo(0)
    , extRi(0)
    , comment()
    , ACSolver(0)
    , dT(0)
    , previousSolutionFile()
    , PrevType(0)
    , DoForceMaxMeshArea(false)
    , DoSmartMesh(true)
    , nodelist()
    , linelist()
    , arclist()
    , labellist()
    , agelist()
    , nodeproplist()
    , lineproplist()
    , blockproplist()
    , circproplist()
    , pathName()
    //, solutionFile()
    , filetype(ftype)
    , blockMap()
    , circuitMap()
    , d_EditMode( EditMode::Invalid )
    , undonodelist()
    , undolinelist()
    , undoarclist()
    , undolabellist()
{}
