
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
#include "PostProcessor.h"

#include "femmcomplex.h"
#include "femmconstants.h"
#include "fparse.h"
#include "linsolve/backend_factory.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

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

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//using namespace std;
using namespace femm;

PProcIface::~PProcIface()
{
}

PProcIface::PProcIface()
{
}

/**
 * Constructor for the PostProcessor class.
 */
femm::PostProcessor::PostProcessor()
{
    // set some default values for problem definition
    d_LineIntegralPoints = 400;
    Smooth = true;
    NumList = nullptr;
    ConList = nullptr;
    bHasMask = false;
    LengthConv = (double *)calloc(6,sizeof(double));
    LengthConv[0] = 0.0254;   //inches
    LengthConv[1] = 0.001;    //millimeters
    LengthConv[2] = 0.01;     //centimeters
    LengthConv[3] = 1.;       //meters
    LengthConv[4] = 2.54e-05; //mils
    LengthConv[5] = 1.e-06;   //micrometers

    for(int i=0; i<9; i++)
    {
        d_PlotBounds[i][0] = d_PlotBounds[i][1] =
                PlotBounds[i][0] = PlotBounds[i][1] = 0;
    }

    // initialise the warning message function pointer to
    // point to the PrintWarningMsg function
    WarnMessage = &PrintWarningMsg;
}

/**
 * Destructor for the PostProcessor class.
 */
femm::PostProcessor::~PostProcessor()
{
    free(LengthConv);
    if (ConList)
    {
        for(int i=0; i<(int)meshnodes.size(); i++)
            if(ConList[i]!=NULL) free(ConList[i]);
        free(ConList);
    }
    if (NumList)
        free(NumList);
}

int PostProcessor::numElements() const
{
    return (int) meshelems.size();
}

int PostProcessor::numNodes() const
{
    return (int) meshnodes.size();
}

bool PostProcessor::selectBlocklabel(double px, double py)
{
    int idx = InTriangle(px,py);
    if (idx >= 0)
    {

        bHasMask = false;
        problem->labellist[meshelems[idx]->lbl]->ToggleSelect();
        return true;
    }
    return false;
}

void PostProcessor::selectConductor(int idx)
{
    for (auto &node: problem->nodelist)
        if (idx == node->InConductor)
            node->ToggleSelect();
    for (auto &line: problem->linelist)
        if (idx == line->InConductor)
            line->ToggleSelect();
    for (auto &arc: problem->arclist)
        if (idx == arc->InConductor)
            arc->ToggleSelect();
    for (auto &mnode: meshnodes)
    {
        mnode->ToggleSelect();
    }
}

void PostProcessor::setSmoothing(bool value)
{
    Smooth = value;
}

void PostProcessor::setMessageCallback(PostProcessor::MessageCB msgFun)
{
    if (msgFun)
        WarnMessage = msgFun;
}

void PostProcessor::addContourPoint(CComplex p)
{
    if (contour.empty() || p!=contour.back())
        contour.push_back(p);
}

void PostProcessor::addContourPointFromNode(double mx, double my)
{
    if (!problem->nodelist.empty())
    {
        int n0=problem->closestNode(mx,my);

        int lineno=-1;
        int arcno=-1;
        CComplex z(problem->nodelist[n0]->x,problem->nodelist[n0]->y);
        if (contour.empty())
        {
            contour.push_back(z);
            //theView->DrawUserContour(FALSE);
            return;
        }
        //check to see if point is the same as last point in the contour;
        CComplex y = contour.back();

        if (y==z)
            return;

        int n1 = problem->closestNode(y.re,y.im);
        CComplex x(problem->nodelist[n1]->x,problem->nodelist[n1]->y);

        //check to see if this point and the last point are ends of an
        //input segment;
        double d1=1.e08;

        if (abs(x-y)<1.e-08)
        {
            for(int k=0; k<(int)problem->linelist.size(); k++)
            {
                if((problem->linelist[k]->n0==n1) && (problem->linelist[k]->n1==n0))
                {
                    double d2=fabs(problem->shortestDistanceFromSegment(mx,my,k));
                    if(d2<d1)
                    {
                        lineno=k;
                        d1=d2;
                    }
                }
                if((problem->linelist[k]->n0==n0) && (problem->linelist[k]->n1==n1))
                {
                    double d2=fabs(problem->shortestDistanceFromSegment(mx,my,k));
                    if(d2<d1){
                        lineno=k;
                        d1=d2;
                    }
                }
            }
        }
        bool reverseOrder = false;
        //check to see if this point and last point are ends of an
        // arc segment; if so, add entire arc to the contour;
        if (abs(x-y)<1.e-08)
        {
            for(int k=0;k<(int)problem->arclist.size();k++)
            {
                if((problem->arclist[k]->n0==n1) && (problem->arclist[k]->n1==n0))
                {
                    double d2=problem->shortestDistanceFromArc(CComplex(mx,my),
                                                              *problem->arclist[k].get());
                    if(d2<d1){
                        arcno=k;
                        lineno=-1;
                        reverseOrder=true;
                        d1=d2;
                    }
                }
                if((problem->arclist[k]->n0==n0) && (problem->arclist[k]->n1==n1))
                {
                    double d2=problem->shortestDistanceFromArc(CComplex(mx,my),
                                                              *problem->arclist[k].get());
                    if(d2<d1){
                        arcno=k;
                        lineno=-1;
                        reverseOrder=false;
                        d1=d2;
                    }
                }
            }
        }
        if((lineno<0) && (arcno<0))
        {
            contour.push_back(z);
            //theView->DrawUserContour(FALSE);
        }
        if(lineno>=0)
        {
            int size=(int) contour.size();
            if(size>1)
            {
                if(abs(contour[size-2]-z)<1.e-08)
                    return;
            }
            contour.push_back(z);
            //theView->DrawUserContour(FALSE);
        }
        if(arcno>=0){
            int k=arcno;
            double R;
            problem->getCircle(*problem->arclist[k].get(),x,R);
            int arcsegments=(int) ceil(problem->arclist[k]->ArcLength/problem->arclist[k]->MaxSideLength);
            if(reverseOrder)
                z=exp(I*problem->arclist[k]->ArcLength*PI/(180.*((double) arcsegments)) );
            else
                z=exp(-I*problem->arclist[k]->ArcLength*PI/(180.*((double) arcsegments)) );

            for(int i=0; i<arcsegments; i++)
            {
                y=(y-x)*z+x;
                int size=(int) contour.size();
                if(size>1)
                {
                    if(abs(contour[size-2]-y)<1.e-08)
                        return;
                }
                contour.push_back(y);
                //theView->DrawUserContour(FALSE);
            }
        }
    }

}


// identical in EPProc, FPProc and HPProc
int femm::PostProcessor::InTriangle(double x, double y) const
{
    static int k=0;
    int j,hi,lo,sz;
    double z;

    sz = meshelems.size();

    if ((k < 0) || (k >= sz)) k = 0;

    // In most applications, the triangle we're looking
    // for is nearby the last one we found.  Since the elements
    // are ordered in a banded structure, we want to check the
    // elements nearby the last one selected first.
    if (InTriangleTest(x,y,k)) return k;

    // wasn't in the last searched triangle, so look through all the
    // elements
    hi = k;
    lo = k;

    for(j=0; j<sz; j+=2)
    {
        hi++;
        if (hi >= sz) hi = 0;
        lo--;
        if (lo < 0)   lo = sz - 1;

        CComplex hiCtr = meshelems[hi]->ctr;
        z = (hiCtr.re - x) * (hiCtr.re - x) + (hiCtr.im - y) * (hiCtr.im - y);

        if (z <= meshelems[hi]->rsqr)
        {
            if (InTriangleTest(x,y,hi))
            {
                k = hi;
                return k;
            }
        }

        CComplex loCtr = meshelems[lo]->ctr;
        z = (loCtr.re-x)*(loCtr.re-x) + (loCtr.im-y)*(loCtr.im-y);

        if (z <= meshelems[lo]->rsqr)
        {
            if (InTriangleTest(x,y,lo))
            {
                k = lo;
                return k;
            }
        }

    }

    return (-1);
}

// EPProc  and FPProc are identical
// FPProc and HPProc differ, but I'm not sure whether hpproc could just use this version instead
bool femm::PostProcessor::InTriangleTest(double x, double y, int i) const
{

    if ((i < 0) || (i >= int(meshelems.size()))) return false;

    int j,k;
    double z;

    for (j=0; j<3; j++)
    {
        k = j + 1;

        if (k == 3) k = 0;

        int p_k = meshelems[i]->p[k];
        int p_j = meshelems[i]->p[j];
        // Case 1: p[k]>p[j]
        if (p_k > p_j)
        {
            z = (meshnodes[p_k]->x - meshnodes[p_j]->x) *
                    (y - meshnodes[p_j]->y) -
                    (meshnodes[p_k]->y - meshnodes[p_j]->y) *
                    (x - meshnodes[p_j]->x);

            if(z<0) return false;
        }
        //Case 2: p[k]<p[j]
        else
        {
            z = (meshnodes[p_j]->x - meshnodes[p_k]->x) *
                    (y - meshnodes[p_k]->y) -
                    (meshnodes[p_j]->y - meshnodes[p_k]->y) *
                    (x - meshnodes[p_k]->x);

            if (z > 0) return false;
        }
    }

    return true;
}

bool PostProcessor::isKosher(int k) const
{
    // If:
    //    1) this is an axisymmetric problem;
    //    2) the selected geometry lies along the r=0 axis, and
    //    3) we have a node on the r=0 axis that we are trying to determine
    //     if we should set to zero.
    // This routine determines whether the node is at the extents of
    // the r=0 domain (or lies at a break in some sub-interval).
    //
    // Returns TRUE if it is OK to define the node as zero;

    if((problem->problemType==PLANAR) || (meshnodes[k]->x>1.e-6)) return true;

    int score=0;
    for(int i=0;i<NumList[k];i++)
    {
        for(int j=0;j<3;j++)
        {
            int n=meshelems[ConList[k][i]]->p[j];
            if((n!=k) && (meshnodes[n]->x<1.e-6))
            {
                score++;
                if(score>1)
                    return false;
            }
        }
    }

    return true;
}

bool PostProcessor::isSameMaterial(const femmsolver::CElement &e1, const femmsolver::CElement &e2) const
{
    return (problem->blockproplist[e1.blk]->isSameMaterialAs(problem->blockproplist[e2.blk].get()));
}

// identical in FPProc and HPProc
CComplex femm::PostProcessor::Ctr(int i)
{
    CComplex c = 0;
    for(int j=0; j<3; j++)
    {
        int p_j = meshelems[i]->p[j];
        CComplex p(meshnodes[ p_j ]->x/3., meshnodes[ p_j ]->y/3.);
        c+=p;
    }

    return c;
}

// identical in FPProc and HPProc
double femm::PostProcessor::ElmArea(int i) const
{
    int n[3];
    for(int j=0; j<3; j++) n[j]=meshelems[i]->p[j];

    double b0=meshnodes[n[1]]->y - meshnodes[n[2]]->y;
    double b1=meshnodes[n[2]]->y - meshnodes[n[0]]->y;
    double c0=meshnodes[n[2]]->x - meshnodes[n[1]]->x;
    double c1=meshnodes[n[0]]->x - meshnodes[n[2]]->x;
    return (b0*c1-b1*c0)/2.;
}

const std::vector<std::unique_ptr<femmsolver::CMeshNode> > &PostProcessor::getMeshNodes() const
{
    return meshnodes;
}

const FemmProblem *PostProcessor::getProblem() const
{
    return problem.get();
}

bool PostProcessor::hasMultiplyDefinedLabels() const
{
    return bMultiplyDefinedLabels;
}

bool PostProcessor::isSelectionOnAxis() const
{
    if (problem->problemType!=AXISYMMETRIC)
        return false;

    for (const auto &elem: meshelems)
    {
        if(problem->labellist[elem->lbl]->IsSelected)
        {
            for(int j=0;j<3;j++)
                if(meshnodes[elem->p[j]]->x<1.e-6)
                    return true;
        }
    }
    return false;
}

// femm42 also has a #ifdef'ed simple version of the algorithm, that we don't bother to replicate here
bool PostProcessor::makeMask()
{
    if(bHasMask) return true;

    std::unique_ptr<femm::LinearSystemBackend<double>> L =
        femm::create_backend<double>(femm::default_backend_kind());
    double Me[3][3],be[3];     // element matrix;
    double p[3],q[3];       // element shape parameters;
    int n[3];               // numbers of nodes for a particular element;


    const static int plus1mod3[3] = {1, 2, 0};
    const static int minus1mod3[3] = {2, 0, 1};

    // figure out bandwidth--helps speed somethings up;
    int NumEls=(int) meshelems.size();
    int bw=0;
    for(int i=0;i<NumEls;i++)
    {
        for(int j=0;j<3;j++)
        {
            int k=(j+1) % 3;
            int d=abs(meshelems[i]->p[j]-meshelems[i]->p[k]);
            if (d>bw) bw=d;
        }
    }
    bw++;

    int NumNodes=(int) meshnodes.size();
    L->create(NumNodes,bw);

    // Sort through materials to see if they denote air;
    int *matflag=(int*)calloc(problem->blockproplist.size(),sizeof(int));
    int *lblflag=(int*)calloc(problem->labellist.size(),sizeof(int));
    for(int i=0;i<(int)problem->blockproplist.size();i++)
    {
        // k==0 for air, k==1 for other than air
        matflag[i] = (problem->blockproplist[i]->isAir()) ? 0 : 1;
    }

    // Now, sort through the labels to see which ones correspond to air blocks.
    for(int i=0;i<(int)problem->labellist.size();i++)
    {
        lblflag[i]=matflag[problem->labellist[i]->BlockType];
        // FIXME(ZaJ): for magnetics, we need this here, too:
        // if(blocklist[i].InCircuit>=0) lblflag[i]=1;
    }
    free(matflag);

    // Determine which nodal values should be fixed
    // and what values they should be fixed at;
    for(int i=0;i<NumNodes;i++){
        // Note(ZaJ): I have added the field Q to CMeshNode, and set it to -2 for CMMeshNode
        //            this makes the code here equivalent to the fpproc implementation
        if (meshnodes[i]->Q!=-2) L->solution()[i]=0;
        else L->solution()[i]=-1;
    }

    // if the problem is axisymmetric, does the selection lie along r=0?
    bool bOnAxis=isSelectionOnAxis();

    // Figure out which nodes are exterior edges and set them to zero;
    for(int i=0;i<NumEls;i++)
    {
        for(int j=0;j<3;j++)
        {
            if (meshelems[i]->n[j] == 1)
            {
                int k;
                k=meshelems[i]->p[plus1mod3[j]];
                if((!bOnAxis) || (isKosher(k))) L->solution()[k]=0;
                k=meshelems[i]->p[minus1mod3[j]];
                if((!bOnAxis) || (isKosher(k))) L->solution()[k]=0;
            }
        }
    }

    // Set all nodes in a selected block equal to 1;
    for(int i=0;i<NumEls;i++)
    {
        if(problem->labellist[meshelems[i]->lbl]->IsSelected)
        {
            for(int j=0;j<3;j++){
                L->solution()[meshelems[i]->p[j]]=1;
            }
        }
        else if(lblflag[meshelems[i]->lbl]!=0)
        {
            for(int j=0;j<3;j++) L->solution()[meshelems[i]->p[j]]=0;
        }
    }
    // Any nodes that have point currents applied to them but are not in the
    // selected region should also be set to zero so that they don't mess up
    // the force calculation
    if(problem->nodeproplist.size()>0)
    {
        CComplex *p=(CComplex *)calloc(problem->nodelist.size(),sizeof(CComplex));
        int npts = 0;
        for(int i=0;i<(int)problem->nodelist.size();i++)
            if(problem->nodelist[i]->BoundaryMarker>=0)
            {
                p[npts]=problem->nodelist[i]->CC();
                npts++;
            }

        if(npts>0)
            for(int i=0;i<NumNodes;i++)
                for(int j=0;j<npts;j++)
                    if(abs(p[j]-meshnodes[i]->CC())<1.e-8)
                    {
                        if (L->solution()[i]<0) L->solution()[i]=0.;
                        npts--;
                        if(npts>0){
                            p[j]=p[npts];
                            j=npts;
                        }
                        else{
                            j=npts;
                            i=NumNodes;
                        }
                    }
        free(p);
    }

    // ugly filetype check, but I couldn't think of a nice way to isolate this better...
    // currentflow files need the same treatment...
    if (problem->filetype == femm::FileType::ElectrostaticsFile)
    {
        for(int i=0;i<NumNodes;i++)
        {
            const femmsolver::CSMeshNode *node = reinterpret_cast<femmsolver::CSMeshNode*>(meshnodes[i].get());
            if (node->IsSelected==true) L->solution()[i]=1;
        }
    }
    // build up element matrices;
    for(int i=0;i<NumEls;i++)
    {
        // zero out Me;
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++) Me[j][k]=0;
            be[j]=0;
        }

        // Determine shape parameters.
        // l == element side lengths;
        // p corresponds to the `b' parameter in Allaire
        // q corresponds to the `c' parameter in Allaire

        for(int k=0;k<3;k++) n[k] = meshelems[i]->p[k];
        p[0]=meshnodes[n[1]]->y - meshnodes[n[2]]->y;
        p[1]=meshnodes[n[2]]->y - meshnodes[n[0]]->y;
        p[2]=meshnodes[n[0]]->y - meshnodes[n[1]]->y;
        q[0]=meshnodes[n[2]]->x - meshnodes[n[1]]->x;
        q[1]=meshnodes[n[0]]->x - meshnodes[n[2]]->x;
        q[2]=meshnodes[n[1]]->x - meshnodes[n[0]]->x;

        double area = (p[0]*q[1]-p[1]*q[0])/2.; //element area

        // quick check for consistency--
        // if the block is not air and is not selected,
        // all of the nodes in the block better be defined
        // to be zero;  Otherwise, the region for force
        // integration has been selected in an invalid way;
        if ((!problem->labellist[meshelems[i]->lbl]->IsSelected) && (lblflag[meshelems[i]->lbl]))
        {
            int k=0;
            for(int j=0;j<3;j++) if (L->solution()[n[j]]==0) k++;
            if(k<3){
                std::string outmsg = "The selected region is invalid. A valid selection\n"
                                     "cannot abut a region which is not free space.";
                WarnMessage(outmsg.c_str());
                free(lblflag);
                return false;
            }
        }

        // Each element weighted by its region's
        // mesh size specification;
        double v=problem->labellist[meshelems[i]->lbl]->MaxArea;
        if (v<=0) v=sqrt(area); else v=sqrt(v);

        // build element matrix;
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                Me[j][k]+=v*(p[j]*p[k]+q[j]*q[k])/area;

        // process any prescribed nodal values;
        // doing it here saves a lot of time.
        for(int j=0;j<3;j++)
        {
            if(L->solution()[n[j]]>=0)
            {
                for(int k=0;k<3;k++)
                {
                    if(j!=k){
                        be[k]-=Me[k][j]*L->solution()[n[j]];
                        Me[k][j]=0;
                        Me[j][k]=0;
                    }
                }
                be[j]=L->solution()[n[j]]*Me[j][j];
            }
        }

        // combine block matrices into global matrices;
        for (int j=0;j<3;j++)
        {
            for (int k=j;k<3;k++)
                if(Me[j][k]!=0)
                    L->put(L->get(n[j],n[k])-Me[j][k],n[j],n[k]);
            L->rhs()[n[j]]-=be[j];
        }
    }

    free(lblflag);

    // solve the problem;
    femm::SolveOptions opts;
    opts.warm_start = false;
    if (!L->solve(opts).converged)
        return false;

    // Process the results to get one row of elements
    // that runs down the center of the gap away from boundaries.
    for(int i=0;i<NumNodes;i++)
        if (L->solution()[i]>0.5) meshnodes[i]->msk = 1;
        else meshnodes[i]->msk = 0;

    bHasMask=true;
    return true;
}

// identical in FPProc and HPProc
double femm::PostProcessor::ElmArea(femmsolver::CElement *elm)
{
    int n[3];
    for(int j=0; j<3; j++) n[j]=elm->p[j];

    double b0=meshnodes[n[1]]->y - meshnodes[n[2]]->y;
    double b1=meshnodes[n[2]]->y - meshnodes[n[0]]->y;
    double c0=meshnodes[n[2]]->x - meshnodes[n[1]]->x;
    double c1=meshnodes[n[0]]->x - meshnodes[n[2]]->x;
    return (b0*c1-b1*c0)/2.;
}


// identical in FPProc and HPProc
CComplex femm::PostProcessor::HenrotteVector(int k) const
{
    int n[3];
    double b[3],c[3];

    for(int i=0; i<3; i++)
    {
        n[i] = meshelems[k]->p[i];
    }

    b[0]=meshnodes[n[1]]->y - meshnodes[n[2]]->y;
    b[1]=meshnodes[n[2]]->y - meshnodes[n[0]]->y;
    b[2]=meshnodes[n[0]]->y - meshnodes[n[1]]->y;
    c[0]=meshnodes[n[2]]->x - meshnodes[n[1]]->x;
    c[1]=meshnodes[n[0]]->x - meshnodes[n[2]]->x;
    c[2]=meshnodes[n[1]]->x - meshnodes[n[0]]->x;

    double da = (b[0] * c[1] - b[1] * c[0]);

    CComplex v = 0;
    for(int i=0; i<3; i++)
    {
        v -= meshnodes[n[i]]->msk * (b[i] + I * c[i]) / (da * LengthConv[problem->LengthUnits]);  // grad
    }

    return v;
}


// identical in FPProc and HPProc
void femm::PostProcessor::bendContour(double angle, double anglestep)
{
    if (angle==0) return;
    if (anglestep==0) anglestep=1;

    int k,n;
    double d,tta,dtta,R;
    CComplex c,a0,a1;

    // check to see if there are at least enough
    // points to have made one line;
    k = contour.size()-1;
    if (k<1) return;

    // restrict the angle of the contour to 180 degrees;
    if ((angle<-180.) || (angle>180.))
    {
        return;
    }
    n = (int) ceil(fabs(angle/anglestep));
    tta = angle*PI/180.;
    dtta = tta/((double) n);

    // pop last point off of the contour;
    a1 = contour[k];
    contour.erase(contour.begin()+k);
    a0 = contour[k-1];

    // compute location of arc center;
    // and radius of the circle that the
    // arc lives on.
    d = abs(a1-a0);
    R = d / ( 2. * sin(fabs(tta/2.)) );

    if(tta>0)
    {
        c = a0 + (R/d) * (a1-a0) * exp(I*(PI-tta)/2.);
    }
    else
    {
        c = a0 + (R/d) * (a1-a0) * exp(-I*(PI+tta)/2.);
    }

    // add the points on the contour
    for(k=1; k<=n; k++)
    {
        contour.push_back( c + (a0 - c) * exp(k * I * dtta) );
    }
}

void PostProcessor::clearContour()
{
    contour.clear();
}

void PostProcessor::clearSelection()
{
    problem->unselectAll();
    bHasMask = false;
}

void PostProcessor::clearBlockSelection()
{
    bHasMask = false;
    for(auto &block: problem->labellist)
    {
        block->IsSelected = false;
    }
}

void PostProcessor::toggleSelectionForGroup(int group)
{
    for (auto &block : problem->labellist)
    {
        if (group==0 || block->InGroup == group)
        {
            block->ToggleSelect();
        }
    }
    bHasMask = false;
}

const std::vector<CComplex> &PostProcessor::getContour() const
{
    return contour;
}

const std::vector<std::unique_ptr<femmsolver::CElement> > &PostProcessor::getMeshElements() const
{
    return meshelems;
}

// identical in epproc and hpproc
double PostProcessor::AECF(const femmsolver::CElement *elem) const
{
    // Computes the permeability correction factor for axisymmetric
    // external regions.  This is sort of a kludge, but it's the best
    // way I could fit it in.  The structure of the code wasn't really
    // designed to have a permeability that varies with position in a
    // continuous way.

    if (problem->problemType == PLANAR) return 1.; // no correction for planar problems
    if (!problem->labellist[elem->lbl]->IsExternal) return 1; // only need to correct for external regions

    double r=abs(elem->ctr-I*problem->extZo);
    return (r*r)/(problem->extRo*problem->extRi); // permeability gets divided by this factor;
}

// identical in epproc and hpproc
double PostProcessor::AECF(const femmsolver::CElement *elem, CComplex p) const
{
    // Correction factor for a point within the element, rather than
    // for the center of the element.
    if (problem->problemType == PLANAR) return 1.; // no correction for planar problems
    if (!problem->labellist[elem->lbl]->IsExternal) return 1; // only need to correct for external regions
    double r=abs(p-I*problem->extZo);
    if (r==0)
        return AECF(elem);
    return (r*r)/(problem->extRo*problem->extRi); // permeability gets divided by this factor;
}

// almost the same in epproc and hpproc; differences noted by comments
void PostProcessor::getNodalD(CComplex *d, int N) const
{
    // this method is only valid for heatflow and electrostatics problems; otherwise punt
    if (problem->filetype != FileType::HeatFlowFile && problem->filetype != FileType::ElectrostaticsFile )
        return;

    int i,j,k,n,m,p,eos,nos,qn;
    int lf,rt;
    double xi,yi,ii,xx,xy,yy,iv,xv,yv,dx,dy,dv,Ex,Ey,det;
    static int q[21];
    bool flag;

    const auto *elem = reinterpret_cast<const femmsolver::CHSElement*>(meshelems[N].get());
    for(i=0;i<3;i++)
    {
        j=elem->p[i];
        lf=rt=-1;
        flag=false;
        for(eos=0;eos<NumList[j];eos++) if(ConList[j][eos]==N) break;

        // scan ccw
        for(k=0,m=eos,qn=0;k<NumList[j];k++)
        {
            n=ConList[j][m];
            const auto &conElem = getMeshElement(n);
            if(!isSameMaterial(*elem,*conElem)) break;

            // figure out which node is the next one in the ccw direction
            // Note that the ConList has been sorted in ccw order,
            // and the nodes in each element are sorted in ccw order,
            // making this task a bit easier.  The next node
            // ends up in the variable p
            for(nos=0;nos<3;nos++) if(conElem->p[nos]==j) break;
            if(nos==3) break;
            nos--;
            if(nos<0) nos=2;
            p=conElem->p[nos];

            // add this node to the list.  We can have a max of 20 nodes,
            // which should never actually occur (usually about 6 to 8)
            if (qn<20) q[qn++]=p;

            // if this is a fixed boundary, get out of the loop;
            if ((meshnodes[j]->Q!=-2) && (meshnodes[p]->Q!=-2)){
                rt=p;
                break;
            }

            m++; if(m==NumList[j]) m=0;
        }

        // scan cw
        for(k=0,m=eos;k<NumList[j];k++)
        {
            n=ConList[j][m];
            const auto &conElem = getMeshElement(n);
            if(!isSameMaterial(*elem,*conElem)) break;

            // figure out which node is the next one in the cw direction
            // The next node ends up in the variable p
            for(nos=0;nos<3;nos++) if(conElem->p[nos]==j) break;
            if(nos==3) break;
            nos++;
            if(nos>2) nos=0;
            p=conElem->p[nos];

            // add this node to the list.  We can have a max of 20 nodes,
            // which should never actually occur (usually about 6 to 8)
            if (qn<20) q[qn++]=p;

            // if this node has a fixed definition, get out of the loop;
            if((meshnodes[j]->Q!=-2) &&(meshnodes[p]->Q!=-2)){
                lf=p;
                break;
            }

            m--; if(m<0) m=NumList[j]-1;
        }

        // catch some annoying special cases;
        if ((lf==rt) && (rt!=-1) && (meshnodes[j]->Q!=-2))
        {
            // The node of interest is at the end of a conductor; not much to
            // do but punt;
            d[i]=elem->D;
            flag=true;
        }
        else if ((rt!=-1) && (meshnodes[j]->Q!=-2) && (lf==-1))
        {
            // Another instance of a node at the
            // end of a conductor; punt!
            d[i]=elem->D;
            flag=true;
        }
        else if ((lf!=-1) && (meshnodes[j]->Q!=-2) && (rt==-1))
        {
            // Another instance of a node at the
            // end of a conductor; punt!
            d[i]=elem->D;
            flag=true;
        }
        else if((lf==-1) && (rt==-1) && (meshnodes[j]->Q!=-2))
        {
            // The node of interest is an isolated charge. Again, not much to
            // do but punt;
            d[i]=elem->D;
            flag=true;
        }
        else if((lf!=-1) && (rt!=-1) && (meshnodes[j]->Q!=-2))
        {

            // The node of interest is on some boundary where the charge is fixed.
            // if the angle is shallow enough, we can just do the regular thing;
            // Otherwise, we punt.
            CComplex x,y;
            x=meshnodes[lf]->CC()-meshnodes[j]->CC(); x/=abs(x);
            y=meshnodes[j]->CC()-meshnodes[rt]->CC(); y/=abs(y);
            if(std::abs(arg(x/y))>10.0001*PI/180.)
            {
                // if the angle is greater than 10 degrees, punt;
                d[i]=elem->D;
                flag=true;
            }
        }

        if(flag==false)
        {
            // The nominal case.
            // Fit a plane through the nodes in the list to solve for E.
            // Then, multiply by permittivity to get D.
            xi=yi=ii=xx=xy=yy=iv=xv=yv=0;

            q[qn++]=j;

            for(k=0;k<qn;k++)
            {
                dx=meshnodes[q[k]]->x-meshnodes[j]->x;
                dy=meshnodes[q[k]]->y-meshnodes[j]->y;
                switch (problem->filetype) {
                case FileType::ElectrostaticsFile:
                {
                    const auto nodej = reinterpret_cast<femmsolver::CSMeshNode*>(meshnodes[j].get());
                    const auto nodek = reinterpret_cast<femmsolver::CSMeshNode*>(meshnodes[q[k]].get());
                    dv=nodej->V - nodek->V;
                }
                    break;
                case FileType::HeatFlowFile:
                {
                    const auto nodej = reinterpret_cast<femmsolver::CHMeshNode*>(meshnodes[j].get());
                    const auto nodek = reinterpret_cast<femmsolver::CHMeshNode*>(meshnodes[q[k]].get());
                    dv=nodej->T - nodek->T;
                }
                    break;
                default:
                    assert(false && "Impossible file type at this point.");
                    return;
                    break;
                }

                ii+=1.;
                xi+=dx;
                yi+=dy;
                xx+=dx*dx;
                xy+=dx*dy;
                yy+=dy*dy;
                iv+=dv;
                xv+=dx*dv;
                yv+=dy*dv;
            }
            det=(-(ii*xy*xy) + 2*xi*xy*yi - xx*yi*yi - xi*xi*yy + ii*xx*yy)*LengthConv[problem->LengthUnits];

            if (det==0) d[i]=elem->D;
            else{
                Ex=(iv*xy*yi - xv*yi*yi - ii*xy*yv + xi*yi*yv - iv*xi*yy + ii*xv*yy)/det;
                Ey=(iv*xi*xy - ii*xv*xy + xi*xv*yi - iv*xx*yi - xi*xi*yv + ii*xx*yv)/det;
                switch (problem->filetype) {
                case FileType::ElectrostaticsFile:
                {
                    const auto bprop = reinterpret_cast<CSMaterialProp*>(problem->blockproplist[elem->blk].get());
                    d[i] = bprop->ex * Ex * eo + I * bprop->ey * Ey * eo;
                    d[i]/=AECF(elem,meshnodes[j]->CC());
                }
                    break;
                case FileType::HeatFlowFile:
                {
                    const auto nodej = reinterpret_cast<femmsolver::CHMeshNode*>(meshnodes[j].get());
                    const auto bprop = reinterpret_cast<CHMaterialProp*>(problem->blockproplist[elem->blk].get());
                    CComplex kn=bprop->GetK(nodej->T);
                    d[i]= Re(kn)*Ex + I*Im(kn)*Ey;
                }
                    break;
                default:
                    assert(false && "Impossible file type at this point.");
                    break;
                }
            }
        }
    }
}



// identical in FPProc and HPProc
void femm::PostProcessor::FindBoundaryEdges()
{
    int i, j;
    const static int plus1mod3[3] = {1, 2, 0};
    const static int minus1mod3[3] = {2, 0, 1};

    // Init all elements' neigh to be unfinished.
    for(i = 0; i < (int)meshelems.size(); i ++)
    {
        for(j = 0; j < 3; j ++)
            meshelems[i]->n[j] = 0;
    }

    int orgi, desti;
    int ei, ni;
    bool done;

    // Loop all elements, to find and set there neighs.
    for(i = 0; i < (int)meshelems.size(); i ++)
    {
        for(j = 0; j < 3; j ++)
        {
            if(meshelems[i]->n[j] == 0)
            {
                // Get this edge's org and dest node index,
                orgi = meshelems[i]->p[plus1mod3[j]];
                desti = meshelems[i]->p[minus1mod3[j]];
                done = false;
                // Find this edge's neigh from the org node's list
                for(ni = 0; ni < NumList[orgi]; ni ++)
                {
                    // Find a Element around org node contained dest node of this edge.
                    ei = ConList[orgi][ni];
                    if (ei == i) continue; // Skip myself.
                    // Check this Element's 3 vert to see if there exist dest node.
                    if(meshelems[ei]->p[0] == desti) {
                        done = true;
                        break;
                    } else if(meshelems[ei]->p[1] == desti) {
                        done = true;
                        break;
                    } else if(meshelems[ei]->p[2] == desti) {
                        done = true;
                        break;
                    }
                }
                if (!done) {
                    // This edge must be a Boundary Edge.
                    meshelems[i]->n[j] = 1;
                }
            } // Finish One Edge
        } // End of One Element Loop
    } // End of Main Loop
}

// identical in hpproc and epproc
void PostProcessor::getPointD(double x, double y, CComplex &D, const femmsolver::CElement &element) const
{
    // this method is only valid for heatflow and electrostatics problems; otherwise punt
    if (problem->filetype != FileType::HeatFlowFile && problem->filetype != FileType::ElectrostaticsFile )
        return;

    // now we know that we can safely use a reinterpret_cast
    const auto elm = reinterpret_cast<const femmsolver::CHSElement&>(element);

    // elm is a reference to the element that contains the point of interest.
    if(!Smooth){
        D=elm.D;
        return;
    }

    const auto &n0 = meshnodes[elm.p[0]];
    const auto &n1 = meshnodes[elm.p[1]];
    const auto &n2 = meshnodes[elm.p[2]];
    double a[3],b[3],c[3];
    a[0]=n1->x * n2->y - n2->x * n1->y;
    a[1]=n2->x * n0->y - n0->x * n2->y;
    a[2]=n0->x * n1->y - n1->x * n0->y;
    b[0]=n1->y - n2->y;
    b[1]=n2->y - n0->y;
    b[2]=n0->y - n1->y;
    c[0]=n2->x - n1->x;
    c[1]=n0->x - n2->x;
    c[2]=n1->x - n0->x;
    double da=(b[0]*c[1]-b[1]*c[0]);

    D=0;
    for(int i=0;i<3;i++)
    {
        D+=(elm.d[i]*(a[i]+b[i]*x+c[i]*y)/da);
    }
}
