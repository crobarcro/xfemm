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
        sztibi82@gmail.com
        richard.crozier@yahoo.co.uk
        johannes.zarl-zierl@jku.at

   Contributions by Johannes Zarl-Zierl were funded by
   Linz Center of Mechatronics GmbH (LCM)
*/

// hsolver.cpp : implementation of the HSolver class
//

#include "femmcomplex.h"
#include "femmconstants.h"
#include "linsolve/backend_factory.h"
#include "fparse.h"
#include "hsolver.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <malloc.h>

// template instantiation:
#include "../libfemm/feasolver.cpp"
#include "../libfemm/cuthill.cpp"
template class FEASolver<
        femm::CHPointProp
        , femm::CHBoundaryProp
        , femm::CHMaterialProp
        , femm::CHConductor
        , femm::CHBlockLabel
        , femmsolver::CElement
        >;

#ifndef _WIN32
#define _strnicmp strncasecmp
#endif

namespace {
//conversion to internal working units of m
constexpr double units[]={0.0254,0.001,0.01,1,2.54e-5,1.e-6};
constexpr double sqr(double x)
{
    return x*x;
}
} // anon namespace

using namespace std;
using namespace femm;


/////////////////////////////////////////////////////////////////////////////
// HSolver construction/destruction

HSolver::HSolver()
    : meshnode(nullptr)
    , Tprev(nullptr)
{

    // initialise the warning message box function pointer to
    // point to the PrintWarningMsg function
    WarnMessage = &PrintWarningMsg;

    bMultiplyDefinedLabels = false;
}

HSolver::~HSolver()
{
    CleanUp();
}

void HSolver::CleanUp()
{
    FEASolver_type::CleanUp();
    if (meshnode)
    {
        delete[] meshnode;
        meshnode = nullptr;
    }
    if (Tprev)
    {
        delete[] Tprev;
        Tprev = nullptr;
    }
}

/////////////////////////////////////////////////////////////////////////////
// HSolver commands

void HSolver::MsgBox(const char* message)
{
    printf("%s\n", message);
}

bool HSolver::LoadProblemFile ()
{
    std::string fehFile = PathName+".feh";

    bool ret = FEASolver_type::LoadProblemFile(fehFile);
    return ret;
}

int HSolver::LoadPrev()
{
    if (previousSolutionFile.empty())
    {
        if (dT!=0)
        {
            dT = 0;
            WarnMessage("Warning: no previous solution set even though dT is non-zero. Resetting dT...\n");
        }
        return true;
    }

	FILE *fp;
    double x,y;
    int k;
	char s[1024],q[256];

    if ((fp=fopen(previousSolutionFile.c_str(),"rt"))==NULL)
	{
		return BADELEMENTFILE;
	}

	// parse the file
	k=0;
	while (fgets(s,1024,fp)!=NULL)
	{
		sscanf(s,"%s",q);
		if( _strnicmp(q,"[solution]",11)==0){
			k=1;
			break;
		}
	}

	// case where the solution is never found.
	if (k==0)
	{
		fclose(fp);
		return BADELEMENTFILE;
	}

	// read in the solution
	fgets(s,1024,fp);
	sscanf(s,"%i",&k);
	if(k!=NumNodes)
	{
		fclose(fp);
		return BADELEMENTFILE;
	}

    Tprev=new double[NumNodes];

	for(k=0;k<NumNodes;k++)
	{
		fgets(s,1024,fp);
		sscanf(s,"%lf	%lf	%lf",&x,&y,&Tprev[k]);
	}

	return 0;
}

LoadMeshErr HSolver::LoadMesh(bool deleteFiles)
{
	int i,j,k,q,n0,n1,n;
	char infile[256];
	FILE *fp;
	char s[1024];
    double c[]={0.0254,0.001,0.01,1,2.54e-5,1.e-6};


	//read meshnodes;
	sprintf(infile,"%s.node",PathName.c_str());
	if((fp=fopen(infile,"rt"))==NULL){
		return BADELEMENTFILE;
	}
	fgets(s,1024,fp);
	sscanf(s,"%i",&k);
	NumNodes = k;

    meshnode = new CNode[k];
    CNode node;
	for(i = 0; i < k; i++)
	{
		fscanf(fp,"%i",&j);
		fscanf(fp,"%lf",&node.x);
		fscanf(fp,"%lf",&node.y);
		fscanf(fp,"%i",&n);

		if (n > 1)
		{
			// strip off point BC number;
			j = n & 0xffff;
			j = j - 2;
			if (j<0)
			{
			    j = -1;
			}

			// strip off Conductor number
			n = (n - (n & 0xffff))/0x10000 - 1;
		}
		else
		{
			j = -1;
			n = -1;
		}
        node.BoundaryMarker = j;
		node.InConductor = n;

		// convert all lengths to internal working units of mm
		node.x *= c[LengthUnits];
		node.y *= c[LengthUnits];

		meshnode[i] = node;
	}
	fclose(fp);

	//read in periodic boundary conditions;
	sprintf(infile,"%s.pbc",PathName.c_str());
	if((fp=fopen(infile,"rt"))==NULL){
		return BADPBCFILE;
	}
	fgets(s,1024,fp);
	sscanf(s,"%i",&k);
	NumPBCs=k;

    if (k!=0) pbclist.reserve(k);
	CCommonPoint pbc;
	for(i=0;i<k;i++){
		fscanf(fp,"%i",&j);
		fscanf(fp,"%i",&pbc.x);
		fscanf(fp,"%i",&pbc.y);
		fscanf(fp,"%i",&pbc.t);
        pbclist.push_back(pbc);
	}
	fclose(fp);

	// read in elements;
	sprintf(infile,"%s.ele",PathName.c_str());
	if((fp=fopen(infile,"rt"))==NULL){
		return BADELEMENTFILE;
	}
	fgets(s,1024,fp);
	sscanf(s,"%i",&k); NumEls=k;

    meshele.reserve(k);
    femmsolver::CElement elm;

	int defaultLabel;
	for(i=0,defaultLabel=-1;i<NumBlockLabels;i++)
		if (labellist[i].IsDefault) defaultLabel=i;

	for(i=0;i<k;i++){
		fscanf(fp,"%i",&j);
		fscanf(fp,"%i",&elm.p[0]);
		fscanf(fp,"%i",&elm.p[1]);
		fscanf(fp,"%i",&elm.p[2]);
		fscanf(fp,"%i",&elm.lbl);
		elm.lbl--;
		if(elm.lbl<0) elm.lbl=defaultLabel;
		if(elm.lbl<0){
		    string msg = "Material properties have not been defined for\n";
            msg += "all regions. Press the \"Run Mesh Generator\"\n";
            msg += "button to highlight the problem regions.";
            WarnMessage(msg.c_str());

			fclose(fp);
            if (deleteFiles)
            {
                sprintf(infile,"%s.ele",PathName.c_str());
                remove(infile);
                sprintf(infile,"%s.node",PathName.c_str());
                remove(infile);
                sprintf(infile,"%s.pbc",PathName.c_str());
                remove(infile);
                sprintf(infile,"%s.poly",PathName.c_str());
                remove(infile);
                sprintf(infile,"%s.edge",PathName.c_str());
                remove(infile);
            }
            return MISSINGMATPROPS;
		}
		// look up block type out of the list of block labels
		elm.blk=labellist[elm.lbl].BlockType;

        meshele.push_back(elm);
	}
	fclose(fp);

	// initialize edge bc's and element permeabilities;
	for(i=0;i<NumEls;i++)
		for(j=0;j<3;j++)
			meshele[i].e[j] = -1;

	// read in edges to which boundary conditions are applied;

		// first, do a little bookkeeping so that element
		// associated with a give edge can be identified fast
		int *nmbr;
		int **mbr;

		nmbr=(int *)calloc(NumNodes,sizeof(int));

		// Make a list of how many elements that tells how
		// many elements to which each node belongs.
		for(i=0;i<NumEls;i++)
			for(j=0;j<3;j++)
				nmbr[meshele[i].p[j]]++;

		// mete out some memory to build a list of the
		// connectivity...
		mbr=(int **)calloc(NumNodes,sizeof(int *));
		for(i=0;i<NumNodes;i++){
			mbr[i]=(int *)calloc(nmbr[i],sizeof(int));
			nmbr[i]=0;
		}

		// fill up the connectivity information;
		for(i=0;i<NumEls;i++)
			for(j=0;j<3;j++)
			{
				k=meshele[i].p[j];
				mbr[k][nmbr[k]]=i;
				nmbr[k]++;
			}

	sprintf(infile,"%s.edge",PathName.c_str());
	if((fp=fopen(infile,"rt"))==NULL)
	{
		return BADEDGEFILE;
	}
	fscanf(fp,"%i",&k);	// read in number of lines

	fscanf(fp,"%i",&j);	// read in boundarymarker flag;
	for(i=0;i<k;i++)
	{
		fscanf(fp,"%i",&j);
		fscanf(fp,"%i",&n0);
		fscanf(fp,"%i",&n1);
		fscanf(fp,"%i",&n);

		// BC number;
		if (n<0)
		{
			n=(-n);
			j = (n & 0xffff) - 2;
			if (j<0) j = -1;

			// Conductor number;
			n= (n - (n & 0xffff))/0x10000 - 1;
			if (n>=0)
			{
				meshnode[n0].InConductor=n;
				meshnode[n1].InConductor=n;
			}
		}
		else j=-1;

		if (j>=0)
		{
			// search through elements to find one containing the line;
			// set corresponding edge equal to the bc number
            for(q=0,n=false;q<nmbr[n0];q++)
			{
				elm=meshele[mbr[n0][q]];

                if ((elm.p[0] == n0) && (elm.p[1] == n1)) {elm.e[0]=j; n=true;}
                if ((elm.p[0] == n1) && (elm.p[1] == n0)) {elm.e[0]=j; n=true;}

                if ((elm.p[1] == n0) && (elm.p[2] == n1)) {elm.e[1]=j; n=true;}
                if ((elm.p[1] == n1) && (elm.p[2] == n0)) {elm.e[1]=j; n=true;}

                if ((elm.p[2] == n0) && (elm.p[0] == n1)) {elm.e[2]=j; n=true;}
                if ((elm.p[2] == n1) && (elm.p[0] == n0)) {elm.e[2]=j; n=true;}

				meshele[mbr[n0][q]]=elm;

				//this is a little hack: line charge distributions should be applied
				//to at most one element;
				if((lineproplist[j].BdryFormat==2) && (n)) q=nmbr[n0];
			}
		}

	}
	fclose(fp);

	// free up the connectivity information
	free(nmbr);
	for(i=0;i<NumNodes;i++) free(mbr[i]);
	free(mbr);

    if (deleteFiles)
    {
        // clear out temporary files
        sprintf(infile,"%s.ele",PathName.c_str());
        remove(infile);
        sprintf(infile,"%s.node",PathName.c_str());
        remove(infile);
        sprintf(infile,"%s.pbc",PathName.c_str());
        remove(infile);
        sprintf(infile,"%s.poly",PathName.c_str());
        remove(infile);
    }

    return NOERROR;
}

//CComplex CMaterialProp::GetK(double t)
//{
//	int i,j;
//
//	// Kx returned as real part;
//	// Ky returned as imag part
//
//	if (npts==0) return (Kx+I*Ky);
//	if (npts==1) return (Im(Kn[0])*(1+I));
//	if (t<=Re(Kn[0])) return (Im(Kn[0])*(1+I));
//	if (t>=Re(Kn[npts-1])) return (Im(Kn[npts-1])*(1+I));
//
//	for(i=0,j=1;j<npts;i++,j++)
//	{
//		if((t>=Re(Kn[i])) && (t<=Re(Kn[j])))
//		{
//			return (1+I)*(Im(Kn[i])+Im(Kn[j]-Kn[i])*Re(t-Kn[i])/Re(Kn[j]-Kn[i]));
//		}
//	}
//
//	return (Kx+I*Ky);
//}

//////////////////////////////////////////////////////////////////////////////////////////////

int HSolver::AnalyzeProblem(femm::LinearSystemBackend<double> &L)
{
	int i,j,k,bf,pctr=0;
	double Me[3][3],be[3];		// element matrices;
	double l[3],p[3],q[3];		// element shape parameters;
	int n[3],ne[3];				// numbers of nodes for a particular element;
	double a,K,r,z,kludge;
	double bta,Tinf,Tlast,*Vo;
    int IsNonlinear=false;
    femmsolver::CElement *El;
	CComplex kn;
	int iter=0;

	Depth*=units[LengthUnits];
	extRo*=units[LengthUnits];
	extRi*=units[LengthUnits];
	extZo*=units[LengthUnits];
	kludge=1;

	//TheView->SetDlgItemText(IDC_FRAME1,"Matrix Construction");

	Vo=(double *) calloc(NumNodes,sizeof(double));

	// scan through the problem to see if there are any elements
	// with a nonlinear conductivity
	for(i=0;i<NumNodes;i++)
	{
		if (blockproplist[meshele[i].blk].npts>0){
            IsNonlinear=true;
            break;
		}
	}

	do{
		// copy old solution
		for(i=0;i<NumNodes;i++) Vo[i]=L.solution()[i];
		L.wipe();

		// do some book-keeping related to fixed boundary conditions;
		// The P vector denotes which nodes have an assigned value
		// The V vector denotes the assigned value
		for(i=0;i<NumNodes;i++)
		{
			L.node_flag()[i] = -2;
            if(meshnode[i].BoundaryMarker >= 0)
                if(nodeproplist[meshnode[i].BoundaryMarker].qp == 0)
				{
                    L.solution()[i] = nodeproplist[meshnode[i].BoundaryMarker].V;
					L.node_flag()[i] = -1;
				}

			if(meshnode[i].InConductor >= 0)
				if(circproplist[meshnode[i].InConductor].CircType == 1)
				{
					L.solution()[i] = circproplist[meshnode[i].InConductor].V;
					L.node_flag()[i] = meshnode[i].InConductor;
				}
		}

		// account for fixed boundary conditions along segments;
		for(i=0;i<NumEls;i++)
		{
			for(j=0;j<3;j++){
				k=j+1; if(k==3) k=0;
				if(meshele[i].e[j]>=0)
				{
					if(lineproplist[ meshele[i].e[j] ].BdryFormat==0)
					{
						L.solution()[meshele[i].p[j]]=lineproplist[meshele[i].e[j]].Tset;
						L.solution()[meshele[i].p[k]]=lineproplist[meshele[i].e[j]].Tset;
						L.node_flag()[meshele[i].p[j]]=-1;
						L.node_flag()[meshele[i].p[k]]=-1;
					}
				}
			}
		}


		// build element matrices using the matrices derived in Allaire's book.
		for(i=0;i<NumEls;i++)
		{
			// update progress bar
			j=5*((i*20)/NumEls);
			if(j!=pctr){ pctr=j;
               // TheView->m_prog1.SetPos(pctr);
            }

			// zero out Me, be;
			for(j=0;j<3;j++){
				for(k=0;k<3;k++) Me[j][k]=0;
				be[j]=0;
			}

			// Determine shape parameters.
			// l's are element side lengths;
			// p's corresponds to the `b' parameter in Allaire
			// q's corresponds to the `c' parameter in Allaire
			El=&meshele[i];

			for(k=0;k<3;k++) n[k]=El->p[k];
			p[0]=meshnode[n[1]].y - meshnode[n[2]].y;
			p[1]=meshnode[n[2]].y - meshnode[n[0]].y;
			p[2]=meshnode[n[0]].y - meshnode[n[1]].y;
			q[0]=meshnode[n[2]].x - meshnode[n[1]].x;
			q[1]=meshnode[n[0]].x - meshnode[n[2]].x;
			q[2]=meshnode[n[1]].x - meshnode[n[0]].x;
			for(j=0,k=1;j<3;k++,j++){
				if (k==3) k=0;
				l[j]=sqrt( sqr(meshnode[n[k]].x-meshnode[n[j]].x) +
						   sqr(meshnode[n[k]].y-meshnode[n[j]].y) );
			}
			a=(p[0]*q[1]-p[1]*q[0])/2.;
			r=(meshnode[n[0]].x+meshnode[n[1]].x+meshnode[n[2]].x)/3.;

			// get the thermal conductivites to use for this element;
			kn = (blockproplist[El->blk].GetK(Vo[n[0]]) +
				  blockproplist[El->blk].GetK(Vo[n[1]]) +
				  blockproplist[El->blk].GetK(Vo[n[2]]))/3.;

			if (ProblemType==AXISYMMETRIC){
				Depth=2.*PI*r;

				// "Warp" the permeability of this element is part of
				// the conformally mapped external region
				if(labellist[meshele[i].lbl].IsExternal)
				{
					z=(meshnode[n[0]].y+meshnode[n[1]].y+meshnode[n[2]].y)/3. - extZo;
					kludge=(r*r+z*z)/(extRi*extRo);
				}
				else kludge=1;
			}


			// x-contribution;
			K = -Depth*Re(kn)/(4.*a)/kludge;
			for(j=0;j<3;j++)
				for(k=j;k<3;k++)
				{
					Me[j][k] += K*p[j]*p[k];
					if (j!=k) Me[k][j]+=K*p[j]*p[k];
				}

			// y-contribution;
			K = -Depth*Im(kn)/(4.*a)/kludge;
			for(j=0;j<3;j++)
				for(k=j;k<3;k++)
				{
					Me[j][k] +=K*q[j]*q[k];
					if (j!=k) Me[k][j]+=K*q[j]*q[k];
				}

			// contribution to Me and be from time-transient term
/*			if (dT!=0)
			{
				K = -Depth*blockproplist[El->blk].Kt*a/(12.*dT);

				Me[0][0]+=2.*K;
				Me[1][1]+=2.*K;
				Me[2][2]+=2.*K;
				Me[0][1]+=K; Me[1][0]+=K;
				Me[0][2]+=K; Me[2][0]+=K;
				Me[1][2]+=K; Me[2][1]+=K;

				be[0]+=K*(2.*Tprev[n[0]] +    Tprev[n[1]] +    Tprev[n[2]]);
				be[1]+=K*(   Tprev[n[0]] + 2.*Tprev[n[1]] +    Tprev[n[2]]);
				be[2]+=K*(   Tprev[n[0]] +    Tprev[n[1]] + 2.*Tprev[n[2]]);
			} */

			if (dT!=0)
			{
				K = -Depth*blockproplist[El->blk].Kt*a/(3.*dT);

				Me[0][0]+=K;
				Me[1][1]+=K;
				Me[2][2]+=K;

				be[0]+=K*Tprev[n[0]];
				be[1]+=K*Tprev[n[1]];
				be[2]+=K*Tprev[n[2]];
			}

			// contribution to be[] from volume charge density
			for(j = 0;j<3;j++){
				K = -Depth*(blockproplist[El->blk].qv)*a/3.;
				be[j]+=K;
			}


			for(j=0;j<3;j++)
			{
				if (El->e[j] >= 0)
				{
					k=j+1; if(k==3) k=0;

					if (ProblemType==AXISYMMETRIC)
						Depth=PI*(meshnode[n[j]].x + meshnode[n[k]].x);

					// contributions to Me, be from derivative boundary conditions;
					// !!! need to put in contribution here for radiation....
					bf=lineproplist[El->e[j]].BdryFormat;
					if ((bf==1) || (bf==2) || (bf==3))
					{
						double c0,c1;

						switch(bf)
						{
							case 1:
								c1=lineproplist[El->e[j]].qs;
								c0=0;
								break;
							case 2:
								c0=lineproplist[El->e[j]].h;
								c1=-c0*lineproplist[El->e[j]].Tinf;
								break;
							case 3:
                                IsNonlinear=true;
								bta =lineproplist[El->e[j]].beta;
								Tinf=lineproplist[El->e[j]].Tinf;
								Tlast=(Vo[n[j]]+Vo[n[k]])/2.;

								c0 = 4.*bta*Ksb*pow(Tlast,3.);
								c1 = -(bta*Ksb*(pow(Tinf,4.) + 3.*pow(Tlast,4.)));

								break;
                            default:
                                assert(false); // can't happen according to if statement above
								break;
						}

						if (ProblemType==AXISYMMETRIC)
						{
							K =-2.*PI*c0*l[j]/6.;
							Me[j][j]+=K*2. *(3.*meshnode[n[j]].x + meshnode[n[k]].x)/4.;
							Me[k][k]+=K*2. *(meshnode[n[j]].x + 3.*meshnode[n[k]].x)/4.;
							Me[j][k]+=K    *(meshnode[n[j]].x + meshnode[n[k]].x)/2.;
							Me[k][j]+=K    *(meshnode[n[j]].x + meshnode[n[k]].x)/2.;

							K = 2.*PI*c1*l[j]/2.;
							be[j]+=K*(2.*meshnode[n[j]].x + meshnode[n[k]].x)/3.;
							be[k]+=K*(meshnode[n[j]].x + 2.*meshnode[n[k]].x)/3.;
						}
						else
						{
							K =-Depth*c0*l[j]/6.;
							Me[j][j]+=K*2.;
							Me[k][k]+=K*2.;
							Me[j][k]+=K;
							Me[k][j]+=K;

							K = Depth*c1*l[j]/2.;
							be[j]+=K;
							be[k]+=K;
						}
					}
				/*
					// contribution to be[] from surface heating
					if (lineproplist[El->e[j]].BdryFormat==2)
					{
						K =-Depth*lineproplist[El->e[j]].qs*l[j]/2.;
						be[j]+=K;
						be[k]+=K;
					}
				*/
				}
			}

			// process any prescribed nodal values;
			for(j=0;j<3;j++)
			{
				if(L.node_flag()[n[j]]!=-2)
				{
					for(k=0;k<3;k++)
					{
						if(j!=k){
							be[k]-=Me[k][j]*L.solution()[n[j]];
							Me[k][j]=0;
							Me[j][k]=0;
						}
					}
					be[j]=L.solution()[n[j]]*Me[j][j];
				}
			}

			// combine block matrices into global matrices;
			for (j=0;j<3;j++)
			{
				ne[j]=n[j];
				if(meshnode[n[j]].InConductor>=0)
					if(circproplist[meshnode[n[j]].InConductor].CircType==0)
						ne[j]=meshnode[n[j]].InConductor+NumNodes;
			}
			for (j=0;j<3;j++){
				for (k=j;k<3;k++)
                L.put(L.get(ne[j],ne[k])-Me[j][k],ne[j],ne[k]);
				L.rhs()[ne[j]]-=be[j];

				if(ne[j]!=n[j])
				{
					L.put(L.get(n[j],n[j])-Me[j][j],n[j],n[j]);
					L.put(L.get(n[j],ne[j])+Me[j][j],n[j],ne[j]);
				}
			}

		} // end of loop that builds element matrices

		// add in contribution from point charge density;
		for(i=0;i<NumNodes;i++)
		{
            if((meshnode[i].BoundaryMarker>=0) && (L.node_flag()[i]==-2))
			{
				if (ProblemType==AXISYMMETRIC) Depth=2.*PI*meshnode[i].x;
                L.rhs()[i]+=(Depth*nodeproplist[meshnode[i].BoundaryMarker].qp);
				L.node_flag()[i]=-1;
			}

			// some bookkeeping to denote which nodes we can smooth over
			if(meshnode[i].InConductor>=0) L.node_flag()[i]=meshnode[i].InConductor;
		}

		// Apply any periodicity/antiperiodicity boundary conditions that we have
		for(k=0,pctr=0;k<NumPBCs;k++)
		{
			if (pbclist[k].t==0) L.constrain_periodic(pbclist[k].x,pbclist[k].y,false);
			if (pbclist[k].t==1) L.constrain_periodic(pbclist[k].x,pbclist[k].y,true);
		}

		// Finish building the equations that assign conductor voltage;
		for(i=0;i<NumCircProps;i++)
		{
			// put a placeholder on the main diagonal;
			k=NumNodes+i;

			if (circproplist[i].CircType==1)
			{
				K=L.get(0,0);
				L.put(K,k,k);
				L.rhs()[k]=K*circproplist[i].V;
			}

			if(circproplist[i].CircType==0)
			{
				for(j=0,K=0;j<L.dimension();j++) if(j!=k) K+=L.get(k,j);
				if(K!=0){
					L.put(-K,k,k);
					L.rhs()[k]=circproplist[i].q;
				}
				else L.put(L.get(0,0),k,k);


			}
		}

		// solve the problem;
        femm::SolveOptions opts;
        opts.warm_start = (iter > 0);
        if (L.solve(opts).converged==false){
			free(Vo);
            return false;
		}
        iter++;

        if (IsNonlinear == true)
		{
			double e1=0;
			double e2=0;
			int prog;
			char fmsg[256];

			sprintf(fmsg,"Iteration(%i) ",iter);
            printf("%s", fmsg);
			//TheView->SetDlgItemText(IDC_FRAME2,fmsg);

			for(i=0;i<NumNodes;i++){
				e1+=(L.solution()[i]-Vo[i])*(L.solution()[i]-Vo[i]);
				e2+=(Vo[i]*Vo[i]);
			}
			if(e2!=0)
			{
				// test to see if we have converged.
                if(sqrt(e1/e2) < Precision*100.) IsNonlinear=false;
				prog=(int)  (100.*log10(e1/e2)/(log10(Precision)+2.));
				if (prog>100) prog=100;
			//	TheView->m_prog2.SetPos(prog);

			}
		}

    }while(IsNonlinear == true);

	// compute total charge on conductors
	// with a specified voltage
	for(i=0;i<NumCircProps;i++)
		if(circproplist[i].CircType==1)
			circproplist[i].q=ChargeOnConductor(i,L);

	free(Vo);
    return true;
}

bool HSolver::runSolver(bool verbose)
{
    // load mesh
    LoadMeshErr err = LoadMesh();
    if (err != NOERROR)
    {
        WarnMessage("problem loading mesh:\n");
        WarnMessage(getErrorString(err).c_str());
        return false;
    }

    if (!LoadPrev() && verbose)
    {
        PrintMessage("Loading previous solution\n");
    }

    // renumber using Cuthill-McKee
    if (verbose)
        PrintMessage("renumbering nodes\n");
    if (!Cuthill())
    {
        WarnMessage("problem renumbering node points\n");
        return false;
    }

    if (verbose)
    {
        PrintMessage("solving...");
        std::string stats = "Problem Statistics:\n";
        stats += to_string(NumNodes) + " nodes\n";
        stats += to_string(NumEls) + " elements\n";
        PrintMessage(stats.c_str());
        std::cout << "Precision: " << Precision << "\n";
    }

    std::unique_ptr<femm::LinearSystemBackend<double>> L =
        femm::create_backend<double>(femm::default_backend_kind());
    if (!L)
    {
        WarnMessage("couldn't create linear system backend\n");
        return false;
    }
    L->set_precision(Precision);
    if (!L->create(NumNodes+NumCircProps,BandWidth))
    {
        WarnMessage("couldn't allocate enough space for matrices\n");
        return false;
    }

    if (!AnalyzeProblem(*L))
    {
        WarnMessage("Couldn't solve the problem\n");
        return false;
    }

    if (verbose)
        PrintMessage("Problem solved\n");

    if (!WriteResults(*L))
    {
       WarnMessage("couldn't write results to disk\n");
       return 6;
    }
    if (verbose)
        PrintMessage("results written to disk\n");

    return true;
}

//=========================================================================
//=========================================================================

int HSolver::WriteResults(femm::LinearSystemBackend<double> &L)
{
	// write solution to disk;

	char c[1024];
	FILE *fp,*fz;
	int i;
	double cf;
	// first, echo input .feh file to the .anh file;
	sprintf(c,"%s.feh",PathName.c_str());

	fz=fopen(c,"rt");
	if(fz==NULL)
    {
		printf("Couldn't open %s.feh\n", PathName.c_str());
        return false;
	}

    sprintf(c,"%s.anh",PathName.c_str());
    fp=fopen(c,"wt");
	if(fp==NULL)
    {
		printf("Couldn't write to %s.anh",PathName.c_str());
        return false;
	}

	while(fgets(c,1024,fz)!=NULL)
    {
        fputs(c,fp);
    }
	fclose(fz);

	// then print out node, line, and element information
	fprintf(fp,"[Solution]\n");
    // get conversion factor for conversion from internal working units of
    // mm to the specified length units
	cf = units[LengthUnits];
	fprintf(fp,"%i\n",NumNodes);
	for(i=0;i<NumNodes;i++)
    {
		fprintf(fp,"%.17g	%.17g	%.17g	%i\n",meshnode[i].x/cf,meshnode[i].y/cf,L.solution()[i],L.node_flag()[i]);
    }

	fprintf(fp,"%i\n",NumEls);

	for(i=0;i<NumEls;i++)
    {
		fprintf(fp,"%i	%i	%i	%i\n",
			meshele[i].p[0],meshele[i].p[1],meshele[i].p[2],meshele[i].lbl);
    }

	// print out circuit info
	fprintf(fp,"%i\n",NumCircProps);
	for(i=0;i<NumCircProps;i++)
    {
		fprintf(fp,"%.17g	%.17g\n",L.solution()[NumNodes+i],circproplist[i].q);
    }

	fclose(fp);
    return true;
}

//=========================================================================
//=========================================================================


double HSolver::ChargeOnConductor(int u, femm::LinearSystemBackend<double> &L)
{
	int i,k;
	double b[3],c[3];		// element shape parameters;
	int n[3];				// numbers of nodes for a particular element;
	double a,da,Dx,Dy,vx,vy,Z;
	CComplex kn;

	for(i=0;i<NumNodes;i++)
		if(meshnode[i].InConductor==u) L.scratch()[i]=1;
		else L.scratch()[i]=0;

	// build element matrices using the matrices derived in Allaire's book.
	for(i=0,Z=0;i<NumEls;i++)
	{
		for(k=0;k<3;k++) n[k]=meshele[i].p[k];

		if((L.scratch()[n[0]]!=0) || (L.scratch()[n[1]]!=0) || (L.scratch()[n[2]]!=0))
		{
			// Determine shape parameters.
			b[0]=meshnode[n[1]].y - meshnode[n[2]].y;
			b[1]=meshnode[n[2]].y - meshnode[n[0]].y;
			b[2]=meshnode[n[0]].y - meshnode[n[1]].y;
			c[0]=meshnode[n[2]].x - meshnode[n[1]].x;
			c[1]=meshnode[n[0]].x - meshnode[n[2]].x;
			c[2]=meshnode[n[1]].x - meshnode[n[0]].x;
			da=(b[0]*c[1]-b[1]*c[0]);
			a=da/2.;
			if (ProblemType==AXISYMMETRIC)
				a*=(2.*PI*(meshnode[n[0]].x+meshnode[n[1]].x+meshnode[n[2]].x)/3.);
			else a*=Depth;
			// get normal vector and element flux density;
			for(k=0,kn=0,vx=0,vy=0,Dx=0,Dy=0;k<3;k++)
			{
				vx-=(L.scratch()[n[k]]*b[k])/da;
				vy-=(L.scratch()[n[k]]*c[k])/da;
				Dx-=(L.solution()[n[k]]*b[k])/da;
				Dy-=(L.solution()[n[k]]*c[k])/da;
				kn+=blockproplist[meshele[i].blk].GetK(L.solution()[n[k]])/3.;
			}
			Dx*=Re(kn);
			Dy*=Im(kn);

			Z+=a*(Dx*vx+Dy*vy);
		}
	}

	return Z;
}


// SortNodes: sorts mesh nodes based on a new numbering
void HSolver::SortNodes (std::vector<int> newnum)
{
    // sort mesh nodes based on newnum;
    for(int i = 0; i < NumNodes; i++)
    {
        while(newnum[i] != i)
        {
            int j = newnum[i];
            swap(newnum[i],newnum[j]);
            swap(meshnode[i],meshnode[j]);
        }
    }
}

bool HSolver::handleToken(const string &token, istream &input, ostream &err)
{
    if( token == "[dt]" )
    {
        expectChar(input, '=', err);
        parseValue(input, dT, err);
        return true;
    }
    if( token == "[frequency]")
    {
        err << "Warning: [frequency] is not an allowed parameter for heat flow problems!\n";
        // ignore line
        return true;
    }
    return false;
}
