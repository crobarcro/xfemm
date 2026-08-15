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

// fmesher.cpp : implementation of FMesher Class
//

#include "femmconstants.h"
#include "fmesher.h"
#include "fparse.h"
#include "IntPoint.h"
#include "make_unique.h"

#include "triangle_version.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace femm;
using namespace fmesher;

namespace {
// set up some default behaviors
constexpr double DEFAULT_MINANGLE=30.;
}

FMesher::FMesher()
    : problem(nullptr)
    , Verbose(true)
    , WarnMessage(&PrintWarningMsg)
    , TriMessage(nullptr)
{
    // initialize the problem data structures
    // and default behaviour etc.
    Initialize(femm::FileType::Unknown);
}


FMesher::FMesher(std::shared_ptr<FemmProblem> p)
    : problem(p)
    , Verbose(true)
    , WarnMessage(&PrintWarningMsg)
    , TriMessage(nullptr)
{
}

bool FMesher::Initialize(femm::FileType t)
{
    // clear out all current lines, nodes, and block labels
    meshnode.clear ();
    meshline.clear ();
    greymeshline.clear ();
    probdescstrings.clear ();

    problem = std::make_shared<femm::FemmProblem>(t);
    // set problem attributes to generic ones;
    problem->MinAngle = DEFAULT_MINANGLE;

    return true;
}

femm::FileType FMesher::GetFileType (string PathName)
{
    // find the position of the last '.' in the string
    size_t dotpos = PathName.rfind ('.');

    if (dotpos == string::npos)
    {
        // no '.' found
        return femm::FileType::Unknown;
    }

    // compare different file extensions and return the appropriate string
    if ( PathName.compare (dotpos, string::npos, ".fem") == 0 )
    {
        return femm::FileType::MagneticsFile;
    }
    else if ( PathName.compare (dotpos, string::npos, ".feh") == 0 )
    {
        return femm::FileType::HeatFlowFile;
    }
    else if ( PathName.compare(dotpos, string::npos, ".fee") == 0)
    {
        return femm::FileType::ElectrostaticsFile;
    }
    else
    {
        return femm::FileType::Unknown;
    }

}

const char *FMesher::backendName(Backend selected)
{
    switch (selected)
    {
    case Backend::Triangle:
        return "Triangle";
    case Backend::Tangle:
        return "Tangle";
    }
    return "Unknown";
}

void FMesher::addFileStr (char *s)
{
    probdescstrings.push_back (string (s));
}

bool FMesher::SaveFEMFile(string PathName)
{
    FILE *fp;
    unsigned int i, j;
    int k,t;
    string s;

    // check to see if we are ready to write a datafile;
    if ((fp = fopen(PathName.c_str(),"wt"))==NULL)
    {
        WarnMessage("Couldn't write to specified file.\nPerhaps the file is write-protected?");
        return false;
    }

    // echo the start of the input file
    for (i = 0; i < probdescstrings.size (); i++)
    {
        fprintf(fp,"%s", probdescstrings[i].c_str () );
    }

    // write out node list
    fprintf(fp,"[NumPoints] = %i\n", (int) problem->nodelist.size());
    for(i=0; i<problem->nodelist.size(); i++)
    {
        for(j=0,t=0; j<problem->nodeproplist.size(); j++)
            if(problem->nodeproplist[j]->PointName==problem->nodelist[i]->BoundaryMarkerName) t=j+1;
        fprintf(fp,"%.17g\t%.17g\t%i\t%i",problem->nodelist[i]->x,problem->nodelist[i]->y,t,
                problem->nodelist[i]->InGroup);

        if (problem->filetype == femm::FileType::HeatFlowFile
                || problem->filetype == femm::FileType::ElectrostaticsFile )
        {
            // find and write number of conductor property group
            for (j=0,t=0; j<problem->circproplist.size (); j++)
                if (problem->circproplist[j]->CircName==problem->nodelist[i]->InConductorName) t=j+1;

            fprintf(fp,"\t%i",t);
        }

        fprintf(fp,"\n");
    }

    // write out segment list
    fprintf(fp,"[NumSegments] = %i\n", (int) problem->linelist.size());
    for(i=0; i<problem->linelist.size(); i++)
    {
        for(j=0,t=0; j<problem->lineproplist.size(); j++)
            if(problem->lineproplist[j]->BdryName==problem->linelist[i]->BoundaryMarkerName) t=j+1;

        fprintf(fp,"%i\t%i\t",problem->linelist[i]->n0,problem->linelist[i]->n1);

        if(problem->linelist[i]->MaxSideLength<0)
        {
            fprintf(fp,"-1\t");
        }
        else
        {
            fprintf(fp,"%.17g\t",problem->linelist[i]->MaxSideLength);
        }

        fprintf(fp,"%i\t%i\t%i",t,problem->linelist[i]->Hidden,problem->linelist[i]->InGroup);

        if (problem->filetype == femm::FileType::HeatFlowFile
                || problem->filetype == femm::FileType::ElectrostaticsFile )
        {
            // find and write number of conductor property group
            for(j=0,t=0;j<problem->circproplist.size ();j++)
            {
                if(problem->circproplist[j]->CircName==problem->linelist[i]->InConductorName) t = j + 1;
            }
            fprintf(fp,"\t%i",t);
        }

        fprintf(fp,"\n");
    }

    // write out arc segment list
    fprintf(fp,"[NumArcSegments] = %i\n", (int) problem->arclist.size());
    for(i=0; i<problem->arclist.size(); i++)
    {
        for(j=0,t=0; j<problem->lineproplist.size(); j++)
            if(problem->lineproplist[j]->BdryName==problem->arclist[i]->BoundaryMarkerName) t=j+1;

        fprintf( fp,"%i\t%i\t%.17g\t%.17g\t%i\t%i\t%i",
                 problem->arclist[i]->n0,
                 problem->arclist[i]->n1,
                 problem->arclist[i]->ArcLength,
                 problem->arclist[i]->MaxSideLength,
                 t,
                 problem->arclist[i]->Hidden,
                 problem->arclist[i]->InGroup );

        if (problem->filetype == femm::FileType::HeatFlowFile
                || problem->filetype == femm::FileType::ElectrostaticsFile )
        {
            // find and write number of conductor property group
            for(j=0,t=0;j<problem->circproplist.size ();j++)
                if(problem->circproplist[j]->CircName==problem->arclist[i]->InConductorName) t=j+1;
            fprintf(fp,"\t%i",t);
        }
        else if (problem->filetype == femm::FileType::MagneticsFile)
        {
            std::cout << "fmesher.cpp SaveFEMFile, mySideLength: " << problem->arclist[i]->mySideLength << std::endl;
            fprintf(fp,"\t%.17g",problem->arclist[i]->mySideLength);
        }
        fprintf(fp,"\n");
    }

    // write out list of holes;
    for(i=0,j=0; i<problem->labellist.size(); i++)
    {
        if(problem->labellist[i]->BlockTypeName=="<No Mesh>")
        {
            j++;
        }
    }

    fprintf(fp,"[NumHoles] = %i\n",j);
    for(i=0,k=0; i<problem->labellist.size(); i++)
    {
        if(problem->labellist[i]->BlockTypeName=="<No Mesh>")
        {
            fprintf(fp,"%.17g\t%.17g\t%i\n",problem->labellist[i]->x,problem->labellist[i]->y,
                    problem->labellist[i]->InGroup);
            k++;
        }
    }

    fclose(fp);

    return true;
}

bool FMesher::LoadMesh(string PathName)
{
    int i,j,k,q,nl;
    string pathname,rootname,infile;
    FILE *fp;
    char s[1024];

    // clear out the old mesh...
    meshnode.clear();
    meshline.clear();
    greymeshline.clear();

    pathname = PathName;
    if (pathname.length()==0)
    {
        WarnMessage("No mesh to display");
        return false;
    }

    rootname = pathname.substr(0,pathname.find_last_of('.'));

    //read meshnodes;
    infile = rootname + ".node";
    if((fp=fopen(infile.c_str(),"rt"))==NULL)
    {
        WarnMessage("No mesh to display");
        return false;
    }
    fgets(s,1024,fp);
    sscanf(s,"%i",&k);
    meshnode.resize(k);
    CNode node;
    for(i=0; i<k; i++)
    {
        fgets(s,1024,fp);
        sscanf(s,"%i\t%lf\t%lf",&j,&node.x,&node.y);
        meshnode[i] = node.clone();
    }
    fclose(fp);

    //read meshlines;
    infile = rootname + ".edge";
    if((fp=fopen(infile.c_str(),"rt"))==NULL)
    {
        WarnMessage("No mesh to display");
        return false;
    }
    fgets(s,1024,fp);
    sscanf(s,"%i",&k);
    meshline.resize(k);
    fclose(fp);

    infile = rootname + ".ele";
    if((fp=fopen(infile.c_str(),"rt"))==NULL)
    {
        WarnMessage("No mesh to display");
        return false;
    }
    fgets(s,1024,fp);
    sscanf(s,"%i",&k);

    IntPoint segm;
    int n[3],p;
    for(i=0,nl=0; i<k; i++)
    {
        fgets(s,1024,fp);
        sscanf(s,"%i	%i	%i	%i	%i",&q,&n[0],&n[1],&n[2],&j);
        for(q=0; q<3; q++)
        {
            p=q+1;
            if(p==3) p=0;
            if (n[p]>n[q])
            {
                segm.x = n[p];
                segm.y = n[q];

                if (j != 0)
                {
                    meshline[nl++] = MAKE_UNIQUE<femm::IntPoint>(segm);
                }
                else
                {
                    greymeshline.push_back(MAKE_UNIQUE<femm::IntPoint>(segm));
                }
            }
        }
    }
    meshline.resize(nl);
    fclose(fp);

    // clear out temporary files
    infile = rootname + ".ele";
    remove(infile.c_str());
    infile = rootname + ".node";
    remove(infile.c_str());
    infile = rootname + ".edge";
    remove(infile.c_str());
    infile = rootname + ".pbc";
    remove(infile.c_str());
    infile = rootname + ".poly";
    remove(infile.c_str());

    return true;
}

//bool FMesher::ScanPreferences()
//{
//	FILE *fp;
//	CStdString fname;
//
//	fname = BinDir + "femme.cfg";
//
//	fp=fopen(fname,"rt");
//	if (fp!=NULL)
//	{
//		bool flag=false;
//		char s[1024];
//		char q[1024];
//		char *v;
//
//		// parse the file
//		while (fgets(s,1024,fp)!=NULL)
//		{
//			sscanf(s,"%s",q);
//
//			if( _strnicmp(q,"<Precision>",11)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%lf",&d_prec);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<MinAngle>",10)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%lf",&d_minangle);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<Frequency>",11)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%lf",&d_freq);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<Depth>",7)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%lf",&d_depth);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<Coordinates>",13)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%i",&d_coord);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<LengthUnits>",13)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%i",&d_length);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<ProblemType>",13)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%i",&d_type);
//			  q[0] = '\0';
//			}
//
//			if( _strnicmp(q,"<ACSolver>",8)==0)
//			{
//			  v = StripKey(s);
//			  sscanf(v,"%i",&d_solver);
//			  q[0] = '\0';
//			}
//		}
//		fclose(fp);
//			return true;
//	}
//
//	return false;
//}


//--------------------------------------------------------------


//void CFemmeView::lnuMakeMesh()
//{
//	OnMakeMesh();
//}



//void CFemmeView::lnu_purge_mesh()
//{
//	OnPurgemesh();
//}
//
//void CFemmeView::lnu_show_mesh()
//{
//	OnShowMesh();
//}
//
//void CFemmeView::lnu_analyze(int bShow)
//{
//	if (bShow) bLinehook=HiddenLua;
//	else bLinehook=NormalLua;
//	OnMenuAnalyze();
//}
//
//void CFemmeView::lua_zoomnatural()
//{
//	OnZoomNatural();
//}
//
//void CFemmeView::lua_zoomout()
//{
//	OnZoomOut();
//}
//
//void CFemmeView::lua_zoomin()
//{
//	OnZoomIn();
//}

//bool FMesher::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
//{
//	// TODO: Add your specialized code here and/or call the base class
//	if (bLinehook!=false) return true;
//	return CDocument::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
//}

//void FMesher::OnEditExterior()
//{
//	// TODO: Add your command handler code here
//	CExteriorProps dlg;
//
//	dlg.m_Ro=extRo;
//	dlg.m_Ri=extRi;
//	dlg.m_Zo=extZo;
//	if(dlg.DoModal()==IDOK)
//	{
//		extRo=dlg.m_Ro;
//		extRi=dlg.m_Ri;
//		extZo=dlg.m_Zo;
//	}
//}



//void FMesher::OnEditMatprops()
//{
//	// TODO: Add your command handler code here
//	CPtProp pProp;
//
//	pProp.pblockproplist=&blockproplist;
//	pProp.PropType=2;
//	pProp.ProblemType=ProblemType;
//
//	pProp.DoModal();
//}

//void FMesher::OnEditPtprops()
//{
//	// TODO: Add your command handler code here
//	CPtProp pProp;
//
//	pProp.pnodeproplist=&nodeproplist;
//	pProp.PropType=0;
//
//	pProp.DoModal();
//}

//void FMesher::OnEditSegprops()
//{
//	// TODO: Add your command handler code here
//	CPtProp pProp;
//
//	pProp.plineproplist=&lineproplist;
//	pProp.PropType=1;
//	pProp.ProblemType=ProblemType;
//	pProp.DoModal();
//}

//void FMesher::OnEditCircprops()
//{
//	CPtProp pProp;
//
//	pProp.pcircproplist=&circproplist;
//	pProp.ProblemType=ProblemType;
//	pProp.PropType=3;
//
//	pProp.DoModal();
//}

//bool FMesher::OpBlkDlg()
//{
//	int i,j,k,nselected,cselected;
//	double a;
//	COpBlkDlg zDlg;
//
//	zDlg.ProblemType=ProblemType;
//
//	// check to see how many (if any) blocks are selected.
//	for(i=0,k=0,nselected=0,cselected=0;i<blocklist.size();i++)
//	{
//		if (blocklist[i].IsSelected!=0){
//			if(nselected==0){
//				nselected++;
//				zDlg.m_ingroup=blocklist[i].InGroup;
//			}
//			else if(blocklist[i].BlockTypeName!=blocklist[k].BlockTypeName)
//				nselected++;
//			if(cselected==0) cselected++;
//			else if(blocklist[i].InCircuit!=blocklist[k].InCircuit)
//				cselected++;
//			if(blocklist[i].InGroup!=zDlg.m_ingroup) zDlg.m_ingroup=0;
//			k=i;
//			if(blocklist[i].IsExternal) zDlg.m_isexternal=true;
//		}
//	}
//	if (nselected==0) return false;
//
//	// find average block size;
//	for(i=0,j=0,a=0.;i<blocklist.size();i++)
//		if(blocklist[i].IsSelected!=0)
//			if (blocklist[i].MaxArea>a) a=blocklist[i].MaxArea;
//
//	zDlg.m_sidelength=floor(2.e07*sqrt(a/PI)+0.5)/1.e07;
//
//	zDlg.pblockproplist=&blockproplist;
//	zDlg.pcircproplist=&circproplist;
//
//	if (nselected==1){
//		if(blocklist[k].BlockTypeName=="<No Mesh>") zDlg.cursel=1;
//		else for(i=0,zDlg.cursel=0;i<blockproplist.size();i++)
//			if (blockproplist[i].BlockName==blocklist[k].BlockTypeName)
//				zDlg.cursel=i+2;
//	}
//	else zDlg.cursel=0;
//
//	if (cselected==1){
//		for(i=0,zDlg.circsel=0;i<circproplist.size();i++)
//			if (circproplist[i].CircName==blocklist[k].InCircuit)
//				zDlg.circsel=i+1;
//	}
//	else zDlg.circsel=0;
//
//	zDlg.m_magdir=blocklist[k].MagDir;
//	zDlg.m_magdirfctn=blocklist[k].MagDirFctn;
//	zDlg.m_turns=blocklist[k].Turns;
//
//	if (zDlg.DoModal()==IDOK){
//		for(i=0;i<blocklist.size();i++)
//		{
//			if(blocklist[i].IsSelected!=0){
//				blocklist[i].MaxArea=PI*zDlg.m_sidelength*zDlg.m_sidelength/4.;
//				blocklist[i].MagDir=zDlg.m_magdir;
//				blocklist[i].MagDirFctn=zDlg.m_magdirfctn;
//				blocklist[i].Turns=zDlg.m_turns;
//				if (blocklist[i].Turns==0) blocklist[i].Turns++;
//				if (zDlg.cursel==0) blocklist[i].BlockTypeName="<None>";
//				else if(zDlg.cursel==1) blocklist[i].BlockTypeName="<No Mesh>";
//				else blocklist[i].BlockTypeName=blockproplist[zDlg.cursel-2].BlockName;
//				if (zDlg.circsel==0) blocklist[i].InCircuit="<None>";
//				else blocklist[i].InCircuit=circproplist[zDlg.circsel-1].CircName;
//				blocklist[i].InGroup=zDlg.m_ingroup;
//				blocklist[i].IsExternal=zDlg.m_isexternal;
//			}
//		}
//		return true;
//	}
//
//	return false;
//}

//void FMesher::OpNodeDlg()
//{
//	int i,k,nselected;
//	COpNodeDlg zDlg;
//
//	// check to see how many (if any) nodes are selected.
//	for(i=0,k=0,nselected=0;i<nodelist.size();i++)
//	{
//		if (nodelist[i].IsSelected!=0){
//			if(nselected==0){
//				nselected++;
//				zDlg.m_ingroup=nodelist[i].InGroup;
//			}
//			else if(nodelist[i].BoundaryMarkerName!=nodelist[k].BoundaryMarkerName)
//				nselected++;
//			if(nodelist[i].InGroup!=zDlg.m_ingroup) zDlg.m_ingroup=0;
//			k=i;
//		}
//	}
//	if (nselected==0) return;
//
//	zDlg.pnodeproplist=&nodeproplist;
//	if (nselected==1){
//		for(i=0,zDlg.cursel=0;i<nodeproplist.size();i++)
//			if (nodeproplist[i].PointName==nodelist[k].BoundaryMarkerName)
//				zDlg.cursel=i+1;
//	}
//	else zDlg.cursel=0;
//
//	if (zDlg.DoModal()==IDOK){
//		for(i=0;i<nodelist.size();i++)
//		{
//			if(nodelist[i].IsSelected!=0){
//				if (zDlg.cursel==0) nodelist[i].BoundaryMarkerName="<None>";
//				else nodelist[i].BoundaryMarkerName=nodeproplist[zDlg.cursel-1].PointName;
//				nodelist[i].InGroup=zDlg.m_ingroup;
//			}
//		}
//	}
//}

//void FMesher::OpSegDlg()
//{
//	int i,j,k,nselected;
//	COpSegDlg zDlg;
//
//	// check to see how many (if any) nodes are selected.
//	for(i=0,k=0,nselected=0;i<linelist.size();i++)
//	{
//		if (linelist[i].IsSelected!=0){
//			if(nselected==0){
//				nselected++;
//				zDlg.m_ingroup=linelist[i].InGroup;
//			}
//			else if(linelist[i].BoundaryMarkerName!=linelist[k].BoundaryMarkerName)
//				nselected++;
//			if(linelist[i].InGroup!=zDlg.m_ingroup) zDlg.m_ingroup=0;
//			k=i;
//		}
//	}
//	if (nselected==0) return;
//
//	// find properties block size;
//	double ms=0;
//	zDlg.m_hide=false;
//	for(i=0,j=0;i<linelist.size();i++)
//		if(linelist[i].IsSelected!=0){
//			if(linelist[i].MaxSideLength<0) ms=-1;
//			if (ms>=0) ms+=linelist[k].MaxSideLength;
//			j++;
//			if (linelist[i].Hidden==true) zDlg.m_hide=true;
//		}
//	ms/=(double) j;
//
//	zDlg.plineproplist=&lineproplist;
//	if (nselected==1){
//		for(i=0,zDlg.cursel=0;i<lineproplist.size();i++)
//			if (lineproplist[i].BdryName==linelist[k].BoundaryMarkerName)
//				zDlg.cursel=i+1;
//	}
//	else zDlg.cursel=0;
//
//	if (ms<0){
//		zDlg.m_automesh=true;
//		zDlg.m_linemeshsize=0;
//	}
//	else{
//		zDlg.m_automesh=false;
//		zDlg.m_linemeshsize=ms;
//	}
//
//	if (zDlg.DoModal()==IDOK){
//		for(i=0;i<linelist.size();i++)
//		{
//			if(linelist[i].IsSelected!=0){
//
//				if (zDlg.m_automesh==true) linelist[i].MaxSideLength=-1;
//				else{
//					if (zDlg.m_linemeshsize>0)
//						linelist[i].MaxSideLength=zDlg.m_linemeshsize;
//					else zDlg.m_linemeshsize=-1;
//				}
//				if (zDlg.cursel==0) linelist[i].BoundaryMarkerName="<None>";
//				else linelist[i].BoundaryMarkerName=lineproplist[zDlg.cursel-1].BdryName;
//
//				linelist[i].Hidden=zDlg.m_hide;
//				linelist[i].InGroup=zDlg.m_ingroup;
//			}
//		}
//	}
//}

//void FMesher::OpArcSegDlg()
//{
//	int i,j,k,nselected;
//	COpArcSegDlg zDlg;
//
//	// check to see how many (if any) arcs are selected.
//
//
//	for(i=0,k=0,nselected=0;i<arclist.size();i++)
//	{
//		if (arclist[i].IsSelected!=0){
//			if(nselected==0){
//				nselected++;
//				zDlg.m_ingroup=arclist[i].InGroup;
//			}
//			else if(arclist[i].BoundaryMarkerName!=arclist[k].BoundaryMarkerName)
//				nselected++;
//			if(arclist[i].InGroup!=zDlg.m_ingroup) zDlg.m_ingroup=0;
//			k=i;
//		}
//	}
//	if (nselected==0) return;
//
//	// find properties block size;
//	double ms=0;
//	zDlg.m_hide=false;
//	for(i=0,j=0;i<arclist.size();i++)
//		if(arclist[i].IsSelected!=0){
//			ms+=arclist[k].MaxSideLength;
//			j++;
//			if(arclist[i].Hidden==true) zDlg.m_hide=true;
//		}
//
//	ms/=(double) j;
//
//	zDlg.plineproplist=&lineproplist;
//	if (nselected==1){
//		zDlg.m_MaxSeg=ms;
//		for(i=0,zDlg.cursel=0;i<lineproplist.size();i++)
//			if (lineproplist[i].BdryName==arclist[k].BoundaryMarkerName)
//				zDlg.cursel=i+1;
//	}
//	else{
//		zDlg.cursel=0;
//		zDlg.m_MaxSeg=ms;
//	}
//
//	if (zDlg.DoModal()==IDOK){
//		for(i=0;i<arclist.size();i++)
//		{
//			if(arclist[i].IsSelected!=0){
//				if (zDlg.cursel==0) arclist[i].BoundaryMarkerName="<None>";
//				else arclist[i].BoundaryMarkerName=lineproplist[zDlg.cursel-1].BdryName;
//				arclist[i].MaxSideLength=zDlg.m_MaxSeg;
//				arclist[i].Hidden=zDlg.m_hide;
//				arclist[i].InGroup=zDlg.m_ingroup;
//			}
//		}
//	}
//}

//void FMesher::OpGrpDlg()
//{
//	COpGrp dlg;
//	bool bFlag=false;
//	int nsel=0;
//	int i;
//
//	dlg.m_ingroup=0;
//
//	for(i=0;i<arclist.size();i++)
//		if (arclist[i].IsSelected!=0){
//			nsel++;
//			if(dlg.m_ingroup!=arclist[i].InGroup)
//			{
//				if(bFlag==false){
//					dlg.m_ingroup=arclist[i].InGroup;
//					bFlag=true;
//				}
//				else dlg.m_ingroup=0;
//			}
//		}
//
//	for(i=0;i<linelist.size();i++)
//		if (linelist[i].IsSelected!=0)
//		{
//			nsel++;
//			if(dlg.m_ingroup!=linelist[i].InGroup)
//			{
//				if(bFlag==false){
//					dlg.m_ingroup=linelist[i].InGroup;
//					bFlag=true;
//				}
//				else dlg.m_ingroup=0;
//			}
//		}
//
//	for(i=0;i<blocklist.size();i++)
//		if (blocklist[i].IsSelected!=0)
//		{
//			nsel++;
//			if(dlg.m_ingroup!=blocklist[i].InGroup)
//			{
//				if(bFlag==false){
//					dlg.m_ingroup=blocklist[i].InGroup;
//					bFlag=true;
//				}
//				else dlg.m_ingroup=0;
//			}
//		}
//
//	for(i=0;i<nodelist.size();i++)
//		if (nodelist[i].IsSelected!=0)
//		{
//			nsel++;
//			if(dlg.m_ingroup!=nodelist[i].InGroup)
//			{
//				if(bFlag==false){
//					dlg.m_ingroup=nodelist[i].InGroup;
//					bFlag=true;
//				}
//				else dlg.m_ingroup=0;
//			}
//		}
//	if (nsel==0) return;
//
//	if(dlg.DoModal()==IDOK){
//
//		for(i=0;i<nodelist.size();i++)
//			if(nodelist[i].IsSelected!=0)
//				nodelist[i].InGroup=dlg.m_ingroup;
//
//		for(i=0;i<linelist.size();i++)
//			if(linelist[i].IsSelected!=0)
//				linelist[i].InGroup=dlg.m_ingroup;
//
//		for(i=0;i<arclist.size();i++)
//			if(arclist[i].IsSelected!=0)
//				arclist[i].InGroup=dlg.m_ingroup;
//
//		for(i=0;i<blocklist.size();i++)
//			if(blocklist[i].IsSelected!=0)
//				blocklist[i].InGroup=dlg.m_ingroup;
//
//	}
//
//}

string fmesher::triangleVersionString()
{
    return TRIANGLE_VERSION;
}
