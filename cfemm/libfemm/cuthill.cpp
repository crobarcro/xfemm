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
*/

// does Cuthill-McKee algorithm as described in Hoole;

#include<stdio.h>
#include<math.h>
#include "malloc.h"
#include "femmcomplex.h"
#include "femmconstants.h"
#include "femmenums.h"
//#include "spars.h"
#include "feasolver.h"

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
int FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::SortElements()
{
    // Comb Sort -- see http://en.wikipedia.org/wiki/Comb_sort
    int *Score;
    int i,j,k,gap;

    Score=(int*)calloc(NumEls,sizeof(int));

    for(k=0; k<NumEls; k++)
    {
        Score[k]=meshele[k].p[0]+meshele[k].p[1]+meshele[k].p[2];
    }

    gap = NumEls;

    do
    {
        //update the gap value for a next comb
        if (gap > 1)
        {
            gap=(gap*10)/13;
            if ((gap==10) || (gap==9)) gap=11;

        }

        //a single "comb" over the input list
        for(j=0,i=0; (j+gap)<NumEls; j++)
        {
            if (Score[j]>Score[j+gap])
            {
                k=j+gap;
                i=Score[k];
                Score[k]=Score[j];
                Score[j]=i;
                std::swap(meshele[k],meshele[j]);
                i=1;
            }
        }
    }
    while((gap>1)&&(i>0));


    free(Score);
    return true;
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
int FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::Cuthill(bool deletefiles)
{
    char infile[256];
    sprintf(infile, "%s.edge", PathName.c_str());
    FILE *fp = fopen(infile, "rt");
    if (fp == nullptr) {
        printf("Couldn't open %s", infile);
        return false;
    }
    long int nLines, marker;
    if (fscanf(fp, "%li%li", &nLines, &marker) != 2) {
        fclose(fp);
        return false;
    }
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(static_cast<std::size_t>(nLines));
    for (long int i = 0; i < nLines; ++i) {
        long int id, boundary;
        int first, second;
        if (fscanf(fp, "%li%i%i%li", &id, &first, &second, &boundary) != 4 ||
            first < 0 || second < 0) {
            fclose(fp);
            return false;
        }
        edges.emplace_back(static_cast<std::size_t>(first),
                           static_cast<std::size_t>(second));
    }
    fclose(fp);
    if (deletefiles)
        remove(infile);
    return Cuthill(edges);
}

template< class PointPropT
          , class BoundaryPropT
          , class BlockPropT
          , class CircuitPropT
          , class BlockLabelT
          , class MeshElementT
          >
int FEASolver<PointPropT,BoundaryPropT,BlockPropT,CircuitPropT,BlockLabelT,MeshElementT>
::Cuthill(const std::vector<std::pair<std::size_t, std::size_t>> &edges)
{
    int i, n0, n1, n, newwide;
    long int j;
    const long int n_lines = static_cast<long int>(edges.size());
    std::vector<std::vector<int>> ocon(NumNodes);
    std::vector<int> newnum(NumNodes, -1), numcon(NumNodes, 0), nxtnum(NumNodes, 0);

    for (const auto &edge : edges) {
        if (edge.first >= static_cast<std::size_t>(NumNodes) ||
            edge.second >= static_cast<std::size_t>(NumNodes))
            return false;
        ++numcon[edge.first];
        ++numcon[edge.second];
    }
    for (i = 0; i < NumNodes; ++i)
        ocon[i].reserve(numcon[i]);
    for (const auto &edge : edges) {
        ocon[edge.first].push_back(static_cast<int>(edge.second));
        ocon[edge.second].push_back(static_cast<int>(edge.first));
    }

    // sort connections in order of increasing connectivity;
    // I'm lazy, so I'm doing a bubble sort;
    for(n0=0; n0<NumNodes; n0++)
    {
        for(i=1; i<numcon[n0]; i++)
            for(j=1; j<numcon[n0]; j++)
                if(numcon[ocon[n0][j]]<numcon[ocon[n0][j-1]])
                {
                    n1=ocon[n0][j];
                    ocon[n0][j]=ocon[n0][j-1];
                    ocon[n0][j-1]=n1;
                }
    }


    // search for a node to start with;
    j = numcon[0];
    n0 = 0;
    for(i=1; i<NumNodes; i++)
    {
        if(numcon[i]<j)
        {
            j=numcon[i];
            n0=i;
        }
        if(j==2) i=n_lines;	// break out if j==2,
        // because this is the best we can do
    }

    // do renumbering algorithm;
    for(i=0; i<NumNodes; i++) nxtnum[i]=-1;
    newnum[n0]=0;
    n=1;
    nxtnum[0]=n0;

    do
    {
        // renumber in order of increasing number of connections;

        for(i=0; i<numcon[n0]; i++)
        {
            if (newnum[ocon[n0][i]]<0)
            {
                newnum[ocon[n0][i]]=n;
                nxtnum[n]=ocon[n0][i];
                n++;
            }
        }

        // need to catch case in which problem is multiply
        // connected and still renumber right.
        if(nxtnum[newnum[n0]+1]<0)
        {
            //	WarnMessage("Multiply Connected!");
            //	exit(0);

            // first, get a node that hasn't been visited yet;
            for(i=0; i<NumNodes; i++)
                if(newnum[i]<0)
                {
                    j=numcon[i];
                    n0=i;
                    break;
                }


            // now, get a new starting node;
            for(i=0; i<NumNodes; i++)
            {
                if((newnum[i]<0) && (numcon[i]<j))
                {
                    j=numcon[i];
                    n0=i;
                }
                if(j==2) break;	// break out if j==2,
                // because this is the
                // best we can do
            }

            // now, set things to restart;
            newnum[n0]=n;
            nxtnum[n]=n0;
            n++;
        }
        else n0=nxtnum[newnum[n0]+1];


    }
    while(n<NumNodes);

    // remap connectivities;
    for(i=0; i<NumNodes; i++)
        for(j=0; j<numcon[i]; j++)
            ocon[i][j]=newnum[ocon[i][j]];

    // remap (anti)periodic boundary points
    for(i=0; i<NumPBCs; i++)
    {
        pbclist[i].x=newnum[pbclist[i].x];
        pbclist[i].y=newnum[pbclist[i].y];
    }

	// remap air gap element information
	for(i=0; i<NumAirGapElems; i++)
	{
		for(auto &point : agelist[i].innerRingTopology) point.n0=newnum[point.n0];
		for(auto &point : agelist[i].outerRingTopology) point.n0=newnum[point.n0];
		for(auto &node : agelist[i].nodeNums) node=newnum[node];
		for(int k=0; k<=agelist[i].totalArcElements; k++)
		{
			agelist[i].quadNode[k].n0=newnum[agelist[i].quadNode[k].n0];
			agelist[i].quadNode[k].n1=newnum[agelist[i].quadNode[k].n1];
			agelist[i].quadNode[k].n2=newnum[agelist[i].quadNode[k].n2];
			agelist[i].quadNode[k].n3=newnum[agelist[i].quadNode[k].n3];
		}
	}

    // find new bandwidth;

    // PBCs fuck up the banding, som could have to do
    // something like:
    // if(NumPBCs!=0) BandWidth=0;
    // else{
    // but if we apply the PCBs the last thing before the
    // solver is called, we can take advantage of banding
    // speed optimizations without messing things up.
    for(n0=0,newwide=0; n0<NumNodes; n0++)
    {
        for(i=0; i<numcon[n0]; i++)
            if(abs(newnum[n0]-ocon[n0][i])>newwide)
            {
                newwide=abs(newnum[n0]-ocon[n0][i]);
            }
    }

    BandWidth=newwide+1;
    // }

    // free up the variables that we needed during the routine....
    //free(numcon);
    //free(nxtnum);
    //free(ocon[0]);
    //free(ocon);

    // new mapping remains in newnum;
    // apply this mapping to elements first.
    for(i=0; i<NumEls; i++)
        for(j=0; j<3; j++)
            meshele[i].p[j]=newnum[meshele[i].p[j]];

//    // now, sort nodes based on newnum;
//    for(i=0; i<NumNodes; i++)
//    {
//        while(newnum[i]!=i)
//        {
//            CNode swap;
//
//            j=newnum[i];
//            n=newnum[j];
//            newnum[j]=newnum[i];
//            newnum[i]=n;
//            swap=meshnode[j];
//            meshnode[j]=meshnode[i];
//            meshnode[i]=swap;
//        }
//    }

    // virtual method that must be overridden by child classes
    // as the mesh nodes class type varies
    SortNodes (newnum);

    //free(newnum);

    SortElements();

    return true;
}
