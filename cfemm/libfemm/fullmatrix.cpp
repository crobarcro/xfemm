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

#include <cstdlib>
#include <math.h>
#if defined(__APPLE__)
#include <stdlib.h>  // macOS has no malloc.h
#else
#include "malloc.h"
#endif
#include "femmcomplex.h"
#include "fullmatrix.h"


CFullMatrix::CFullMatrix()
{
    n=0;
    M=NULL;
    b=NULL;
}

CFullMatrix::CFullMatrix(int d)
{
    n=0;
    M=NULL;
    b=NULL;
    Create(d);
}

CFullMatrix::~CFullMatrix()
{
    if(n==0) return;

    int i;
    for(i=0; i<n; i++) free(M[i]);
    free(M);
    free(b);
    n=0;
}

void CFullMatrix::Wipe()
{
    /* fills up a square matrix with zeros in every entry */
    int i,j;

    for(i=0; i<n; i++)
    {
        b[i]=0;
        for(j=0; j<n; j++)
            M[i][j]=0;
    }
}


int CFullMatrix::Create(int d)
{
    int i;

    M=(double **)calloc(d,sizeof(double *));
    b=(double *)calloc(d,sizeof(double));
    if ((M==NULL) || (b==NULL)) return false;
    for(i=0; i<d; i++)
    {
        M[i]=(double *)calloc(d,sizeof(double));
        if (M[i]==NULL) return false;
    }
    n=d;
    return true;
}

int CFullMatrix::GaussSolve()
{
    /* solves the linear system m x = b for x.  The result is returned
       in b, m is destroyed in the process */

    int i,j,k,q = 0;
    double *z;
    double max,f;

    for(i=0; i<n; i++)
    {
        for(j=i,max=0; j<n; j++)
            if (fabs(M[j][i])>fabs(max))
            {
                max=M[j][i];
                q=j;
            }
        if(max==0) return false;
        z=M[i];
        M[i]=M[q];
        M[q]=z;
        f=b[i];
        b[i]=b[q];
        b[q]=f;
        for(j=i+1; j<n; j++)
        {
            f=M[j][i]/M[i][i];
            b[j]=b[j]-f*b[i];
            for (k=i; k<n; k++)
                M[j][k]-=(f*M[i][k]);
        }
    }

    for(i=n-1; i>=0; i--)
    {
        for(j=n-1,f=0; j>i; j--)
            f+=M[i][j]*b[j];
        b[i]=(b[i]-f)/M[i][i];
    }
    return true;
}

CComplexFullMatrix::CComplexFullMatrix()
{
    n=0;
    M=NULL;
    b=NULL;
}

CComplexFullMatrix::CComplexFullMatrix(int d)
{
    n=0;
    M=NULL;
    b=NULL;
    Create(d);
}

CComplexFullMatrix::~CComplexFullMatrix()
{
    if(n==0) return;

    int i;
    for(i=0; i<n; i++) free(M[i]);
    free(M);
    free(b);
    n=0;
}

void CComplexFullMatrix::Wipe()
{
    /* fills up a square matrix with zeros in every entry */
    int i,j;

    for(i=0; i<n; i++)
    {
        b[i]=0;
        for(j=0; j<n; j++)
            M[i][j]=0;
    }
}


int CComplexFullMatrix::Create(int d)
{
    int i;

    M=(CComplex **)calloc(d,sizeof(CComplex *));
    b=(CComplex *)calloc(d,sizeof(CComplex));
    if ((M==NULL) || (b==NULL)) return false;
    for(i=0; i<d; i++)
    {
        M[i]=(CComplex *)calloc(d,sizeof(CComplex));
        if (M[i]==NULL) return false;
    }
    n=d;
    return true;
}

int CComplexFullMatrix::GaussSolve()
{
    /* solves the linear system m x = b for x.  The result is returned
       in b, m is destroyed in the process */

    int i,j,k,q = 0;
    CComplex max,f;

    for(i=0; i<n; i++)
    {
        for(j=i,max=0; j<n; j++)
            if (abs(M[j][i])>abs(max))
            {
                max=M[j][i];
                q=j;
            }
        if(max==0) return false;
        std::swap(M[i],M[q]);
        std::swap(b[i],b[q]);
        for(j=i+1; j<n; j++)
        {
            f=M[j][i]/M[i][i];
            b[j]=b[j]-f*b[i];
            for (k=i; k<n; k++)
                M[j][k]-=(f*M[i][k]);
        }
    }

    for(i=n-1; i>=0; i--)
    {
        for(j=n-1,f=0; j>i; j--)
            f+=M[i][j]*b[j];
        b[i]=(b[i]-f)/M[i][i];
    }
    return true;
}



