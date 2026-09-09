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

#ifndef SPARS_H
#define SPARS_H

class CEntry
{
public:

    double x;				// value stored in the entry
    int c;					// column that the entry lives in
    CEntry *next;			// pointer to next entry in row;
    CEntry();

private:
};


class CBigLinProb
{
public:

    // data members

    double *V;				// solution
    double *P;				// search direction;
    double *R;				// residual;
    double *U;				// A * P;
    double *Z;
    double *b;				// RHS of linear equation
    CEntry **M;				// pointer to list of matrix entries;
    int n;					// dimensions of the matrix;
    int bdw;				// Optional matrix bandwidth parameter;
    double Precision;		// error tolerance for solution
    double Lambda;			// relaxation factor;

    int *Q; ///< Used by esolver and hsolver.

    // member functions

    // constructor
    CBigLinProb();
    // destructor
    ~CBigLinProb();
    virtual int Create(int d, int bw);	// initialize the problem
    void Put(double v, int p, int q);
    // use to create/set entries in the matrix
    double Get(int p, int q);
    /** Preconditioned conjugate gradients.
     *  @param flag non-zero if a starting guess for V is present.
     *  @param maxIterations iteration ceiling; 0 selects the default.
     *  @return false if the iteration broke down, or ran out of
     *          iterations without reaching Precision.
     */
    bool PCGSolve(int flag, int maxIterations = 0);
    void MultPC(const double *X, double *Y);
    void AddTo(double v, int p, int q);
    void MultA(double *X, double *Y);
    void SetValue(int i, double x);
    void Periodicity(int i, int j);
    void AntiPeriodicity(int i, int j);
    void Wipe();
    double Dot(double *X, double *Y);
    void ComputeBandwidth();

//		CFknDlg *TheView;

private:

};

#endif
