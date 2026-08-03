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
#include "CMaterialProp.h"

#include "fullmatrix.h"
#include "femmcomplex.h"
#include "femmconstants.h"
#include "fparse.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctype.h>
#include <istream>
#include <stdexcept>

#ifdef DEBUG_MEX
#include "mex.h"
#endif // DEBUG_MEX

#define ElementsPerSkinDepth 10

#ifdef DEBUG_CMATERIALPROP
#define debug std::cerr << __func__ << "(): "
#else
#define debug while(false) std::cerr
#endif

using namespace std;
using namespace femm;

CMaterialProp::CMaterialProp()
    : BlockName("New Material")
{
}

CMaterialProp::~CMaterialProp()
{
}

CMMaterialProp::CMMaterialProp()
    : CMaterialProp()
    , mu_x(1.)
    , mu_y(1.)
    , BHpoints(0)
    , Bdata()
    , Hdata()
    , slope()
    , LamType(0)
    , LamFill(1.)
    , H_c(0.)
    , Nrg(0.)
    , J(0.)
    , Cduct(0.)
    , Lam_d(0.)
    , Theta_hn(0.)
    , Theta_hx(0.)
    , Theta_hy(0.)
    , NStrands(0)
    , WireD(0)
    , mu_fdx()
    , mu_fdy()
    , Frequency(0.)
{
}

CMMaterialProp::~CMMaterialProp()
{
}

CMMaterialProp::CMMaterialProp( const CMMaterialProp& other )
{
    BlockName = other.BlockName;
    mu_x = other.mu_x;
    mu_y = other.mu_y;            // permeabilities, relative
    BHpoints = other.BHpoints;

    if (other.BHpoints !=(int) other.Bdata.size())
        std::cerr << "Size difference:" << other.BHpoints << " vs " << other.Bdata.size() << std::endl;
    Bdata = other.Bdata;
    Hdata = other.Hdata;
    slope = other.slope;

    H_c = other.H_c;                // magnetization, A/m
    Nrg = other.Nrg;
    J = other.J;
    Cduct = other.Cduct;            // conductivity of the material, MS/m
    Lam_d = other.Lam_d;            // lamination thickness, mm
    Theta_hn = other.Theta_hn;        // hysteresis angle, degrees
    Theta_hx = other.Theta_hx;        // hysteresis angle, degrees
    Theta_hy = other.Theta_hy;        // hysteresis angle, degrees
    NStrands = other.NStrands;            // number of strands per wire
    WireD = other.WireD;
    LamFill = other.LamFill;            // lamination fill factor;
    LamType = other.LamType;            // type of lamination;
}

void CMMaterialProp::clearSlopes()
{
    slope.clear();
}

void CMMaterialProp::GetSlopes(double omega)
{
    if (BHpoints < 0 || BHpoints == 1)
        throw std::invalid_argument("B-H curve must contain either zero points or at least two points");
    if (Bdata.size() != static_cast<size_t>(BHpoints)
            || Hdata.size() != static_cast<size_t>(BHpoints))
        throw std::invalid_argument("B-H curve point count does not match its B and H data sizes");
    for (int point = 0; point < BHpoints; ++point) {
        if (!std::isfinite(Bdata[point]) || !std::isfinite(Hdata[point].re)
                || !std::isfinite(Hdata[point].im))
            throw std::invalid_argument("B-H curve contains a non-finite B or H component at point "
                                        + std::to_string(point));
        if (point > 0 && !(Bdata[point] > Bdata[point - 1]))
            throw std::invalid_argument("B-H curve B values must be strictly increasing (invalid interval at point "
                                        + std::to_string(point) + ")");
    }
    if (BHpoints==0) return; // catch trivial case;
    if (!slope.empty()) return; // already have computed the slopes;

    int i,k;
    bool CurveOK=false;
    bool ProcessedLams=false;
    CComplexFullMatrix L;
    double l1,l2;
    CComplex *hn;
    double *bn;
    CComplex mu;

    L.Create(BHpoints);
    bn   =(double *)  calloc(BHpoints,sizeof(double));
    hn   =(CComplex *)calloc(BHpoints,sizeof(CComplex));
    slope.reserve(BHpoints);


    // strip off some info that we can use during the first
    // nonlinear iteration;
    mu_x = Bdata[1] / (muo*abs(Hdata[1]));
    mu_y = mu_x;
    Theta_hx = Theta_hn;
    Theta_hy = Theta_hn;

    // first, we need to doctor the curve if the problem is
    // being evaluated at a nonzero frequency.
    if(omega!=0)
    {
        debug << "fixing the curve for a harmonic problem.\n";
        // Make an effective B-H curve for harmonic problems.
        // this one convolves B(H) where H is sinusoidal with
        // a sine at the same frequency to get the effective
        // amplitude of B
        double munow,mumax=0;
        for(i=1;i<BHpoints;i++)
        {
            hn[i]=Hdata[i];
            for(k=1,bn[i]=0;k<=i;k++)
            {
                bn[i]+=Re((4.*(Hdata[k]*Bdata[k-1] -
                           Hdata[k-1]*Bdata[k])*(-cos((Hdata[k-1]*PI)/(2.*Hdata[i])) +
                        cos((Hdata[k]*PI)/(2.*Hdata[i]))) + (-Bdata[k-1] +
                        Bdata[k])*((Hdata[k-1] - Hdata[k])*PI +
                        Hdata[i]*(-sin((Hdata[k-1]*PI)/Hdata[i]) +
                        sin((Hdata[k]*PI)/Hdata[i]))))/ ((Hdata[k-1] - Hdata[k])*PI));
            }
        }

        for(i=1;i<BHpoints;i++)
        {
            Bdata[i]=bn[i];
            Hdata[i]=hn[i];
            munow=Re(Bdata[i]/Hdata[i]);
            if (munow>mumax) mumax=munow;
        }

        // apply complex permeability to approximate the
        // effects of hysteresis.  We will follow a kludge
        // suggested by O'Kelly where hysteresis angle
        // is proportional to permeability.  This implies
        // that loss goes with B^2
        for(i=1;i<BHpoints;i++)
        {
            Hdata[i]*=exp(I*Bdata[i]*Theta_hn*DEG/(Hdata[i]*mumax));
        }
        MuMax = mumax / muo;
    }

    const int maxCurveRepairIterations = 100;
    int curveRepairIterations = 0;
    while(CurveOK!=true)
    {
        if (++curveRepairIterations > maxCurveRepairIterations) {
            free(bn);
            free(hn);
            throw std::runtime_error("B-H curve repair did not converge after "
                                     + std::to_string(maxCurveRepairIterations) + " iterations");
        }
        for (i = 1; i < BHpoints; ++i) {
            if (!std::isfinite(Bdata[i]) || !std::isfinite(Hdata[i].re)
                    || !std::isfinite(Hdata[i].im) || !(Bdata[i] > Bdata[i - 1])) {
                free(bn);
                free(hn);
                throw std::runtime_error("B-H curve processing produced an invalid interpolation interval at point "
                                         + std::to_string(i));
            }
        }
        slope.clear();
        debug << "curve not ok yet.\n";
        // make sure that the space for computing slopes is cleared out
        L.Wipe();

        // impose natural BC on the `left'
        l1=Bdata[1]-Bdata[0];
        L.M[0][0]=4./l1;
        L.M[0][1]=2./l1;
        L.b[0]=6.*(Hdata[1]-Hdata[0])/(l1*l1);

        // impose natural BC on the `right'
        int n=BHpoints;
        l1=Bdata[n-1]-Bdata[n-2];
        L.M[n-1][n-1]=4./l1;
        L.M[n-1][n-2]=2./l1;
        L.b[n-1]=6.*(Hdata[n-1]-Hdata[n-2])/(l1*l1);

        for(i=1;i<BHpoints-1;i++)
        {
            l1=Bdata[i]-Bdata[i-1];
            l2=Bdata[i+1]-Bdata[i];

            L.M[i][i-1]=2./l1;
            L.M[i][i]=4.*(l1+l2)/(l1*l2);
            L.M[i][i+1]=2./l2;

            L.b[i]=6.*(Hdata[i]-Hdata[i-1])/(l1*l1) +
                    6.*(Hdata[i+1]-Hdata[i])/(l2*l2);
        }

        L.GaussSolve();
        for(i=0;i<BHpoints;i++) slope.push_back(L.b[i]);

        // now, test to see if there are any "bad" segments in there.
        // it is probably sufficient to do this test just on the
        // real part of the BH curve...
        for(i=1,CurveOK=true;i<BHpoints;i++)
        {
            double L,c0,c1,c2,d0,d1,u0,u1,X0,X1;

            // it is probably sufficient to do this test just on the
            // real part of the BH curve.  We do the test on just the
            // real part of the curve by taking the real parts of
            // slope and Hdata
            d0=slope[i-1].re;
            d1=slope[i].re;
            u0=Hdata[i-1].re;
            u1=Hdata[i].re;
            L=Bdata[i]-Bdata[i-1];

            c0=d0;
            c1= -(2.*(2.*d0*L + d1*L + 3.*u0 - 3.*u1))/(L*L);
            c2= (3.*(d0*L + d1*L + 2.*u0 - 2.*u1))/(L*L*L);
            X0=-1.;
            X1=-1.;

            u0=c1*c1-4.*c0*c2;
            // check out degenerate cases
            if(c2==0)
            {
                if(c1!=0) X0=-c0/c1;
            }
            else if(u0>0)
            {
                u0=sqrt(u0);
                X0=-(c1 + u0)/(2.*c2);
                X1=(-c1 + u0)/(2.*c2);
            }

            //now, see if we've struck gold!
            if (((X0>=0.)&&(X0<=L))||((X1>=0.)&&(X1<=L)))
                CurveOK=false;
        }

        if(CurveOK!=true)  //remedial action
        {
            // Smooth out input points
            // to get rid of rapid transitions;
            // Uses a 3-point moving average
            for(i=1;i<BHpoints-1;i++)
            {
                bn[i]=(Bdata[i-1]+Bdata[i]+Bdata[i+1])/3.;
                hn[i]=(Hdata[i-1]+Hdata[i]+Hdata[i+1])/3.;
            }

            bool madeProgress = false;
            for(i=1;i<BHpoints-1;i++)
                madeProgress = madeProgress || bn[i] != Bdata[i] || hn[i] != Hdata[i];
            if (!madeProgress) {
                free(bn);
                free(hn);
                throw std::runtime_error("B-H curve repair cannot make progress");
            }

            for(i=1;i<BHpoints-1;i++){
                Hdata[i]=hn[i];
                Bdata[i]=bn[i];
            }
        }


        if((CurveOK==true) && (ProcessedLams==false))
        {
            // if the material is laminated and has a non-zero conductivity,
            // we have to do some work to find the "right" apparent BH
            // curve for the material for AC problems.  Essentially, we have
            // to solve a 1-D finite element problem at each B-H point to
            // figure out how much flux the lamination is really carrying.
            if((omega>0) && (Lam_d!=0) && (Cduct!=0))
            {
                // Calculate a new apparent b and h for each point on the
                // curve to account for the laminations.
                for(i=1;i<BHpoints;i++)
                {
                    mu=LaminatedBH(omega,i);
                    bn[i]=abs(mu*Hdata[i]);
                    hn[i]=bn[i]/mu;
                }

                // Replace the original BH points with the apparent ones
                for(i=1;i<BHpoints;i++)
                {
                    Bdata[i]=bn[i];
                    Hdata[i]=hn[i];
                }

                // Make the program check the consistency and recompute the
                // proper slopes of the new BH curve one more time;
                CurveOK=false;
            }

            // take care of LamType=0 situation by changing the apparent B-H curve.
            if((LamType==0) && (LamFill!=1))
            {
                // this is changed from the previous version so that
                // a complex-valued H can be properly accommodated
                // in the algorithm.
                for(i=1;i<BHpoints;i++)
                {
                    mu=LamFill*Bdata[i]/Hdata[i]+(1.-LamFill)*muo;
                    Bdata[i]=abs(mu*Hdata[i]);
                    Hdata[i]=Bdata[i]/mu;
                }
                // Make the program check the consistency and recompute the
                // proper slopes of the new BH curve one more time;
                CurveOK=false;
            }

            ProcessedLams=true;
        }

    }

    free(bn);
    free(hn);
    return;
}


CComplex CMMaterialProp::LaminatedBH(double w, int i)
{
    debug << "CMMaterialProp::LaminatedBH("<<w<<", "<<i<<")\n";
    int k,n,iter=0;
    CComplex *m0,*m1,*b,*x;
    double L,o,d,ds,B,lastres;
    double res=0;
    double Relax=1;
    CComplex mu,vo,vi,c,H;
    CComplex Md,Mo;
    bool Converged=false;

    // Base the required element spacing on the skin depth
    // at the surface of the material
    mu=Bdata[i]/Hdata[i];
    o=Cduct*1.e6;
    d=(Lam_d*0.001)/2.;
    ds=sqrt(2/(w*o*abs(mu)));
    n= ElementsPerSkinDepth * ((int) ceil(d/ds));
    L=d/((double) n);

    x =(CComplex *)calloc(n+1,sizeof(CComplex));
    b =(CComplex *)calloc(n+1,sizeof(CComplex));
    m0=(CComplex *)calloc(n+1,sizeof(CComplex));
    m1=(CComplex *)calloc(n+1,sizeof(CComplex));

    do{
        // make sure that the old stuff is wiped out;
        for(k=0;k<=n;k++)
        {
            m0[k]=0;
            m1[k]=0;
            b[k] =0;
        }

        // build matrix
        for(k=0;k<n;k++)
        {
            if(iter!=0)
            {
                B=abs(x[k+1]-x[k])/L;
                vi=GetdHdB(B);
                vo=GetH(CComplex(B))/B;
            }
            else
            {
                vo=1./mu;
                vi=1./mu;
            }

            Md = ( (vi+vo)/(2.*L) + I*w*o*L/4.);
            Mo = (-(vi+vo)/(2.*L) + I*w*o*L/4.);
            m0[k]   += Md;
            m1[k]   += Mo;
            m0[k+1] += Md;

            Md = ( (vi-vo)/(2.*L));
            Mo = (-(vi-vo)/(2.*L));
            b[k]    += (Md*x[k] + Mo*x[k+1]);
            b[k+1]  += (Mo*x[k] + Md*x[k+1]);
        }

        // set boundary conditions
        m1[0] = 0;
        b[0]  = 0;
        b[n] += Hdata[i];

        // solve tridiagonal equation
        // solution ends up in b;
        for(k = 0; k < n; k++)
        {
            c = m1[k]/m0[k];
            m0[k+1] -= m1[k]*c;
            b[k+1]  -= b[k]*c;
        }
        b[n] = b[n]/m0[n];
        for(k = n-1; k >= 0; k--)
            b[k] = (b[k] - m1[k]*b[k + 1])/m0[k];

        iter++;

        lastres=res;
        res=abs(b[n]-x[n])/d;

        if (res<1.e-8) Converged=true;

        // Do the same relaxation scheme as is implemented
        // in the solver to make sure that this effective
        // lamination permeability calculation converges
        if(iter>5)
        {
            if ((res>lastres) && (Relax>0.1)) Relax/=2.;
            else Relax+= 0.1 * (1. - Relax);
        }

        for(k=0;k<=n;k++) x[k]=Relax*b[k]+(1.0-Relax)*x[k];

    }while(Converged!=true);


    mu = x[n]/(Hdata[i]*d);

    free(x );
    free(b );
    free(m0);
    free(m1);

    return mu;
}

CComplex CMMaterialProp::GetdHdB(const double B) const
{
    double b,z,l;
    CComplex h;
    int i;

    b=fabs(B);

    if(BHpoints==0)    return CComplex(b/(mu_x*muo));

    if(b>Bdata[BHpoints-1])
        return slope[BHpoints-1];

    for(i=0;i<BHpoints-1;i++)
        if((b>=Bdata[i]) && (b<=Bdata[i+1])){
            l=(Bdata[i+1]-Bdata[i]);
            z=(b-Bdata[i])/l;
            h=6.*z*(z-1.)*Hdata[i]/l +
                    (1.-4.*z+3.*z*z)*slope[i] +
                    6.*z*(1.-z)*Hdata[i+1]/l +
                    z*(3.*z-2.)*slope[i+1];
            return h;
        }

    return CComplex(0);
}

double CMMaterialProp::GetH(const double x) const
{
    return Re(GetH(CComplex(x)));
}

CComplex CMMaterialProp::GetH(const CComplex x) const
{
    double b,z,z2,l;
    CComplex p,h;
    int i;

    b=abs(x);
    if((BHpoints==0) || (b==0))    return 0;
    p=x/b;

    if(b>Bdata[BHpoints-1])
        return p*(Hdata[BHpoints-1] + slope[BHpoints-1]*(b-Bdata[BHpoints-1]));

    for(i=0;i<BHpoints-1;i++)
        if((b>=Bdata[i]) && (b<=Bdata[i+1])){
            l=Bdata[i+1]-Bdata[i];
            z=(b-Bdata[i])/l;
            z2=z*z;
            h=(1.-3.*z2+2.*z2*z)*Hdata[i] +
                    z*(1.-2.*z+z2)*l*slope[i] +
                    z2*(3.-2.*z)*Hdata[i+1] +
                    z2*(z-1.)*l*slope[i+1];
            return p*h;
        }

    return 0;
}

double CMMaterialProp::GetB(const double hc) const
{
    if (BHpoints==0) return muo*mu_x*hc;

    double b,bo;

    b=0;
    do{
        bo = b;
        b  = bo + (hc-GetH(bo))/Re(GetdHdB(bo));
    }while (fabs(b-bo)>1.e-8);

    return b;
}

// GetEnergy for the magnetostatic case
double CMMaterialProp::GetEnergy(const double x) const
{
    double b,z,z2,l,nrg;
    double b0,b1,h0,h1,dh0,dh1;
    int i;

    b=fabs(x);
    nrg=0.;

    if(BHpoints==0)    return 0;

    for(i=0;i<BHpoints-1;i++){

        b0=Bdata[i];    h0=Re(Hdata[i]);
        b1=Bdata[i+1];    h1=Re(Hdata[i+1]);
        dh0=Re(slope[i]);
        dh1=Re(slope[i+1]);

        if((b>=b0) && (b<=b1)){
            l=b1-b0;
            z=(b-b0)/l;
            z2=z*z;

            nrg += (dh0*l*l*(6. + z*(-8. + 3.*z))*z2)/12. +
                    (h0*l*z*(2. + (-2. + z)*z2))/2. -
                    (h1*l*(-2. + z)*z2*z)/2. +
                    (dh1*l*l*(-4. + 3.*z)*z2*z)/12;

            return nrg;
        }
        else{
            // point isn't in this section, but add in the
            // energy required to pass through it.
            b0=Bdata[i];    h0=Re(Hdata[i]);
            b1=Bdata[i+1];    h1=Re(Hdata[i+1]);
            dh0=Re(slope[i]);
            dh1=Re(slope[i+1]);
            nrg += ((b0 - b1)*((b0 - b1)*(dh0 - dh1) -
                               6.*(h0 + h1)))/12.;
        }
    }

    // if we've gotten to this point, the point is off the scale,
    // so we have to extrapolate the rest of the way...
    h0=Re(Hdata[BHpoints-1]);
    dh0=Re(slope[BHpoints-1]);
    b0=Bdata[BHpoints-1];

    nrg += ((b - b0)*(b*dh0 - b0*dh0 + 2*h0))/2.;

    return nrg;
}

double CMMaterialProp::GetCoEnergy(const double b) const
{
    return (fabs(b)*GetH(b) - GetEnergy(b));
}

double CMMaterialProp::DoEnergy(const double b1, const double b2)
{
    // calls the raw routine to get point energy,
    // but deals with the load of special cases that
    // arise because I insisted on trying to deal with
    // laminations on a continuum basis.

    double nrg,biron,bair;
    nrg=biron=bair = 0;

    // easiest case: the material is linear!
    if (BHpoints==0)
    {
        double h1 = 0;
        double h2 = 0;

        if(LamType==0){        // laminated in-plane
            h1=b1/((1.+LamFill*(mu_x-1.))*muo);
            h2=b2/((1.+LamFill*(mu_y-1.))*muo);
        }

        if(LamType==1){        // laminated parallel to x;
            h1=b1/((1.+LamFill*(mu_x-1.))*muo);
            h2=b1*(LamFill/(mu_y*muo) + (1. - LamFill)/muo);
        }

        if(LamType==2){        // laminated parallel to x;
            h2=b1/((1.+LamFill*(mu_y-1.))*muo);
            h1=b1*(LamFill/(mu_x*muo) + (1. - LamFill)/muo);
        }

        if(LamType>2){
            h1=b1/muo;
            h2=b2/muo;
        }

        return ((h1*b1+h2*b2)/2.);
    }

    // Rats! The material is nonlinear.  Now, we have to do
    // a bit of work to get the energy.
    if(LamType==0) nrg = GetEnergy(sqrt(b1*b1+b2*b2));

    if(LamType==1){
        biron=sqrt((b1/LamFill)*(b1/LamFill) + b2*b2);
        bair=b2;
        nrg = LamFill*GetEnergy(biron)+(1-LamFill)*bair*bair/(2.*muo);
    }

    if(LamType==2){
        biron=sqrt((b2/LamFill)*(b2/LamFill) + b1*b1);
        bair=b1;
        nrg = LamFill*GetEnergy(biron)+(1-LamFill)*bair*bair/(2.*muo);
    }

    return nrg;
}

double CMMaterialProp::DoCoEnergy(const double b1, const double b2)
{
    double nrg,biron,bair;
    nrg=biron=bair = 0;

    // easiest case: the material is linear!
    // in this case, energy and coenergy are the same!
    if (BHpoints==0) return DoEnergy(b1,b2);

    if(LamType==0) nrg = GetCoEnergy(sqrt(b1*b1+b2*b2));

    if(LamType==1){
        biron=sqrt((b1/LamFill)*(b1/LamFill) + b2*b2);
        bair=b2;
        nrg = LamFill*GetCoEnergy(biron)+(1-LamFill)*bair*bair/(2.*muo);
    }

    if(LamType==2){
        biron=sqrt((b2/LamFill)*(b2/LamFill) + b1*b1);
        bair=b1;
        nrg = LamFill*GetCoEnergy(biron)+(1-LamFill)*bair*bair/(2.*muo);
    }

    return nrg;
}


double CMMaterialProp::DoEnergy(const CComplex b1, const CComplex b2)
{
    // This one is meant for the frequency!=0 case.
    // Fortunately, there's not so much effort in this case.
    CComplex mu1,mu2,h1,h2;

    GetMu(b1,b2,mu1,mu2);
    h1=b1/(mu1*muo);
    h2=b2/(mu2*muo);
    return (Re(h1*conj(b1)+h2*conj(b2))/4.);

}

double CMMaterialProp::DoCoEnergy(const CComplex b1, const CComplex b2)
{
    return DoEnergy(b1,b2);
}

bool CMMaterialProp::isAir() const
{
    // from femm/makemask.cpp
    if ((mu_x!=1) || (mu_y!=1)) return false;
    if (BHpoints!=0) return false;
    if (LamType!=0) return false;
    if (H_c!=0) return false;
    if((J.re!=0) || (J.im!=0)) return false;
    if (Cduct!=0) return false;
    if ((Theta_hn!=0) || (Theta_hx!=0) || (Theta_hy!=0)) return false;
    return true;
}

bool CMMaterialProp::isSameMaterialAs(const CMaterialProp *) const
{
    // stub
    return false;
}

void CMMaterialProp::toStream(ostream &) const
{
    assert(false && "CMMaterialProp::toStream() should never be called. Did you mean to call CMSolverMaterialProp::toStream()?");
}

void CMMaterialProp::GetMu(const CComplex b1, const CComplex b2,
                           CComplex &mu1, CComplex &mu2)
{
    // gets the permeability, given a flux density
    // version for frequency!=0

    CComplex biron;

    // easiest case: the material is linear!
    if (BHpoints==0)
    {
        mu1=mu_fdx;
        mu2=mu_fdy;
        return;
    }

    // Rats! The material is nonlinear.
    else{
        CComplex muiron;

        if(LamType==0){
            biron=sqrt(b1*conj(b1)+b2*conj(b2));
            if(abs(biron)<1.e-08) mu1=1./slope[0]; //catch degenerate case
            else mu1=biron/GetH(biron);
            mu2=mu1;
        }

        if(LamType==1){
            biron=sqrt((b1/LamFill)*(b1/LamFill) + b2*b2);
            if(abs(biron)<1.e-08) muiron=1./slope[0];
            else muiron=biron/GetH(biron);
            mu1=muiron*LamFill;
            mu2=1./(LamFill/muiron + (1. - LamFill)/muo);

        }

        if(LamType==2){
            biron=sqrt((b2/LamFill)*(b2/LamFill) + b1*b1);
            if(abs(biron)<1.e-08) muiron=1./slope[0];
            else muiron=biron/GetH(biron);
            mu2=muiron*LamFill;
            mu1=1./(LamFill/muiron + (1. - LamFill)/muo);
        }
    }

    // convert to relative permeability
    mu1/=muo;
    mu2/=muo;

    return;
}


void CMMaterialProp::GetMu(const double b1, const double b2, double &mu1, double &mu2)
{
    // gets the permeability, given a flux density
    //

    double biron;

    mu1=mu2=muo;            // default

    // catch the incremental permeability case
	if (MuMax>0) // for DC problems, MuMax is sleazily used as a flag to tell that it's incremental
	{
		IncrementalPermeability(std::sqrt(b1*b1 + b2*b2), mu1, mu2);
		return;
	}

    // easiest case: the material is linear!
    if (BHpoints==0)
    {
        if(LamType==0){        // laminated in-plane
            mu1=((1.+LamFill*(mu_x-1.))*muo);
            mu2=((1.+LamFill*(mu_y-1.))*muo);
        }

        if(LamType==1){        // laminated parallel to x;
            mu1=((1.+LamFill*(mu_x-1.))*muo);
            mu2=1./(LamFill/(mu_y*muo) + (1. - LamFill)/muo);
        }

        if(LamType==2){        // laminated parallel to x;
            mu2=((1.+LamFill*(mu_y-1.))*muo);
            mu1=1./(LamFill/(mu_x*muo) + (1. - LamFill)/muo);
        }
    }

    // Rats! The material is nonlinear.
    else{
        double muiron;

        if(LamType==0){
            biron=sqrt(b1*b1+b2*b2);
            if(biron<1.e-08) mu1=1./Re(slope[0]); //catch degenerate case
            else mu1=biron/GetH(biron);
            mu2=mu1;
        }

        if(LamType==1){
            biron=sqrt((b1/LamFill)*(b1/LamFill) + b2*b2);
            if(biron<1.e-08) muiron=1./Re(slope[0]);
            else muiron=biron/GetH(biron);
            mu1=muiron*LamFill;
            mu2=1./(LamFill/muiron + (1. - LamFill)/muo);
        }

        if(LamType==2){
            biron=sqrt((b2/LamFill)*(b2/LamFill) + b1*b1);
            if(biron<1.e-08) muiron=1./Re(slope[0]);
            else muiron=biron/GetH(biron);
            mu2=muiron*LamFill;
            mu1=1./(LamFill/muiron + (1. - LamFill)/muo);
        }

    }

    // convert to relative permeability
    mu1 /= muo;
    mu2 /= muo;

    return;
}

void CMMaterialProp::incrementalPermeability(const double B, const double w, CComplex &mu1, CComplex &mu2)
{
    // B == flux density in Tesla
    // w == frequency in rad/s

    // get incremental permeability of the DC material
    // (i.e. incremental permeability at the offset)
    double muinc=1./(muo*Re(GetdHdB(B)));
    double murel=1./(muo*Re(Get_v(B)));

    // if material is not laminated, just apply hysteresis lag...
    if ((Lam_d==0) || (LamFill==0))
    {
        mu1=muinc*exp(-I*Theta_hn*DEG*muinc/MuMax);
        mu2=murel*exp(-I*Theta_hn*DEG*murel/MuMax);
        return;
    }

    // crap.  Need to make an equivalent permeability that rolls in the effects of laminated
    // eddy currents, using the incremental permeability as the basis for creating the impedance.
    // this can get annoying because we need to back out the iron portion of the permeability
    // in the lamfill<1 case...


    if (Cduct!=0)
    {
        const CComplex deg45=1+I;

        // incremental permeability direction
        double mu = (muinc - (1.-LamFill))/LamFill;
        CComplex halflag=exp(-I*Theta_hn*DEG*mu/(2.*MuMax));
        double ds=sqrt(2./(0.4*PI*w*Cduct*mu));
        CComplex K=halflag*deg45*Lam_d*0.001/(2.*ds);
        mu1=(LamFill*mu*tanh(K)/K + (1.- LamFill));

        // normal permeability direction
        mu = (murel - (1.-LamFill))/LamFill;
        halflag=exp(-I*Theta_hn*DEG*mu/(2.*MuMax));
        ds=sqrt(2./(0.4*PI*w*Cduct*mu));
        K=halflag*deg45*Lam_d*0.001/(2.*ds);
        mu2=(LamFill*mu*tanh(K)/K + (1.- LamFill));
    }
    else{
        // incremental permeability direction
        double mu = (muinc - (1.-LamFill))/LamFill;
        mu1=(mu*exp(-I*Theta_hn*DEG*mu/MuMax)*LamFill + (1.-LamFill));

        // normal permeability direction
        mu = (murel - (1.-LamFill))/LamFill;
        mu2=(mu*exp(-I*Theta_hn*DEG*mu/MuMax)*LamFill + (1.-LamFill));
    }
}

CComplex CMMaterialProp::Get_v(double B)
{
    if (B==0) return slope[0];

    return (GetH(B)/B);
}

// get incremental permeability of a nonlinear material for use in
// incremental permeability formulation about DC offset
void CMMaterialProp::IncrementalPermeability(const double B, double &mu1, double &mu2)
{
	// B == flux density in Tesla

	double muinc, murel;

	// get incremental permeability of the DC material
	// (i.e. incremental permeability at the offset)
	muinc = 1. / (muo*Re(GetdHdB(B)));
	murel = 1. / (muo*Re(Get_v(B)));

	// if material is not laminated, just return
	if ((Lam_d == 0) || (LamFill == 0)){
		mu1 = muinc;
		mu2 = murel;
		return;
	}

	// incremental permeability direction
	mu1 = (muinc*LamFill + (1. - LamFill));

	// normal permeability direction
	mu2 = (murel*LamFill + (1. - LamFill));

	return;
}

ostream &operator<<(ostream &os, const CMMaterialProp &prop)
{
    prop.toStream(os);
    return os;
}


CMSolverMaterialProp::~CMSolverMaterialProp()
{
}

// Constructor
CMSolverMaterialProp::CMSolverMaterialProp()
    : CMMaterialProp()
{
    Theta_m=0.;            // magnetization direction, degrees;
}

// Copy Constructor
CMSolverMaterialProp::CMSolverMaterialProp( const CMSolverMaterialProp &other )
    : femm::CMMaterialProp( other )
{
    Theta_m = other.Theta_m;    // magnetization direction, degrees;
}

CComplex CMSolverMaterialProp::GetH(double B)
{
    double b,z,z2,l;
    CComplex h;
    int i;

    b=fabs(B);

    if(BHpoints==0)	return CComplex(b/(mu_x*muo));

    if(b>Bdata[BHpoints-1])
        return (Hdata[BHpoints-1] + slope[BHpoints-1]*(b-Bdata[BHpoints-1]));

    for(i=0; i<BHpoints-1; i++)
        if((b>=Bdata[i]) && (b<=Bdata[i+1]))
        {
            l=(Bdata[i+1]-Bdata[i]);
            z=(b-Bdata[i])/l;
            z2=z*z;
            h=(1.-3.*z2+2.*z2*z)*Hdata[i] +
              z*(1.-2.*z+z2)*l*slope[i] +
              z2*(3.-2.*z)*Hdata[i+1] +
              z2*(z-1.)*l*slope[i+1];
            return h;
        }

    return CComplex(0);
}


CComplex CMSolverMaterialProp::Get_dvB2(double B)
{
    if (B==0) return 0;

    return 0.5*(GetdHdB(B)/(B*B) - GetH(B)/(B*B*B));
}

void CMSolverMaterialProp::GetBHProps(double B, double &v, double &dv)
{
    // version to use in the magnetostatic case in
    // which we know that v and dv ought to be real-valued.
    CComplex vc,dvc;

    GetBHProps(B,vc,dvc);
    v =Re(vc);
    dv=Re(dvc);
}

void CMSolverMaterialProp::GetBHProps(double B, CComplex &v, CComplex &dv)
{
    double b,z,z2,l;
    CComplex h,dh;
    int i;

    b=fabs(B);

    if(BHpoints==0)
    {
        v=mu_x;
        dv=0;
        return;
    }

    if(b==0)
    {
        v=slope[0];
        dv=0;
        return;
    }

    if(b>Bdata[BHpoints-1])
    {
        h=(Hdata[BHpoints-1] + slope[BHpoints-1]*(b-Bdata[BHpoints-1]));
        dh=slope[BHpoints-1];
        v=h/b;
        dv=0.5*(dh/(b*b) - h/(b*b*b));
        return;
    }

    for(i=0; i<BHpoints-1; i++)
        if((b>=Bdata[i]) && (b<=Bdata[i+1]))
        {
            l=(Bdata[i+1]-Bdata[i]);
            z=(b-Bdata[i])/l;
            z2=z*z;
            h=(1.-3.*z2+2.*z2*z)*Hdata[i] +
              z*(1.-2.*z+z2)*l*slope[i] +
              z2*(3.-2.*z)*Hdata[i+1] +
              z2*(z-1.)*l*slope[i+1];
            dh=6.*z*(z-1.)*Hdata[i]/l +
               (1.-4.*z+3.*z*z)*slope[i] +
               6.*z*(1.-z)*Hdata[i+1]/l +
               z*(3.*z-2.)*slope[i+1];
            v=h/b;
            dv=0.5*(dh/(b*b) - h/(b*b*b));
            return;
        }
}

// this can't be immediately merged with femm::CMaterialProp,
// because this uses the slightly different GetH implementation
// from fsolver
CComplex CMSolverMaterialProp::LaminatedBH(double w, int i)
{
    debug << "CMSolverMaterialProp::LaminatedBH("<<w<<", "<<i<<")\n";
    int k,n,iter=0;
    CComplex *m0,*m1,*b,*x;
    double L,o,d,ds,B,lastres,res;
    res=0;
    double Relax=1;
    CComplex mu,vo,vi,c,H;
    CComplex Md,Mo;
    int Converged=false;

    // Base the required element spacing on the skin depth
    // at the surface of the material
    mu=Bdata[i]/Hdata[i];
    o=Cduct*1.e6;
    d=(Lam_d*0.001)/2.;
    ds=sqrt(2/(w*o*abs(mu)));
    n= ElementsPerSkinDepth * ((int) ceil(d/ds));
    L=d/((double) n);

    x =(CComplex *)calloc(n+1,sizeof(CComplex));
    b =(CComplex *)calloc(n+1,sizeof(CComplex));
    m0=(CComplex *)calloc(n+1,sizeof(CComplex));
    m1=(CComplex *)calloc(n+1,sizeof(CComplex));

    do
    {
        // make sure that the old stuff is wiped out;
        for(k=0; k<=n; k++)
        {
            m0[k]=0;
            m1[k]=0;
            b[k] =0;
        }

        // build matrix
        for(k=0; k<n; k++)
        {
            if(iter!=0)
            {
                B=abs(x[k+1]-x[k])/L;
                vi=GetdHdB(B);
                vo=GetH(B)/B; // <- this is the GetH from fsolver
            }
            else
            {
                vo=1./mu;
                vi=1./mu;
            }

            Md = ( (vi+vo)/(2.*L) + I*w*o*L/4.);
            Mo = (-(vi+vo)/(2.*L) + I*w*o*L/4.);
            m0[k]   += Md;
            m1[k]   += Mo;
            m0[k+1] += Md;

            Md = ( (vi-vo)/(2.*L));
            Mo = (-(vi-vo)/(2.*L));
            b[k]    += (Md*x[k] + Mo*x[k+1]);
            b[k+1]  += (Mo*x[k] + Md*x[k+1]);
        }

        // set boundary conditions
        m1[0] = 0;
        b[0]  = 0;
        b[n] += Hdata[i];

        // solve tridiagonal equation
        // solution ends up in b;
        for(k = 0; k < n; k++)
        {
            c = m1[k]/m0[k];
            m0[k+1] -= m1[k]*c;
            b[k+1]  -= b[k]*c;
        }
        b[n] = b[n]/m0[n];
        for(k = n-1; k >= 0; k--)
            b[k] = (b[k] - m1[k]*b[k + 1])/m0[k];

        iter++;

        lastres=res;
        res=abs(b[n]-x[n])/d;

        if (res<1.e-8) Converged=true;

        // Do the same relaxation scheme as is implemented
        // in the solver to make sure that this effective
        // lamination permeability calculation converges
        if(iter>5)
        {
            if ((res>lastres) && (Relax>0.1)) Relax/=2.;
            else Relax+= 0.1 * (1. - Relax);
        }

        for(k=0; k<=n; k++) x[k]=Relax*b[k]+(1.0-Relax)*x[k];

    }
    while(Converged!=true);

    mu = x[n]/(Hdata[i]*d);

    free(x );
    free(b );
    free(m0);
    free(m1);

    return mu;
}

CMSolverMaterialProp CMSolverMaterialProp::fromStream(std::istream &input, std::ostream &err, PropertyParseMode mode)
{
    using namespace femm;
    CMSolverMaterialProp prop;

    if( mode == PropertyParseMode::NoBeginBlock
            || expectToken(input, "<beginblock>", err) )
    {
        string token;
        while (input.good() && token != "<endblock>")
        {
            nextToken(input,&token);

            if( token == "<blockname>" )
            {
                expectChar(input, '=', err);
                parseString(input, &prop.BlockName, err);
                continue;
            }

            if( token == "<mu_x>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.mu_x, err);
                continue;
            }

            if( token == "<mu_y>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.mu_y, err);
                continue;
            }

            if( token == "<h_c>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.H_c, err);
                continue;
            }

            if( token == "<h_cangle>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Theta_m, err);
                continue;
            }

            if( token == "<j_re>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.J.re, err);
                continue;
            }

            if( token == "<j_im>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.J.im, err);
                continue;
            }

            if( token == "<sigma>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Cduct, err);
                continue;
            }

            if( token == "<phi_h>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Theta_hn, err);
                continue;
            }


            if( token == "<phi_hx>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Theta_hx, err);
                continue;
            }

            if( token == "<phi_hy>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Theta_hy, err);
                continue;
            }

            if( token == "<d_lam>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Lam_d, err);
                continue;
            }

            if( token == "<lamfill>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.LamFill, err);
                continue;
            }

            if( token == "<wired>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.WireD, err);
                continue;
            }

            if( token == "<lamtype>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.LamType, err);
                continue;
            }

            if( token == "<nstrands>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.NStrands, err);
                continue;
            }

            if( token == "<bhpoints>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.BHpoints, err);
                if (prop.BHpoints < 0 || prop.BHpoints == 1)
                {
                    err << "CMSolverMaterialProp: B-H curve must contain either zero points or at least two points\n";
                    input.setstate(std::ios::failbit);
                    break;
                }
                if (prop.BHpoints > 0)
                {
                    prop.Hdata.reserve(prop.BHpoints);
                    prop.Bdata.reserve(prop.BHpoints);
                    for(int i=0; i<prop.BHpoints; i++)
                    {
                        double b;
                        CComplex h;
                        input >> b >> h.re;
                        if (!input)
                        {
                            err << "CMSolverMaterialProp: could not read B-H point " << i << "\n";
                            break;
                        }
                        prop.Bdata.push_back(b);
                        prop.Hdata.push_back(h);
                    }
                }
                if (input && (prop.Bdata.size() != static_cast<size_t>(prop.BHpoints)
                              || prop.Hdata.size() != static_cast<size_t>(prop.BHpoints)))
                {
                    err << "CMSolverMaterialProp: B-H curve point count does not match parsed data\n";
                    input.setstate(std::ios::failbit);
                }
                for (int i = 0; input && i < prop.BHpoints; ++i)
                {
                    if (!std::isfinite(prop.Bdata[i]) || !std::isfinite(prop.Hdata[i].re)
                            || !std::isfinite(prop.Hdata[i].im))
                    {
                        err << "CMSolverMaterialProp: non-finite B or H component at B-H point " << i << "\n";
                        input.setstate(std::ios::failbit);
                    }
                    else if (i > 0 && !(prop.Bdata[i] > prop.Bdata[i - 1]))
                    {
                        err << "CMSolverMaterialProp: B values must be strictly increasing; invalid interval at B-H point " << i << "\n";
                        input.setstate(std::ios::failbit);
                    }
                }
                continue;
            }

            if (token != "<endblock>")
                err << "CMSolverMaterialProp: unexpected token: "<<token << "\n";
        }

        // Note: this code is shared between all problem stages - processing needs to be done in the stage-specific code
        // (e.g. in the solver)
    }

    return prop;
}

void CMSolverMaterialProp::toStream(ostream &out) const
{
    out << "  <BeginBlock>\n";
    out << "    <BlockName> = \"" << BlockName << "\"\n";
    out << "    <Mu_x> = " << mu_x << "\n";
    out << "    <Mu_y> = " << mu_y << "\n";
    out << "    <H_c> = " << H_c << "\n";
    out << "    <H_cAngle> = " << Theta_m << "\n";
    out << "    <J_re> = " << J.re << "\n";
    out << "    <J_im> = " << J.im << "\n";
    out << "    <Sigma> = " << Cduct << "\n";
    out << "    <d_lam> = " << Lam_d << "\n";
    out << "    <Phi_h> = " << Theta_hn << "\n";
    out << "    <Phi_hx> = " << Theta_hx << "\n";
    out << "    <Phi_hy> = " << Theta_hy << "\n";
    out << "    <LamType> = " << LamType << "\n";
    out << "    <LamFill> = " << LamFill << "\n";
    out << "    <NStrands> = " << NStrands << "\n";
    out << "    <WireD> = " << WireD << "\n";
    out << "    <BHPoints> = " << BHpoints << "\n";
    for(int i=0; i<BHpoints; i++)
    {
        out << "      " << Bdata.at(i) << "\t" << Hdata.at(i).re << "\n";
    }
    out << "  <EndBlock>\n";
}

CHMaterialProp::CHMaterialProp()
    : CMaterialProp()
    , Kx(1)
    , Ky(1)
    , Kt(0)
    , qv(0)
    , npts(0)
{
}

CHMaterialProp::~CHMaterialProp()
{
}

CHMaterialProp::CHMaterialProp( const CHMaterialProp & other)
{
    Kx = other.Kx;
    Ky = other.Ky;
    Kt = other.Kt;  // volumetric heat capacity
    qv = other.qv;  // volume heat generation

    // properties for nonlinear conductivity
    npts = other.npts;			// number of points in the nonlinear conductivity curve

    for (int i = 0; i < npts; i++)
    {
        // copy the thermal conductivity data points.
        Kn[i] = other.Kn[i];
    }
}

CComplex CHMaterialProp::GetK(double t) const
{
    int i,j;

    // Kx returned as real part;
    // Ky returned as imag part

    if (npts==0) return (Kx+I*Ky);
    if (npts==1) return (Im(Kn[0])*(1+I));
    if (t<=Re(Kn[0])) return (Im(Kn[0])*(1+I));
    if (t>=Re(Kn[npts-1])) return (Im(Kn[npts-1])*(1+I));

    for(i=0,j=1;j<npts;i++,j++)
    {
        if((t>=Re(Kn[i])) && (t<=Re(Kn[j])))
        {
            return (1+I)*(Im(Kn[i])+Im(Kn[j]-Kn[i])*Re(t-Kn[i])/Re(Kn[j]-Kn[i]));
        }
    }

    return (Kx+I*Ky);
}

CHMaterialProp CHMaterialProp::fromStream(std::istream &input, std::ostream &err, PropertyParseMode mode)
{
    CHMaterialProp prop;

    if( mode == PropertyParseMode::NoBeginBlock
            || expectToken(input, "<beginblock>", err) )
    {
        string token;
        while (input.good() && token != "<endblock>")
        {
            nextToken(input,&token);

            if( token == "<kx>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Kx, err);
                continue;
            }

            if( token == "<ky>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Ky, err);
                continue;
            }

            if( token == "<kt>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.Kt, err);
                continue;
            }

            if( token == "<qv>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.qv, err);
                continue;
            }

            if( token == "<blockname>" )
            {
                expectChar(input, '=', err);
                parseString(input, &(prop.BlockName), err);
                continue;
            }

            if( token == "<tkpoints>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.npts, err);
                if (prop.npts > 0)
                {
                    // FIXME: make Kn variable size
                    if (prop.npts>128)
                    {
                        err << "MaterialProp exceeds maximum size! File a bug report...";
                        prop.npts = 128;
                    }
                    for(int i=0; i<prop.npts; i++)
                    {
                        input >> prop.Kn[i].re >> prop.Kn[i].im;
                    }
                }
                continue;
            }
            if( token != "<endblock>")
                err << "CHMaterialProp: unexpected token: "<<token << "\n";
        }
    }

    return prop;
}

bool CHMaterialProp::isAir() const
{
    return false;
}

bool CHMaterialProp::isSameMaterialAs(const CMaterialProp *other) const
{
    // Are the same material trivially if they are the same block type
    if (other == this) return true;

    const CHMaterialProp *m2 = dynamic_cast<const CHMaterialProp*>(other);

    // If the materials are linear and have the same Kx and Ky, we
    // can say that they are the same material;
    if ((Kx==m2->Kx) &&
        (Ky==m2->Ky) &&
        (npts==0) &&
        (m2->npts==0)) return true;

    // If the materials are nonlinear and have all of the same T-k points,
    // they are the same material;
    if (npts>0){
        if (npts==m2->npts)
        {
            for(int k=0;k<npts;k++)
            {
                if ((Kn[k].re!=m2->Kn[k].re) ||
                    (Kn[k].im!=m2->Kn[k].im))
                    return false;
            }
            return true;
        }
    }

    return false;
}

void CHMaterialProp::toStream(std::ostream &out) const
{
    out << "  <BeginBlock>\n";
    out << "    <Kx> = " << Kx << "\n";
    out << "    <Ky> = " << Ky << "\n";
    out << "    <Kt> = " << Kt << "\n";
    out << "    <qv> = " << qv << "\n";
    if (!BlockName.empty())
        out << "    <BlockName> = \"" << BlockName << "\"\n";
    out << "    <TKPoints> = " << npts << "\n";
    if (npts > 0)
    {
        for(int i=0; i<npts; i++)
        {
            out << "        " << Kn[i].re << "\t" << Kn[i].im << "\n";
        }
    }
    out << "  <EndBlock>\n";
}

CSMaterialProp::CSMaterialProp()
    : CMaterialProp()
    , ex(1.)
    , ey(1.)
    , qv(0)
{
}

CSMaterialProp::~CSMaterialProp()
{
}

CSMaterialProp CSMaterialProp::fromStream(istream &input, ostream &err, PropertyParseMode mode)
{
    using namespace femm;
    CSMaterialProp prop;

    if( mode == PropertyParseMode::NoBeginBlock
            || expectToken(input, "<beginblock>", err) )
    {
        string token;
        while (input.good() && token != "<endblock>")
        {
            nextToken(input,&token);

            if( token == "<blockname>" )
            {
                expectChar(input, '=', err);
                parseString(input, &prop.BlockName, err);
                continue;
            }

            if( token == "<ex>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.ex, err);
                continue;
            }

            if( token == "<ey>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.ey, err);
                continue;
            }

            if( token == "<qv>" )
            {
                expectChar(input, '=', err);
                parseValue(input, prop.qv, err);
                continue;
            }

            if (token != "<endblock>")
                err << "CSMaterialProp: unexpected token: "<<token << "\n";
        }
    }

    return prop;
}

bool CSMaterialProp::isAir() const
{
    // from femm/bv_makemask.cpp
    if ((ex!=1) || (ey!=1)) return false;
    if (qv!=0) return false;
    return true;
}

bool CSMaterialProp::isSameMaterialAs(const CMaterialProp *other) const
{
    // Are the same material trivially if they are the same block type
    if (other == this) return true;

    const CSMaterialProp *m2 = dynamic_cast<const CSMaterialProp*>(other);
    return (m2 != nullptr && (ex==m2->ex) &&  (ey==m2->ey));
}

void CSMaterialProp::toStream(ostream &out) const
{
    out << "  <BeginBlock>\n";
    out << "    <BlockName> = \"" << BlockName << "\"\n";
    out << "    <ex> = " << ex << "\n";
    out << "    <ey> = " << ey << "\n";
    out << "    <qv> = " << qv << "\n";
    out << "  <EndBlock>\n";
}
