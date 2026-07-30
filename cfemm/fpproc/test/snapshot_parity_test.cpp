#include "MagneticSolutionSnapshot.h"
#include "fpproc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#define REQUIRE(condition) do { if (!(condition)) { std::fprintf(stderr, "requirement failed at line %d: %s\n", __LINE__, #condition); std::abort(); } } while (false)

namespace {
bool close(CComplex a, CComplex b)
{
    const auto scalarClose = [](double x, double y) {
        return (std::isnan(x) && std::isnan(y)) ||
               std::abs(x-y) <= 1e-11*(1+std::abs(x));
    };
    return scalarClose(a.re,b.re) && scalarClose(a.im,b.im);
}
}

int main(int argc, char **argv)
{
    REQUIRE(argc == 2);
    FPProc persisted;
    REQUIRE(persisted.OpenDocument(argv[1]));
    const auto snapshot = persisted.solutionSnapshot();
    FPProc memory(snapshot);

    REQUIRE(memory.numNodes() == persisted.numNodes());
    REQUIRE(memory.numElements() == persisted.numElements());
    CMPointVals a, b;
    const auto &element = persisted.meshelem.front();
    const double x = (persisted.meshnode[element.p[0]].x + persisted.meshnode[element.p[1]].x + persisted.meshnode[element.p[2]].x) / 3.;
    const double y = (persisted.meshnode[element.p[0]].y + persisted.meshnode[element.p[1]].y + persisted.meshnode[element.p[2]].y) / 3.;
    REQUIRE(persisted.GetPointValues(x, y, a) == memory.GetPointValues(x, y, b));
    REQUIRE(close(a.A, b.A)); REQUIRE(close(a.B1, b.B1)); REQUIRE(close(a.B2, b.B2));

    for (std::size_t i=0;i<persisted.blocklist.size();++i) {
        persisted.blocklist[i].IsSelected=true; memory.blocklist[i].IsSelected=true;
    }
    // Includes energy/current/volume, Lorentz and weighted-stress force/torque families.
    for (int integral=0;integral<=30;++integral) {
        const CComplex fromFile=persisted.BlockIntegral(integral);
        const CComplex fromMemory=memory.BlockIntegral(integral);
        if (!close(fromFile,fromMemory)) {
            std::fprintf(stderr,"block integral %d differs: (%g,%g) != (%g,%g)\n",
                         integral,fromFile.re,fromFile.im,fromMemory.re,fromMemory.im);
            std::abort();
        }
    }
    // The non-AGE fixture supplies a simple interior contour for all line integrals.
    if (persisted.agelist.empty()) {
        const auto &n0=persisted.meshnode[element.p[0]];
        const auto &n1=persisted.meshnode[element.p[1]];
        persisted.contour={CComplex(.9*x+.1*n0.x,.9*y+.1*n0.y),
                           CComplex(.9*x+.1*n1.x,.9*y+.1*n1.y)};
        memory.contour=persisted.contour;
        for (int integral=0;integral<=5;++integral) {
            CComplex fromFile[4]{},fromMemory[4]{};
            persisted.LineIntegral(integral,fromFile);
            memory.LineIntegral(integral,fromMemory);
            for (int i=0;i<4;++i)
                REQUIRE(close(fromFile[i],fromMemory[i]));
        }
    }
    for (std::size_t i=0;i<persisted.circproplist.size();++i) {
        REQUIRE(close(persisted.GetVoltageDrop(static_cast<int>(i)), memory.GetVoltageDrop(static_cast<int>(i))));
        REQUIRE(close(persisted.GetFluxLinkage(static_cast<int>(i)), memory.GetFluxLinkage(static_cast<int>(i))));
    }
    for (const auto &gap : persisted.agelist) {
        int na=0, nb=0;
        REQUIRE(persisted.numGapHarmonics(gap.BdryName,na)==memory.numGapHarmonics(gap.BdryName,nb));
        REQUIRE(na==nb);
        for (int n=0;n<na;++n) {
            CComplex ac1,as1,br1,bs1,bt1,bts1, ac2,as2,br2,bs2,bt2,bts2;
            REQUIRE(persisted.getGapHarmonics(gap.BdryName,n,ac1,as1,br1,bs1,bt1,bts1)==
                   memory.getGapHarmonics(gap.BdryName,n,ac2,as2,br2,bs2,bt2,bts2));
            REQUIRE(close(ac1,ac2)&&close(as1,as2)&&close(br1,br2)&&close(bs1,bs2)&&close(bt1,bt2)&&close(bts1,bts2));
        }
        double tq1=0,tq2=0;
        REQUIRE(persisted.gapDCTorqueIntegral(gap.BdryName,tq1)==memory.gapDCTorqueIntegral(gap.BdryName,tq2));
        REQUIRE(std::abs(tq1-tq2)<=1e-11*(1+std::abs(tq1)));
        CComplex fx1,fy1,fx2,fy2;
        REQUIRE(persisted.gapDCForceIntegral(gap.BdryName,fx1,fy1)==memory.gapDCForceIntegral(gap.BdryName,fx2,fy2));
        REQUIRE(close(fx1,fx2)&&close(fy1,fy2));
    }

    const std::string copy = std::string(argv[1])+".snapshot-copy";
    const auto diskSnapshot = MagneticSolutionSnapshot::fromAnsFile(argv[1]);
    REQUIRE(diskSnapshot.writeAnsFile(copy));
    FPProc roundTrip(MagneticSolutionSnapshot::fromAnsFile(copy));
    REQUIRE(roundTrip.numNodes()==persisted.numNodes());
    std::remove(copy.c_str());
}
