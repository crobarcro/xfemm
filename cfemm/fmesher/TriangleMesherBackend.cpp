#include "TriangleMesherBackend.h"
#include "TangleMesherBackend.h"
#include "SolverMeshFileWriter.h"
#include "fmesher.h"
#include "femmconstants.h"
#include <fstream>
#include <iomanip>

namespace fmesher {
namespace {
std::string root(const std::string &p) { auto n=p.find_last_of('.'); return n==std::string::npos?p:p.substr(0,n); }
void diagnostic(femm::mesh::MeshResult &r, const std::string &s, int code) {
 r.status=femm::mesh::MeshStatus::BackendFailure; r.diagnostics.push_back({femm::mesh::MeshDiagnosticSeverity::Error,s,"Triangle",code});
}
}
femm::mesh::MeshResult TriangleMesherBackend::mesh(femm::FemmProblem &problem,bool periodic,const femm::mesh::MeshingOptions &options){
 femm::mesh::MeshResult result; FMesher legacy(std::shared_ptr<femm::FemmProblem>(&problem,[](femm::FemmProblem*){}));legacy.Verbose=options.verbose;legacy.writePolyFiles=writePolyFiles;if(WarnMessage)legacy.WarnMessage=WarnMessage;if(TriMessage)legacy.TriMessage=TriMessage;
 std::string path=compatibilityPath.empty()?"xfemm-backend.fem":compatibilityPath;int status=periodic?legacy.doPeriodicTriangleWorkflow(path,result.mesh):legacy.doNonPeriodicTriangleWorkflow(path,result.mesh);if(status){diagnostic(result,"Triangle meshing workflow failed",status);return result;}result.status=femm::mesh::MeshStatus::Success;return result;
}

bool SolverMeshFileWriter::write(const femm::mesh::SolverMesh &m,const femm::FemmProblem &p,const std::string &path,int(*warn)(const char*,...)){
 auto base=root(path);const double inv=1.0/femm::LengthConvMeters[p.LengthUnits];std::ofstream f(base+".node");if(!f)return false;f<<m.nodes.size()<<"\t2\t0\t1\n"<<std::setprecision(17);for(size_t i=0;i<m.nodes.size();i++)f<<i<<'\t'<<m.nodes[i].x*inv<<'\t'<<m.nodes[i].y*inv<<'\t'<<m.nodes[i].boundaryMarker<<'\n';
 f.close();f.open(base+".edge");if(!f)return false;f<<m.edges.size()<<"\t1\n";for(size_t i=0;i<m.edges.size();i++)f<<i<<'\t'<<m.edges[i].first<<'\t'<<m.edges[i].second<<'\t'<<m.edges[i].boundaryMarker<<'\n';
 f.close();f.open(base+".ele");if(!f)return false;f<<m.elements.size()<<"\t3\t1\n";for(size_t i=0;i<m.elements.size();i++)f<<i<<'\t'<<m.elements[i].nodes[0]<<'\t'<<m.elements[i].nodes[1]<<'\t'<<m.elements[i].nodes[2]<<'\t'<<m.elements[i].regionAttribute<<'\n';
 f.close();f.open(base+".pbc");if(!f)return false;f<<m.periodicConstraints.size()<<'\n';for(size_t i=0;i<m.periodicConstraints.size();i++)f<<i<<' '<<m.periodicConstraints[i].first<<' '<<m.periodicConstraints[i].second<<' '<<(m.periodicConstraints[i].periodicity==femm::mesh::SolverMesh::Periodicity::Antiperiodic)<<'\n';f<<m.airGaps.size()<<'\n';for(const auto&a:m.airGaps){f<<std::quoted(a.boundaryName)<<'\n'<<(a.periodicity==femm::mesh::SolverMesh::Periodicity::Antiperiodic)<<' '<<a.innerAngleDegrees<<' '<<a.outerAngleDegrees<<' '<<a.innerRadius*inv<<' '<<a.outerRadius*inv<<' '<<a.totalArcLengthDegrees<<' '<<a.centerX*inv<<' '<<a.centerY*inv<<' '<<a.totalArcElements<<' '<<a.innerShift<<' '<<a.outerShift<<'\n';for(const auto&q:a.quadraturePoints){for(int k=0;k<4;k++)f<<q.nodes[k]<<' '<<q.weights[k]<<(k==3?'\n':' ');}}if(!f&&warn)warn("Couldn't write mesh files");return bool(f);
}
}

namespace fmesher {
int FMesher::DoNonPeriodicBCTriangulation(std::string path) {
    femm::mesh::MeshingOptions options; options.verbose=Verbose;
    femm::mesh::MeshResult result;
    if (backend == Backend::Tangle) {
        TangleMesherBackend selected;
        result=selected.mesh(*problem,false,options);
    } else {
        TriangleMesherBackend selected; selected.WarnMessage=WarnMessage; selected.TriMessage=TriMessage;
        selected.writePolyFiles=writePolyFiles; selected.compatibilityPath=path;
        result=selected.mesh(*problem,false,options);
    }
    return result.succeeded() && SolverMeshFileWriter::write(result.mesh,*problem,path,WarnMessage) ? 0 : -1;
}
int FMesher::DoPeriodicBCTriangulation(std::string path) {
    femm::mesh::MeshingOptions options; options.verbose=Verbose;
    femm::mesh::MeshResult result;
    if (backend == Backend::Tangle) {
        TangleMesherBackend selected;
        result=selected.mesh(*problem,true,options);
    } else {
        TriangleMesherBackend selected; selected.WarnMessage=WarnMessage; selected.TriMessage=TriMessage;
        selected.writePolyFiles=writePolyFiles; selected.compatibilityPath=path;
        result=selected.mesh(*problem,true,options);
    }
    return result.succeeded() && SolverMeshFileWriter::write(result.mesh,*problem,path,WarnMessage) ? 0 : -1;
}
}
