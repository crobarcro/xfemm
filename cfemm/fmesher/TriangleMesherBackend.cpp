#include "TriangleMesherBackend.h"
#include "SolverMeshFileWriter.h"
#include "fmesher.h"
#include "femmconstants.h"
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace fmesher {
namespace {
std::string root(const std::string &p) { auto n=p.find_last_of('.'); return n==std::string::npos?p:p.substr(0,n); }
void diagnostic(femm::mesh::MeshResult &r, const std::string &s, int code) {
 r.status=femm::mesh::MeshStatus::BackendFailure; r.diagnostics.push_back({femm::mesh::MeshDiagnosticSeverity::Error,s,"Triangle",code});
}
bool readLegacy(const std::string &path, const femm::FemmProblem &problem, femm::mesh::SolverMesh &out) {
 const auto base=root(path); std::ifstream f(base+".node"); size_t count; int dim,attrs,markers;
 if(!(f>>count>>dim>>attrs>>markers)) return false;
 const double scale=femm::LengthConvMeters[problem.LengthUnits];
 out.nodes.resize(count); for(size_t q=0;q<count;q++){size_t id;int marker; f>>id>>out.nodes[id].x>>out.nodes[id].y>>marker; out.nodes[id].x*=scale;out.nodes[id].y*=scale;out.nodes[id].boundaryMarker=marker;}
 f.close(); f.open(base+".ele"); int corners; if(!(f>>count>>corners>>attrs))return false; out.elements.resize(count);
 for(size_t q=0;q<count;q++){size_t id;int region;f>>id>>out.elements[id].nodes[0]>>out.elements[id].nodes[1]>>out.elements[id].nodes[2]>>region;out.elements[id].regionAttribute=region;}
 f.close();f.open(base+".edge");if(!(f>>count>>markers))return false;out.edges.resize(count);
 for(size_t q=0;q<count;q++){size_t id;int marker;f>>id>>out.edges[id].first>>out.edges[id].second>>marker;out.edges[id].boundaryMarker=marker;}
 f.close();f.open(base+".pbc");if(!(f>>count))return false;out.periodicConstraints.resize(count);
 for(size_t q=0;q<count;q++){size_t id;int anti;f>>id>>out.periodicConstraints[id].first>>out.periodicConstraints[id].second>>anti;out.periodicConstraints[id].periodicity=anti?femm::mesh::SolverMesh::Periodicity::Antiperiodic:femm::mesh::SolverMesh::Periodicity::Periodic;}
 size_t gaps=0;if(!(f>>gaps))return false;out.airGaps.reserve(gaps);
 for(size_t g=0;g<gaps;g++){femm::mesh::SolverMesh::AirGap a; f>>std::quoted(a.boundaryName);int anti;size_t n;f>>anti>>a.innerAngleDegrees>>a.outerAngleDegrees>>a.innerRadius>>a.outerRadius>>a.totalArcLengthDegrees>>a.centerX>>a.centerY>>n>>a.innerShift>>a.outerShift;a.periodicity=anti?femm::mesh::SolverMesh::Periodicity::Antiperiodic:femm::mesh::SolverMesh::Periodicity::Periodic;a.totalArcElements=n;a.innerRadius*=scale;a.outerRadius*=scale;a.centerX*=scale;a.centerY*=scale;a.quadraturePoints.resize(n+1);for(auto &qp:a.quadraturePoints)for(int k=0;k<4;k++)f>>qp.nodes[k]>>qp.weights[k];out.airGaps.push_back(std::move(a));}
 return bool(f);
}
}
femm::mesh::MeshResult TriangleMesherBackend::mesh(femm::FemmProblem &problem,bool periodic,const femm::mesh::MeshingOptions &options){
 femm::mesh::MeshResult result; FMesher legacy(std::shared_ptr<femm::FemmProblem>(&problem,[](femm::FemmProblem*){}));legacy.Verbose=options.verbose;legacy.writePolyFiles=writePolyFiles;if(WarnMessage)legacy.WarnMessage=WarnMessage;if(TriMessage)legacy.TriMessage=TriMessage;
 std::string path=compatibilityPath.empty()?"xfemm-backend.fem":compatibilityPath;int status=periodic?legacy.doPeriodicTriangleWorkflow(path):legacy.doNonPeriodicTriangleWorkflow(path);if(status){diagnostic(result,"Triangle meshing workflow failed",status);return result;}if(!readLegacy(path,problem,result.mesh)){diagnostic(result,"Could not convert Triangle output to SolverMesh",-1);return result;}result.status=femm::mesh::MeshStatus::Success;return result;
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
    TriangleMesherBackend backend; backend.WarnMessage=WarnMessage; backend.TriMessage=TriMessage;
    backend.writePolyFiles=writePolyFiles; backend.compatibilityPath=path;
    femm::mesh::MeshingOptions options; options.verbose=Verbose;
    auto result=backend.mesh(*problem,false,options);
    return result.succeeded() && SolverMeshFileWriter::write(result.mesh,*problem,path,WarnMessage) ? 0 : -1;
}
int FMesher::DoPeriodicBCTriangulation(std::string path) {
    TriangleMesherBackend backend; backend.WarnMessage=WarnMessage; backend.TriMessage=TriMessage;
    backend.writePolyFiles=writePolyFiles; backend.compatibilityPath=path;
    femm::mesh::MeshingOptions options; options.verbose=Verbose;
    auto result=backend.mesh(*problem,true,options);
    return result.succeeded() && SolverMeshFileWriter::write(result.mesh,*problem,path,WarnMessage) ? 0 : -1;
}
}
