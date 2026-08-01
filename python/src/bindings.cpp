#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"
#include "fpproc.h"

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

class PythonSession {
public:
    explicit PythonSession(const std::string &filename)
    {
        std::unique_ptr<femm::FemmProblem> owned(
            new femm::FemmProblem(femm::FileType::MagneticsFile));
        std::shared_ptr<femm::FemmProblem> view(owned.get(), [](femm::FemmProblem *) {});
        std::ostringstream errors;
        femm::MagneticsReader reader(view, errors);
        if (reader.parse(filename) != femm::F_FILE_OK)
            throw std::runtime_error("could not load magnetic problem: " + errors.str());
        m_solver = std::make_shared<femm::FSolverAnalysisBackend>();
        m_session.reset(new femm::AnalysisSession(
            femm::ModelDefinition(std::move(owned)), m_solver));
    }

    void setBackend(const std::string &name)
    {
        if (name == "triangle")
            m_session->setMesher(std::make_shared<fmesher::TriangleMesherBackend>());
        else if (name == "tangle")
            m_session->setMesher(std::make_shared<fmesher::TangleMesherBackend>());
        else
            throw std::invalid_argument("backend must be 'triangle' or 'tangle'");
    }

    std::size_t mesh() { return m_session->ensureMesh()->elements.size(); }

    void setCircuit(const std::string &name, const std::string &kind,
                    std::complex<double> value)
    {
        const auto id = m_session->model().circuit(name);
        const CComplex native(value.real(), value.imag());
        if (kind == "current") m_session->setCircuitCurrent(id, native);
        else if (kind == "voltage") m_session->setCircuitVoltage(id, native);
        else if (kind == "open") m_session->setCircuitOpen(id);
        else if (kind == "coupled") m_session->setCircuitCoupled(id);
        else throw std::invalid_argument(
            "constraint must be current, voltage, open, or coupled");
    }

    void setAgePosition(const std::string &name, double inner, double outer)
    {
        m_session->setAirGapAngle(m_session->model().airGap(name), inner, outer);
    }

    void setFrequency(double value) { m_session->setFrequency(value); }
    void setTime(double value) { m_session->setTime(value); }

    py::dict solve()
    {
        m_trial.reset(new femm::TrialSolution(m_session->solve()));
        m_post.reset(new FPProc);
        if (!m_post->OpenDocument(m_session->model().problem(), m_solver->solvedSolver(),
                                  m_solver->solvedSystem()))
            throw std::runtime_error("could not initialize in-memory post-processor");
        const auto mesh = m_session->mesh();
        py::dict out;
        out["success"] = true;
        out["id"] = m_trial->id;
        out["time"] = m_trial->time;
        out["nodeCount"] = mesh ? mesh->nodes.size() : 0;
        out["elementCount"] = mesh ? mesh->elements.size() : 0;
        out["meshGenerationCount"] = m_session->meshGenerationCount();
        out["solveCount"] = m_solver->solveCount();
        out["operatorAssemblyCount"] = m_solver->operatorAssemblyCount();
        out["rightHandSideAssemblyCount"] = m_solver->rightHandSideAssemblyCount();
        return out;
    }

    py::dict result() const
    {
        const auto &trial = requireTrial();
        py::dict out;
        out["id"] = trial.id;
        out["time"] = trial.time;
        py::array_t<std::complex<double>> current(trial.circuits.size());
        py::array_t<std::complex<double>> flux(trial.circuits.size());
        py::array_t<std::complex<double>> voltage(trial.circuits.size());
        auto cu = current.mutable_unchecked<1>();
        auto fl = flux.mutable_unchecked<1>();
        auto vo = voltage.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(trial.circuits.size()); ++i) {
            const auto &c = trial.circuits[static_cast<std::size_t>(i)];
            cu(i) = {c.current.re, c.current.im};
            fl(i) = {c.fluxLinkage.re, c.fluxLinkage.im};
            vo(i) = c.terminalVoltage
                  ? std::complex<double>(c.terminalVoltage->re, c.terminalVoltage->im)
                  : std::complex<double>(NAN, NAN);
        }
        out["current"] = current;
        out["fluxLinkage"] = flux;
        out["terminalVoltage"] = voltage;
        const std::vector<double> empty;
        const auto &real = trial.real;
        out["A"] = vectorArray(real ? real->nodal.magneticVectorPotential : empty);
        out["x"] = vectorArray(real ? real->nodal.x : empty);
        out["y"] = vectorArray(real ? real->nodal.y : empty);
        return out;
    }

    void accept()
    {
        m_accepted = m_session->acceptSolution(requireTrial());
        m_session->setInitialState(m_accepted);
    }

    void reject() { m_trial.reset(); m_post.reset(); }

    py::array_t<std::complex<double>> pointValues(const py::array_t<double> &points)
    {
        requirePost();
        auto in = points.unchecked<2>();
        py::array_t<std::complex<double>> out(
            std::vector<py::ssize_t>{in.shape(0), 14});
        auto values = out.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < in.shape(0); ++i) {
            CMPointVals v;
            if (!m_post->GetPointValues(in(i, 0), in(i, 1), v)) {
                for (int j = 0; j < 14; ++j) values(i, j) = {NAN, NAN};
                continue;
            }
            values(i,0)={v.A.re,v.A.im}; values(i,1)={v.B1.re,v.B1.im};
            values(i,2)={v.B2.re,v.B2.im}; values(i,3)=v.c; values(i,4)=v.E;
            values(i,5)={v.H1.re,v.H1.im}; values(i,6)={v.H2.re,v.H2.im};
            values(i,7)={v.Je.re,v.Je.im}; values(i,8)={v.Js.re,v.Js.im};
            values(i,9)={v.mu1.re,v.mu1.im}; values(i,10)={v.mu2.re,v.mu2.im};
            values(i,11)=v.Pe; values(i,12)=v.Ph; values(i,13)=v.ff;
        }
        return out;
    }

    void smooth(bool enabled) { requirePost(); m_post->Smooth = enabled; }
    void clearContour() { requirePost(); m_post->contour.clear(); }
    void addContour(const py::array_t<double> &points) {
        requirePost(); auto in=points.unchecked<2>();
        for(py::ssize_t i=0;i<in.shape(0);++i) m_post->contour.emplace_back(in(i,0),in(i,1));
    }
    py::array_t<std::complex<double>> lineIntegral(int type) {
        requirePost(); if(type<0||type>5) throw std::invalid_argument("line integral type must be 0..5");
        CComplex z[4]{}; m_post->LineIntegral(type,z); int count=(type==3&&m_post->Frequency!=0)?4:(type==2?2:(type==0||type==1||type==5?2:1));
        py::array_t<std::complex<double>> out(count); auto a=out.mutable_unchecked<1>();
        if(type==2){a(0)=z[0].re;a(1)=z[0].im;} else for(int i=0;i<count;++i)a(i)={z[i].re,z[i].im}; return out;
    }
    void selectBlock(double x,double y){requirePost();int e=m_post->InTriangle(x,y);if(e>=0)m_post->blocklist[m_post->meshelem[e].lbl].IsSelected=true;}
    void groupSelectBlock(int group){requirePost();for(auto &label:m_post->blocklist)if(label.InGroup==group)label.IsSelected=true;}
    void selectAllBlocks(){requirePost();for(auto &label:m_post->blocklist)label.IsSelected=true;}
    void clearBlock(){requirePost();for(auto &label:m_post->blocklist)label.IsSelected=false;}
    std::complex<double> blockIntegral(int type){requirePost();if(type>=18&&type<=23)m_post->MakeMask();auto z=m_post->BlockIntegral(type);return {z.re,z.im};}
    py::array_t<double> problemInfo(){requirePost();py::array_t<double> o(4);auto a=o.mutable_unchecked<1>();a(0)=m_post->problemType;a(1)=m_post->Frequency;a(2)=m_post->Depth;a(3)=m_post->LengthConv[m_post->LengthUnits];return o;}
    py::array_t<std::complex<double>> circuitProps(const std::string &name){requirePost();int k=-1;for(std::size_t i=0;i<m_post->circproplist.size();++i)if(m_post->circproplist[i].CircName==name)k=i;if(k<0)throw std::invalid_argument("unknown circuit");py::array_t<std::complex<double>>o(3);auto a=o.mutable_unchecked<1>();auto c=m_post->circproplist[k].Amps,v=m_post->GetVoltageDrop(k),f=m_post->GetFluxLinkage(k);a(0)={c.re,c.im};a(1)={v.re,v.im};a(2)={f.re,f.im};return o;}
    int numNodes(){requirePost();return m_post->numNodes();}
    int numElements(){requirePost();return m_post->numElements();}

private:
    static py::array_t<double> vectorArray(const std::vector<double> &source)
    {
        py::array_t<double> result(source.size());
        auto values = result.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(source.size()); ++i)
            values(i) = source[static_cast<std::size_t>(i)];
        return result;
    }

    const femm::TrialSolution &requireTrial() const
    {
        if (!m_trial) throw std::logic_error("there is no trial solution; call solve first");
        return *m_trial;
    }
    void requirePost() const { if(!m_post) throw std::logic_error("there is no post-processor; call solve first"); }

    std::unique_ptr<femm::AnalysisSession> m_session;
    std::shared_ptr<femm::FSolverAnalysisBackend> m_solver;
    std::unique_ptr<femm::TrialSolution> m_trial;
    std::shared_ptr<const femm::AcceptedState> m_accepted;
    std::unique_ptr<FPProc> m_post;
};

} // namespace

PYBIND11_MODULE(_xfemm, module)
{
    module.doc() = "Native xfemm analysis session";
    py::class_<PythonSession>(module, "NativeSession")
        .def(py::init<const std::string &>())
        .def("set_backend", &PythonSession::setBackend)
        .def("mesh", &PythonSession::mesh)
        .def("set_circuit", &PythonSession::setCircuit)
        .def("set_age_position", &PythonSession::setAgePosition)
        .def("set_frequency", &PythonSession::setFrequency)
        .def("set_time", &PythonSession::setTime)
        .def("solve", &PythonSession::solve)
        .def("result", &PythonSession::result)
        .def("accept", &PythonSession::accept)
        .def("reject", &PythonSession::reject)
        .def("getpointvalues", &PythonSession::pointValues)
        .def("smooth", &PythonSession::smooth)
        .def("clearcontour", &PythonSession::clearContour)
        .def("addcontour", &PythonSession::addContour)
        .def("lineintegral", &PythonSession::lineIntegral)
        .def("selectblock", &PythonSession::selectBlock)
        .def("groupselectblock", &PythonSession::groupSelectBlock)
        .def("selectallblocks", &PythonSession::selectAllBlocks)
        .def("clearblock", &PythonSession::clearBlock)
        .def("blockintegral", &PythonSession::blockIntegral)
        .def("getprobleminfo", &PythonSession::problemInfo)
        .def("getcircuitprops", &PythonSession::circuitProps)
        .def("nummeshnodes", &PythonSession::numNodes)
        .def("numelements", &PythonSession::numElements);
}
