#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"

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

    void reject() { m_trial.reset(); }

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

    std::unique_ptr<femm::AnalysisSession> m_session;
    std::shared_ptr<femm::FSolverAnalysisBackend> m_solver;
    std::unique_ptr<femm::TrialSolution> m_trial;
    std::shared_ptr<const femm::AcceptedState> m_accepted;
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
        .def("reject", &PythonSession::reject);
}
