#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"
#include "femmenums.h"

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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
        m_lengthScale = femm::LengthConvMeters[owned->LengthUnits];
        m_isPlanar = owned->problemType == femm::PLANAR;
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

    std::size_t mesh()
    {
        return m_session->ensureMesh()->elements.size();
    }

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
        buildMeshOrderedPotential();
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

    void reject()
    {
        m_trial.reset();
        m_meshPotential.clear();
    }

    std::size_t numNodes() const { return requireMesh()->nodes.size(); }
    std::size_t numElements() const { return requireMesh()->elements.size(); }

    py::array_t<double> vertices(const py::array_t<std::int64_t> &requested) const
    {
        const auto mesh = requireMesh();
        const auto ids = elementIds(requested, mesh->elements.size());
        py::array_t<double> out(std::vector<py::ssize_t>{static_cast<py::ssize_t>(ids.size()), 6});
        auto values = out.mutable_unchecked<2>();
        for (py::ssize_t row = 0; row < static_cast<py::ssize_t>(ids.size()); ++row) {
            const auto &e = mesh->elements[ids[static_cast<std::size_t>(row)]];
            for (std::size_t j = 0; j < 3; ++j) {
                const auto &node = mesh->nodes[e.nodes[j]];
                values(row, 2 * j) = node.x / m_lengthScale;
                values(row, 2 * j + 1) = node.y / m_lengthScale;
            }
        }
        return out;
    }

    py::array_t<double> elements(const py::array_t<std::int64_t> &requested) const
    {
        const auto mesh = requireMesh();
        const auto ids = elementIds(requested, mesh->elements.size());
        py::array_t<double> out(std::vector<py::ssize_t>{static_cast<py::ssize_t>(ids.size()), 7});
        auto values = out.mutable_unchecked<2>();
        for (py::ssize_t row = 0; row < static_cast<py::ssize_t>(ids.size()); ++row) {
            const auto &e = mesh->elements[ids[static_cast<std::size_t>(row)]];
            for (std::size_t j = 0; j < 3; ++j) values(row, j) = e.nodes[j] + 1;
            const auto geometry = elementGeometry(*mesh, e);
            values(row, 3) = geometry.cx / m_lengthScale;
            values(row, 4) = geometry.cy / m_lengthScale;
            values(row, 5) = geometry.area / (m_lengthScale * m_lengthScale);
            const auto label = e.regionAttribute > 0
                ? static_cast<std::size_t>(e.regionAttribute - 1) : 0;
            const auto &labels = m_session->model().problem().labellist;
            values(row, 6) = label < labels.size() ? labels[label]->InGroup : 0;
        }
        return out;
    }

    py::array_t<double> centroids(const py::array_t<std::int64_t> &requested) const
    {
        const auto mesh = requireMesh();
        const auto ids = elementIds(requested, mesh->elements.size());
        py::array_t<double> out(std::vector<py::ssize_t>{static_cast<py::ssize_t>(ids.size()), 2});
        auto values = out.mutable_unchecked<2>();
        for (py::ssize_t row = 0; row < static_cast<py::ssize_t>(ids.size()); ++row) {
            const auto g = elementGeometry(*mesh, mesh->elements[ids[row]]);
            values(row, 0) = g.cx / m_lengthScale;
            values(row, 1) = g.cy / m_lengthScale;
        }
        return out;
    }

    py::array_t<double> areas(const py::array_t<std::int64_t> &requested) const
    {
        const auto mesh = requireMesh();
        const auto ids = elementIds(requested, mesh->elements.size());
        py::array_t<double> out(ids.size());
        auto values = out.mutable_unchecked<1>();
        for (py::ssize_t row = 0; row < static_cast<py::ssize_t>(ids.size()); ++row)
            values(row) = elementGeometry(*mesh, mesh->elements[ids[row]]).area
                        / (m_lengthScale * m_lengthScale);
        return out;
    }

    py::array_t<double> getA(const py::array_t<double> &points) const
    {
        requireTrial();
        const auto input = points.unchecked<2>();
        py::array_t<double> out(input.shape(0));
        auto values = out.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < input.shape(0); ++i) {
            const auto sample = pointValue(input(i, 0) * m_lengthScale,
                                           input(i, 1) * m_lengthScale);
            values(i) = sample.a;
        }
        return out;
    }

    py::array_t<double> getB(const py::array_t<double> &points) const
    {
        requireTrial();
        if (!m_isPlanar)
            throw std::runtime_error(
                "getb for axisymmetric solutions requires the full native post-processor");
        const auto input = points.unchecked<2>();
        py::array_t<double> out(std::vector<py::ssize_t>{input.shape(0), 2});
        auto values = out.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < input.shape(0); ++i) {
            const auto sample = pointValue(input(i, 0) * m_lengthScale,
                                           input(i, 1) * m_lengthScale);
            values(i, 0) = sample.bx;
            values(i, 1) = sample.by;
        }
        return out;
    }

private:
    struct Geometry { double cx, cy, area, determinant; };
    struct PointResult { double a, bx, by; };

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

    std::shared_ptr<const femm::mesh::SolverMesh> requireMesh() const
    {
        const auto mesh = m_session->mesh();
        if (!mesh) throw std::logic_error("there is no mesh; call mesh or solve first");
        return mesh;
    }

    static std::vector<std::size_t> elementIds(
        const py::array_t<std::int64_t> &requested, std::size_t count)
    {
        const auto input = requested.unchecked<1>();
        std::vector<std::size_t> result;
        if (input.shape(0) == 0) {
            result.resize(count);
            for (std::size_t i = 0; i < count; ++i) result[i] = i;
            return result;
        }
        result.reserve(input.shape(0));
        for (py::ssize_t i = 0; i < input.shape(0); ++i) {
            if (input(i) < 1 || static_cast<std::size_t>(input(i)) > count)
                throw py::index_error("element index is out of range");
            result.push_back(static_cast<std::size_t>(input(i) - 1));
        }
        return result;
    }

    static Geometry elementGeometry(const femm::mesh::SolverMesh &mesh,
                                    const femm::mesh::SolverMesh::Element &e)
    {
        const auto &a = mesh.nodes[e.nodes[0]];
        const auto &b = mesh.nodes[e.nodes[1]];
        const auto &c = mesh.nodes[e.nodes[2]];
        const double det = (b.y - c.y) * (a.x - c.x)
                         + (c.x - b.x) * (a.y - c.y);
        return {(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0,
                std::abs(det) / 2.0, det};
    }

    void buildMeshOrderedPotential()
    {
        const auto &trial = requireTrial();
        if (!trial.real) throw std::runtime_error("solution has no real nodal field");
        const auto mesh = requireMesh();
        m_meshPotential.assign(mesh->nodes.size(), NAN);
        const auto &nodal = trial.real->nodal;
        if (nodal.x.size() != mesh->nodes.size()
                || nodal.y.size() != mesh->nodes.size()
                || nodal.magneticVectorPotential.size() != mesh->nodes.size())
            throw std::runtime_error("solved node arrays do not match the session mesh");
        using CoordinateKey = std::pair<long long, long long>;
        constexpr double keyScale = 1e10;
        std::map<CoordinateKey, std::vector<std::size_t>> byCoordinate;
        for (std::size_t j = 0; j < nodal.x.size(); ++j)
            byCoordinate[{std::llround(nodal.x[j] * keyScale),
                          std::llround(nodal.y[j] * keyScale)}].push_back(j);
        for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
            const CoordinateKey key{std::llround(mesh->nodes[i].x * keyScale),
                                    std::llround(mesh->nodes[i].y * keyScale)};
            auto found = byCoordinate.find(key);
            if (found == byCoordinate.end() || found->second.empty())
                throw std::runtime_error("could not map solved nodes to the session mesh");
            const std::size_t j = found->second.back();
            found->second.pop_back();
            m_meshPotential[i] = nodal.magneticVectorPotential[j];
        }
    }

    PointResult pointValue(double x, double y) const
    {
        const auto mesh = requireMesh();
        constexpr double tolerance = 1e-12;
        for (const auto &e : mesh->elements) {
            const auto &a = mesh->nodes[e.nodes[0]];
            const auto &b = mesh->nodes[e.nodes[1]];
            const auto &c = mesh->nodes[e.nodes[2]];
            const double det = (b.y - c.y) * (a.x - c.x)
                             + (c.x - b.x) * (a.y - c.y);
            if (std::abs(det) < 1e-30) continue;
            const double l0 = ((b.y - c.y) * (x - c.x)
                             + (c.x - b.x) * (y - c.y)) / det;
            const double l1 = ((c.y - a.y) * (x - c.x)
                             + (a.x - c.x) * (y - c.y)) / det;
            const double l2 = 1.0 - l0 - l1;
            if (l0 < -tolerance || l1 < -tolerance || l2 < -tolerance) continue;
            const double av[3] = {m_meshPotential[e.nodes[0]],
                                  m_meshPotential[e.nodes[1]],
                                  m_meshPotential[e.nodes[2]]};
            const double potential = l0 * av[0] + l1 * av[1] + l2 * av[2];
            const double dadx = (av[0] * (b.y - c.y) + av[1] * (c.y - a.y)
                               + av[2] * (a.y - b.y)) / det;
            const double dady = (av[0] * (c.x - b.x) + av[1] * (a.x - c.x)
                               + av[2] * (b.x - a.x)) / det;
            return {potential, dady, -dadx};
        }
        return {NAN, NAN, NAN};
    }

    std::unique_ptr<femm::AnalysisSession> m_session;
    std::shared_ptr<femm::FSolverAnalysisBackend> m_solver;
    std::unique_ptr<femm::TrialSolution> m_trial;
    std::shared_ptr<const femm::AcceptedState> m_accepted;
    std::vector<double> m_meshPotential;
    double m_lengthScale = 1.0;
    bool m_isPlanar = true;
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
        .def("nummeshnodes", &PythonSession::numNodes)
        .def("numelements", &PythonSession::numElements)
        .def("getvertices", &PythonSession::vertices)
        .def("getelements", &PythonSession::elements)
        .def("getcentroids", &PythonSession::centroids)
        .def("getareas", &PythonSession::areas)
        .def("geta", &PythonSession::getA)
        .def("getb", &PythonSession::getB);
}
