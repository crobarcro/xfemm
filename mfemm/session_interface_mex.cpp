#include "mex.h"

#define CLASS_HANDLE_SIGNATURE 0xA51E5510
#include "postproc/class_handle.hpp"

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string stringValue(const mxArray *value, const char *name)
{
    if (!mxIsChar(value))
        throw std::invalid_argument(std::string(name) + " must be a character vector");
    char *text = mxArrayToString(value);
    if (!text) throw std::runtime_error("could not convert character vector");
    std::string result(text);
    mxFree(text);
    return result;
}

double scalarValue(const mxArray *value, const char *name)
{
    if (!mxIsNumeric(value) || mxGetNumberOfElements(value) != 1)
        throw std::invalid_argument(std::string(name) + " must be a numeric scalar");
    return mxGetScalar(value);
}

CComplex inputComplex(const mxArray *value, const char *name)
{
    const double real = scalarValue(value, name);
    const double imaginary = mxIsComplex(value) ? *mxGetPi(value) : 0.0;
    return CComplex(real, imaginary);
}

class SessionGateway {
public:
    void backend(const std::string &name)
    {
        if (name == "triangle")
            m_session->setMesher(std::make_shared<fmesher::TriangleMesherBackend>());
        else if (name == "tangle")
            m_session->setMesher(std::make_shared<fmesher::TangleMesherBackend>());
        else
            throw std::invalid_argument("backend must be 'triangle' or 'tangle'");
        m_backendName = name;
    }

    static std::unique_ptr<SessionGateway> create(const std::string &filename)
    {
        std::unique_ptr<SessionGateway> result(new SessionGateway);
        std::unique_ptr<femm::FemmProblem> owned(new femm::FemmProblem(femm::FileType::MagneticsFile));
        // The reader's shared_ptr is a non-owning view; ownership moves directly
        // from this gateway factory into ModelDefinition after parsing.
        std::shared_ptr<femm::FemmProblem> problem(owned.get(), [](femm::FemmProblem *) {});
        std::ostringstream errors;
        femm::MagneticsReader reader(problem, errors);
        if (reader.parse(filename) != femm::F_FILE_OK)
            throw std::runtime_error("could not load magnetic problem: " + errors.str());
        result->m_solver = std::make_shared<femm::FSolverAnalysisBackend>();
        result->m_session.reset(new femm::AnalysisSession(
            femm::ModelDefinition(std::move(owned)), result->m_solver));
        result->m_filename = filename;
        return result;
    }

    void mesh() { m_session->ensureMesh(); }
    femm::AnalysisSession &session() { return *m_session; }
    const std::string &backendName() const { return m_backendName; }
    femm::TrialSolution &trial() {
        if (!m_trial) throw std::logic_error("there is no trial solution; call solve first");
        return *m_trial;
    }
    void solve() { m_trial.reset(new femm::TrialSolution(m_session->solve())); }
    void writeSolution(const std::string &path) { trial(); m_solver->writeSolution(path); }
    void accept() { m_accepted = m_session->acceptSolution(trial()); m_session->setInitialState(m_accepted); }
    void reject() { m_trial.reset(); }
    void restore(std::shared_ptr<const femm::AcceptedState> state) {
        m_accepted = std::move(state); m_session->setInitialState(m_accepted);
    }
    std::shared_ptr<const femm::AcceptedState> accepted() const { return m_accepted; }

private:
    SessionGateway() = default;
    std::unique_ptr<femm::AnalysisSession> m_session;
    std::shared_ptr<femm::FSolverAnalysisBackend> m_solver;
    std::unique_ptr<femm::TrialSolution> m_trial;
    std::shared_ptr<const femm::AcceptedState> m_accepted;
    std::string m_filename, m_backendName = "triangle";
};

mxArray *trialStruct(const femm::TrialSolution &trial)
{
    const char *fields[] = {"id", "time", "current", "fluxLinkage",
                            "terminalVoltage", "A", "x", "y"};
    mxArray *out = mxCreateStructMatrix(1, 1, 8, fields);
    mxSetField(out, 0, "id", mxCreateDoubleScalar(static_cast<double>(trial.id)));
    mxSetField(out, 0, "time", mxCreateDoubleScalar(trial.time));
    const mwSize count = trial.circuits.size();
    mxArray *current = mxCreateDoubleMatrix(count, 1, mxCOMPLEX);
    mxArray *flux = mxCreateDoubleMatrix(count, 1, mxCOMPLEX);
    mxArray *voltage = mxCreateDoubleMatrix(count, 1, mxCOMPLEX);
    for (mwSize i = 0; i < count; ++i) {
        mxGetPr(current)[i] = trial.circuits[i].current.re; mxGetPi(current)[i] = trial.circuits[i].current.im;
        mxGetPr(flux)[i] = trial.circuits[i].fluxLinkage.re; mxGetPi(flux)[i] = trial.circuits[i].fluxLinkage.im;
        if (trial.circuits[i].terminalVoltage) {
            mxGetPr(voltage)[i] = trial.circuits[i].terminalVoltage->re;
            mxGetPi(voltage)[i] = trial.circuits[i].terminalVoltage->im;
        } else mxGetPr(voltage)[i] = mxGetNaN();
    }
    mxSetField(out, 0, "current", current); mxSetField(out, 0, "fluxLinkage", flux);
    mxSetField(out, 0, "terminalVoltage", voltage);
    const std::vector<double> empty;
    const auto &a = trial.real ? trial.real->nodal.magneticVectorPotential : empty;
    mxArray *nodes = mxCreateDoubleMatrix(a.size(), 1, mxREAL);
    for (mwSize i = 0; i < a.size(); ++i) mxGetPr(nodes)[i] = a[i];
    mxSetField(out, 0, "A", nodes);
    const auto &xValues = trial.real ? trial.real->nodal.x : empty;
    const auto &yValues = trial.real ? trial.real->nodal.y : empty;
    mxArray *x = mxCreateDoubleMatrix(xValues.size(), 1, mxREAL);
    mxArray *y = mxCreateDoubleMatrix(yValues.size(), 1, mxREAL);
    for (mwSize i = 0; i < xValues.size(); ++i) {
        mxGetPr(x)[i] = xValues[i];
        mxGetPr(y)[i] = yValues[i];
    }
    mxSetField(out, 0, "x", x); mxSetField(out, 0, "y", y);
    return out;
}

} // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
try {
    if (nrhs < 1) throw std::invalid_argument("first input must be a command");
    const std::string command = stringValue(prhs[0], "command");
    if (command == "new") {
        if (nrhs != 2 || nlhs != 1) throw std::invalid_argument("new requires one filename and one output");
        plhs[0] = convertPtr2Mat(SessionGateway::create(stringValue(prhs[1], "filename")).release());
        return;
    }
    if (nrhs < 2) throw std::invalid_argument("second input must be a session handle");
    if (command == "delete") { destroyObject<SessionGateway>(prhs[1]); return; }
    SessionGateway *gateway = convertMat2Ptr<SessionGateway>(prhs[1]);
    auto &session = gateway->session();
    if (command == "backend") gateway->backend(stringValue(prhs[2], "backend"));
    else if (command == "mesh") { gateway->mesh(); if (nlhs) plhs[0] = mxCreateDoubleScalar(session.mesh()->elements.size()); }
    else if (command == "frequency") session.setFrequency(scalarValue(prhs[2], "frequency"));
    else if (command == "time") session.setTime(scalarValue(prhs[2], "time"));
    else if (command == "circuit") {
        const auto id = session.model().circuit(stringValue(prhs[2], "circuit name"));
        const std::string kind = stringValue(prhs[3], "constraint");
        if (kind == "current") session.setCircuitCurrent(id, inputComplex(prhs[4], "value"));
        else if (kind == "voltage") session.setCircuitVoltage(id, inputComplex(prhs[4], "value"));
        else if (kind == "open") session.setCircuitOpen(id);
        else if (kind == "coupled") session.setCircuitCoupled(id);
        else throw std::invalid_argument("constraint must be current, voltage, open, or coupled");
    } else if (command == "age") session.setAirGapAngle(session.model().airGap(stringValue(prhs[2], "AGE name")), scalarValue(prhs[3], "inner angle"), scalarValue(prhs[4], "outer angle"));
    else if (command == "solve") { gateway->solve(); if (nlhs) plhs[0] = trialStruct(gateway->trial()); }
    else if (command == "result") plhs[0] = trialStruct(gateway->trial());
    else if (command == "export") gateway->writeSolution(stringValue(prhs[2], "solution path"));
    else if (command == "accept") gateway->accept();
    else if (command == "reject") gateway->reject();
    else if (command == "state") {
        const auto state = gateway->accepted();
        if (!state) { plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL); return; }
        const char *fields[] = {"modelRevision", "parameterRevision", "time", "backendState"};
        plhs[0] = mxCreateStructMatrix(1, 1, 4, fields);
        mxSetField(plhs[0],0,"modelRevision",mxCreateDoubleScalar(state->modelRevision));
        mxSetField(plhs[0],0,"parameterRevision",mxCreateDoubleScalar(state->parameterRevision));
        mxSetField(plhs[0],0,"time",mxCreateDoubleScalar(state->time));
        mxSetField(plhs[0],0,"backendState",mxCreateString(state->backendState.c_str()));
    } else if (command == "restore") {
        const mxArray *value = prhs[2];
        if (!mxIsStruct(value)) throw std::invalid_argument("saved state must be a struct");
        std::shared_ptr<femm::AcceptedState> state(new femm::AcceptedState);
        const mxArray *field = mxGetField(value, 0, "modelRevision");
        if (!field) throw std::invalid_argument("saved state has no modelRevision");
        state->modelRevision = static_cast<std::uint64_t>(scalarValue(field, "modelRevision"));
        field = mxGetField(value, 0, "parameterRevision");
        if (!field) throw std::invalid_argument("saved state has no parameterRevision");
        state->parameterRevision = static_cast<std::uint64_t>(scalarValue(field, "parameterRevision"));
        field = mxGetField(value, 0, "time");
        if (!field) throw std::invalid_argument("saved state has no time");
        state->time = scalarValue(field, "time");
        field = mxGetField(value, 0, "backendState");
        if (!field) throw std::invalid_argument("saved state has no backendState");
        state->backendState = stringValue(field, "backendState");
        gateway->restore(state);
    } else throw std::invalid_argument("unknown session command: " + command);
} catch (const std::exception &error) {
    mexErrMsgIdAndTxt("MFEMM:session:error", "%s", error.what());
}
