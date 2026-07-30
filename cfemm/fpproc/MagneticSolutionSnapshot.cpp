#include "MagneticSolutionSnapshot.h"

#include "fpproc.h"

#include <fstream>
#include <iterator>
#include <stdexcept>

MagneticSolutionSnapshot::MagneticSolutionSnapshot(
        std::shared_ptr<const FPProc> state,
        std::shared_ptr<const std::string> ansRepresentation)
    : m_state(std::move(state)), m_ansRepresentation(std::move(ansRepresentation))
{
    if (!m_state)
        throw std::invalid_argument("a magnetic solution snapshot requires solution state");
}

MagneticSolutionSnapshot MagneticSolutionSnapshot::fromAnsFile(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
        throw std::runtime_error("could not read magnetic solution file: " + path);
    auto representation = std::make_shared<const std::string>(
        std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());

    auto processor = std::make_shared<FPProc>();
    if (!processor->OpenDocument(path))
        throw std::runtime_error("invalid magnetic solution file: " + path);
    return MagneticSolutionSnapshot(std::move(processor), std::move(representation));
}

bool MagneticSolutionSnapshot::writeAnsFile(const std::string &path) const
{
    if (!m_ansRepresentation)
        return false;
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output)
        return false;
    output.write(m_ansRepresentation->data(),
                 static_cast<std::streamsize>(m_ansRepresentation->size()));
    return output.good();
}
