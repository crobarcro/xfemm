#ifndef MAGNETICSOLUTIONSNAPSHOT_H
#define MAGNETICSOLUTIONSNAPSHOT_H

#include <memory>
#include <string>

class FPProc;

/**
 * Immutable, in-memory result of a magnetic analysis.
 *
 * The snapshot deliberately has no mutating accessors.  FPProc is a view over
 * a snapshot; .ans files are only persistence adapters used to create or save
 * one.
 */
class MagneticSolutionSnapshot
{
public:
    /** Read the legacy FEMM .ans representation into an immutable snapshot. */
    static MagneticSolutionSnapshot fromAnsFile(const std::string &path);

    /** Persist the original .ans representation, when this snapshot came from one. */
    bool writeAnsFile(const std::string &path) const;

private:
    friend class FPProc;
    MagneticSolutionSnapshot(std::shared_ptr<const FPProc> state,
                             std::shared_ptr<const std::string> ansRepresentation = {});

    std::shared_ptr<const FPProc> m_state;
    std::shared_ptr<const std::string> m_ansRepresentation;
};

#endif
