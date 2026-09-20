#ifndef FEMM_JSON_H
#define FEMM_JSON_H

#include <string>
#include <utility>
#include <vector>

namespace femm {
namespace json {

/**
 * A minimal JSON value used by the tiled-model file format.
 *
 * This is deliberately self-contained (no third-party dependency) so it works
 * in every build configuration, including the MATLAB/Octave MEX build. It
 * supports the JSON subset the format needs: objects, arrays, strings,
 * numbers, booleans, and null.
 */
class Value
{
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Value() = default;

    static Value boolean(bool value);
    static Value number(double value);
    static Value string(std::string value);
    static Value array();
    static Value object();

    Type type() const { return m_type; }
    bool isNull() const { return m_type == Type::Null; }
    bool isBool() const { return m_type == Type::Bool; }
    bool isNumber() const { return m_type == Type::Number; }
    bool isString() const { return m_type == Type::String; }
    bool isArray() const { return m_type == Type::Array; }
    bool isObject() const { return m_type == Type::Object; }

    bool boolValue() const { return m_bool; }
    double numberValue() const { return m_number; }
    const std::string &stringValue() const { return m_string; }

    std::vector<Value> &arrayItems() { return m_array; }
    const std::vector<Value> &arrayItems() const { return m_array; }
    void push(Value value);

    bool has(const std::string &key) const { return find(key) != nullptr; }
    const Value *find(const std::string &key) const;
    void set(const std::string &key, Value value);
    const std::vector<std::pair<std::string, Value>> &members() const { return m_object; }

    /** Serialize as indented JSON. */
    std::string dump(int indent = 2) const;

    /** Parse \p text. On failure returns false and sets \p error. */
    static bool parse(const std::string &text, Value &out, std::string &error);

private:
    void dumpTo(std::string &out, int indent, int depth) const;

    Type m_type = Type::Null;
    bool m_bool = false;
    double m_number = 0.0;
    std::string m_string;
    std::vector<Value> m_array;
    std::vector<std::pair<std::string, Value>> m_object;
};

} // namespace json
} // namespace femm

#endif
