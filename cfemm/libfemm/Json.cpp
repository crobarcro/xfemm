#include "Json.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace femm {
namespace json {

Value Value::boolean(bool value)
{
    Value v;
    v.m_type = Type::Bool;
    v.m_bool = value;
    return v;
}

Value Value::number(double value)
{
    Value v;
    v.m_type = Type::Number;
    v.m_number = value;
    return v;
}

Value Value::string(std::string value)
{
    Value v;
    v.m_type = Type::String;
    v.m_string = std::move(value);
    return v;
}

Value Value::array()
{
    Value v;
    v.m_type = Type::Array;
    return v;
}

Value Value::object()
{
    Value v;
    v.m_type = Type::Object;
    return v;
}

void Value::push(Value value)
{
    m_type = Type::Array;
    m_array.push_back(std::move(value));
}

const Value *Value::find(const std::string &key) const
{
    for (const auto &member : m_object)
        if (member.first == key)
            return &member.second;
    return nullptr;
}

void Value::set(const std::string &key, Value value)
{
    m_type = Type::Object;
    for (auto &member : m_object)
        if (member.first == key) {
            member.second = std::move(value);
            return;
        }
    m_object.emplace_back(key, std::move(value));
}

namespace {

void escapeString(const std::string &value, std::string &out)
{
    out.push_back('"');
    for (char c : value) {
        switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b"; break;
        case '\f': out += "\\f"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
            if (static_cast<unsigned char>(c) < 0x20) {
                char buffer[8];
                std::snprintf(buffer, sizeof buffer, "\\u%04x", c);
                out += buffer;
            } else {
                out.push_back(c);
            }
        }
    }
    out.push_back('"');
}

std::string numberToString(double value)
{
    if (!std::isfinite(value))
        return "0";
    char buffer[40];
    std::snprintf(buffer, sizeof buffer, "%.17g", value);
    return buffer;
}

struct Parser {
    const std::string &text;
    std::size_t pos = 0;
    std::string error;

    explicit Parser(const std::string &input) : text(input) {}

    bool fail(const std::string &message)
    {
        if (error.empty()) {
            std::ostringstream stream;
            stream << message << " at offset " << pos;
            error = stream.str();
        }
        return false;
    }

    void skipWhitespace()
    {
        while (pos < text.size()) {
            const char c = text[pos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
                ++pos;
            else
                break;
        }
    }

    bool parseValue(Value &out)
    {
        skipWhitespace();
        if (pos >= text.size())
            return fail("unexpected end of input");
        const char c = text[pos];
        if (c == '{')
            return parseObject(out);
        if (c == '[')
            return parseArray(out);
        if (c == '"') {
            std::string value;
            if (!parseString(value))
                return false;
            out = Value::string(std::move(value));
            return true;
        }
        if (c == 't' || c == 'f')
            return parseBool(out);
        if (c == 'n')
            return parseNull(out);
        return parseNumber(out);
    }

    bool parseObject(Value &out)
    {
        ++pos; // '{'
        out = Value::object();
        skipWhitespace();
        if (pos < text.size() && text[pos] == '}') {
            ++pos;
            return true;
        }
        while (true) {
            skipWhitespace();
            if (pos >= text.size() || text[pos] != '"')
                return fail("expected object key string");
            std::string key;
            if (!parseString(key))
                return false;
            skipWhitespace();
            if (pos >= text.size() || text[pos] != ':')
                return fail("expected ':' after object key");
            ++pos;
            Value value;
            if (!parseValue(value))
                return false;
            out.set(key, std::move(value));
            skipWhitespace();
            if (pos >= text.size())
                return fail("unterminated object");
            if (text[pos] == ',') {
                ++pos;
                continue;
            }
            if (text[pos] == '}') {
                ++pos;
                return true;
            }
            return fail("expected ',' or '}' in object");
        }
    }

    bool parseArray(Value &out)
    {
        ++pos; // '['
        out = Value::array();
        skipWhitespace();
        if (pos < text.size() && text[pos] == ']') {
            ++pos;
            return true;
        }
        while (true) {
            Value value;
            if (!parseValue(value))
                return false;
            out.push(std::move(value));
            skipWhitespace();
            if (pos >= text.size())
                return fail("unterminated array");
            if (text[pos] == ',') {
                ++pos;
                continue;
            }
            if (text[pos] == ']') {
                ++pos;
                return true;
            }
            return fail("expected ',' or ']' in array");
        }
    }

    bool parseString(std::string &out)
    {
        ++pos; // opening quote
        out.clear();
        while (pos < text.size()) {
            const char c = text[pos++];
            if (c == '"')
                return true;
            if (c == '\\') {
                if (pos >= text.size())
                    return fail("unterminated escape");
                const char e = text[pos++];
                switch (e) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    if (pos + 4 > text.size())
                        return fail("truncated unicode escape");
                    unsigned int code = 0;
                    for (int i = 0; i < 4; ++i) {
                        const char h = text[pos++];
                        code <<= 4;
                        if (h >= '0' && h <= '9')
                            code |= static_cast<unsigned int>(h - '0');
                        else if (h >= 'a' && h <= 'f')
                            code |= static_cast<unsigned int>(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F')
                            code |= static_cast<unsigned int>(h - 'A' + 10);
                        else
                            return fail("invalid unicode escape");
                    }
                    // Encode as UTF-8 (BMP only; surrogate pairs are not needed
                    // by this format).
                    if (code < 0x80) {
                        out.push_back(static_cast<char>(code));
                    } else if (code < 0x800) {
                        out.push_back(static_cast<char>(0xC0 | (code >> 6)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    } else {
                        out.push_back(static_cast<char>(0xE0 | (code >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    break;
                }
                default:
                    return fail("invalid escape sequence");
                }
                continue;
            }
            out.push_back(c);
        }
        return fail("unterminated string");
    }

    bool parseBool(Value &out)
    {
        if (text.compare(pos, 4, "true") == 0) {
            pos += 4;
            out = Value::boolean(true);
            return true;
        }
        if (text.compare(pos, 5, "false") == 0) {
            pos += 5;
            out = Value::boolean(false);
            return true;
        }
        return fail("invalid literal");
    }

    bool parseNull(Value &out)
    {
        if (text.compare(pos, 4, "null") == 0) {
            pos += 4;
            out = Value();
            return true;
        }
        return fail("invalid literal");
    }

    bool parseNumber(Value &out)
    {
        const char *start = text.c_str() + pos;
        char *end = nullptr;
        const double value = std::strtod(start, &end);
        if (end == start)
            return fail("invalid number");
        pos += static_cast<std::size_t>(end - start);
        out = Value::number(value);
        return true;
    }
};

} // namespace

void Value::dumpTo(std::string &out, int indent, int depth) const
{
    const auto newline = [&](int d) {
        if (indent > 0) {
            out.push_back('\n');
            out.append(static_cast<std::size_t>(indent * d), ' ');
        }
    };

    switch (m_type) {
    case Type::Null: out += "null"; break;
    case Type::Bool: out += m_bool ? "true" : "false"; break;
    case Type::Number: out += numberToString(m_number); break;
    case Type::String: escapeString(m_string, out); break;
    case Type::Array:
        if (m_array.empty()) {
            out += "[]";
            break;
        }
        out.push_back('[');
        for (std::size_t i = 0; i < m_array.size(); ++i) {
            if (i)
                out.push_back(',');
            newline(depth + 1);
            m_array[i].dumpTo(out, indent, depth + 1);
        }
        newline(depth);
        out.push_back(']');
        break;
    case Type::Object:
        if (m_object.empty()) {
            out += "{}";
            break;
        }
        out.push_back('{');
        for (std::size_t i = 0; i < m_object.size(); ++i) {
            if (i)
                out.push_back(',');
            newline(depth + 1);
            escapeString(m_object[i].first, out);
            out.push_back(':');
            if (indent > 0)
                out.push_back(' ');
            m_object[i].second.dumpTo(out, indent, depth + 1);
        }
        newline(depth);
        out.push_back('}');
        break;
    }
}

std::string Value::dump(int indent) const
{
    std::string out;
    dumpTo(out, indent, 0);
    return out;
}

bool Value::parse(const std::string &text, Value &out, std::string &error)
{
    Parser parser(text);
    if (!parser.parseValue(out)) {
        error = parser.error;
        return false;
    }
    parser.skipWhitespace();
    if (parser.pos != text.size()) {
        error = "trailing characters after JSON value";
        return false;
    }
    return true;
}

} // namespace json
} // namespace femm
