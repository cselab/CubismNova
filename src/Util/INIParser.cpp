// File       : INIParser.cpp
// Created    : Mon Dec 23 2019 10:20:51 PM (+0100)
// Author     : Fabian Wermelinger
// Description: .ini file parser implementation.
// Copyright 2019 ETH Zurich. All Rights Reserved.
#include "Util/INIParser.h"
#include "third/inih/ini.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <stdexcept>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

INIParser::INIParser(const std::string &filename)
{
    recursiveIncludes_(filename);
}

INIParser::INIParser(const char *buffer, size_t buffer_size)
{
    std::string content(buffer, buffer_size);
    ferror_["fbuffer"] = ini_parse_string(content.c_str(), valueHandler, this);
}

bool INIParser::parseError() const
{
    bool err = false;
    for (const auto &f : ferror_) {
        err = f.second;
    }
    return err;
}

std::map<std::string, int> INIParser::fileErrors() const { return ferror_; }

std::string INIParser::get(const std::string &section,
                           const std::string &name) const
{
    std::string key = makeKey(section, name);
    if (!values_.count(key)) {
        throw std::runtime_error("get: key=" + name + " in section=" + section +
                                 " does not exist");
    }

    std::string val = values_.find(key)->second;
    if (val.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
        throw std::runtime_error("get: key=" + name + " in section=" + section +
                                 " has no value");
    }
    return val;
}

std::string INIParser::getString(const std::string &section,
                                 const std::string &name) const
{
    return get(section, name);
}

std::vector<std::string>
INIParser::getStringArray(const std::string &section,
                          const std::string &name) const
{
    std::vector<std::string> ret;
    std::istringstream iss(get(section, name));
    for (std::string val; iss >> val;) {
        ret.push_back(val);
    }
    if (ret.empty()) {
        // this code should never be executed since get() does not allow empty
        // values
        throw std::runtime_error("getStringArray: empty container");
    }
    return ret;
}

long INIParser::getInteger(const std::string &section,
                           const std::string &name) const
{
    const std::string valstr = get(section, name);
    const char *value = valstr.c_str();
    char *end;
    // This parses "1234" (decimal) and also "0x4D2" (hex)
    long n = strtol(value, &end, 0);
    if (end <= value) {
        throw std::runtime_error("getInteger: can not convert '" + valstr +
                                 "' to integer for key=" + name +
                                 " in section=" + section);
    }
    return n;
}

std::vector<long> INIParser::getIntegerArray(const std::string &section,
                                             const std::string &name) const
{
    std::vector<long> ret;
    std::istringstream iss(get(section, name));
    for (std::string valstr; iss >> valstr;) {
        const char *value = valstr.c_str();
        char *end;
        // This parses "1234" (decimal) and also "0x4D2" (hex)
        long n = strtol(value, &end, 0);
        if (end <= value) {
            throw std::runtime_error("getIntegerArray: can not convert '" +
                                     valstr + "' to integer for key=" + name +
                                     " in section=" + section);
        }
        ret.push_back(n);
    }
    if (ret.empty()) {
        // this code should never be executed since get() does not allow empty
        // values
        throw std::runtime_error("getIntegerArray: empty container");
    }
    return ret;
}

double INIParser::getReal(const std::string &section,
                          const std::string &name) const
{
    const std::string valstr = get(section, name);
    const char *value = valstr.c_str();
    char *end;
    double n = strtod(value, &end);
    if (end <= value) {
        throw std::runtime_error("getReal: can not convert '" + valstr +
                                 "' to floating point for key=" + name +
                                 " in section=" + section);
    }
    return n;
}

std::vector<double> INIParser::getRealArray(const std::string &section,
                                            const std::string &name) const
{
    std::vector<double> ret;
    std::istringstream iss(get(section, name));
    for (std::string valstr; iss >> valstr;) {
        const char *value = valstr.c_str();
        char *end;
        double n = strtod(value, &end);
        if (end <= value) {
            throw std::runtime_error("getRealArray: can not convert '" +
                                     valstr + "' to floating point for key=" +
                                     name + " in section=" + section);
        }
        ret.push_back(n);
    }
    if (ret.empty()) {
        // this code should never be executed since get() does not allow empty
        // values
        throw std::runtime_error("getRealArray: empty container");
    }
    return ret;
}

bool INIParser::getBoolean(const std::string &section,
                           const std::string &name) const
{
    std::string valstr = get(section, name);
    // Convert to lower case to make std::string comparisons case-insensitive
    std::transform(valstr.begin(), valstr.end(), valstr.begin(), ::tolower);
    if (valstr == "true" || valstr == "yes" || valstr == "on" ||
        valstr == "1") {
        return true;
    } else if (valstr == "false" || valstr == "no" || valstr == "off" ||
               valstr == "0") {
        return false;
    } else {
        throw std::runtime_error("getBoolean: can not convert '" + valstr +
                                 "' to boolean for key=" + name +
                                 " in section=" + section);
    }
}

std::vector<bool> INIParser::getBooleanArray(const std::string &section,
                                             const std::string &name) const
{
    std::string valstr = get(section, name);
    // Convert to lower case to make std::string comparisons case-insensitive
    std::transform(valstr.begin(), valstr.end(), valstr.begin(), ::tolower);
    std::istringstream iss(valstr);
    std::vector<bool> ret;
    for (std::string val; iss >> val;) {
        bool bval;
        if (val == "true" || val == "yes" || val == "on" || val == "1") {
            bval = true;
        } else if (val == "false" || val == "no" || val == "off" ||
                   val == "0") {
            bval = false;
        } else {
            throw std::runtime_error("getBooleanArray: can not convert '" +
                                     val + "' to boolean for key=" + name +
                                     " in section=" + section);
        }
        ret.push_back(bval);
    }
    if (ret.empty()) {
        // this code should never be executed since get() does not allow empty
        // values
        throw std::runtime_error("getBooleanArray: empty container");
    }
    return ret;
}

bool INIParser::hasSection(const std::string &section) const
{
    const std::string key = makeKey(section, "");
    std::map<std::string, std::string>::const_iterator pos =
        values_.lower_bound(key);
    if (pos == values_.end()) {
        return false;
    }
    // Does the key at the lower_bound pos start with "section"?
    return pos->first.compare(0, key.length(), key) == 0;
}

bool INIParser::hasValue(const std::string &section,
                         const std::string &name) const
{
    const std::string key = makeKey(section, name);
    return values_.count(key);
}

std::string INIParser::makeKey(const std::string &section,
                               const std::string &name)
{
    return section + "=" + name;
}

int INIParser::valueHandler(void *user,
                            const char *section,
                            const char *name,
                            const char *value)
{
    if (!name) // Happens when INI_CALL_HANDLER_ON_NEW_SECTION enabled
        return 1;
    INIParser *parser = static_cast<INIParser *>(user);

    // append [include] section to file queue
    // XXX: [fabianw@mavt.ethz.ch; 2019-12-24] THIS IS NOT THREAD SAFE!
    std::string sec(section);
    std::transform(sec.begin(), sec.end(), sec.begin(), ::tolower);
    if (sec == "include" && value) {
        parser->fqueue_.push(value);
        // append existing key.
        const std::string key = makeKey(sec, name);
        if (parser->values_.count(key)) {
            parser->values_[key] += " | ";
        }
        parser->values_[key] += value;
    } else {
        // overwrite existing key.
        const std::string key = makeKey(section, name);
        parser->values_[key] = value ? value : "";
    }
    return 1;
}

void INIParser::recursiveIncludes_(const std::string &filename)
{
    const auto it = ferror_.find(filename);
    if (it != ferror_.end()) {
        throw std::runtime_error(
            "recursiveIncludes_: Cyclic inclusion of file '" + filename + "'");
    }
    ferror_[filename] = ini_parse(filename.c_str(), valueHandler, this);
    if (fqueue_.empty()) {
        return;
    } else {
        const std::string fnext = fqueue_.front();
        fqueue_.pop();
        recursiveIncludes_(fnext);
    }
}

std::ostream &operator<<(std::ostream &os, const INIParser &p)
{
    for (const auto &item : p.values_) {
        std::istringstream iss(item.first);
        std::string section, name;
        std::getline(iss, section, '=');
        std::getline(iss, name, '=');
        os << "[" << section << "]: " << name << " = " << item.second << '\n';
    }
    return os;
}

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)
