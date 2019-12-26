// File       : INIParser.h
// Created    : Mon Dec 23 2019 09:39:13 PM (+0100)
// Author     : Fabian Wermelinger
// Description: .ini file parser.  This is a modified class based on
//              third/inih/cpp/INIReader.h
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef INIPARSER_H_1WLRUC2B
#define INIPARSER_H_1WLRUC2B

#include "Core/Common.h"

#include <map>
#include <ostream>
#include <queue>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

class INIParser
{
public:
    /// @brief Construct INIParser and parse given filename. See ini.h for more
    /// info about the parsing.
    ///
    /// @param filename
    explicit INIParser(const std::string &filename);

    /// @brief Construct INIParser and parse given buffer. See ini.h for more
    /// info about the parsing.
    ///
    /// @param buffer
    /// @param buffer_size
    explicit INIParser(const char *buffer, size_t buffer_size);

    /// @brief Return the global result of ini_parse() for all files.
    ///
    /// @return false if no errors occurred during parsing
    bool parseError() const;

    /// @brief Return the result of ini_parse().
    ///
    /// @return 0 on success, line number of first error or -1 on file open
    /// error for each file
    std::map<std::string, int> fileErrors() const;

    /// @brief Get a string value from INI file.  Throws runtime error if not
    /// found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return string value
    std::string get(const std::string &section, const std::string &name) const;

    /// @brief Get a string value from INI file.  Throws runtime error if not
    /// found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed string
    std::string getString(const std::string &section,
                          const std::string &name) const;

    /// @brief Get an array of string values from INI file.  Throws runtime
    /// error if not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed array of strings
    std::vector<std::string> getStringArray(const std::string &section,
                                            const std::string &name) const;

    /// @brief Get an integer (long) value from INI file.  Valid integers are
    /// decimal "1234", "-1234", or hex "0x4d2".  Throws runtime error if not
    /// found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed integer
    long getInteger(const std::string &section, const std::string &name) const;

    /// @brief Get an array of integer (long) values from INI file.  Valid
    /// integers are decimal "1234", "-1234", or hex "0x4d2".  Throws runtime
    /// error if not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed array of integer
    std::vector<long> getIntegerArray(const std::string &section,
                                      const std::string &name) const;

    /// @brief Get a real (floating point double) value from INI file. Valid
    /// floating point values according to strtod().  Throws runtime error if
    /// not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed floating point value
    double getReal(const std::string &section, const std::string &name) const;

    /// @brief Get an array of real (floating point double) values from INI
    /// file. Valid floating point values according to strtod().  Throws runtime
    /// error if not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed array of floating point values
    std::vector<double> getRealArray(const std::string &section,
                                     const std::string &name) const;

    /// @brief Get a boolean value from INI file.  Valid true values are "true",
    /// "yes", "on", "1", and valid false values are "false", "no", "off", "0"
    /// (not case sensitive).  Throws runtime error if not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed boolean value
    bool getBoolean(const std::string &section, const std::string &name) const;

    /// @brief Get an array of boolean values from INI file.  Valid true values
    /// are "true", "yes", "on", "1", and valid false values are "false", "no",
    /// "off", "0" (not case sensitive).  Throws runtime error if not found.
    ///
    /// @param section
    /// @param name
    ///
    /// @return parsed array of boolean values
    std::vector<bool> getBooleanArray(const std::string &section,
                                      const std::string &name) const;

    /// @brief Return true if the given section exists (section must contain at
    /// least one name=value pair).
    ///
    /// @param section
    bool hasSection(const std::string &section) const;

    /// @brief Return true if a value exists with the given section and field
    /// names.
    ///
    /// @param section
    /// @param name
    bool hasValue(const std::string &section, const std::string &name) const;

    /// @brief Write currently loaded configuration to INI file format
    ///
    /// @param filename
    void write(const std::string &filename = "runtime.ini") const;

    /// @brief Redirection to ostream
    ///
    /// @param os
    /// @param p
    friend std::ostream &operator<<(std::ostream &os, const INIParser &p);

private:
    std::queue<std::string> fqueue_;    // file queue for [include] files
    std::map<std::string, int> ferror_; // file error map
    std::map<std::string, std::string> values_; // raw values

    /// @brief Key generator for values_ map
    ///
    /// @param section
    /// @param name
    static std::string makeKey(const std::string &section,
                               const std::string &name);

    /// @brief This implementation of value handler does not append multiple
    /// section/value pairs but overwrite existing ones.
    ///
    /// @param user
    /// @param section
    /// @param name
    /// @param value
    static int valueHandler(void *user,
                            const char *section,
                            const char *name,
                            const char *value);

    /// @brief Handle include sections in config files in a depth-first manner.
    /// The included files are processed in the order they appear in the file,
    /// overwriting existing keys.
    void recursiveIncludes_(const std::string &filename);
};

std::ostream &operator<<(std::ostream &os, const INIParser &p);

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* INIPARSER_H_1WLRUC2B */
