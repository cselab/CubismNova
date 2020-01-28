// File       : INIParser.h
// Created    : Mon Dec 23 2019 09:39:13 PM (+0100)
// Author     : Fabian Wermelinger
// Description: .ini file parser.  This is a modified class based on
//              third/inih/cpp/INIReader.h
// Copyright 2019 ETH Zurich. All Rights Reserved.
#ifndef INIPARSER_H_1WLRUC2B
#define INIPARSER_H_1WLRUC2B

#include "Cubism/Common.h"
#include <map>
#include <ostream>
#include <queue>
#include <string>
#include <vector>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(Util)

/**
 * @ingroup Util
 * @brief INI config file parser
 *
 * @rst
 * Simple configuration file parser for use on the application level.
 * CubismNova does not depend on this parser.  The .ini file format is explained
 * in more detail on Wikipedia_ for example.  The python ``configparser`` module
 * supports this format for convenient simulation case setup using python.
 *
 * .. _Wikipedia: https://en.wikipedia.org/wiki/INI_file
 * @endrst
 * */
class INIParser
{
public:
    /**
     * @brief Construct INIParser and parse given filename.
     * @param filename Configuration file path
     */
    explicit INIParser(const std::string &filename);

    /**
     * @brief Construct INIParser and parse given buffer.
     * @param buffer Configuration character buffer
     * @param buffer_size Number of bytes in buffer
     */
    explicit INIParser(const char *buffer, size_t buffer_size);

    /**
     * @brief Return the global result of ``ini_parse()`` for all files.
     * @return False if no errors occurred during parsing
     */
    bool parseError() const;

    /**
     * @brief Return the result of ``ini_parse()``
     * @return Map for each INI file found in the main configuration
     *
     * @rst
     * For each key in the map: ``0`` on success, line number of first error or
     * ``-1`` on file open error.
     * @endrst
     */
    std::map<std::string, int> fileErrors() const;

    /**
     * @brief Get a string value from INI file.
     * @param section Section name
     * @param name Key name
     * @return String key value
     *
     * Throws a runtime error if not found.
     */
    std::string get(const std::string &section, const std::string &name) const;

    /**
     * @brief Get a string value from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed string key value
     *
     * Throws a runtime error if not found.
     */
    std::string getString(const std::string &section,
                          const std::string &name) const;

    /**
     * @brief Get an array of string values from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed array of strings
     *
     * Throws a runtime error if not found.
     */
    std::vector<std::string> getStringArray(const std::string &section,
                                            const std::string &name) const;

    /**
     * @brief Get an integer (``long``) value from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed integer key value
     *
     * @rst
     * Valid integers are decimal ``1234``, ``-1234``, or hex ``0x4d2``.  Throws
     * a runtime error if not found.
     * @endrst
     */
    long getInteger(const std::string &section, const std::string &name) const;

    /**
     * @brief Get an array of integer values (``long``) from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed array of integers
     *
     * @rst
     * Valid integers are decimal ``1234``, ``-1234``, or hex ``0x4d2``.  Throws
     * a runtime error if not found.
     * @endrst
     */
    std::vector<long> getIntegerArray(const std::string &section,
                                      const std::string &name) const;

    /**
     * @brief Get a real (``double``) value from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed floating point key value
     *
     * @rst
     * Valid floating point values according to ``strtod()``.  Throws a runtime
     * error if not found.
     * @endrst
     */
    double getReal(const std::string &section, const std::string &name) const;

    /**
     * @brief Get an array of real (``double``) values from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed array of floating point values
     *
     * @rst
     * Valid floating point values according to ``strtod()``.  Throws a runtime
     * error if not found.
     * @endrst
     */
    std::vector<double> getRealArray(const std::string &section,
                                     const std::string &name) const;

    /**
     * @brief Get a boolean (``bool``) value from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed boolean value
     *
     * @rst
     * Valid true values are ``true``, ``yes``, ``on``, ``1``, and valid false
     * values are ``false``, ``no``, ``off``, ``0`` (not case sensitive).
     * Throws a runtime error if not found.
     * @endrst
     */
    bool getBoolean(const std::string &section, const std::string &name) const;

    /**
     * @brief Get an array of boolean (``bool``) values from INI file.
     * @param section Section name
     * @param name Key name
     * @return Parsed array of boolean values
     *
     * @rst
     * Valid true values are ``true``, ``yes``, ``on``, ``1``, and valid false
     * values are ``false``, ``no``, ``off``, ``0`` (not case sensitive).
     * Throws a runtime error if not found.
     * @endrst
     */
    std::vector<bool> getBooleanArray(const std::string &section,
                                      const std::string &name) const;

    /**
     * @brief Return ``true`` if the given section exists
     * @param section Section name
     *
     * @rst
     * The section must contain at least one ``name=value`` pair.
     * @endrst
     */
    bool hasSection(const std::string &section) const;

    /**
     * @brief Return true if a value exists
     * @param section Section name
     * @param name Key name
     */
    bool hasValue(const std::string &section, const std::string &name) const;

    /**
     * @brief Write currently loaded configuration to INI file format
     * @param filename Output file path
     */
    void write(const std::string &filename = "runtime.ini") const;

    /**
     * @brief Redirection to ``std::ostream``
     * @param os Output stream
     * @param p INI parser instance
     */
    friend std::ostream &operator<<(std::ostream &os, const INIParser &p);

private:
    std::queue<std::string> fqueue_;    // file queue for [include] files
    std::map<std::string, int> ferror_; // file error map
    std::map<std::string, std::string> values_; // raw values

    /**
     * @brief Key generator for values_ map
     * @param section
     * @param name
     */
    static std::string makeKey(const std::string &section,
                               const std::string &name);

    /**
     * @brief Value handler
     * @param user
     * @param section
     * @param name
     * @param value
     *
     * This implementation of value handler does not append multiple
     * section/value pairs but overwrite existing ones.
     */
    static int valueHandler(void *user,
                            const char *section,
                            const char *name,
                            const char *value);

    /**
     * @brief Handle include sections in config files in a depth-first manner.
     *
     * The included files are processed in the order they appear in the file,
     * overwriting existing keys.
     */
    void recursiveIncludes_(const std::string &filename);
};

std::ostream &operator<<(std::ostream &os, const INIParser &p);

NAMESPACE_END(Util)
NAMESPACE_END(Cubism)

#endif /* INIPARSER_H_1WLRUC2B */
