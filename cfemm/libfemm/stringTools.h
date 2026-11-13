/* Copyright 2016-2017 Johannes Zarl-Zierl <johannes.zarl-zierl@jku.at>
 * Contributions by Johannes Zarl-Zierl were funded by Linz Center of
 * Mechatronics GmbH (LCM)
 *
 * The source code in this file is heavily derived from
 * FEMM by David Meeker <dmeeker@ieee.org>.
 * For more information on FEMM see http://www.femm.info
 * This modified version is not endorsed in any way by the original
 * authors of FEMM.
 *
 * License:
 * This software is subject to the Aladdin Free Public Licence
 * version 8, November 18, 1999.
 * The full license text is available in the file LICENSE.txt supplied
 * along with the source code.
 */
#ifndef FEMM_STRINGTOOLS_H
#define FEMM_STRINGTOOLS_H

#include <string>
#include <algorithm>
#include <locale>

/**
 * \file stringTools.h
 * \brief String functions for use within femm
 * I've added these functions here so that no external libraries like boost are needed.
 */

namespace femm
{
/**
 * @brief to_lower transforms a string to lower case.
 * @param s the string to be transformed in place
 */
inline void to_lower(std::string &s)
{
    // transform s to lower case
    transform(s.begin(), s.end(), s.begin(), ::tolower);
}

/**
 * @brief ltrim trims whitespace from the beginning of a string
 * \note taken from https://stackoverflow.com/posts/217605/revisions
 * @param s the string to be transformed in place
 */
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [](unsigned char c) { return !std::isspace(c); }));
}

/**
 * @brief rtrim trims whitespace from the end of a string
 * \note taken from https://stackoverflow.com/posts/217605/revisions
 * @param s the string to be transformed in place
 */
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](unsigned char c) { return !std::isspace(c); }).base(), s.end());
}

/**
 * @brief trim whitespace from the beginning and end of a string
 * @param s the string to be transformed in place
 */
inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

/**
 * @brief Check whether the string \c value ends with the string \c ending.
 * @param value
 * @param ending
 * @return
 */
inline bool ends_with(std::string const & value, std::string const & ending)
{
    // taken from https://stackoverflow.com/questions/874134/find-if-string-ends-with-another-string-in-c
    if (ending.size() > value.size())
         return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/**
 * @brief Check whether the string \c value starts with the string \c beginning.
 * @param value
 * @param beginning
 * @return
 */
inline bool begins_with(std::string const & value, std::string const & beginning)
{
    // taken from https://stackoverflow.com/questions/874134/find-if-string-ends-with-another-string-in-c
    if (beginning.size() > value.size())
         return false;
    return std::equal(beginning.begin(), beginning.end(), value.begin());
}
}
#endif
