// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef RPERM_TOOLS_H_
#define RPERM_TOOLS_H_


#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <vector>

#include <boost/range/irange.hpp>


namespace pdb {

/// open filename and check for errors
void open_input(std::ifstream& ifs, std::string const& filename);


/// split with regexp delm
std::vector<std::string> split(std::string const& str ,std::string const& delm);


/// split with char delm
std::vector<std::string> split(std::string const& str, char const delm);


/// Generate a random string of specified length consists of alnum ([a-zA-Z0-9])
/// @return random_string
/// @param len the number of characters returned
std::string rand_str(unsigned const len);



/// @brief  Return the basename of given path
/// @return basename of given path
/// @retval filename.txt if path = "/home/hoge/dir/filename.txt"
/// @param  path a string which includes the basename we want
std::string basename(std::string const& path);



/// @brief            Check if the filename is a file and exists
/// @return whether   the file exists
/// @retval true      if exists
/// @retval fales     if not exists
/// @param  filename  a filepath to check
bool is_file_exist(std::string const& filename);


/// @brief  Read from is into a string stream ss.
void is2ss(std::istream & is, std::stringstream & ss);


/// boost::irange wrapper
template <typename T>
decltype(auto) range(T n) {
  return boost::irange(static_cast<T>(0), n);
}



#ifdef LOGGING
/// @beirf            Output log with unified format
void log(std::string const& msg);
#else
void log(std::string);
#endif // ifdef LOGGING



#ifdef DEBUG
void debug_log(std::string const& msg);
#else
void debug_log(std::string);
#endif // ifdef DEBUG


#ifndef QUIET
void warning(std::string const& msg);
#else
void warning(std::string);
#endif // ifndef QUIET



/// for dryrunmode
#ifdef NATIVE_DRYRUN
std::string gen_seq_file(std::string const& pdb_file);
#endif // ifdef NATIVE_DRYRUN

#if defined(DRYRUN) || defined(NATIVE_DRYRUN)
unsigned long long factorial(unsigned long long const n);
#endif // ifdef DRYRUN




} // namespace rperm


#endif // ifndef RPERM_TOOLS_H_
