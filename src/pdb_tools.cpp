// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "pdb/exceptions.h"
#include "pdb/tools.h"

namespace pdb {

// ****************************************************************************************
// Function open_input()
// ****************************************************************************************

void open_input(std::ifstream& ifs, std::string const& filename) {
  ifs.open(filename);
  if (ifs.is_open() == false) {
    throw open_file_error(filename);
  }
}



// ****************************************************************************************
// Function split()
// ****************************************************************************************

std::vector<std::string> split(std::string const& str ,std::string const& delm) {
    std::vector<std::string> ret;

    std::regex delim(delm);
    auto first = std::sregex_token_iterator( str.cbegin(), str.cend(), delim, -1 );
    auto last  = std::sregex_token_iterator();

    for ( auto itr = first; itr != last; itr++ ) {
      if (itr->str().size() != 0) {
        ret.push_back( itr->str() );
      }
    }
    return ret;
}



// ****************************************************************************************
// Function split()
// ****************************************************************************************

std::vector<std::string> split(std::string const& str, char const delm) {
  std::istringstream is(str);

  std::vector<std::string> ret;
  for (std::string buff; std::getline(is, buff, delm);) {
    ret.push_back(buff);
  }
  return ret;
}



// ****************************************************************************************
// Function rand_str()
// ****************************************************************************************

std::string rand_str(unsigned const len) {
  static std::string const alnum = 
          "abcdefghijklmnopqrstuvwxyz"
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
          "0123456789";

  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<> dist(0, alnum.size() - 1);

  std::string ret;
  ret.reserve(len);

  std::generate_n(std::back_inserter(ret), len, 
                  [&]() { return alnum[dist(rng)]; });

  return ret;
}




// ****************************************************************************************
// Function basename()
// ****************************************************************************************

std::string basename(std::string const& path) {
  return path.substr(path.find_last_of('/') + 1);
} // function basename()




// ****************************************************************************************
// Function is_file_exist()
// ****************************************************************************************

bool is_file_exist(std::string const& filename) {
  return std::ifstream(filename).is_open();
} // function is_file_exists()



// ****************************************************************************************
// Function is2ss()
// ****************************************************************************************

void is2ss(std::istream & is, std::stringstream & ss) {
  ss << is.rdbuf();
}



// ****************************************************************************************
// Function log()
// ****************************************************************************************

#ifdef LOGGING
void log(std::string const& msg) {
  constexpr char const start[] = "\033[1;38;5;250m";
  constexpr char const end[] = "\033[00m";
  std::clog << start << "[LOG]   " << end
            << msg << std::endl;
}
#else
void log(std::string) {}
#endif // ifdef LOGGING



// ****************************************************************************************
// Function debug_log()
// ****************************************************************************************

#ifdef DEBUG
void debug_log(std::string const& msg) {
  constexpr char const start[] = "\033[1;38;5;107m";
  constexpr char const end[] = "\033[00m";
  std::clog << start << "[DEBUG] " << end
            << msg << std::endl;
}
#else
void debug_log(std::string) {}
#endif // ifdef DEBUG




// ****************************************************************************************
// Function warning()
// ****************************************************************************************

#ifndef QUIET
void warning(std::string const& msg) {
  constexpr char const start[] = "\033[1;38;5;185m";
  constexpr char const end[] = "\033[00m";
  std::clog << start << "[WARNING] " << end
            << msg << std::endl;
}
#else
void warning(std::string) {}
#endif // ifndef QUIET




// ****************************************************************************************
// Function gen_seq_file()
// ****************************************************************************************
#ifdef NATIVE_DRYRUN
std::string gen_seq_file(std::string const& pdb_file) {
  std::string const tmp_seq_file = "./stride_seq_" + rperm::basename(pdb_file) + "-" +
                                   rperm::rand_str(5) + ".seq";

  std::string const cmd = "(echo -n 'SSseq:' && stride " + pdb_file +
                          "|grep --color=never '^ASG'|cut -b 25|sed 's/[^HE]/-/g'|tr -d '\\n'|"
                          "sed 's/HE/--/g'|sed 's/EH/--/g') >"+
                          tmp_seq_file;
  auto const ret = std::system(cmd.c_str());
  if (ret != 0) {
    throw std::runtime_error("stride returned non-zero exit status '" + std::to_string(ret) + "'");
  }
  return tmp_seq_file;
}

#endif // ifdef NATIVE_DRYRUN


#if defined(DRYRUN) || defined(NATIVE_DRYRUN)
// ****************************************************************************************
// Function gen_seq_file()
// ****************************************************************************************
unsigned long long factorial(unsigned long long const n) {
  if (n == 0) {
    return 1;
  }
  return n * factorial(n - 1);
}
#endif //ifdef DRYRUN or NATIVE_DRYRUN



} // namespace rperm

