// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PDB_EXCEPTIONS_H_
#define PDB_EXCEPTIONS_H_

#include <iostream>
#include <stdexcept>
#include <string>

#include "pdb/constants.h"
#include "pdb/tools.h"

namespace pdb {

// ****************************************************************************************
// Base Exception Class Definitions
// ****************************************************************************************

class exception_base: public std::exception {
public:
  explicit exception_base(std::string const& msg_): std::exception{}, msg{msg_} {}

  char const * what() const noexcept override {
    return msg.c_str();
  }

  void say() const {
    std::clog << "\033[1;38;5;250m"  << msg << "\033[00m" << std::endl;
  }

protected:
  std::string const msg;
};

// *******************************************************************
// Base Exception Class fatal_error_base()
// *******************************************************************
class fatal_error_base: public exception_base {
public:
  explicit fatal_error_base(std::string const& msg_): exception_base{"\033[1;38;5;160m[ERROR]: \033[00m" + msg_} {}
}; // Exception class fatal_error_base




// *******************************************************************
// Base Exception Class warning_base()
// *******************************************************************
class warning_base : public exception_base {
public:
  explicit warning_base(std::string const& msg_): exception_base{"\033[1;38;5;3m[WARNING]: \033[00m" + msg_} {}
}; // Exception class warning_base




// ****************************************************************************************
// Derived Exception Class Definitions
// ****************************************************************************************

// *******************************************************************
// Exception class invalid_sse_range
// *******************************************************************
//
class invalid_sse_range: public fatal_error_base {
public:
  invalid_sse_range(int const init, int const end, std::string const& msg_=""):
    fatal_error_base("INVALID SSE RANGE: "
                     "init = '" + std::to_string(init) + "', "
                     "end = '" + std::to_string(end) + "' : " + msg_) {}
}; // Exception class invalid_sse_range




// *******************************************************************
// Exception class unknown_sse_type
// *******************************************************************
//
class unknown_sse_type: public fatal_error_base {
public:
  explicit unknown_sse_type(char const type, std::string const& msg_=""):
    fatal_error_base{"Unknown SSE type '" + std::string{type} + "'. " + msg_} {}
}; // Exception class unknown_sse_type



// *******************************************************************
// Exception class open_file_error
// *******************************************************************
//
class open_file_error: public fatal_error_base {
public:
  explicit open_file_error(std::string const& filename, std::string const& msg_=""):
    fatal_error_base("CANNOT OPEN FILE '" + filename + "': " + msg_){}
}; // Exception class opne_file_error



// *******************************************************************
// Exception class non_sse_resnum
// *******************************************************************
//
class non_sse_resnum: public warning_base {
public:
  explicit non_sse_resnum(int const resnum, std::string const& msg_="")
    : warning_base("Resnum '" + std::to_string(resnum) + "' is not in any SSEs. : " + msg_) {}
}; // Exception class non_sse_resnum



// *******************************************************************
// Exception class padding_atom_found
// *******************************************************************
//
class padding_atom_found: public warning_base {
public:
  explicit padding_atom_found(std::string const& msg_=""):
    warning_base{"Padding ATOM found during generating representative ATOMs. "
                 "Represendative ATOMs will be turned off.: " + msg_} {}
};


// *******************************************************************
// Exception class stride_failed
// *******************************************************************
//
class stride_failed: public fatal_error_base {
public:
  stride_failed(std::string const& filename, int const ret_code):
    fatal_error_base{"COMMAND 'stride " + filename +
                     "' RETURNED NON-ZERO EXIT STATUS '" + std::to_string(ret_code) + "'"} {}
}; // Exception class stride_failed



// *******************************************************************
// Exception class resnum_out_of_range
// *******************************************************************
//
class resnum_out_of_range: public fatal_error_base {
public:
  resnum_out_of_range(int const resnum):
    fatal_error_base{"RESNUM '" + std::to_string(resnum) + "' is out of range."} {}
}; // Exception class resnum_out_of_range



// *******************************************************************
// Exception class loop_access_out_of_range
// *******************************************************************
//
class loop_access_out_of_range: public fatal_error_base {
public:
  loop_access_out_of_range(std::size_t const loop_id, std::size_t const max):
    fatal_error_base{"ACCESSING OUT OF RANGE LOOP[" + std::to_string(loop_id) + "]. "
                     "(MAX = " + std::to_string(max) + ")"} {}
}; // Exception class loop_access_out_of_range



} // namespace pdb
#endif // ifndef PDB_EXCEPTIONS_H_
