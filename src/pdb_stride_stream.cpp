// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "pdb/stride_stream.h"

namespace pdb {


// *********************************************************************************
// Function pdb2stride_stream()
// *********************************************************************************

stride_stream pdb2stride_stream(std::string const& pdb_file) {
  // check pdb_file existence
  if (pdb::is_file_exist(pdb_file) == false) {
    throw open_file_error(pdb_file, "In Pairs::run_stride()");
  }

  // set the filename of the stride output
  std::string const outfile = "./stride_" + pdb::basename(pdb_file) + "_" +
                              pdb::rand_str(10) + ".stride";

  // run stride
  std::string const cmd = "stride -h " + pdb_file + " >" + outfile;

  // check the return code of stride command
  auto const ret_code = std::system(cmd.c_str());
  if (ret_code) {
    // remove tmp file (because the destructor will not be called)
    std::remove(outfile.c_str());

    throw stride_failed(pdb::basename(pdb_file), ret_code);
  }

  // Generate stride_stream object
  stride_stream ss{outfile};

  /// Delete temporaly stride_file. (See also http://en.cppreference.com/w/cpp/io/c/remove)
  std::remove(outfile.c_str());

  return ss;
}




} // namespace pdb

