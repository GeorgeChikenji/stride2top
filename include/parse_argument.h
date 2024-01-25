// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PARSE_ARGUMENT_H_
#define PARSE_ARGUMENT_H_

#include <iostream>

#include <boost/program_options.hpp>

#include "pdb/exceptions.h"
#include "pdb/stride_stream.h"
#include "pdb/tools.h"

namespace bpo = boost::program_options;


namespace arg {

/// An exception class for any errors in parsing arguments.
class argument_error: public pdb::fatal_error_base {
public:
  explicit argument_error(std::string const& msg_=""):
    pdb::fatal_error_base{"Error while parsing arguments.: " + msg_} {}
};


/// An exception class thrown in the help mode, indicates the help mode has just finished.
class help_mode: public argument_error {
public:
  help_mode():
    argument_error("") {}
};



/// A helper function to parse_arguments().
/// Define the acceptable options and return boost::program_options::options_description.
/// !! pdb-file OPTION MUST BE ADDED AFTER THIS FUNCTION !!
void define_options(bpo::options_description & opts,
                    bpo::options_description & advanced_opts,
                    bpo::positional_options_description& p_opts);



/// Parse the given command line options.
bpo::variables_map parse_cmd_arg(int const argc, char* const * argv,
                                 bpo::options_description& opts,
                                 bpo::positional_options_description& p_opts);


/// An interface function.
bpo::variables_map parse_arguments(int const argc, char* const * argv);


/// A function object to output help message.
class Help {
public:
  Help(std::string const& exec_path, bpo::options_description const& opts,
       bpo::options_description const& adv_ots):
    exec_name{pdb::basename(exec_path)},
    normal_opts{opts}, advanced_opts{adv_ots} {}

  std::string const exec_name{""};
  bpo::options_description const normal_opts{};
  bpo::options_description const advanced_opts{};

  /// Show help message using a boost option description object
  void operator()(std::ostream & os=std::clog, std::string const& msg="") const;

  /// Show help messages including the advanced options.
  void all(std::ostream & os=std::clog, std::string const& msg="") const;

private:
  void print_head(std::ostream & os) const;
  void print_opts(std::ostream & os) const;
  void print_adv_opts(std::ostream & os) const;
  void print_examples(std::ostream & os) const;
};


/// Get stride_stream based on the pattern of arguements given.
pdb::stride_stream stride_from_argument(bpo::variables_map const& vm);


} // namespace arg

#endif // ifndef PARSE_ARGUMENT_H_
