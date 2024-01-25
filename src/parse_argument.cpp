// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "parse_argument.h"

namespace arg {

// *********************************************************************************
// Function define_options()
// *********************************************************************************
void define_options(bpo::options_description & opts,
                    bpo::options_description & advanced_opts,
                    bpo::positional_options_description& p_opts) {
  // Normal Options
  opts.add_options()
    ("help,h", "Show this help message and exit.")

    ("help-all", "Show the help message with the advanced options and exit.")

    ("extract-sheets,e", bpo::value<std::size_t>(),
     "Extract specified number of adjacent (forming hydrogen bonding each other) "
     "beta-strands from the beta-sheets in any possible patterns. And output the "
     "topology string for them.")

    ("pdb-file,f", bpo::value<std::string>()->required(),
     "Positional option 'PDB_FILE' can also be specified by this option.")

    ("graphviz,g", bpo::value<std::string>(),
     "Output the arrangement of strands inside the sheets in graphviz dot file format. "
     "If '-' is specified, dot file will be printed to standard output and "
     "other outputs (such as substrands or sheets information) will be turned off.")

    ("no-stride-sse,n", bpo::bool_switch()->default_value(false),
     "If specified, use the SSE (Secondary Structure Element) assignments in the PDB_FILE "
     "instead of ones from the STRIDE_FILE. Hydrogen bonding information will be read "
     "from the STRIDE_FILE either way.")

    ("output-file,o", bpo::value<std::string>(),
     "Output file to write the results. If not specified, output to standard output.")

    ("stride-file,s", bpo::value<std::string>(),
     "Positional option 'STRIDE_FILE' can also be specified by this option. "
     "Ignored if '--with-stride' (or just '-a') option is given.")

    ("format-type,t", bpo::value<std::size_t>()->default_value(0),
     "The type of output format. 0: like PDB format, 1: like mmcif format.")

    ("with-stride,w", bpo::bool_switch()->default_value(false),
     "Invoke stride command for the given pdb-file inside this program. "
     "Useful if you have 'stride' command in your $PATH. "
     "When specified, 'stride-file' (or '-s') option will be ignored.")
    ;

  advanced_opts.add_options()
    ("max-mid-residues", bpo::value<unsigned>()->default_value(60),
     "The max number of residues between the first and the last strand.")

    ("max-mid-strands", bpo::value<unsigned>()->default_value(1),
     "The max number of strands on the same sheet between the first and the last strand.")

    ("cutoff-left-score", bpo::value<double>()->default_value(0.6, "0.6"),
     "The cut off value that divides right and left handedness. If the score is greater "
     "than this value, judged as left-handed. Otherwise, right-handed. This value must "
     "be in range [0.0, 1.0].")

    ("min-side-dist", bpo::value<double>()->default_value(1.0, "1.0"),
     "The minimum distance between a triangle between strands and a CA atom to judge the "
     "handedness. CA atoms nearer than this distance to a triangle will NOT be counted.")

    ("apj-max-allowed-jump", bpo::value<unsigned>()->default_value(1),
     "Anti-Parallel strands with larger jumps than this value will be output "
     "as 'rare topology' when '-a' option is specified.")

    ("pcc-min-allowed-jump", bpo::value<unsigned>()->default_value(1),
     "Parallel Crossover Connections with smaller jumps than this value will be output "
     "as 'rare topology' when '-p' option is specified.")
    ;

  // Positional Options
  p_opts.add("pdb-file", 1);
  p_opts.add("stride-file", 1);
}



// *********************************************************************************
// Function parse_arguments()
// *********************************************************************************

bpo::variables_map parse_cmd_arg(int const argc, char* const * argv,
                                   bpo::options_description& opts,
                                   bpo::positional_options_description& p_opts) {
  bpo::variables_map vm;
  bpo::store(bpo::command_line_parser(argc, argv).options(opts).positional(p_opts).run(), vm);
  return vm;
}



// *********************************************************************************
// Function parse_arguments()
// *********************************************************************************

bpo::variables_map parse_arguments(int const argc, char* const * argv) {

  bpo::options_description normal_opts("Option Descriptions");
  bpo::options_description advanced_opts("Advanced Options");
  bpo::positional_options_description p_opts;
  define_options(normal_opts, advanced_opts, p_opts);

  // Prepare the help function object
  Help const help{argv[0], normal_opts, advanced_opts};

  try {

    // If no arguments, just showt the help and exit.
    if (argc == 1) {
      help(std::cout);
      throw help_mode{};
    }

    // merge all option groups
    normal_opts.add(advanced_opts);
    auto vm = parse_cmd_arg(argc, argv, normal_opts, p_opts);

    if (vm.count("help")) {
      help(std::cout);
      throw help_mode{};
    } else if (vm.count("help-all")) {
      help.all(std::cout);
      throw help_mode{};
    }

    bpo::notify(vm);
    return vm;

  } catch (bpo::error const& e) {
    std::cerr << e.what() << "\n"
              << "Try '" << help.exec_name
              << " --help' for more information." << std::endl;
    throw argument_error{e.what()};
  }
} // function parse_arguments()


// *********************************************************************************
// Class Help
// *********************************************************************************

// ********************************************************
// Public Member operator ()
// ********************************************************
void Help::operator()(std::ostream & os, std::string const& msg) const {
  os << (msg != "" ? msg + "\n\n" : "");
  print_head(os);
  print_opts(os);
  print_examples(os);
} // public member operator ()



// ********************************************************
// Public Member Function all()
// ********************************************************

void Help::all(std::ostream & os, std::string const& msg) const {
  os << (msg != "" ? msg + "\n\n" : "");
  print_head(os);
  print_opts(os);
  print_adv_opts(os);
  print_examples(os);
} // public member function all()



// ********************************************************
// Private Member Function print_head()
// ********************************************************

void Help::print_head(std::ostream & os) const {
  os << "Usage: " << exec_name << " [OPTIONS] [PDB_FILE] [STRIDE_FILE]\n"
     << "   or: " << exec_name << " [OPTIONS] -f pdb_file -s stride_file\n"
     << "   or: " << exec_name << " [OPTIONS] -w [PDB_FILE]\n\n"
     << "Output the arrangements of the beta-strands in a sheet, and detect the rare\n"
     << "topologies included in the PDB_FILE.\n\n"
     << "  'PDB_FILE' and 'STRIDE_FILE' are required. You can also specify these paths\n"
     << "  using '--pdb-file' and '--stride-file' options.\n"
     << "  STRIDE_FILE must contain the output of 'stride -h' command.\n"
     << "  If STRIDE_FILE is not given, the input stride file should be\n"
     << "  given from standard input.\n"
     << "\n";
} // private member function print_head()



// ********************************************************
// Private Member Function print_opts()
// ********************************************************

void Help::print_opts(std::ostream & os) const {
  os << normal_opts
     << "\n";
} // private member function print_opts()




// ********************************************************
// Private Member Function print_adv_opts()
// ********************************************************

void Help::print_adv_opts(std::ostream & os) const {
  os << advanced_opts
     << "\n";
} // private member function print_adv_opts()



// ********************************************************
// Private Member Function print_examples()
// ********************************************************
void Help::print_examples(std::ostream & os) const {
  os << "Examples:\n"
     << "  % " << exec_name << " -R example.pdb example.stride\n"
     << "      Enable all the available rare topology detection and output information\n"
     << "      about the found rare topologies in addition to the default output.\n\n"
     << "  % stride -h example.pdb | " << exec_name << " example.pdb\n"
     << "      Give the output of 'stride -h' command through the standard input.\n\n"
     << "  % " << exec_name << " -wg- example.pdb | dot -Tpng > example.png \n"
     << "      Generate a PNG file of the graph, using the graphviz 'dot' command. \n"
     << "      Default output will be discarded. Graphviz needs to be installed on your system."
     << "\n\n"
     << std::endl;
} // private member function print_examples()



// *********************************************************************************
// Function stride_from_argument()
// *********************************************************************************

pdb::stride_stream stride_from_argument(bpo::variables_map const& vm) {
  // Run stride with system command.
  if (vm["with-stride"].as<bool>()) {
    return pdb::pdb2stride_stream(vm["pdb-file"].as<std::string>());
  } else {

    // Open the specified stride file.
    if (vm.count("stride-file")) {
      return pdb::stride_stream{vm["stride-file"].as<std::string>()};

    // Listen to standard input.
    } else {
      return pdb::stride_stream{std::cin};
    }
  }
} // function stride_stream_from_argument()


} // namespace arg
