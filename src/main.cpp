// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <iostream>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "adj_out.h"
#include "data_store.h"
#include "functions.h"
#include "handedness.h"
#include "parse_argument.h"
#include "sheets_out.h"
#include "substrands.h"

#include "pdb/sses.h"
#include "pdb/tools.h"
#include "bab/filter.h"
#include "sheet/directed_adjacency_list.h"
#include "sheet/filter.h"


#ifdef DEBUG

class opt_to_stream {
public:
  opt_to_stream(std::ostream & _os, boost::program_options::variables_map const& _vm):
    os{_os}, vm{_vm} {os << std::boolalpha;}

template<typename T>
void print(std::string const& key) {
  if (has(key)) {
    os << key << ":\t" << vm[key].as<T>() << "\n";
  } else {
    print_no_type(key);
  }
}

void print_no_type(std::string const& key) {
  os << key << ":\t" << has(key) << "\n";
}

private:
  bool has(std::string const& key) const { return static_cast<bool>(vm.count(key)); }

  std::ostream & os;
  boost::program_options::variables_map const& vm;
};


void print_arguments(boost::program_options::variables_map const& vm) {
  std::clog << "Debug: Options given\n";
  opt_to_stream opt_to_clog{std::clog, vm};

  opt_to_clog.print_no_type("help");
  opt_to_clog.print_no_type("help-all");
  opt_to_clog.print<std::size_t>("extract-sheets");
  opt_to_clog.print<std::string>("pdb-file");
  opt_to_clog.print<std::string>("graphviz");
  opt_to_clog.print<bool>("no-stride-sse");
  opt_to_clog.print<std::string>("output-file");
  opt_to_clog.print<std::string>("stride-file");
  opt_to_clog.print<bool>("with-stride");

  opt_to_clog.print<unsigned>("max-mid-residues");
  opt_to_clog.print<unsigned>("max-mid-strands");
  opt_to_clog.print<double>("cutoff-left-score");
  opt_to_clog.print<double>("min-side-dist");
  opt_to_clog.print<unsigned>("apj-max-allowed-jump");
  opt_to_clog.print<unsigned>("pcc-min-allowed-jump");
}
#endif


int main(int const argc, char* const * argv) {

  try {
    // Boost Variables Map
    auto const vm = arg::parse_arguments(argc, argv);

    #ifdef DEBUG
    print_arguments(vm);
    #endif

    // Prepare the Directed Adjacency List object
    auto stride = arg::stride_from_argument(vm);

    // If no-stride-sse, use the SSE headers in PDB_FILE
    auto const sses = vm["no-stride-sse"].as<bool>() ?
                        pdb::SSES{vm["pdb-file"].as<std::string>()} :
                        pdb::SSES{vm["pdb-file"].as<std::string>(), stride};

    sheet::DirectedAdjacencyList const dir_adj_list{sses, stride};
    auto const sheet_id_map = substrands::gen_sheet_id_map(dir_adj_list);


    // ***************
    // Graphviz
    // ***************
    if (vm.count("graphviz")) {
      // graphviz-only mode
      if (vm["graphviz"].as<std::string>() == "-") {
        graphviz::adj_list_to_dot(std::cout, dir_adj_list);
        return 0;

      // continue to other output
      } else {
        std::ofstream ofs_graphviz{vm["graphviz"].as<std::string>()};
        graphviz::adj_list_to_dot(ofs_graphviz, dir_adj_list);
      }
    }


    // Prepare an object to store the all output data
    data_store::Data<table::Set> output_data{std::make_tuple(
        table::TBLSubStrand{"substrand",
                            {"SubStrand_ID", "Sheet_ID", "Ini", "End"}},

        table::TBLHelix{"helix", {"SSE_ID", "Ini", "End"}},

        table::TBLSheet{"sheet",
                        {"Sheet_ID", "N_strands", "Cycle", "Undirected",
                         "With_branch", "Consecutive", "All_para", "All_anti",
                         "Member", "Nomenclature_R", "Nomenclature_C"}},

        table::TBLExtractedSheet{"EXT_Sheet",
                                 {"Sheet_ID", "N_strands",
                                  "Same_as_Original", "Member",
                                  "Nomenclature_C"}},

        table::TBLCycle{"cycle", {"Sheet_ID", "N_strands", "Member"}},

        table::TBLSubStrandsPair{"substrands_pair",
                                 {"B1", "B2", "Sheet", "Dir", "PorA",
                                  "Jump", "D1", "D2", "Bridge",
                                  "Score", "SSEs_LBTS", "NumRes_LBTS"}},
        table::TBLResiduePair{"residue_pair", {"ResNum1", "ResNum2",
                                               "PorA", "Pair-type", "ForB"}}
    )};

    // set default out (if out_stdout)

    // for output into a file
    std::ofstream ofs;
    bool const out_stdout = vm.count("output-file") == 0;
    if (not out_stdout) {
      ofs.open(vm["output-file"].as<std::string>());
    }
    std::ostream & out_stream = out_stdout ? std::cout : ofs;


    // ***************
    // Cycles
    // ***************
    cycles::output_cycles(output_data.table<table::Cycle>(), dir_adj_list);


    // ***************
    // Default Output
    // ***************
    substrands::substrands_out(output_data.table<table::SubStrand>(), dir_adj_list,
                               sheet_id_map);
    substrands::helices_out(output_data.table<table::Helix>(), dir_adj_list);
    sheets_out::print_sheet(output_data.table<table::Sheet>(), dir_adj_list);

    auto bab_filter = bab::BabFilter{dir_adj_list, std::greater<double>(),
                                     vm["max-mid-residues"].as<unsigned>(),
                                     vm["max-mid-strands"].as<unsigned>(),
                                     vm["cutoff-left-score"].as<double>(),
                                     vm["min-side-dist"].as<double>()};

    substrands::substrands_pair_out(output_data.table<table::SubStrandsPair>(),
                                    dir_adj_list, sheet_id_map, bab_filter);


    // ***************
    // Extract Sheet
    // ***************
    if (vm.count("extract-sheets")) {
      sheets_out::extracted_adjacent_substr_out(output_data.table<table::ExtractedSheet>(),
                                                vm["extract-sheets"].as<std::size_t>(),
                                                dir_adj_list);
    }


    // ***************
    // Residue Pairs
    // ***************
    rpo::residue_pair_out(output_data.table<table::ResiduePair>(), dir_adj_list);


    // Actually output the results
    output_data.format_out(out_stream, vm["format-type"].as<std::size_t>());

    // if mmcif output
    if (vm["format-type"].as<std::size_t>() == 1) {
      adj_out::adj_list_out(out_stream, dir_adj_list);
      rare::output_handedness(out_stream, dir_adj_list, vm);
    }


  // if help mode ('--help')
  } catch (arg::help_mode const&) {
    return 0;

  // If exited while processing arguments or 'help' mode.
  } catch (arg::argument_error const&) {
    return 1;

  #ifdef DEBUG

  } catch (...) {
    throw;
  }

  #else // ifdef DEBUG

  // For Other Possible Errors
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return 2;

  // For Unknown Errors
  } catch (...) {
    std::cerr << "Unknown Error" << std::endl;
    return 2;
  }

  #endif // ifdef DEBUG
}

