// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <numeric>
#include <ostream>
#include <vector>
#include <boost/program_options.hpp>

#include "pdb/tools.h"
#include "sheet/directed_adjacency_list.h"
#include "bab/filter.h"
#include "functions.h"
#include "handedness.h"

namespace bpo = boost::program_options;

namespace rare {
 
// **********************************************************************************
// Function print_handedness()
// **********************************************************************************

std::vector<bab::BabFilterResult> get_handedness(sheet::DirectedAdjacencyList const& adj,
                                                 bab::BabFilter & bab_filter) {
  std::vector<bab::BabFilterResult> found_bab{};

  auto const n_sse = adj.sses.size;
  std::vector<pdb::IndexType> pseudo_seq(n_sse);
  std::iota(pseudo_seq.begin(), pseudo_seq.end(), 0);

  auto const b = pseudo_seq.cbegin();

  // b+i: an iterator to the i-th element in pseudo_seq.
  for (auto const i : pdb::range(n_sse)) {
    // b+j: past the end iterator
    for (auto const j : boost::irange(i+1, n_sse+1)) {
      bab_filter(b+i, b+j);
      auto const& result = bab_filter.result();
      // If the filter ran correctly,
      // no matter whether the Beta-Alpha-Beta is right handed or not.
      if (result.success) {
        found_bab.push_back(result);
      }
    }
  }

  return found_bab;
} // function get_handedness().



// **********************************************************************************
// Function output_handedness()
// **********************************************************************************

bool output_handedness(std::ostream & os, sheet::DirectedAdjacencyList const& adj,
                       bpo::variables_map const& vm) {
  bab::BabFilter bab_filter{adj, std::greater<double>(),
                            vm["max-mid-residues"].as<unsigned>(),
                            vm["max-mid-strands"].as<unsigned>(),
                            vm["cutoff-left-score"].as<double>(),
                            vm["min-side-dist"].as<double>()};

  auto const found_bab = get_handedness(adj, bab_filter);

  mmcif::mmcif_like mmcif_out{os, "handedness"};
  mmcif_out.key_value("num", found_bab.size());

  if (found_bab.size()) {
    mmcif_out.loop_head({"SubStrand_0", "SubStrand_1", "score", "mid_residues",
                         "mid_strands", "type", "jump"});

    for (auto const& result : found_bab) {
      assert(result.connection_type);

      std::string type_string = "";
      switch (result.connection_type) {
        case 1:
        case 3:
          type_string = "beta-alpha-beta";
          break;
        case 2:
          type_string = "beta-loop-beta";
          break;
        case 4:
        case 6:
          type_string = "beta-beta-beta";
          break;
        case 5:
        case 7:
          type_string = "beta-alpha(+beta)-beta";
          break;
      }

      os << boost::format("%s %|8t|%s %|16t|%4.2f %3d %2d %22s %d\n")
             % result.sub_first.string()
             % result.sub_last.string()
             % result.left_score
             % result.mid_res_len
             % result.n_mid_str
             % type_string
             % result.jump;
    }
    return true;
  } else {
    return false;
  }
} // function output_handedness().

} // namespace rare

