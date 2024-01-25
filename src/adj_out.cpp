// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <ostream>
#include <string>
#include <unordered_map>

#include "adj_out.h"
#include "functions.h"

#include "sheet/adj_list_with_sub.h"
#include "sheet/directed_adjacency_list.h"
#include "sheet/sheets.h"

namespace adj_out {


// **********************************************************************************
// Function adj_list_out()
// **********************************************************************************

void adj_list_out(std::ostream & os, sheet::DirectedAdjacencyList const& adj) {
  mmcif::mmcif_like const out{os, "adjacency_list"};

  auto const adj_vec = gen_adj_list_vec(adj);
  out.key_value("num", adj_vec.size());

  // if no adj_lists
  if (adj_vec.size() == 0) {
    return;
  }

  out.loop_head({"sheet_id", "direction", "delta_1", "delta_2", "num_bridges",
                 "substr_0", "substr_1"});

  for (auto const str : adj_vec) {
    os << str << "\n";
  }

} // function adj_list_out()



// **********************************************************************************
// Function gen_adj_list_vec()
// **********************************************************************************

std::vector<std::string> gen_adj_list_vec(sheet::DirectedAdjacencyList const& adj) {
  std::vector<std::string> adj_vec;

  auto const n_sheets = adj.sheets.size();
  out::substr2str ss_writer(std::make_shared<sheet::DirectedAdjacencyList>(adj));

  unsigned counter = 0;
  for (unsigned sheet_idx = 0; sheet_idx < n_sheets; ++sheet_idx) {
    auto const& sheet = adj.sheets[sheet_idx];
    for (auto const& pair_key : sheet.substr_keys()) {
      auto const& data = adj.adj_sub().map().at(pair_key);
      auto const str = (boost::format("%3d  %13s %3d %3d %3d  %5s %5s")
                                % sheet_idx
                                % (data.direction ? "Parallel" : "Anti-Parallel")
                                % data.delta_1
                                % data.delta_2
                                % data.residue_pairs
                                % ss_writer(pair_key.sub0())
                                % ss_writer(pair_key.sub1())
                       ).str();
      adj_vec.push_back(str);
      ++counter;
    }
  }

  return adj_vec;
} // function gen_adj_list_vec()


} // namespace adj_out

