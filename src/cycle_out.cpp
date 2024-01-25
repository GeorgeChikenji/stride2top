// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <numeric>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include "table.h"
#include "functions.h"
#include "table.h"

namespace cycles {

void output_cycles(table::Table<table::Cycle> & tbl,
                   sheet::DirectedAdjacencyList const& adj) {

  auto const cycles = gen_cycles_vec(adj);
  auto const n_cycle = cycles.size();

  // if no cycle
  if (n_cycle == 0) {
    return;
  }

  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};

  for (std::size_t i = 0; i < n_cycle; ++i) {
    std::size_t const orig_sheet_id = std::get<0>(cycles[i]);
    auto const& c_vec = std::get<1>(cycles[i]);

    // A comma separated list of member substrands
    auto const member_list = "'" +
                             out::join(c_vec.cbegin(), c_vec.cend(), ",",
                                       [&ss_writer](auto a){return ss_writer(*a);}) +
                             "'";

    tbl.add(std::make_tuple(orig_sheet_id, c_vec.size(), member_list));
  }
}


CyclesVec gen_cycles_vec(sheet::DirectedAdjacencyList const& adj) {

  CyclesVec cycles_vec;
  std::size_t sheet_id = 0;
  for (auto const sheet : adj.sheets) {
    sheet::FindCycle const cycle{sheet.substr_keys()};
    for (auto const& one_cycle : cycle.cycles) {
      cycles_vec.push_back(std::make_tuple(sheet_id, rotate_to_smallest(one_cycle.cbegin(),
                                                                        one_cycle.cend())));
    }
    ++sheet_id;
  }
  return cycles_vec;
}

} // namespace cycles
