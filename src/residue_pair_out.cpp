// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include "functions.h"
#include "table.h"

namespace rpo {

int get_resnum(pdb::SSES const& sses, std::size_t const serial_str_id, std::size_t const serial_res_id) {
  auto const& strand_indices_vec = sses.gen_index_vec('E');
  auto const sse_id = strand_indices_vec[serial_str_id];
  return sses[sse_id].init + serial_res_id;
}


void residue_pair_out(table::TBLResiduePair & tbl,
                      sheet::DirectedAdjacencyList const& adj) {
  std::array<std::string, 5> type = {{"NoBridge", "ParallelHbonds", "ParallelNoHBonds", "SmallRing", "LargeRing"}};

  std::vector<table::ResiduePair> residue_pairs;

  unsigned sse_count = 0;
  for (auto const& strict_strand : adj.get_strict_zone().strict) {
    unsigned residue_count = 0;
    for (auto const& zi : strict_strand) {
      if (not zi.colored) {
        ++residue_count;
        continue;
      }

      auto const res0 = get_resnum(adj.sses, sse_count, residue_count);

      for (int i = 0; i < 2; ++i) {
        if (zi.adj_set[i]) {
          assert(zi.side != sheet::ZoneInfo::SideStatus::Undefined);
          assert(zi.bridge_type[i] != sheet::ZoneInfo::BridgeType::NoBridge);

          auto const res1 =
            get_resnum(adj.sses,
                       zi.adj_residues[i].serial_str_id,
                       zi.adj_residues[i].serial_res_id);

          auto const lambda_switch_PorA = [&res0, &res1](sheet::ZoneInfo::BridgeType const t){
            switch (t) {
              case sheet::ZoneInfo::BridgeType::ParallelNoHbonds:
                  return std::make_tuple("para", "A");
              case sheet::ZoneInfo::BridgeType::ParallelHbonds:
                  return std::make_tuple("para", "B");
              case sheet::ZoneInfo::BridgeType::LargeRing:
                  return std::make_tuple("anti", "Non-H-bonded");
              case sheet::ZoneInfo::BridgeType::SmallRing:
                  return std::make_tuple("anti", "H-bonded");
              default:
                  throw invalid_zone_info_exception{"No Bridge found between residue[" +
                    std::to_string(res0) +
                    "] and residue[" +
                    std::to_string(res1) +
                    "]"};
            }
          };

          auto const PorA_and_type = lambda_switch_PorA(zi.bridge_type[i]);

          residue_pairs.push_back(std::make_tuple(
                res0, res1,
                std::get<0>(PorA_and_type),
                std::get<1>(PorA_and_type),
                zi.side  ? "Front" : "Back"));
        }
      }

      ++residue_count;
    }
    ++sse_count;
  }

  std::sort(residue_pairs.begin(), residue_pairs.end(),
            [](auto const& a, auto const& b){
              if (std::get<0>(a) < std::get<0>(b)) {
                return true;
              } else if (std::get<0>(a) == std::get<0>(b)) {
                if (std::get<1>(a) < std::get<1>(b)) {
                  return true;
                }
              }
              return false;
            });

  for (auto const& rp : residue_pairs) {
    tbl.add(rp);
  }
}

} // namespace rpo

