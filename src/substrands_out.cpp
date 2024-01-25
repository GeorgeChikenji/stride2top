// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/format.hpp>
#include <boost/range/irange.hpp>

#include "functions.h"
#include "substrands.h"
#include "table.h"

#include "sheet/directed_adjacency_list.h"
#include "bab/filter.h"

namespace substrands {

// *******************************************************************************
// Function substrands_out()
// *******************************************************************************

void substrands_out(table::Table<table::SubStrand> & tbl,
                    sheet::DirectedAdjacencyList const& adj,
                    SubStrandStr2SheetIdxMap const& sheet_id_map) {
  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};

  for (auto const& sub : adj.substrs().vec()) {
    auto const substrand_id = ss_writer(sub);
    tbl.add(std::make_tuple(substrand_id, sheet_id_map.at(substrand_id),
                            adj.substrs().n_term_res(sub),
                            adj.substrs().c_term_res(sub)));
  }

} // function substrands_out()



// *******************************************************************************
// Function gen_sheet_id_map()
// *******************************************************************************

SubStrandStr2SheetIdxMap gen_sheet_id_map(sheet::DirectedAdjacencyList const& adj) {
  SubStrandStr2SheetIdxMap map;

  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};
  for (std::size_t sheet_id = 0; sheet_id < adj.sheets.size(); ++sheet_id) {
    for (auto const& ss : adj.sheets[sheet_id].member()) {
      map.insert({ss_writer(ss), sheet_id});
    }
  }
  return map;
} // function gen_sheet_id_map()


// *******************************************************************************
// Function helices_out()
// *******************************************************************************

void helices_out(table::Table<table::Helix> & tbl,
                 sheet::DirectedAdjacencyList const& adj) {
  for (std::size_t i = 0; i < adj.sses.size; ++i) {
    if (adj.sses[i].type == 'H') {
      tbl.add(std::make_tuple(i, adj.sses[i].init, adj.sses[i].end));
    }
  }
} // function helices_out()



// *******************************************************************************
// Function substrands_pair_out()
// *******************************************************************************

void substrands_pair_out(table::TBLSubStrandsPair & tbl,
                         sheet::DirectedAdjacencyList const& adj,
                         SubStrandStr2SheetIdxMap const& sheet_id_map,
                         bab::BabFilter & bab) {

  auto const N_SUBSTR = adj.substrs().vec().size();
  auto const ss_b = adj.substrs().vec().cbegin();

  std::vector<pdb::IndexType> pseudo_seq(adj.sses.size);
  std::iota(pseudo_seq.begin(), pseudo_seq.end(), 0);
  auto const pseudo_b = pseudo_seq.cbegin();

  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};

  for (auto const i : pdb::range(N_SUBSTR)) {
    for (auto const j : boost::irange(i + 1, N_SUBSTR)) {
      auto const& ss0 = *(ss_b + i);
      auto const& ss1 = *(ss_b + j);

      // sses_lbts and numres_lbts can be calculated for all SubStrand pairs.
      auto const sses_lbts = check_connection_type(ss0, ss1, adj, sheet_id_map, ss_writer);
      auto const numres_lbts = static_cast<std::size_t>(adj.substrs().n_term_res(ss1)-
                                                        adj.substrs().c_term_res(ss0)-1);


      sheet::SubStrandsPairKey const seq_key{ss0, ss1};
      auto const  rev_key = seq_key.reverse();

      auto const ss0_str = ss_writer(ss0);
      auto const ss1_str = ss_writer(ss1);

      // if not on the same sheet
      if (sheet_id_map.at(ss0_str) != sheet_id_map.at(ss1_str)) {
        tbl.add(std::make_tuple(ss0_str, ss1_str, "other", "", "",
                                0, -1, -1, 0, -1.0, "", 0));
        continue;
      }

      // Attributes for the path between ss_0 and ss_1
      auto const& seq_attr = adj.adj_attr.at(seq_key);
      auto const& rev_attr = adj.adj_attr.at(rev_key);

      // unreachable
      if (not seq_attr.reachable and not rev_attr.reachable) {
        tbl.add(std::make_tuple(ss0_str, ss1_str, "same", "?", "????",
                                100, -1, -1, 0, -1.0, sses_lbts, numres_lbts));
        continue;
      }

      bool const undirected = adj.sheets[sheet_id_map.at(ss0_str)].undirected();
      bool const in_same_cycle = in_cycle(ss0, ss1,
                                          adj.sheets[sheet_id_map.at(ss0_str)]);

      // If ss0 and ss1 are in the same cycle,
      // and 1 of 2 pathways in opposite directions are not found.
      if (in_same_cycle and not (seq_attr.reachable and rev_attr.reachable)) {
        if (not seq_attr.reachable and rev_attr.reachable) {
          throw one_directional_cycle_exception{ss0_str, ss1_str};
        } else if (seq_attr.reachable and not rev_attr.reachable) {
          throw one_directional_cycle_exception{ss1_str, ss0_str};
        }
      }

      auto const& attr = in_same_cycle ?
                            (seq_attr.jump < rev_attr.jump ? seq_attr: rev_attr) :
                            (seq_attr.reachable ? seq_attr : rev_attr);

      auto const& key = in_same_cycle ?
                            (seq_attr.jump < rev_attr.jump ? seq_key : rev_key) :
                            (seq_attr.reachable ? seq_key : rev_key);

      int d1 = -1;
      int d2 = -1;
      int br = 0;
      if (attr.jump == 0) {
        auto const& data = adj.adj_sub().map(key);
        d1 = data.delta_1;
        d2 = data.delta_2;
        br = data.residue_pairs;
        assert(data.direction == attr.direction);
      }

      std::string sheet = "same";
      std::string dir = seq_attr.reachable ? "-->" : "<--";

      // if reachable in both direction
      if (in_same_cycle and undirected) {
        sheet = "same_undir_cycle";
        dir = "?";
      } else if (in_same_cycle) {
        sheet = "same_in_cycle";
        dir = seq_attr.jump < rev_attr.jump ? "-->" : "<--";
      } else if (undirected) {
        sheet = "same_undirected";
      }

      std::string const PorA = attr.direction ? "para" : "anti";
      std::size_t const jump = attr.jump;


      // ********************************************
      // run bab_filter
      // ********************************************
      // If ss0 and ss1 are anti-parallel, calculate the score as if ss1 were reversed.
      unsigned const reversed = attr.direction ? 0 : (1u << adj.strand_indices[ss1.str]);

      // The fourth argument should be a past the end iterator.
      // (that's why '+1' is appended.)
      // The order of ss_0 and ss_1 must be sequential.
      // So the order of key is ignored and it's OK.
      bab(ss0, ss1,
          pseudo_b + adj.strand_indices[ss0.str],
          pseudo_b + adj.strand_indices[ss1.str] + 1, reversed);

      auto const& filter_result = bab.result();
      double score = 0.0;
      if (filter_result.success) {
        score = filter_result.left_score;
      } else {
        score = -1.0;
      }

      tbl.add(std::make_tuple(ss0_str, ss1_str, sheet, dir, PorA,
                              jump, d1, d2, br,
                              score, sses_lbts, numres_lbts));
    }
  }

} // function substrands_pair_out()



// *******************************************************************************
// Function in_cycle()
// *******************************************************************************

bool in_cycle(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1,
              sheet::Sheet const& sheet) {
  for (auto const& c : sheet.cycles()) {
    if (std::find(c.cbegin(), c.cend(), ss0) != c.cend() and
        std::find(c.cbegin(), c.cend(), ss1) != c.cend()) {
      return true;
    }
  }
  return false;
} // function in_cycle()



// *******************************************************************************
// Function check_connection_type()
// *******************************************************************************

std::string check_connection_type(sheet::SubStrand const& ss0,
                                  sheet::SubStrand const& ss1,
                                  sheet::DirectedAdjacencyList const& adj,
                                  SubStrandStr2SheetIdxMap const& sheet_id_map,
                                  out::substr2str const& ss_writer) {
  auto const& substrs = adj.substrs().vec();
  auto itr0 = std::lower_bound(substrs.cbegin(), substrs.cend(), ss0);
  auto itr1 = std::lower_bound(substrs.cbegin(), substrs.cend(), ss1);

  assert(itr0 != substrs.cend());
  assert(itr1 != substrs.cend());

  std::size_t ctype = 2u;
  auto const last_sse_id = adj.strand_indices[ss1.str];
  for (auto sse_id = adj.strand_indices[ss0.str] + 1; sse_id < last_sse_id; ++sse_id) {
    if (adj.sses[sse_id].type == 'H') {
      ctype |= 1u;
      break;
    }
  }

  bool const other_sheet = check_middle_ss_sheet(itr0, itr1,
                                      [](auto const a, auto const b) {return a != b;},
                                                 sheet_id_map, ss_writer);
  bool const same_sheet = check_middle_ss_sheet(itr0, itr1,
                                      [](auto const a, auto const b) {return a == b;},
                                                 sheet_id_map, ss_writer);
  if (other_sheet) {
    ctype |= 4u;
  }

  return gen_sses_lbts(ctype, same_sheet);
} // function check_connection_type()



// *******************************************************************************
// Function check_middle_ss_sheet()
// *******************************************************************************

bool check_middle_ss_sheet(std::vector<sheet::SubStrand>::const_iterator itr0,
                           std::vector<sheet::SubStrand>::const_iterator itr1,
                           std::function<bool(std::size_t, std::size_t)> cmp,
                           SubStrandStr2SheetIdxMap const& sheet_id_map,
                           out::substr2str const& ss_writer) {
  ++itr0;
  for (;itr0 != itr1; ++itr0) {
    if (cmp(sheet_id_map.at(ss_writer(*itr0)), sheet_id_map.at(ss_writer(*itr1)))) {
      return true;
    }
  }
  return false;
} // function check_middle_ss_sheet()



// *******************************************************************************
// Function substrands_pair_out()
// *******************************************************************************

std::string gen_sses_lbts(std::size_t const c_type, bool const on_same_sheet) {
  std::string type = "b-";

  // if loop only
  if (c_type == 2u and not on_same_sheet) {
    return "b-c-b";
  }

  // check Helix
  if (c_type & 1u) {
    type += "a";
  }

  // check Strand on the same sheet
  if (on_same_sheet) {
    type += "b";
  }

  // check Strand on other sheet
  if (c_type & 4u) {
    type += "b'";
  }

  return type + "-b";
} // function gen_sses_lbts()


} // namespace substrands
