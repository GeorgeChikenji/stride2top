// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/format.hpp>
#include <boost/range/irange.hpp>

#include "functions.h"
#include "sheets_out.h"
#include "table.h"
#include "pdb/sses.h"
#include "pdb/tools.h"
#include "sheet/sheets.h"
#include "sheet/directed_adjacency_list.h"

namespace sheets_out {

// **********************************************************************************
// Function print_sheet()
// **********************************************************************************

void print_sheet(table::Table<table::Sheet> & tbl,
                 sheet::DirectedAdjacencyList const& adj) {
  // if no sheet
  if (adj.sheets.size() == 0) {
    return;
  }

  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};

  std::size_t sheet_id = 0;
  for (auto const& sheet : adj.sheets) {
    bool const with_branch = sheet.size() != sheet.member().size();

    bool const consec_beta = is_all_consec(sheet, adj);;
    topology_string const topo{sheet, adj};
    std::string const topo_str_richardson = topo.str();
    std::string const topo_str_cohen = topo.str(2);

    bool all_p = false, all_ap = false;
    std::tie(all_p, all_ap) = check_all_pap(sheet, adj);

    auto const seq_ss = sort_sheet_members(sheet);
    auto const member_substrs = "'" +
                                out::join(seq_ss.cbegin(), seq_ss.cend(), ",",
                                       [&ss_writer](auto a){return ss_writer(*a);}) +
                                "'";

    tbl.add(std::make_tuple(sheet_id, sheet.member().size(), sheet.cycles().size(),
                            (sheet.undirected() ? 'T' : 'F'),
                            (with_branch ? 'T' : 'F'),
                            (consec_beta ? 'T' : 'F'),
                            (all_p ? 'T' : 'F'),
                            (all_ap ? 'T' : 'F'),
                            member_substrs,
                            topo_str_richardson, topo_str_cohen
                            ));
    ++sheet_id;
  }
} // function print_sheet()



// **********************************************************************************
// Function sort_sheet_members()
// **********************************************************************************

std::vector<sheet::SubStrand> sort_sheet_members(sheet::Sheet const& sheet) {
  std::vector<sheet::SubStrand> sorted_ss;
  sorted_ss.reserve(sheet.member().size());
  sorted_ss.insert(sorted_ss.end(), sheet.member().cbegin(), sheet.member().cend());
  std::sort(sorted_ss.begin(), sorted_ss.end());
  return sorted_ss;
} // function sort_sheet_members()




// **********************************************************************************
// Function is_all_consec()
// **********************************************************************************

bool is_all_consec(sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj) {
  // Generate a sorted vector of Sub-Strands in this sheet.
  auto const seq_ss = sort_sheet_members(sheet);

  bool consec_beta = true;

  // The number of Sub-Strand Pairs.
  assert(seq_ss.size() != 0);
  auto const n_ss_pairs = seq_ss.size() - 1;

  auto const& substrs = adj.substrs().vec();

  // for all Sub-Strand pairs in this sheet.
  for (std::size_t i = 0; i < n_ss_pairs; ++i) {
    auto const iter_0 = std::lower_bound(substrs.cbegin(), substrs.cend(), seq_ss[i]);
    auto const iter_1 = std::lower_bound(substrs.cbegin(), substrs.cend(), seq_ss[i+1]);

    // If you follow the chain sequentially, it goes to other sheet and then
    // return to this sheet afterwards.
    if (std::distance(iter_0, iter_1) != 1) {
      consec_beta = false;
      break;
    }
    auto const sse_id_diff = adj.strand_indices[seq_ss[i+1].str] -
                             adj.strand_indices[seq_ss[i].str];
    if (sse_id_diff != 0 and sse_id_diff != 1) {
      consec_beta = false;
      break;
    }
  }

  return consec_beta;
} // function is_all_consec()





// **********************************************************************************
// Function check_all_pap()
// **********************************************************************************

std::tuple<bool, bool> check_all_pap(sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj) {
  bool all_p = true;
  bool all_ap = true;

  auto const& substrs = sheet.member();
  auto const n_substrs = substrs.size();

  // inefficient !!!
  std::vector<sheet::SubStrand> substr_vec;
  substr_vec.reserve(n_substrs);
  substr_vec.insert(substr_vec.begin(), substrs.cbegin(), substrs.cend());

  for (auto const i : pdb::range(n_substrs)) {
    for (auto const j : boost::irange(i+1, n_substrs)) {
      auto const& attr = adj.attr(substr_vec[i], substr_vec[j]);
      if (attr.jump != 0) {
        continue;
      }

      // if this is parallel
      if (attr.direction) {
        all_ap = false;

      // if this is anti-parallel
      } else {
        all_p = false;
      }
    }
  }
  return std::make_tuple(all_p, all_ap);
} // function check_all_pap



// **********************************************************************************
// Function extracted_adjacent_substr_out()
// **********************************************************************************

void extracted_adjacent_substr_out(table::Table<table::ExtractedSheet> & tbl,
                                   unsigned const n,
                                   sheet::DirectedAdjacencyList const& adj) {

  // Run extraction
  std::vector<std::tuple<std::vector<sheet::SubStrand>, unsigned, bool>> extracted;

  for (unsigned sheet_idx = 0; sheet_idx < adj.sheets.size(); ++sheet_idx) {

    auto const tmp_ext = extract_adjacent_substr(n, adj.sheets[sheet_idx], adj);
    for (auto const& each_tmp_ext : tmp_ext) {
      extracted.push_back(std::make_tuple(each_tmp_ext, sheet_idx,
                          n == adj.sheets[sheet_idx].member().size()));
    }
  }
  
  // if no extract
  if (extracted.size() == 0) {
    return;
  }
  out::substr2str ss_writer{std::make_shared<sheet::DirectedAdjacencyList>(adj)};

  for (unsigned i = 0; i < extracted.size(); ++i) {
    auto const substr_vec = sort_substr_vec(std::get<0>(extracted[i]));
    assert(substr_vec.size() == n);

    auto const sheet_idx = std::get<1>(extracted[i]);
    bool const with_cycle = cycle_checker(substr_vec, adj);
    topology_string const topo{substr_vec, with_cycle, adj};
    std::string const topo_str = topo.str(2);

    auto const member = "'" + out::join(substr_vec.cbegin(), substr_vec.cend(), ",",
                                       [&ss_writer](auto a){return ss_writer(*a);}) +
                        "'";

    tbl.add(std::make_tuple(sheet_idx, n, (std::get<2>(extracted[i]) ? 'T' : 'F'),
                            member, topo_str));
  }

} // function extracted_adjacent_substr_out()



// **********************************************************************************
// Function cycle_checker()
// **********************************************************************************

bool cycle_checker(std::vector<sheet::SubStrand> const& ss_vec,
                   sheet::DirectedAdjacencyList const& adj) {

  std::vector<sheet::SubStrandsPairKey> key_vec;
  assert(ss_vec.size() != 0);
  key_vec.reserve(ss_vec.size() - 1);

  // Get involved pair keys
  for (auto const& key_value: adj.adj_sub().map()) {
    auto const& pair_key = key_value.first;
    auto const ss0 = pair_key.sub0();
    auto const ss1 = pair_key.sub1();
    if (std::find(ss_vec.cbegin(), ss_vec.cend(), ss0) == ss_vec.cend() or
        std::find(ss_vec.cbegin(), ss_vec.cend(), ss1) == ss_vec.cend()) {
      continue;
    }
    key_vec.push_back(pair_key);
  }

  sheet::FindCycle cycles{key_vec};
  return cycles.cycles.size() != 0;
} // function cycle_checker()



// **********************************************************************************
// Function extract_adjacent_substr()
// **********************************************************************************

sheet::SubStrandsPairKeyVec adj_sub_vec2pair_key_vec(sheet::AdjSubVec const& adj_sub_vec) {
  sheet::SubStrandsPairKeyVec pair_key_vec;
  // concatenate vectors
  for (auto const& key_value : adj_sub_vec) {
    auto const& vec = key_value.second;
    pair_key_vec.insert(pair_key_vec.cend(), vec.cbegin(), vec.cend());
  }
  return pair_key_vec;
} // function cycle_checker()



// **********************************************************************************
// Function extract_adjacent_substr()
// **********************************************************************************

std::vector<std::vector<sheet::SubStrand>> extract_adjacent_substr(unsigned const n, sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj) {
  std::vector<std::vector<sheet::SubStrand>> extracted_sheets;

  // for all member SubStrands
  for (auto const& ss : sheet.member()) {
    auto const new_sheets = extract_from_one_substr(ss, n, sheet, adj);
    extracted_sheets.insert(extracted_sheets.cend(), new_sheets.cbegin(), new_sheets.cend());
  }
  return extracted_sheets;
} // function extract_adjacent_substr()



// **********************************************************************************
// Function extract_from_one_substr()
// **********************************************************************************

std::vector<std::vector<sheet::SubStrand>> extract_from_one_substr(sheet::SubStrand const& start_ss, unsigned const n, sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj) {
  std::vector<std::vector<sheet::SubStrand>> new_extract;
  recursive_extract(std::vector<sheet::SubStrand>{start_ss}, new_extract, sheet, n, adj);

  return new_extract;
} // function extract_from_one_substr()



// **********************************************************************************
// Function recursive_extract()
// **********************************************************************************

void recursive_extract(std::vector<sheet::SubStrand> const& current_path,
                       std::vector<std::vector<sheet::SubStrand>> & found_path,
                       sheet::Sheet const& sheet,
                       unsigned const n, sheet::DirectedAdjacencyList const& adj) {
  // if current_path has specified number of SubStrands
  if (current_path.size() == n) {
    auto const tmp = adj.adj_sub().substr_vec2adj_sub_vec(current_path.cbegin(),
                                                          current_path.cend());
    // if not undirected (if directed)
    if (not sheet::check_undirected(adj_sub_vec2pair_key_vec(tmp))) {
      found_path.push_back(current_path);
    }
    return;
  }

  // for each SubStrands pairs in this sheet 
  for (auto const& key : sheet.substr_keys()) {
    auto const& ss0 = key.sub0();
    auto const& ss1 = key.sub1();

    // if the first SubStrand differs from the target one
    if (ss0 != current_path.back()) {
      continue;
    }

    // if ss1 is already in current_path
    if (std::find(current_path.cbegin(), current_path.cend(), ss1)
        != current_path.cend()) {
      continue;
    }

    std::vector<sheet::SubStrand> new_path = current_path;
    new_path.push_back(ss1);
    recursive_extract(new_path, found_path, sheet, n, adj);
  }
} // function recursive_extract()



// **********************************************************************************
// Function sort_substr_vec()
// **********************************************************************************

std::vector<sheet::SubStrand> sort_substr_vec(std::vector<sheet::SubStrand> const& ss_vec) {
  auto seq_ss = ss_vec;
  std::sort(seq_ss.begin(), seq_ss.end());
  return seq_ss;
} // function sort_substr_vec()



// **********************************************************************************
// Public Member Function str()
// **********************************************************************************

std::string topology_string::str(unsigned const style) const {
  switch (style) {
    // Richardson
    case 0:
      return pair_style_str();
    case 1:
      return position_style_str(ss_position_style, true);
    // Cohen
    case 2:
      return position_style_2_str(ss_position_style);
    default:
      throw std::runtime_error{"Unkown style specified in topology_string::str()"};
  }
} // public member function str()



// **********************************************************************************
// Private Member Function init_pair_style()
// **********************************************************************************

std::vector<topology_string::ss_pair_arrangement> topology_string::init_pair_style(std::vector<sheet::SubStrand> const& ss_vec, sheet::AdjSubVec const& adj_sub_vec, bool const with_cycle) const {

  // Sort input Sub-Strand vector
  auto const seq_ss = sort_substr_vec(ss_vec);


  // The number of Sub-Strand Pairs.
  assert(seq_ss.size() != 0);
  auto const n_ss_pairs = seq_ss.size() - 1;

  // A vector to be returned
  std::vector<ss_pair_arrangement> ret{};

  // With cycle
  if (with_cycle) {
    // define the direction for topology string inside the sheet
    std::array<sheet::StrandsPairAttribute, 2> const first_pair = {{
      adj.search(seq_ss[0], seq_ss[1], adj_sub_vec),
      adj.search(seq_ss[1], seq_ss[0], adj_sub_vec)
    }};

    assert(first_pair[0].reachable);
    assert(first_pair[1].reachable);

    unsigned const idx = first_pair[0].jump < first_pair[1].jump ? 0 : 1;

    for (std::size_t i = 0; i < n_ss_pairs; ++i) {
      auto const& attr = idx == 0 ?
        adj.search(seq_ss[i], seq_ss[i+1], adj_sub_vec):
        adj.search(seq_ss[i+1], seq_ss[i], adj_sub_vec);

      assert(attr.reachable);
      ret.push_back(ss_pair_arrangement{static_cast<int>(attr.jump + 1),
                                        attr.direction});
    }


  // no cycle
  } else {
    for (std::size_t i = 0; i < n_ss_pairs; ++i) {
      auto const& attr_f = adj.search(seq_ss[i], seq_ss[i+1], adj_sub_vec);
      auto const& attr_r = adj.search(seq_ss[i+1], seq_ss[i], adj_sub_vec);

      #ifndef NDEBUG
      if (attr_f.reachable) {
        assert(not attr_r.reachable);
      } else {
        assert(attr_r.reachable);
      }
      #endif // ifndef NDEBUG

      auto const& attr = attr_f.reachable ? attr_f : attr_r;
      int const mult = attr_f.reachable ? 1 : -1;
      ret.push_back(ss_pair_arrangement{mult * static_cast<int>(attr.jump + 1),
                                        attr.direction});
    }
  }
  return modify_to_plus(ret);
} // private member function init_pair_style()



// **********************************************************************************
// Private Member Function init_pair_style()
// **********************************************************************************
std::vector<topology_string::ss_pair_arrangement> topology_string::init_pair_style(sheet::Sheet const& sheet) const {
  // If undirected and branches
  if (sheet.undirected() or sheet.size() != sheet.member().size()) {
    return std::vector<ss_pair_arrangement>{};
  } else {
  // Generate a sorted vector of Sub-Strands in this sheet.
    return init_pair_style(sort_sheet_members(sheet), adj.adj_sub().adj_substr_vec(),
                           sheet.has_cycle());
  }
} // private member function init_pair_style()



// **********************************************************************************
// Private Member Function modify_to_plus()
// **********************************************************************************
std::vector<topology_string::ss_pair_arrangement> topology_string::modify_to_plus(std::vector<topology_string::ss_pair_arrangement> const& orig) const {

  if (orig[0].to_next < 0) {

    std::vector<ss_pair_arrangement> ret;
    ret.reserve(orig.size());

    for (auto const& a : orig) {
      ret.push_back(ss_pair_arrangement{-1 * a.to_next, a.direction});
    }
    return ret;

  } else {
    return orig;
  }
} // private member function modify_to_plus()




// **********************************************************************************
// Private Member Function init_position_style()
// **********************************************************************************

std::vector<topology_string::ss_position> topology_string::init_position_style() const {
  // if NA
  if (pair_style.size() == 0) {
    return std::vector<ss_position>{};
  }

  // The number of the strand of the target sheet
  auto const n_str = pair_style.size() + 1;

  std::vector<ss_position> first;
  first.reserve(n_str);

  first.push_back(ss_position{0, true});
  std::size_t const n_pairs = pair_style.size();
  for (std::size_t i = 0; i < n_pairs; ++i) {
    // add i+1 th element (depends on i th elements of first and pair_style) of first.
    first.push_back(ss_position{first[i].pos + pair_style[i].to_next,
                                first[i].direction == pair_style[i].direction});
  }

  std::vector<ss_position> ret;
  ret.reserve(first.size());
  int const min_pos = std::min_element(first.cbegin(), first.cend(),
                                       [](auto const& a, auto const& b) {
                                         return a.pos < b.pos;
                                       })->pos;

  //  set the smallest to 1
  int const diff = - min_pos + 1;
  // sequential counter
  unsigned seq_c = 0;
  for (auto const& f : first) {
    ret.push_back(ss_position{seq_c, f.pos + diff, f.direction});
    ++seq_c;
  }


  // fllip horizontally if the position of the first strand is not in the left side half.
  if (static_cast<int>(n_str)/2 < ret[0].pos) {
    for (auto & strand : ret) {
      strand.pos = n_str - strand.pos + 1;
    }
  }

  return ret;
} // private member function init_position_style()




// **********************************************************************************
// Private Member Function pair_style_str()
// **********************************************************************************
std::string topology_string::pair_style_str() const {
  // if not generated
  if (pair_style.size() == 0) {
    return "NA";
  }

  std::string topo = "";
  for (auto const& pair : pair_style) {
    topo += (boost::format("%+d") % pair.to_next).str() + (pair.direction ? "x" : "");
  }
  return topo;
} // private member function pair_style_str()




// **********************************************************************************
// Private Member Function position_style_str()
// **********************************************************************************

std::string topology_string::position_style_str(std::vector<ss_position> const& pos_vec,
                                                bool const use_pos) const {
  // if not generated
  if (pos_vec.size() == 0) {
    return "NA";
  }

  std::string topo = "";
  for (auto const& ss_pos : pos_vec) {
    topo += (ss_pos.direction ? "+_" : "-_") +
            std::to_string((use_pos ? ss_pos.pos : ss_pos.seq_id)) + ",";
  }
  return topo;
} // private member function position_style_str()




// **********************************************************************************
// Private Member Function position_style_str()
// **********************************************************************************

std::string topology_string::position_style_2_str(std::vector<ss_position> const& pos_vec) const {
  auto ss_pos_cpy = pos_vec;

  // add sequential id
  for (std::size_t i = 0; i < ss_pos_cpy.size(); ++i) {
    ss_pos_cpy[i].seq_id = i + 1;
  }

  std::sort(ss_pos_cpy.begin(), ss_pos_cpy.end(),
            [] (auto const& a, auto const& b) {
              return a.pos < b.pos;
            });
  return position_style_str(ss_pos_cpy);
} // private member function position_style_str()


} // namespace sheets_out
