// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_ADJACENT_SUBSTRAND_H_
#define SHEET_ADJACENT_SUBSTRAND_H_

#include "sheet/adj_list_with_sub.h"
#include "sheet/sheets.h"
#include "sheet/sub_strands_range.h"

namespace sheet {

using SubStrandVectors = std::array<std::vector<SubStrand>, 2>;

class AdjacentSubStrands {
public:

  /// @brief  Add new relation between key and ss.
  /// @param  key   A base Sub-Strand.
  /// @param  ss    A Sub-Strand that is adjacent to key.
  /// @param  side  True if ss is on right of key, false otherwise.
  void add(SubStrand const& key, SubStrand const& ss, bool const side, bool const rel_dir) {
    // if side == true, ss is on the right of key (key is on the left of ss)
    // if Parallel, side_r is same as side, but if Anti-Parallel, side_r is inverted.
    bool const side_r = rel_dir ? !side : side;
    add_helper(key, ss, side);
    add_helper(ss, key, side_r);
  }


  /// Given adj_list_with_sub and sheets (with cycle information)
  /// remove some unneeded paths and make the undirected graph partialy directed.
  void fix_undirected_paths(AdjListWithSub & adj_sub, Sheets const& sheets) const {
    for (auto const& sheet : sheets) {
      if (sheet.undirected()) {
        remove_adj_paths(adj_sub, sheet);
      }
    }
  }



  // ******************************************************************
  // Accessor Member Functions
  // ******************************************************************

  auto const& at(SubStrand const& key, bool const right) const {
    return data.at(key)[right];
  }

  /// Access the data directory.
  auto const& d() const { return data; }


protected:

  /// Helper Function for add(). Actually inserts a new data.
  void add_helper(SubStrand const& key, SubStrand const& ss, IndexType const side) {
    if (data.count(key) == 0) {
      data.insert({key, SubStrandVectors{}});
    }

    auto & v = data.at(key)[side];
    if (std::find(v.cbegin(), v.cend(), ss) == v.cend()) {
      v.push_back(ss);
    }

    #if defined(DEBUG) && defined(LOG_SHEET)
    pdb::debug_log("add_helper: key = " + key.string() + ", value = " + ss.string() + ", index = " + (side ? "1" : "0"));
    #endif

  }



  void remove_adj_paths(AdjListWithSub & adj_sub, Sheet const& sheet) const {
    constexpr std::array<bool, 2> const sides = {{true, false}};
    auto const subs_in_cyels = sheet.subs_in_cycles();
    for (auto const& sub_in_cycle : subs_in_cyels) {
      for (bool const side : sides) {
        for (auto const& sub : data.at(sub_in_cycle)[side]) {
          if (subs_in_cyels.count(sub) == 0) {
            SubStrandSet visited;
            visited.insert(sub_in_cycle);
            remove_one_path(visited, subs_in_cyels, adj_sub, sub_in_cycle, sub, side);
          }
        }
      }
    }
  }



  /// Actually remove the unneeded path.
  /// @param  side  The side of next related to sub.
  void remove_one_path(SubStrandSet & visited, SubStrandSet const& cycle_subs,
                       AdjListWithSub & adj_sub, SubStrand const& start,
                       SubStrand const& next, bool const side) const {
    auto const remove_key = side ? SubStrandsPairKey{next, start} :
                                   SubStrandsPairKey{start, next};
    adj_sub.remove(remove_key);
    visited.insert(next);

    // The side of next_next should be the opposite to start.
    bool const next_next_side = ! side_of(start, next);
    for (auto const& next_next : data.at(next)[next_next_side]) {
      if (visited.count(next_next) == 0 and cycle_subs.count(next_next) == 0) {
        remove_one_path(visited, cycle_subs, adj_sub, next, next_next, side);
      }
    }
  }


  /// Get on which side relative to base Sub-Strand, adj is located.
  bool side_of(SubStrand const& adj, SubStrand const& base) const {
    constexpr std::array<bool, 2> const sides = {{true, false}};
    for (bool const side : sides) {
      if(std::find(data.at(base)[side].cbegin(), data.at(base)[side].cend(), adj)
         != data.at(base)[side].cend()) {
        return side;
      }
    }
    throw AdjacentSubStrandNotFound{"Base = " + base.string() + ", adj = " + adj.string()};
  }


  /// data.at(i_substr): An array of 2 (1 for each h-bonded side of a Sub-Strand) elements,
  ///                 which store the list of adjacent Sub-Strands.
  /// Sub-Strands in data.at(i_substr)[0] and data.at(i_substr)[1] are on the other side.
  std::unordered_map<SubStrand, SubStrandVectors, SubStrandHasher> data{};
};


  
} // namespace sheet {

#endif // ifndef SHEET_ADJACENT_SUBSTRAND_H_
