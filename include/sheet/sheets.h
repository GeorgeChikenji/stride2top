// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_SHEETS_H_
#define SHEET_SHEETS_H_

#include <algorithm>
#include <vector>
#include <boost/range/irange.hpp>

#include "sheet/adj_list_with_sub.h"
#include "sheet/find_cycle.h"
#include "sheet/substr_pair_attr.h"

namespace sheet {


/// @brief  Check if the given vector of keys is undirected (by fallback mode) or not.
/// @param  this_key_vec  A vector of Keys that should be an output of keys().
/// @return True if this sheet is undirected, false otherwise.
bool check_undirected(SubStrandsPairKeyVec const& this_key_vec);



class Sheet {
public:

  // *******************************
  // Public Member Functions
  // *******************************

  Sheet() = default;

  /// Add the Sub-Strands in \c pair to \c member_substr if they don't exist in it.
  /// Update the longest jump +2 value if any longer value is found.
  void add(SubStrandsPairKey const& pair, StrandsPairAttribute const& attr);

  /// Merges other Sheet object into this Sheet.
  void merge(Sheet const& other);

  /// Run the final processes.
  void finish(SubStrandsPairKeyVec const& all_key_vec);

  /// Generates a set of Sub-Strands that are in any cycles.
  SubStrandSet subs_in_cycles() const;

  /// Copy the SubStrandsPairKeys that involves this sheets. (initialize substr_pair_keys)
  void add_key_vec(SubStrandsPairKeyVec const& all_key_vec_after_directed);


  // *******************************
  // Public Accessor Member Functions
  // *******************************

  auto const& member() const { return member_substr; }
  auto const& cycles() const { return cycles_vec; }
  auto size() const { return sheet_size; }
  bool has_cycle() const { return with_cycle; }
  bool undirected() const { return is_undirected; }
  auto const& substr_keys() const { return substr_pair_keys; };

protected:

  // *******************************
  // Protected Member Functions
  // *******************************

  void add_jump(SubStrandsPairKey const& pair, IndexType const jump);

  /// Pull in only the keys whose both Sub-Strands are in this Sheet.
  SubStrandsPairKeyVec keys(SubStrandsPairKeyVec const& all_key_vec) const;


  /// Check if there is any cycles in this sheet.
  void cycle_check(SubStrandsPairKeyVec const& this_key_vec);


  // *******************************
  // Protected Member Variables
  // *******************************

  /// A set of member Sub-Strands.
  SubStrandSet member_substr{};

  /// A vector of cycles that this sheet has.
  std::vector<std::vector<SubStrand>> cycles_vec{};

  /// {Max length between any 2 Sub-Strands in \x member_substr (Jump)} + 2.
  /// If this sheet has any cycles (When has_cycle is true),
  /// this value is the size of the largest cycle.
  std::size_t sheet_size{0};

  /// A pair of Sub-Strands that gives the max sheet_size.
  SubStrandsPairKey max_key{0, 0, 0, 0};

  /// True if this sheet has any cycle.
  bool with_cycle{false};

  bool is_undirected{false};

  /// Keys of the AdjListWithSub inside this sheet.
  SubStrandsPairKeyVec substr_pair_keys{};

};


/// Store information about the sheets. A collection of Sheet class objects.
class Sheets {
public:
  Sheets() = default;

  void add(SubStrandsPairKey const& key, StrandsPairAttribute const& attr) {
    auto sheet_iter = find_sheet(attr.jumped_substrs);
    if (sheet_iter == data.end()) {
      data.push_back(Sheet{});
      sheet_iter = data.end() - 1;
    }
    sheet_iter->add(key, attr);
    reconstruct();
  }

  void cycle_check(AdjListWithSubData const& adj_map) {
    // Prepare a vector of all keys
    SubStrandsPairKeyVec keys;
    keys.reserve(adj_map.size());
    for (auto const& pair : adj_map) {
      keys.push_back(pair.first);
    }

    for (auto & sheet : data) {
      sheet.finish(keys);
    }
  }


  /// Copy the sheet_key_vec from adj_list_with_sub after fix_undirected_paths.
  void add_key_vec(AdjListWithSub const& adj_sub) {
    SubStrandsPairKeyVec keys;
    keys.reserve(adj_sub.map().size());
    for (auto const& pair : adj_sub.map()) {
      keys.push_back(pair.first);
    }

    for (auto& sheet : data) {
      sheet.add_key_vec(keys);
    }
  }



  /// Sort the containts of data variable sequentially.
  void sort_sheets() {
    std::sort(data.begin(), data.end(),
        [](auto const& a, auto const& b) {
          auto const& a_min = *std::min_element(a.member().cbegin(), a.member().cend());
          auto const& b_min = *std::min_element(b.member().cbegin(), b.member().cend());
          return a_min < b_min; 
        });
  }


  auto const& operator[](std::size_t const pos) const { return data.at(pos); }
  auto size() const { return data.size(); }
  auto begin() const { return data.cbegin(); }
  auto end() const { return data.cend(); }
  auto cbegin() const { return data.cbegin(); }
  auto cend() const { return data.cend(); }


  using SheetVecConstIter = std::vector<Sheet>::const_iterator;

  /// @brief  const version of find_sheet() and argument is a 1 Sub-Strand
  SheetVecConstIter find_sheet(SubStrand const& sub) const {
    for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
      if (std::find(iter->member().cbegin(), iter->member().cend(), sub)
          != iter->member().cend()) {
        return iter;
      }
    }
    return data.cend();
  }


protected:

  using SheetVecIter = std::vector<Sheet>::iterator;

  /// @return Iterator to the first sheet that contains common Sub-Strand in \c substrs.
  ///         data.cend() if not found.
  SheetVecIter find_sheet(std::vector<SubStrand> const& substrs) {
    for (SheetVecIter iter = data.begin(); iter != data.end(); ++iter) {
      if (std::find_first_of(iter->member().cbegin(), iter->member().cend(),
                             substrs.cbegin(), substrs.cend()) != iter->member().cend()) {
        return iter;
      }
    }
    return data.end();
  }



  /// Find any pairs of Sheets in \c data that have common Sub-Strand.
  void reconstruct() {
    for (auto iter_i = data.begin(); iter_i != data.end();++iter_i) {
      for (auto iter_j = iter_i + 1; iter_j != data.end();) {
        // Merge 2 Sheets if a common Sub-Strand found.
        if (std::find_first_of(iter_i->member().cbegin(), iter_i->member().cend(),
                               iter_j->member().cbegin(), iter_j->member().cend())
            != iter_i->member().cend()) {
          iter_i->merge(*iter_j);
          iter_j = data.erase(iter_j);
          continue;
        }
        ++iter_j;
      }
    }
  }

  std::vector<Sheet> data{};
};


} // namespace sheet

#endif // ifndef SHEET_SHEETS_H_
