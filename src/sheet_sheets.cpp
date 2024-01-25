// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "sheet/adj_list_with_sub.h"
#include "sheet/find_cycle.h"
#include "sheet/sheets.h"

namespace sheet {




bool check_undirected(SubStrandsPairKeyVec const& this_key_vec) {
  std::unordered_set<SubStrandsPairKey, SubStrandsPairKeyHasher> key_set;
  key_set.insert(this_key_vec.cbegin(), this_key_vec.cend());
  assert(key_set.size() != 0);

  bool erase_occured = false;

  for (auto iter = key_set.cbegin(); iter != key_set.cend();) {
    auto const rev_key = iter->reverse();
    if (key_set.count(rev_key)) {
      key_set.erase(rev_key);
      iter = key_set.erase(iter);

      erase_occured = true;

      continue;
    }
    ++iter;
  }

  if (key_set.size() == 0) {
    assert(erase_occured);
  }

  return erase_occured;
}




void Sheet::add(SubStrandsPairKey const& pair, StrandsPairAttribute const& attr) {
  member_substr.insert(pair.sub0());
  member_substr.insert(pair.sub1());
  member_substr.insert(attr.jumped_substrs.cbegin(), attr.jumped_substrs.cend());
  add_jump(pair, attr.jump);
}



void Sheet::merge(Sheet const& other) {
  member_substr.insert(other.member_substr.cbegin(), other.member_substr.cend());
  add_jump(other.max_key, other.sheet_size-2);
}



void Sheet::finish(SubStrandsPairKeyVec const& all_key_vec) {
  auto const this_key_vec = keys(all_key_vec);
  is_undirected = check_undirected(this_key_vec);
  cycle_check(this_key_vec);

  // if undirected, with_cycle must be true
  assert(not is_undirected or with_cycle);
}



void Sheet::add_key_vec(SubStrandsPairKeyVec const& all_key_vec_after_directed) {
  for (auto const& key: all_key_vec_after_directed) {
    if (member_substr.count(key.sub0()) and member_substr.count(key.sub1())) {
      substr_pair_keys.push_back(key);
    }
  }
  std::sort(substr_pair_keys.begin(), substr_pair_keys.end());
}



SubStrandSet Sheet::subs_in_cycles() const {
  SubStrandSet set;
  for (auto const& one_cycle : cycles_vec) {
    set.insert(one_cycle.cbegin(), one_cycle.cend());
  }
  return set;
}




void Sheet::add_jump(SubStrandsPairKey const& pair, IndexType const jump) {
  if (sheet_size < jump + 2u) {
    sheet_size = jump + 2u;
    max_key = pair;
  }
}



SubStrandsPairKeyVec Sheet::keys(SubStrandsPairKeyVec const& all_key_vec) const {
  SubStrandsPairKeyVec this_key_vec;
  for (auto const& one_key : all_key_vec) {
    if (member_substr.count(one_key.sub0()) and member_substr.count(one_key.sub1())) {
      this_key_vec.push_back(one_key);
    }
  }
  return this_key_vec;
}



/// Check if there is any cycles in this sheet.
void Sheet::cycle_check(SubStrandsPairKeyVec const& this_key_vec) {
  // Get all cycles
  auto const tmp_cycles = sheet::FindCycle{this_key_vec}.cycles;
  cycles_vec.reserve(tmp_cycles.size());


  // Push all cycles into member varible with rotating to the smallest element.
  for (auto const& one_cycle : tmp_cycles) {
    cycles_vec.push_back(rotate_to_smallest(one_cycle.cbegin(), one_cycle.cend()));
  }
  // Find the max cycle.
  auto const itr = std::max_element(cycles_vec.cbegin(), cycles_vec.cend(),
                                    [](auto const& a, auto const& b) {
                                        return a.size() < b.size(); });
  if (itr != cycles_vec.cend()) {
    sheet_size = itr->size();
    with_cycle = true;
  }
}




} // namespace sheet

