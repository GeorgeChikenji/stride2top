// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cassert>

#include "sheet/adj_list_with_sub.h"
#include "sheet/exceptions.h"

namespace sheet {


// ********************************************************************
// Public Member Function add_sheet()
// ********************************************************************

void AdjListWithSub::add_sheet() {
  sheets_key_vec.push_back(SubStrandsPairKeyVec{});
} // public member function add_sheet()



// ********************************************************************
// Public Member Function erase_too_short()
// ********************************************************************

bool AdjListWithSub::erase_too_short(SubStrandSet const& short_subs) {
  bool erased = false;
  auto & key_vec = sheets_key_vec.back();
  for (auto key_itr = key_vec.cbegin(); key_itr != key_vec.cend();) {
    if (short_subs.count(key_itr->sub0()) or short_subs.count(key_itr->sub1())) {
      erased = true;
      data.erase(*key_itr);
      key_itr = key_vec.erase(key_itr);
      continue;
    }
    ++key_itr;
  }
  return erased;
} // public member function erase_too_short()



// ********************************************************************
// Public Member Function ensure_undirected()
// ********************************************************************

void AdjListWithSub::ensure_undirected() {
  std::vector<SubStrandsPairKey> to_add;
  for (auto const key : sheets_key_vec.back()) {
    if (data.count(key.reverse()) == 0) {
      to_add.push_back(key);
    }
  }

  for (auto const& key : to_add) {
    data.insert({key.reverse(), data[key]});
    sheets_key_vec.back().push_back(key.reverse());
  }

} // public member funciton ensure_undirected()



// ********************************************************************
// Public Member Function update_key_substr()
// ********************************************************************

void AdjListWithSub::update_key_substr(std::unordered_map<SubStrand, SubStrand, SubStrandHasher> const& conv) {
  std::unordered_map<SubStrandsPairKey, SubStrandsPairNode,
                     SubStrandsPairKeyHasher> to_add;
  std::vector<SubStrandsPairKey> to_remove;
  for (auto const& key_value : data) {
    auto const& ss0 = key_value.first.sub0();
    auto const& ss1 = key_value.first.sub1();

    auto const& new_ss0 = conv.count(ss0) ? conv.at(ss0) : ss0;
    auto const& new_ss1 = conv.count(ss1) ? conv.at(ss1) : ss1;
    auto const new_key = SubStrandsPairKey{new_ss0, new_ss1};

    if (new_key == key_value.first) {
      continue;
    }
    to_add.insert({new_key, key_value.second});
    replace_erased_itr(key_value.first, new_key);
    to_remove.push_back(key_value.first);
  }

  for (auto const& key : to_remove) {
    data.erase(key);
  }
  for (auto const& key_value : to_add) {
    data.insert(key_value);
  }

} // public member funciton update_key_substr()




// ********************************************************************
// Public Member Function cleanup()
// ********************************************************************


void AdjListWithSub::cleanup(std::vector<IndexType> const& strand_indices) {
  for (IndexType i = 0; i < sheets_key_vec.size(); ++i) {
    if (not sheets_directed[i]) {
      continue;
    }

    for (auto key_itr = sheets_key_vec[i].cbegin(); key_itr != sheets_key_vec[i].cend();) {
      auto const rev_key = key_itr->reverse();

      // if the strands pair in reverse direcion was not found (had already processed)
      if (data.count(rev_key) == 0) {
        ++key_itr;
        continue;
      }

      auto const n_itr = data.at(*key_itr).residue_pairs;
      auto const n_rev = data.at(rev_key).residue_pairs;
      if (0.66 < std::min(n_itr, n_rev) / std::max(n_itr, n_rev)) {

        throw substrand_cleanup_failure{strand_indices[rev_key.str1], rev_key.substr1, n_itr,
                                        strand_indices[rev_key.str0], rev_key.substr0, n_rev};
      }

      if (n_itr < n_rev) {
        data.erase(*key_itr);
        key_itr = sheets_key_vec[i].erase(key_itr);
      } else {
        auto const d = std::find_if(sheets_key_vec[i].cbegin(), sheets_key_vec[i].cend(),
                                [&rev_key](auto const& a) { return a == rev_key; });

        // d must be found
        assert(d != sheets_key_vec[i].cend());

        // backup the value of key_itr becausee erase() may invalidate it.
        auto const key_back = *key_itr;

        // invalidates key_itr
        sheets_key_vec[i].erase(d);
        data.erase(rev_key);

        // find the iterator pointing to the original key value
        key_itr = std::find(sheets_key_vec[i].cbegin(), sheets_key_vec[i].cend(), key_back);

        // itr must be found
        assert(key_itr != sheets_key_vec[i].cend());
        ++key_itr;
      }
    } // for each edge

    // Erase the sheets that became empty after removing too short strands.
    if (sheets_key_vec[i].size() == 0) {
      sheets_key_vec.erase(sheets_key_vec.cbegin() + i);
      sheets_directed.erase(sheets_directed.cbegin() + i);
      --i;
    }

  } // for each sheet
} // public member function cleanup()



// ********************************************************************
// Public Member Function remove()
// ********************************************************************
void AdjListWithSub::remove(SubStrandsPairKey const& key) {
  // register as removed
  removed_keys.push_back(key);

  // erase from data
  data.erase(key);

  // erase from sheets_key_vec
  for (auto & key_vec : sheets_key_vec) {
    auto iter = std::find(key_vec.cbegin(), key_vec.cend(), key);
    if (iter != key_vec.cend()) {
      key_vec.erase(iter);
      break;
    }
  }

  // erase from adj_sub_vec
  auto & adj_vec = adj_sub_vec.at(key.sub0());
  auto iter = std::find(adj_vec.cbegin(), adj_vec.cend(), key);
  if (iter != adj_vec.cend()) {
    adj_vec.erase(iter);
  }

} // public member function remove


// ********************************************************************
// Public Member Function gen_adj_sub_vec()
// ********************************************************************

void AdjListWithSub::gen_adj_sub_vec() {
  for (auto const& pair : data) {
    SubStrand const ss0{pair.first.str0, pair.first.substr0};
    SubStrand const ss1{pair.first.str1, pair.first.substr1};

    if (adj_sub_vec.count(ss0) == 0) {
      adj_sub_vec.insert({ss0, std::vector<SubStrandsPairKey>()});
    }
    adj_sub_vec[ss0].push_back(pair.first);
  }

} // public member function gen_adj_sub_vec




// ********************************************************************
// Public Member Function update_delta()
// ********************************************************************

void AdjListWithSub::update_delta(SubStrandsPairKey const& key,
                                  int const delta_1, int const delta_2) {
  data.at(key).delta_1 = delta_1;
  data.at(key).delta_2 = delta_2;
} // public member function update_delta



// ********************************************************************
// Public Member Function replace_erased_iter()
// ********************************************************************
void AdjListWithSub::replace_erased_itr(SubStrandsPairKey const& old_key,
                                        SubStrandsPairKey const& new_key) {

  for (auto & vec : sheets_key_vec) {
    auto const d_iter = std::find(vec.begin(), vec.end(), old_key);

    if (d_iter == vec.end()) {
      continue;
    }
    vec.erase(d_iter);
    vec.push_back(new_key);
    return;
  }

} // protected member function replace_erased_iter()



} // namespace sheet

