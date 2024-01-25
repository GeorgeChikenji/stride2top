// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_ADJ_LIST_WITH_SUB_H_
#define SHEET_ADJ_LIST_WITH_SUB_H_

#ifdef DEBUG
#include <boost/format.hpp>
#endif

#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "pdb/sses.h"
#include "sheet/common.h"

namespace sheet {

struct SubStrand {
  SubStrand() = default;
  SubStrand(IndexType const s, IndexType const ss) noexcept: str{s}, substr{ss} {}

  std::string string() const {
    return std::to_string(str) + "_" + std::to_string(substr);
  }
  std::string string(std::vector<IndexType> const& strand_indices) const {
    return "SSE[" + std::to_string(strand_indices[str]) + "]-" + std::to_string(substr);
  }


  bool operator==(SubStrand const& other) const noexcept {
    return str == other.str and substr == other.substr;
  }
  bool operator!=(SubStrand const& other) const noexcept {
    return !(*this == other);
  }
  bool operator<(SubStrand const& other) const noexcept {
    return str < other.str or (str == other.str and substr < other.substr);
  }

  IndexType str{0};
  IndexType substr{0};
};


struct SubStrandHasher {
  std::size_t operator()(SubStrand const& k) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, boost::hash_value(k.str));
    boost::hash_combine(seed, boost::hash_value(k.substr));
    return seed;
  }
};



/// Data struct stored in AdjListWithSub
struct SubStrandsPairNode {

  SubStrandsPairNode() = default;
  SubStrandsPairNode(SubStrandsPairNode const&) = default;

  SubStrandsPairNode(bool const dir) noexcept:
    direction{dir} {}

  SubStrandsPairNode(bool const dir, unsigned const pair_count) noexcept:
    direction{dir}, residue_pairs{pair_count} {}


  /// true  Parallel
  /// false Anti-Parallel
  bool direction{false};

  /// The number of residues between the first residue of each strands.
  /// If the first residues of paired strands are adjacent (paired), this value equals to 0.
  /// The differences of the number of residues
  int delta_1{0};
  int delta_2{1};

  /// The number of residue pairs between the 2 strands.
  unsigned residue_pairs{1};

};




/// Key struct to access the stored data in AdjListWithSub
struct SubStrandsPairKey {
  SubStrandsPairKey(IndexType const s0, IndexType const ss0,
                    IndexType const s1, IndexType const ss1) noexcept:
    str0{s0}, substr0{ss0}, str1{s1}, substr1{ss1} {}

  SubStrandsPairKey(SubStrand const& sub0, SubStrand const& sub1) noexcept:
    str0{sub0.str}, substr0{sub0.substr}, str1{sub1.str}, substr1{sub1.substr} {}

  SubStrand sub0() const noexcept { return SubStrand{str0, substr0}; }
  SubStrand sub1() const noexcept { return SubStrand{str1, substr1}; }

  bool operator==(SubStrandsPairKey const& other) const noexcept {
    return str0 == other.str0 and substr0 == other.substr0 and
           str1 == other.str1 and substr1 == other.substr1;
  }

  bool operator<(SubStrandsPairKey const& other) const noexcept {
    if (str0 < other.str0) {
      return true;
    } else if (str0 == other.str0) {
      if (substr0 < other.substr0) {
        return true;
      } else if (substr0 == other.substr0) {
        if (str1 < other.str1) {
          return true;
        } else if (str1 == other.str1) {
          return substr1 < other.substr1;
        }
      }
    }
    return false;
  }

  SubStrandsPairKey reverse() const noexcept {
    return SubStrandsPairKey{str1, substr1, str0, substr0};
  }

  #ifdef DEBUG
  /// @brief Return the string representation of this key for debugging
  std::string str(std::vector<IndexType> const& strand_indices) const {
    return ((boost::format("SSE[%2d]_%2d -> SSE[%2d]_%2d")
                           % +strand_indices[str0]
                           % +substr0
                           % +strand_indices[str1]
                           % +substr1).str());
  }
  std::string str() const {
    return ((boost::format("%2d_%2d -> %2d_%2d")
                           % +str0
                           % +substr0
                           % +str1
                           % +substr1).str());
  }
  #endif // ifdef DEBUG



  IndexType str0{0};
  IndexType substr0{0};

  IndexType str1{0};
  IndexType substr1{0};
};



/// Hasher for SubStrandsPairKey
struct SubStrandsPairKeyHasher {
  std::size_t operator()(SubStrandsPairKey const& k) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, boost::hash_value(k.str0));
    boost::hash_combine(seed, boost::hash_value(k.substr0));
    boost::hash_combine(seed, boost::hash_value(k.str1));
    boost::hash_combine(seed, boost::hash_value(k.substr1));
    return seed;
  }
};




// Data type

using AdjListWithSubData = std::unordered_map<SubStrandsPairKey, SubStrandsPairNode,
                                              SubStrandsPairKeyHasher>;

using AdjSubVec = std::unordered_map<SubStrand, std::vector<SubStrandsPairKey>,
                                     SubStrandHasher>;

/// A vector contains the keys to the elements of AdjListWithSubData in the same sheet.
using SubStrandsPairKeyVec = std::vector<SubStrandsPairKey>;

using SubStrandSet = std::unordered_set<SubStrand, SubStrandHasher>;


/// Adjacency List with sub_strands.
class AdjListWithSub {
public:

  // ***************************************************************************
  // Public Member Functions
  // ***************************************************************************

  void add_sheet();

  /// @brief    Erase the edges that contain too short node.
  /// @return   Wheter the erasing has occured.
  bool erase_too_short(SubStrandSet const& short_subs);

  /// Ensure that all the edges have the counterpart with a reverse direction,
  /// if fallback mode
  void ensure_undirected();


  /// Change the id of SubStrands in keys, after sorting the substrands based on resnums.
  void update_key_substr(std::unordered_map<SubStrand, SubStrand, SubStrandHasher> const& conv);


  /// @brief  Remove the smaller connections
  ///         If there are 2 pairs in both direction,
  ///         remove the pair with less paired residues (weaker connection).
  ///         Ignore the undirected sheets generated through the fallback mode.
  void cleanup(std::vector<IndexType> const& strand_indices);


  /// @brief  Remove the path data specified by key from data and sheets_key_vec.
  void remove(SubStrandsPairKey const& key);


  /// @brief  Translate data into adj_sub_vec
  void gen_adj_sub_vec();


  // Accessor Methods

  /// data
  auto const& map() const { return data; }
  auto const& map(SubStrandsPairKey const& key) const { return data.at(key); }
  auto insert_map(SubStrandsPairKey const& key, bool const dir) {
    return data.insert({key, SubStrandsPairNode{dir}});
  }
  void add_count_pairs_map(SubStrandsPairKey const& key) {++data[key].residue_pairs;}
  void update_delta(SubStrandsPairKey const& key, int const delta_1, int const delta_2);

  /// adj_sub_vec
  auto const& adj_substr_vec() const { return adj_sub_vec; }

  /// Convert a set of substrnd into an AdjSubVec object.
  /// @param  first An iterator to the first element of the Sub-Strand set.
  /// @param  last  Past the end iterator of the Sub-Strand set.
  template <typename Iter>
  AdjSubVec substr_vec2adj_sub_vec(Iter first, Iter last) const {
    SubStrandSet tmp_substr_set;
    tmp_substr_set.insert(first, last);

    AdjSubVec new_adj_sub;

    // for each input vec
    for (; first != last; ++first) {

      // if edge (has no adjacent substrand)
      if (adj_sub_vec.count(*first) == 0) {
        continue;
      }

      new_adj_sub.insert({*first, SubStrandsPairKeyVec{}});

      // for each adjacent candidates
      for (auto const& pair_candidate : adj_sub_vec.at(*first)) {
        #ifndef NDEBUG
        auto const ss0 = pair_candidate.sub0();
        #endif // ifndef NDEBUG
        auto const ss1 = pair_candidate.sub1();

        assert(*first == ss0);

        if (tmp_substr_set.count(ss1)) {
          new_adj_sub[*first].push_back(pair_candidate);
        }
      }

      // if no adjacent found
      if (new_adj_sub[*first].size() == 0) {
        new_adj_sub.erase(*first);
      }
    }
    return new_adj_sub;
  }



  [[deprecated("Use sheet::Sheets::substr_keys instead.")]]
  auto const& key_vec_in_sheet() const { return sheets_key_vec; }
  auto register_key_vec(SubStrandsPairKey const& key) {
    return sheets_key_vec.back().push_back(key);
  }


  [[deprecated("Use sheet::Sheets:: instead.")]]
  auto const& directed_flags_sheet() const { return sheets_directed; }
  auto register_sheet_directed(bool const directed) {
    return sheets_directed.push_back(directed);
  };

protected:

  // ***************************************************************************
  // Protected Member Variables
  // ***************************************************************************

  /// Data variable, same as the accumulation of sheets_key_vec.
  /// For undirected sheets, some adjacent substrands were directized in
  /// DirectedAdjacencyList::init_sheets()
  AdjListWithSubData data{};

  /// Key: SubStrand, Value: a vector of SubStrands that is on the right of key substrand.
  AdjSubVec adj_sub_vec{};

  /// A vector of sheets
  std::vector<SubStrandsPairKeyVec> sheets_key_vec{};

  /// Whether the corresponding sheet in \c sheets_key_vec is directed.
  /// \c true for directeed, \c false otherwise.
  std::vector<bool> sheets_directed{};


  /// A vector of keys that were removed in remove().
  std::vector<SubStrandsPairKey> removed_keys{};


  // ***************************************************************************
  // Protected Member Functions
  // ***************************************************************************

  void replace_erased_itr(SubStrandsPairKey const& old_key,
                          SubStrandsPairKey const& new_key);

};



} // namespace sheet

#endif // ifndef SHEET_ADJ_LIST_WITH_SUB_H_
