// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_SUB_STRANDS_RANGE_H_
#define SHEET_SUB_STRANDS_RANGE_H_

#include <array>
#include <limits>
#include <vector>

#include "sheet/adj_list_with_sub.h"
#include "sheet/common.h"
#include "sheet/zone_residue.h"

namespace sheet {


class SubStrandsRange {
public:
  // *********************************************************************************
  // Public Member Types
  // *********************************************************************************

  struct SubStrandsIters {
    using SubStrandsIdxIter = std::vector<SubStrand>::const_iterator;
    explicit SubStrandsIters(SubStrandsIdxIter const f_, SubStrandsIdxIter const l_):
      f{f_}, l{l_}, erased{f == l} {}
    auto begin() const {return f;}
    auto end() const {return l;}
    SubStrandsIdxIter const f;
    SubStrandsIdxIter const l;
    bool const erased{false};
  };


  // *********************************************************************************
  // Public Member Functions
  // *********************************************************************************

  // ********************************************
  // Construction Related Functions
  // ********************************************

  /// Initialize data with 1 \c default_range sub strand for each strand.
  SubStrandsRange(std::size_t const n_strands):
    data{n_strands, std::vector<std::array<int, 2>>{1, default_range}}
  {}


  /// @brief  Extend the range of known sub_strand if already exists, otherwise create.
  /// @param  res  A residue to check the sub_strand it belongs to.
  void extend_substrand(ZoneResidue const& res);


  /// @brief    * Add new default_ranges for the next iteration of bfs.
  ///           * Remove the ranges that are too short.
  /// @return   A set of too short sub_strands.
  SubStrandSet cleanup_sheet();


  /// * Remove the last default_range from all sub_strands
  /// * Sort the Sub-Strands based on the initial resnum.
  /// * Set up index_vec and sub_strands_iters_vec.
  std::unordered_map<SubStrand, SubStrand, SubStrandHasher> finish();


  // ********************************************
  // Accessor Functions
  // ********************************************

  /// Return the substr id of the last substr in the specified str id.
  /// @throw std::out_of_range if \c s does not exist.
  IndexType last_substr_id(IndexType const str_id) const;


  /// Return the residue number of the N-term edge of the specified Sub-Strand.
  /// @throw std::out_of_range if \c s does not exist.
  int n_term_res(SubStrand const& s) const;


  /// Return the residue number of the C-term edge of the specified Sub-Strand.
  /// @throw std::out_of_range if \c s does not exist.
  int c_term_res(SubStrand const& s) const;


  /// Return the residue number of the specified edge of the specified Sub-Strand.
  /// @param  i    0 for N-term, 1 for C-term.
  /// @throw std::out_of_range if \c s does not exist.
  int term_res(SubStrand const& s, IndexType const i) const;


  /// Return a SubStrand that is at the most N-term side of the Strand.
  /// @param str  A serial strand id of the wanted strand.
  SubStrand n_term_sub(IndexType const str) const;

  /// @brief C-term version of n_term_sub().
  ///        The returned Sub-Strand is the same if the str has only 1 Sub-Strand
  /// @sa n_term_sub
  SubStrand c_term_sub(IndexType const str) const;



  /// @return A sorted vector of SubStrands managed by this class.
  auto const& vec() const { return index_vec; }

  /// Return a set of begin() and end() of the index_vec on a strand str.
  /// @param  str   A serial strand ID.
  SubStrandsIters const& vec(IndexType const str) const {
    return sub_strands_iters_vec.at(str);
  }


  // ********************************************
  // Logging Functions for Debugging
  // ********************************************

  #if defined(LOGGING) && defined(LOG_SHEET)
  void print(std::vector<IndexType> const& strand_indices) const;
  #endif // ifdef LOGGING

protected:
  // *********************************************************************************
  // Protected Member Functions
  // *********************************************************************************

  /// Helper function of finish(). Generates a vector of indices of sorted SubStrands.
  std::vector<IndexType> sorted_indices(std::vector<std::array<int, 2>> const& substrs);


  /// Helper function of finish(). Convert index vector into a map of SubStrands.
  std::unordered_map<SubStrand, SubStrand, SubStrandHasher> convert_to_substrand(std::vector<std::vector<IndexType>> const& changed_index_map) const;

  /// Helper function of finish(). Generates a vector of SubStrands.
  void init_index_vec();


  /// Helper function of finish().
  void init_sub_strands_iters_vec();

  /// @brief  Helper function of init_sub_strands_iters_vec().
  ///         Generates the iterator to \c sub_strands_list used for Range-Based for loops.
  ///         Find the range using binary-search.
  /// @return An object with methods begin() and end().
  SubStrandsIters find_sub_strands(IndexType const serial_str_id) const;


  // *********************************************************************************
  // Protected Member Variables
  // *********************************************************************************


  /// default range for sub_strands, this doesn't mean anything.
  static constexpr std::array<int, 2> default_range = {{std::numeric_limits<int>::min(),
                                                        std::numeric_limits<int>::min()}};

  /// Minimum number of residues for any Sub-Strands.
  static constexpr IndexType min_sub_str_len = 2;


  /// Store the range data
  std::vector<std::vector<std::array<int, 2>>> data{};


  /// A vector of SubStrands. SubStrands are sequentially sorted.
  std::vector<SubStrand> index_vec{};

  /// A vector of SubStrandsIters.
  /// Each element is for accessing index_vec by range-based-for loops.
  std::vector<SubStrandsIters> sub_strands_iters_vec{};
};


} // namespace sheeet

#endif // ifndef SHEET_SUB_STRANDS_RANGE_H_
