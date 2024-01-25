// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_DIRECTED_ADJACENCY_LIST_H_
#define SHEET_DIRECTED_ADJACENCY_LIST_H_

#include <string>

#include "pdb/sses.h"
#include "pdb/stride_stream.h"

#include "sheet/adjacent_substrand.h"
#include "sheet/adj_list_with_sub.h"
#include "sheet/common.h"
#include "sheet/cb_side.h"
#include "sheet/pairs.h"
#include "sheet/sheets.h"
#include "sheet/sub_strands_range.h"
#include "sheet/substr_pair_attr.h"


namespace sheet {



// ***********************************************************************************
// Struct LastPair
// ***********************************************************************************
struct LastPair {
  LastPair(int const i, int const j, int const dir, bool const f) :
    res0{i}, res1{j}, hbond_dir{dir}, first{f} {}

  LastPair() = default;

  int res0{0};
  int res1{0};
  int hbond_dir{-1};
  bool first{true};
};



using ATOM_vec_iter = std::vector<pdb::ATOM>::const_iterator;

// ***********************************************************************************
// Class DirectedAdjacencyList
// ***********************************************************************************

class DirectedAdjacencyList {
public:

  // ***************************************************************************
  // Public Member Types
  // ***************************************************************************

  using AdjIndList = std::unordered_map<IndexType, std::vector<IndexType>>;
  using AdjAttrMap = std::unordered_map<SubStrandsPairKey, StrandsPairAttribute,
                                        SubStrandsPairKeyHasher>;

  // ***************************************************************************
  // Public Member Functions
  // ***************************************************************************

  /// Constructor
  DirectedAdjacencyList(pdb::SSES const& sses_, pdb::stride_stream & stride) :
    sses{sses_},
    strand_indices{sses.gen_index_vec('E')},
    sub_strands_range{strand_indices.size()},
    adj_substrands{},
    strict_zone{sses, sub_strands_range, adj_substrands},
    adj_list_with_sub(init_list(Pairs{sses, stride})),
    adj_list{translate_sub()},
    adj_index_list{init_adj_index_list()},
    // Run fix_undirected_paths inside init_sheets()
    sheets{init_sheets()},
    adj_attr{init_adj_attr()}
  {}



  /// @brief  An accessor method to \c sub_strands_range .
  /// @return A reference to the \c sub_strands with const qualifier.
  SubStrandsRange const& substrs() const { return sub_strands_range; }

  /// @brief  An accessor method to \c adj_substrands .
  /// @return A reference to it with a const qualifier.
  AdjacentSubStrands const& adj_substrs() const { return adj_substrands; }


  /// Access the first ATOM with \c SubStrand class object.
  ATOM_vec_iter atom_cbegin(SubStrand const& ss) const;

  /// Return the past the end iterator of the vector of ATOMs,
  /// given \c SubStrand class object.
  ATOM_vec_iter atom_cend(SubStrand const& ss) const;


  /// @brief  Access the generated adjacency list with sub strands by an SSE ID.
  /// @return A reference to the target data (of type \c SubStrandsPairNode ).
  /// @param  sse_id        A SSE ID of the wanted sse. Not a serial SSE ID.
  /// @param  substrand_id  An index of the target sub-strand. [default: 0]
  /// @exception std::out_of_range if no such hbonded strands pair found.
  SubStrandsPairNode const& edge_info(IndexType const sse_id_0, IndexType const sse_id_1,
                                      IndexType const substrand_id_0 = 0,
                                      IndexType const substrand_id_1 = 0) const;


  /// @brief  Generate \c SubStrand class object from the given \c sse_id .
  SubStrand substr(IndexType const sse_id, IndexType const substr_id = 0) const;


  /// @brief  Run BFS and find the path from ss0 to ss1.
  ///         Use \c adj_attr instead of this method.
  ///         This version is implemented as a wrapper of the second version and uses
  ///         all the adjacent information available.
  /// @return An object containing the information of the found path.
  StrandsPairAttribute search(SubStrand const& ss0, SubStrand const& ss1) const;

  /// @brief  This version only uses adjacency information given through adj_sub_vec.
  StrandsPairAttribute search(SubStrand const& ss0, SubStrand const& ss1,
                              AdjSubVec const& adj_sub_vec) const;

  /// @brief  Run search() from SSE_IDs.
  ///         Use \c adj_attr instead of this method.
  StrandsPairAttribute search(IndexType const sse_id_0, IndexType const sse_id_1) const;


  /// @brief  Look up \c adj_attr also in the reversed direction,
  //          if it failed to reach in the default direction. (ss0 -> ss1 or ss1 -> ss0).
  StrandsPairAttribute const& attr(SubStrand const& ss0, SubStrand const& ss1) const;

  /// Accessor to the AdjListWithSub class object.
  auto const& adj_sub() const { return adj_list_with_sub; }

  /// Accessor method for strict zone
  auto const& get_strict_zone() const { return strict_zone; }

  // *******************************************************************************
  // DEBUG member function
  // *******************************************************************************
  #ifdef DEBUG
  void print_list(AdjList const& adj_list, std::string const& title="") const;
  #endif // ifdef DEBUG
  // *******************************************************************************


  // ***************************************************************************
  // Public Member Variables
  // ***************************************************************************

  pdb::SSES const& sses;
  std::vector<IndexType> const strand_indices{};
protected:

  /// will be initialized inside the StrictZone class in init_list()
  SubStrandsRange sub_strands_range{0};


  /// will be initialized inside the StrictZone class in init_list()
  AdjacentSubStrands adj_substrands{};

  // Store the region of residues whose sides of CB atoms (Upper or Lower) can be determined.
  /// will be initialized in init_list()
  StrictZone strict_zone;

  /// Adjacency Hash List with Sub-Strands
  AdjListWithSub adj_list_with_sub{};

public:

  /// Traditional Adjacency Hash List
  AdjList const adj_list{};
  AdjIndList const adj_index_list{};


  /// Stores the sheets information (Number of Strands in the sheet etc ...).
  Sheets const sheets{};

  /// A map of strands pair attributes.
  AdjAttrMap const adj_attr{};

protected:

  // ***************************************************************************
  // Protected Member Functions
  // ***************************************************************************

  AdjListWithSub init_list(Pairs const& pairs);

  /// create a undirected adjacency list
  std::unordered_set<IndexType> add_undirected_adj_list(AdjList & undirected_adj_list, 
                                                        Pairs::PairsVec const& involved_pairs,
                                                        pdb::SSE const& e) const;


  /// generate a vector of sse_ids adjacent to a target SSE
  std::unordered_set<IndexType> adj_id_set(Pairs::PairsVec const& involved_pairs) const;


  /// Given the number of votes for parallel and anti-parallel, return the direction
  bool direction_by_majority(IndexType const sse_id, pdb::SSE const& e,
                            std::unordered_map<IndexType, std::unordered_map<bool, unsigned>> direction_counter) const;


  /// @brief  Ensure that the undirected_adj_list is symmetric.
  /// @param  undirected_adj_list An undirected adjacency list to check the symmetry.
  /// @exception non_symmetric When no reverse directional strands pair exists.
  /// @exception non_symmetric_dir_opposite When the direction of the strand pair differs.
  void undirected_adj_symmetry_check(AdjList const& undirected_adj_list) const;


  /// Generate directed adjacency list from a undirected adjacency list
  void gen_strict_zone(AdjList const& undirected_adj_list,
                       Pairs const& pairs, std::vector<IndexType>& strands,
                       std::unordered_map<IndexType, std::unordered_set<IndexType>> const& adj_index_map);


  /// Translate the AdjListWithSub into AdjList
  AdjList translate_sub() const;


  /// Initialize adj_index_list
  AdjIndList init_adj_index_list() const;

  /// Initialize sheets.
  Sheets init_sheets();

  /// Initializee adj_attr
  AdjAttrMap init_adj_attr() const;

  // *******************************
  // Helper Functions for search()
  // *******************************

  using PathParents = std::unordered_map<SubStrand, SubStrand, SubStrandHasher>;

  /// @return A vector of reversed path
  void search_backtrace(StrandsPairAttribute & attr, PathParents const& parents) const;


  /// @return  A path from the \c first to the \c last node.
  PathParents search_bfs(SubStrand const& first, SubStrand const& last,
                         AdjSubVec const& adj_sub_vec) const;

};




} // namespace sheet

#endif //ifndef SHEET_DIRECTED_ADJACENCY_LIST_H_

