// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_CB_SIDE_H_
#define SHEET_CB_SIDE_H_

#include <array>
#include <vector>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "pdb/sses.h"
#include "sheet/adjacent_substrand.h"
#include "sheet/adj_list_with_sub.h"
#include "sheet/common.h"
#include "sheet/pairs.h"
#include "sheet/exceptions.h"
#include "sheet/sub_strands_range.h"
#include "sheet/zone_residue.h"

namespace sheet {

/// Structure of the data stored in the StrictZone.strict
struct ZoneInfo {
  ZoneInfo() = default;

  enum SideStatus: char {Undefined = -1, Upper = 1, Lower = 0};
  enum BridgeType: char {NoBridge, ParallelHbonds, ParallelNoHbonds, SmallRing, LargeRing};

  /// @retval true  If adding succeeded.
  /// @retval false Otherwise.
  bool add_pair(ZoneResidue const& residue, bool const hbonded,
                BridgeType const bridge_type);

  bool colored{false};
  SideStatus side{Undefined};

  /// [paired residue on non-hbonded side, hbonded side]
  std::array<ZoneResidue, 2> adj_residues{{ZoneResidue{}, ZoneResidue{}}};

  /// Store if the corresponding adj_residue item has already been initialized.
  std::array<bool, 2> adj_set{{false, false}};

  /// Store the type of the bridge paired residue on [non-hbonded, hbonded] side.
  std::array<BridgeType, 2> bridge_type{{NoBridge, NoBridge}};
};


struct StrictZone {

  // *************************************************************
  // Public Member Types
  // *************************************************************

  /// @brief  An unordered_set of ZoneResidues.
  ///         Used as a set of remainders or a set of queue contents.
  using ZoneResidueSet = std::unordered_set<ZoneResidue, ZoneResidueHasher>;


  /// Store the relative direction of each strand against the base strand
  /// index : serial SSE ID
  /// value : -  1 (Parallel)     the strand is parallel to the base strand
  ///         -  0 (AntiParallel)    anti-parallel
  ///         - -1 (NotSet)  the element is not initialized yet
  enum RelDir: char {NotSet = -1, AntiParallel = 0, Parallel = 1};

  /// A vector of \c RelDir
  using RelativeDirsVec = std::vector<RelDir>;



  // *************************************************************
  // Public Member Functions
  // *************************************************************

  StrictZone(pdb::SSES const& sses_, SubStrandsRange & sub_strands_ref,
             AdjacentSubStrands & adj_substrands_ref):
    sses{sses_}, strand_indices{sses.gen_index_vec('E')}, strict{init_strict()},
    sub_strands{sub_strands_ref}, adj_substrands{adj_substrands_ref}
  {}

  /// @param  sse_id  SSE ID of the target residue
  /// @param  resnum  Residue number of the target residue
  /// @param  paired_sse_id
  /// @param  paired_resnum
  /// @param  hbonded true if the first residue is hbonding to the direction of the second.
  void On(IndexType const sse_id, int const resnum, IndexType const paired_sse_id,
          int const paired_resnum, bool const hbonded,
          ZoneInfo::BridgeType const bridge_type);
  void On_one(ZoneResidue const& res0, ZoneResidue const& res1, bool const hbonded,
              ZoneInfo::BridgeType const bridge_type);


  /// Run BFS (Breadth First Search) and set the side data for all registered residues.
  AdjListWithSub decide_side(AdjList const& undirected_adj_list);


  /// @brief  Accessor to the this->strict using a ZoneResidue class object
  /// @param  zone_res  An object of class ZoneResidue to access a content of strict[][].
  /// @retval A const reference to a corresponding object of class ZoneInfo.
  ZoneInfo const& strict_c_info(ZoneResidue const& zone_res) const;


  /// @brief  Accessor to the this->strict using a ZoneResidue class object
  ///         Non-const version of strict_c_info().
  /// @param  zone_res  An object of class ZoneResidue to access a content of strict[][].
  /// @retval A (non-const) reference to a corresponding object of class ZoneInfo.
  ZoneInfo & strict_info(ZoneResidue const& zone_res);




#ifdef LOGGING
  void PrintZone() const;

  void PrintSide() const;

  void PrintList(AdjListWithSub const& adj_list) const;
#endif // ifdef LOGGING


  // *************************************************************
  // Public Member Variables
  // *************************************************************


  /// Reference to SSES class
  pdb::SSES const& sses;

  /// a vector of sse_id for strands
  std::vector<IndexType> strand_indices;


  /// Data storage for zone info. Convert the SSE_ID to the serial SSE_ID when accessing.
  /// 2-dimensional vector
  /// N_SSE * N_res_SSE
  std::vector<std::vector<ZoneInfo>> strict{};


  /// Just resize the vector by the number of strands
  /// array: resnums (not indices from 0) of the first and the last residue of the range
  SubStrandsRange & sub_strands;

  /// A reference to the protected member variable adj_substrands in DirectedAdjacencyList.
  AdjacentSubStrands & adj_substrands;

protected:

  // *************************************************************
  // Protected Member Functions
  // *************************************************************

  std::vector<std::vector<ZoneInfo>> init_strict() const;


  // Helper functions for decide_side
  /// @brief  Generate a set of unproccessed and colored residues
  void collect_colored(ZoneResidueSet & set) const;


  /// Calculate the gap (the number of residue differences) between paired strands.
  void calc_deltas(AdjListWithSub & adj_list) const;


  /// @brief  Helper function for calc_deltas(). Get the first residue that resides in the
  ///         bridge between the base Sub-Strand and the Sub-Strand adjacent to the base.
  /// @param  zres_base_start  The first residue on the base Sub-Strand to start counting.
  /// @param  toward_c_term    Whether to search from start (zres_base_start) toward
  ///                          C-term side of base.
  /// @param  adj_ss           The Sub-Strand adjacent to the base.
  /// @return A tuple of following 2 values (base_delta, first_adj_res).
  ///         * base_delta     The number of residues inside the strict hydrogen-bonding
  ///                          zone on the base Sub-Strand (this value will be initialized
  ///                          inside this function).
  ///         * first_adj_res  The first residue that is on the Sub-Strand adjacent to
  ///                          the base Sub-Strand.
  std::tuple<int, ZoneResidue> count_delta_base(ZoneResidue const& zres_base_start,
                                                bool const toward_c_term,
                                                SubStrand const& ss_adj) const;


  /// @brief Count the number of colored residues (residues inside the strict hydrogen
  ///        bonding pattern) from \c start to \c last. Both \c start and \c last must be
  ///        belong to the same Sub-Strand.
  /// @param  start  The first residue to start counting.
  /// @param  step   A step value to move at one count. This value should be 1 or -1.
  ///                If 1, count to the C-term side, otherwise N-term side.
  /// @param  last   The final residue where the counting should stop. This residue will
  ///                be included in the count.
  /// @return The number of residues between \c start and \c last .
  int count_delta_adj(ZoneResidue const& start, int const step,
                      ZoneResidue const& last) const;


  int count_delta_1(SubStrand const& ss_base, SubStrand const& ss_adj,
                    bool const dir) const;
  int count_delta_2(SubStrand const& ss_base, SubStrand const& ss_adj,
                    bool const dir) const;


  /// @brief  Run bfs.
  void bfs(AdjListWithSub & adj_list,
           ZoneResidueSet & remainder,
           AdjList const& undirected_adj_list);

  // Helper Functions for bfs()

  /// @brief  Add the count for hbonded residue pair.
  /// @param  key               A key to access the adj_list.
  /// @param  gap               The number of residue difference of parallel strands.
  /// @param  adj_list          Adjacency List
  /// @param  direction_to_adj  Type of strand pairing.
  /// @exception  std::runtime_error  If the direction_to_adj differs from the last value.
  void add_adj_list_count(SubStrandsPairKey const& key,
                          AdjListWithSub & adj_list,
                          bool const direction_to_adj);


  /// @brief  Resolve the relative direction of adj and update the \c relative_directions .
  /// @return If fallback mode should be turned on.
  /// @retval false   Normal result.
  /// @retval true    Indicates that the resolved direction does not coincide with the
  ///                 previous value.
  /// @param target               A residue on the base SSE.
  /// @param adj                  A residue on the paired SSE.
  /// @param direction_to_adj     The direciton of a strands pair between target and adj.
  /// @param relative_directions  A vector of relative directions on the same sheet.
  /// @exception TargetRelativeDirectionNotSet Indicates that the relative direcion of
  ///                                          the target strand is not set.
  bool update_rel_dir(ZoneResidue const& target, ZoneResidue const& adj,
                      bool const direction_to_adj,
                      RelativeDirsVec & relative_directions) const;


  /// Generate the key to add into the adj_list hash and the gap between strand positions.
  SubStrandsPairKey gen_list_key(ZoneResidue const& res0, ZoneResidue const& res1,
                                 bool const hbonded, RelativeDirsVec const& rel_dir) const;


  void add_adj_substrands(ZoneResidue const& res0, ZoneResidue const& res1,
                          bool const hbonded, bool const rel_dir);


  /// @brief  Push the paired or the sequentially adjacent residue into the search queue
  ///         if it's not already in the queue.
  /// @param  queue     Queue to push.
  /// @param  contents  A set of queue contents.
  /// @param  new_res   A residue to be pushed into queue.
  void push_into_queue(std::deque<ZoneResidue>& queue,
                       ZoneResidueSet & contents,
                       ZoneResidue const& new_res) const;



  /// @brief  Check if the residue before or after the target residue is valid.
  /// @return A ZoneResidue class object to be pushed back into a vector.
  ///         \c has_value of the returned ZoneResidue will be \c false if it should
  ///         be ignored, otherwise true.
  /// @param  target  Current target ZoneResidue. (queue.front())
  /// @param  diff    A difference in residue number between the target and the ba.
  ///                 1 for the next residue, -1 for the previous residue.
  /// @param  remainder  A reference to the list of remaining residue to examine.
  ZoneResidue ba_check(ZoneResidue const& target, int const diff,
                       ZoneResidueSet const& remainder) const;

};

} // namespace sheet

#endif // ifndef SHEET_CB_SIDE_H_
