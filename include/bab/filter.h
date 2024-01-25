// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef BAB_FILTER_H_
#define BAB_FILTER_H_

#include <functional>
#include <vector>

#include <Eigen/Core>

#include "pdb/sses.h"
#include "sheet/directed_adjacency_list.h"
#include "sheet/adj_list_with_sub.h"
#include "bab/side.h"

namespace bab {

#define DEFAULT_MAX_RES_LEN 60
#define DEFAULT_MAX_MID_STR 1
#define DEFAULT_MAX_SCORE 0.6
#define DEFAULT_SIDE_MIN_DIST 1.0


/// Store the left-handedness score as a result of filtering.
struct BabFilterResult {
  BabFilterResult() = default;

  BabFilterResult(sheet::SubStrand const& sub_f, sheet::SubStrand const& sub_l):
    sub_first{sub_f}, sub_last{sub_l}
  {}

  #ifdef WITH_STAT
  /// Return a string representing this object.
  /// Output format is
  /// #substr0  #substr1  #score  #mid_res_len
  std::string str() const {
    return sub_first.string() + "\t" + sub_last.string() + "\t" + str_short();
  }


  std::string str_short() const {
    return (success ? std::to_string(left_score) : std::to_string(0.0)) + "\t" +
           std::to_string(mid_res_len) + "\t" +
           std::to_string(n_mid_str) +  "\t" +
           std::to_string(jump) +  "\t" +
           (success ? "T" : "F") + "\t" + std::to_string(non_bab_reason);
  }
  #endif // WITH_STAT


  /// Whether the filtering was processed correctly, false if not.
  bool success{false};


  /// If success == true, this value states the type of connection. Otherwise, no meaning.
  /// The value is [0, 7]. Only 3 bits from the LSB (Least Significant Bit) are used.
  /// These bits has following meanints.
  /// * 1st bit: 1 if the connection has Helices.
  /// * 2nd bit: 1 if the connection has Loops.
  /// * 3rd bit: 1 if the connection has Strands.
  /// The value is one of the followings.
  /// * 0: No meaning default value.
  /// * 1: Beta - Alpha - Beta
  /// * 2: Beta - Loop  - Beta
  /// * 3: beta - Alpha - Beta with loops.
  /// * 4: Beta - Beta  - Beta. Mid-Beta is on other sheets.
  /// * 5: Beta - Alpha - Beta with mid-Beta-strands (on other sheets).
  /// * 6: Beta - Beta  - Beta with loops. Mid-Beta is on other sheets.
  /// * 7: Beta - Alpha - Beta with mid-Beta-strandds (on other sheets) and loops.
  unsigned connection_type{0};


  /// If success == false, this variable states the cause of failure.
  /// Has NO meaning if success == true.
  /// * 0:  No information.
  /// * 1:  One of the edge SSE is not a strand.
  /// * 2:  No mid-SSEs to perform bab filter.
  /// * 3:  ss0 and ss1 is not reachable or not Parallel.
  /// * 4:  The total number of residues in the mid-part exceeded the limit.
  /// * 5:  The total number of mid-strands in the same sheet exceeded the limit.
  unsigned non_bab_reason{0};


  /// A left handed ness score.
  double left_score{0.0};

  /// The number of residues between the last residue of the first Sub-Strand
  /// and the first residue of the last Sub-Strand. (Strands in the same sheet included).
  /// NOT COUNTED FOR LOOP REGIONS UNLESS WITH_LOOP IS DEFINED (counted if WITH_LOOP).
  unsigned mid_res_len{0};

  /// The total number of triangle atom decision.
  unsigned tri_atom_count{0};

  /// The number of residues actually processed by the filtering.
  /// (The number of real, non-padding atoms.)
  unsigned n_pdb_atoms{0};

  /// The number of strands in the same sheet.
  unsigned n_mid_str{0};

  unsigned jump{0};


  sheet::SubStrand sub_first{0, 0};
  sheet::SubStrand sub_last{0, 0};
};



/// key  : A pair of sub-strands.
/// value: A vector of Sides between the sub-strands pair.
using SidesMap = std::unordered_map<sheet::SubStrandsPairKey, std::vector<Side>,
                                    sheet::SubStrandsPairKeyHasher>;

using SeqIter = std::vector<pdb::IndexType>::const_iterator;



/// Function object to run filter
using FilterFunc = std::function<bool(bab::SeqIter const, bab::SeqIter const,
                                      unsigned const)>;



class BabFilter {
public:

  // *********************************************************************************
  // Public Member Functions
  // *********************************************************************************

  explicit BabFilter(sheet::DirectedAdjacencyList const& adj_,
                     std::function<bool(double, double)> const comp_=std::greater<double>(),
                     unsigned const c_res_len=DEFAULT_MAX_RES_LEN,
                     unsigned const c_mid_str=DEFAULT_MAX_MID_STR,
                     double const c_score=DEFAULT_MAX_SCORE,
                     double const c_side_min_dist=DEFAULT_SIDE_MIN_DIST):
    sses{adj_.sses},
    adj{adj_},
    sides_map{init_sides_map()},
    comp{comp_},
    cut_off_res_len{c_res_len},
    cut_off_mid_str{c_mid_str},
    cut_off_left_score{c_score},
    cut_off_side_min_dist{c_side_min_dist}
  {}


  /// @brief  Run filter for a set of SSEs. If the first and/or the last strands have
  ///         multiple sub-strands, try filtering with all pairs of sub-strands in the
  ///         same sheet.
  ///         * Iterators must be bidirectional.
  ///         * The set of SSEs to which this filtering will be performed is [first, last).
  /// @param first   An iterator pointing to the first element of the set of SSEs.
  /// @param last    An iterator pointing to the end of the range (Past the end iterator).
  /// @param reverse A sequence of bits 0 for non-reversed, 1 or reversed.
  /// @return Whether the given range is a left-handed beta-alpha-beta unit.
  /// @retval true   If left-handed.
  /// @retval false  If right-handed or strictly not a beta-alpha-beta unit.
  bool operator()(SeqIter const first, SeqIter const last, unsigned const reverse=0);


  /// @brief  In this version, the SubStrands had already been specified by user.
  ///         * There is no SSE type check. If the first or (last - 1) don't point to
  ///           strands, out of range access to a vector will occur.
  ///           MAKE SURE THAT THE EDGE SSES ARE STRANEDS !!!
  bool operator()(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1,
              SeqIter const first, SeqIter const last,
              unsigned const reverse=0);

  #ifdef WITH_STAT
  BabFilterResult const& result() const { return last_result; }
  #endif // WITH_STAT


protected:

  // *********************************************************************************
  // Protected Member Functions
  // *********************************************************************************

  /// @brief  Initialize sides map.
  SidesMap init_sides_map() const;

  std::vector<Side> gen_sides_vec(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1) const;


  /// given sides for each strand, return the left ratio of the given atom point
  /// @param b0   SSE_ID of the first strand
  /// @param b1   SSE_ID of the second strand
  std::tuple<unsigned, unsigned> count_left_tri(sheet::SubStrand const& b0,
                                                sheet::SubStrand const& b1,
                                                bool const b0_reverse,
                                                bool const b1_reverse,
                                                ATOM_vec_iter a_begin,
                                                ATOM_vec_iter const a_end) const;


  /// @brief    Run filter for one beta-alpha-beta unit
  /// @return   An object with filter results.
  /// @param first    Iterator to the SSE_ID of the first beta of the unit
  /// @param last     Iterator to the SSE_ID of the last beta of the unit
  /// @param reverse  reverse flag (if i-th bit is 1, an SSE of SSE_ID i is reversed)
  BabFilterResult filter_one_unit(sheet::SubStrand const& b0, sheet::SubStrand const& b1,
                                  bool const b0_reverse, bool const b1_reverse,
                                  SeqIter const first, SeqIter const last);


  /// Check if the given range [first, last) meets the bab unit condition.
  /// @return  Whether the given range meets the condition, and if not, the reason.
  ///          The reason number is consistent with the variable non_bab_reason of
  ///          BabFilterResult class. See non_bab_reason for detailed possible values.
  /// @retval 0 if the range [first, last) meets the condition.
  unsigned non_bab_condition(SeqIter first, SeqIter last) const;



  // *********************************************************************************
  // Protected Member Variables
  // *********************************************************************************


  /// A reference to an SSES class object.
  pdb::SSES const& sses;

  sheet::DirectedAdjacencyList const& adj;

  /// A map from a pair of sub-strands to a vector of Sides between them.
  SidesMap const sides_map;


  /// Comparison function for the result score and cut_off value.
  std::function<bool(double, double)> const comp;

  /// max residue length to run the filter
  unsigned const cut_off_res_len{DEFAULT_MAX_RES_LEN};

  /// Max number of mid-strands in the same sheet.
  unsigned const cut_off_mid_str{DEFAULT_MAX_MID_STR};

  /// max left_score
  /// if calculated left_score is larger than this value, the filtering will return true
  double const cut_off_left_score{DEFAULT_MAX_SCORE};

  /// Minimum distance between the triangle and the CA atom to use in the filtering.
  /// CA atoms nearer than this value will be ignored.
  double const cut_off_side_min_dist{DEFAULT_SIDE_MIN_DIST};

  #ifdef WITH_STAT
  BabFilterResult last_result{};
  #endif // WITH_STAT
};


} // namespace bab

#endif //ifndef BAB_FILTER_H_

