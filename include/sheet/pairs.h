// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_PAIRS_H_
#define SHEET_PAIRS_H_

#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>

#include <boost/functional/hash.hpp>

#include "pdb/sses.h"
#include "pdb/stride_stream.h"
#include "pdb/tools.h"

#include "sheet/common.h"

namespace sheet {


// **************************************************************************************
// Class definition
// **************************************************************************************
//
// ************************************************************************
// Class Pairs
// ************************************************************************


/// Store the information about the hbonds between strands given a stride output file.
class Pairs {
public:

  // **********************************************************
  // Public Member Types
  // **********************************************************

  /// store the hbond information from stride
  /// array : {res0, res1, reverse_flag}
  ///   res0 : residue number on the target SSE
  ///   res1 : residue number on the paired SSE
  ///   reverse_flag :
  ///         0 for non-reversed (hbond from res0[N] to res1[O])
  ///         1 for reversed (hbond from res1[N] to res0[O])
  using PairsVec = std::vector<std::array<int, 3>>;


  // **********************************************************
  // Public Member Functions
  // **********************************************************


  Pairs(pdb::SSES const& sses, pdb::stride_stream & stride) :
    dnr{read_stride_stream(stride.ss)},
    involved_pairs{init_involved_pairs(dnr, sses)} {}

  /// Sort the involved_pairs vector based on the direcion
  /// (Parallel or Anti-Parallel) of the paired Strand.
  PairsVec resort_involved_pairs(IndexType const serial_str_id, AdjList const& adj_list,
                                 pdb::SSES const& sses) const;


  // **********************************************************
  // Public Member Variables
  // **********************************************************
  PairsVec const dnr{};
  std::vector<PairsVec> const involved_pairs{};

protected:

  // **********************************************************
  // Protected Member Functions
  // **********************************************************

  /// @brief  Given the stream of stride file and record, generate a vector of hbond-pairs
  ///         This function reads only DNR records.
  /// @return hbond-pair_vector each pair consists of 2 fields [N_atom_resnum, O_atom_resnum]
  /// @retval {{5, 13}, {7, 11}, ...} if N in the resnum 5 and O in the resnum 13 are hbonded
  /// @param  ifs     in-file-stream of the stride output file in /tmp directory
  PairsVec read_stride_stream(std::istream & is) const;




  /// @brief            Initialize involved_pairs
  /// @return hashmap   map the SSE_ID to PairsVec generated by involve_with()
  /// @param  sses      A vector of SSEs
  std::vector<PairsVec> init_involved_pairs(PairsVec const& dnr, pdb::SSES const& sses) const;



  /// @brief            find pairs whcih includes the residues of an SSE e
  /// @retval PairsVec  A vector of paired resnum data either of them is in the SSE e.
  ///                   Sorted by the residue number on the SSE e, and DNR to ACC.
  /// @param  e         An SSE class object to find pair
  /// @param  offset    Number of offset resnums to be passed to member function SSE::in_range()
  PairsVec involve_with(PairsVec const& dnr, pdb::SSE const& e, int const offset=1) const;



}; // class Pairs


} // namespace sheet
#endif // ifndef SHEET_PAIRS_H_

