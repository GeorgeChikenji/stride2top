// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PDB_LOOP_H_
#define PDB_LOOP_H_

#include <vector>
#include <unordered_map>

#include "pdb/atom.h"
#include "pdb/constants.h"
#include "pdb/sse.h"


namespace pdb {
#ifdef WITH_LOOP

/// Stores the loop information
class LOOP: public SSE_Base {
public:

  // *********************************************************************************
  // Public Member Functions
  // *********************************************************************************

  LOOP(int const i, int const e, IndexType const ind,
       std::unordered_map<int, std::string> const& ca_lines):
    SSE_Base{i, e, ind, ca_lines}, zero_size{false}
   {}

  LOOP(int const i, int const e, IndexType const ind):
    SSE_Base{i, e, ind}, zero_size{true}
   {}

  // ************************************************************************
  // Public Member Variables
  // ************************************************************************

  /// True if the length of this loop is 0.
  /// If this value is true, init is the end resnum of the previous SSE and
  /// end is the first resnum of the next SSE.
  bool const zero_size{false};


};

#endif // WITH_LOOP


} // namespace pdb

#endif // ifndef PDB_LOOP_H_
