// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <string>

#include "sheet/adj_list_with_sub.h"
#include "sheet/directed_adjacency_list.h"

#include "functions.h"


namespace out {

// **********************************************************************
// Member Function of substr2str
// **********************************************************************

std::string substr2str::operator()(sheet::SubStrand const& ss) const {
  return str(ss, *ptr);
} // public member operator()

/// @brief  Output a sheet::SubStrand into a string, converting Strand ID into SSE ID.
std::string substr2str::str(sheet::SubStrand const& ss,
                            sheet::DirectedAdjacencyList const& adj) {
  return std::to_string(adj.strand_indices[ss.str]) + "_" + std::to_string(ss.substr);
} // static member function substr2str::str()

} // namespace out


namespace mmcif {


// **********************************************************************
// Public Member Function loop_head()
// **********************************************************************

void mmcif_like::loop_head(std::vector<std::string> const& keys) const {
  // Header
  os << "#\n"
     << "loop_\n";

  // field names
  for (auto const& key : keys) {
    os << key_head << key << "\n";
  }
} // public member function loop_head()




} // namespace mmcif

