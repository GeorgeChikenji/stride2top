// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef ADJ_OUT_H_
#define ADJ_OUT_H_

#include <ostream>
#include <string>
#include <vector>

#include "sheet/directed_adjacency_list.h"

namespace adj_out {

void adj_list_out(std::ostream & os, sheet::DirectedAdjacencyList const& adj);


std::vector<std::string> gen_adj_list_vec(sheet::DirectedAdjacencyList const& adj);

} // namespace adj_out

#endif // ifndef ADJ_OUT_H_
