// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef HANDEDNESS_H_
#define HANDEDNESS_H_

#include <ostream>
#include <vector>
#include <boost/program_options.hpp>

#include "sheet/directed_adjacency_list.h"
#include "bab/filter.h"

namespace bpo = boost::program_options;

namespace rare {

/// @brief  A helper function of output_handedness. Actually print all left-handed
///         beta-alpha-beta, beta-loop-beta, or beta-beta(on other sheet)-beta.
std::vector<bab::BabFilterResult> get_handedness(sheet::DirectedAdjacencyList const& adj,
                                                 bab::BabFilter & bab_filter);



/// @brief  Read the arguments and generate bab::BabFilter object, and then,
///         call print_handedness() to output the results.
/// @return If there is at least 1 connection to output.
bool output_handedness(std::ostream & os, sheet::DirectedAdjacencyList const& adj,
                       bpo::variables_map const& vm);

} // namespace rare

#endif // ifndef HANDEDNESS_H_
