// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef BAB_SEQ_BAB_PATTERN_H_
#define BAB_SEQ_BAB_PATTERN_H_

#include <string>
#include <vector>

#include <boost/range/irange.hpp>

#include "pdb/constants.h"

namespace bab {

/// @brief  Generate a vector of counterpart pos id.
/// @param  ss_seq must be rr_core::Seq.ss_seq

std::vector<std::vector<pdb::IndexType>> bab_id_pairs(std::string const& ss_seq);


} // namespace bab

#endif // ifndef BAB_SEQ_BAB_PATTERN_H_
