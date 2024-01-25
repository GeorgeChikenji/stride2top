// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "pdb/tools.h"
#include "bab/seq_bab_pattern.h"



namespace bab {

std::vector<std::vector<pdb::IndexType>> bab_id_pairs(std::string const& ss_seq) {
  std::vector<std::vector<pdb::IndexType>> ret;
  ret.reserve(ss_seq.size());

  for (auto const i : pdb::range(ss_seq.size())) {
    if (ss_seq[i] == 'H') {
      ret.push_back(std::vector<pdb::IndexType>{});
    } else if (ss_seq[i] == 'E') {
      
      std::vector<pdb::IndexType> v{};
      for (auto const j : pdb::range(i)) {
        if (ss_seq[j] == 'E') {
          bool found_helix = false;
          for (auto k = j; k < i; ++k) {
            if (ss_seq[k] == 'H') {
              found_helix = true;
            }
          }
          if (found_helix) {
            v.push_back(static_cast<pdb::IndexType>(j));
          }
        }
      }
      ret.push_back(v);
    }
  }
  return ret;
}

} // namespace bab
