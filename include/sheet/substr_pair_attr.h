// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_SUBSTR_PAIR_ATTR_H_
#define SHEET_SUBSTR_PAIR_ATTR_H_

#include "sheet/adj_list_with_sub.h"

namespace sheet {

// ***********************************************************************************
// Struct StrandsPairAttribute
// ***********************************************************************************

struct StrandsPairAttribute {

  StrandsPairAttribute(SubStrand const ss0_, SubStrand const ss1_):
    ss0{ss0_}, ss1{ss1_} {}

  std::string string() const {
    return reachable ? ss0.string() + "\t" + ss1.string() + "\t" +
      std::to_string(jump) + "\t" + (direction ? "P\t" : "AP\t") : "";
  }

  SubStrand ss0;
  SubStrand ss1;

  /// True if the 2 sub-strands are on the same sheet. (Reachable on the graph)
  /// If false, any other attributes are nonsense.
  bool reachable{false};

  /// Number of other sub-strands between ss0 and ss1.
  /// 0 for adjacent sub-strands.
  unsigned jump{0};

  /// Relative direction of sub-strands,
  /// true for Parallel, false for Anti-Parallel.
  bool direction{true};

  /// A vector of jumped Sub-Strands.
  std::vector<SubStrand> jumped_substrs{};
};



} // namespace sheet

#endif // ifndef SHEET_SUBSTR_PAIR_ATTR_H_
