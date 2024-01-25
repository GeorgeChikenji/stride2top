// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_FILTER_H_
#define SHEET_FILTER_H_

#define PCC_MIN_ALLOWED_JUMP 1
#define APJ_MAX_ALLOWED_JUMP 1

#include <tuple>
#include <vector>

#include "sheet/directed_adjacency_list.h"


namespace sheet {

using SeqIter = std::vector<IndexType>::const_iterator;



/// @brief  Determine if the given 2 SSEs are Parallel or Anti-Parallel
///         considering reversing. BE CAREFUL WITH THE RETURN VALUE !!!!!!
/// @retval true  if Anti-Parallel
/// @retval false if Parallel
bool direction_with_reverse(bool const default_direction,
                            bool const rev_first, bool const rev_last);



/// @brief  Get a SubStrand of given iterator considering reverse.
///         Note that sequentially consecutive Sub-Strands must be on differenct sheets.
SubStrand get_substr(DirectedAdjacencyList const& adj, SeqIter const iter,
                     bool const first, bool const reversed);


using FilterFunc = std::function<bool(SeqIter const, bool const, bool const)>;


/// A Dummy Filtering Function Object
struct FalseFilter {
  bool operator()(SeqIter const, bool const, bool const) const {
    return false;
  }
};



/// A tuple that stores a pair of Strands and the number of jump.
using SheetFilterTarget = std::tuple<SubStrand, SubStrand, unsigned>;



class SheetFilterBase {
public:
  SheetFilterBase(DirectedAdjacencyList const& adj_): adj{adj_} {}

  virtual ~SheetFilterBase() = default;

  /// @brief  Return true if the consecutive 2 strands in the same sheet are parallel.
  /// @param  first  An iterator to the first position in the perm vector.
  bool operator()(SeqIter const first, bool const reversed_0, bool const reversed_1);

  DirectedAdjacencyList const& adj;

  /// Stores the last target information.
  SheetFilterTarget last_target{};

private:
  virtual bool actual(bool const direction, unsigned const jump) const = 0;

};




/// A struct to run a filter that eliminates
/// the Parallel Crossover Connections Between Sequential Strands
class PCCFilter: public SheetFilterBase {
public:
  // **************************************************************
  // Public Member Funtions
  // **************************************************************

  PCCFilter(DirectedAdjacencyList const& adj_,
            unsigned const min_allowed_jump_=PCC_MIN_ALLOWED_JUMP):
    SheetFilterBase{adj_}, min_allowed_jump{min_allowed_jump_} {}

  // **************************************************************
  // Public Member Variables
  // **************************************************************

  unsigned const min_allowed_jump{PCC_MIN_ALLOWED_JUMP};

private:

  // **************************************************************
  // Private Member Funtions
  // **************************************************************
  bool actual(bool const direction, unsigned const jump) const override;
};




class APJumpFilter : public SheetFilterBase {
public:
  // **************************************************************
  // Public Member Funtions
  // **************************************************************
  APJumpFilter(DirectedAdjacencyList const& adj_,
               unsigned const max_allowed_jump_=APJ_MAX_ALLOWED_JUMP):
    SheetFilterBase{adj_}, max_allowed_jump{max_allowed_jump_} {}

  // **************************************************************
  // Public Member Variables
  // **************************************************************
  unsigned const max_allowed_jump{APJ_MAX_ALLOWED_JUMP};

private:

  // **************************************************************
  // Private Member Funtions
  // **************************************************************
  bool actual(bool const direction, unsigned const jump) const override;
};


} // namespace sheet

#endif // #ifndef SHEET_FILTER_H_
