// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <tuple>
#include <vector>

#include "sheet/directed_adjacency_list.h"
#include "sheet/filter.h"
#include "sheet/sub_strands_range.h"

namespace sheet {


// ************************************************************************************
// Function direction_with_reverse()
// ************************************************************************************

/// BE CAREFUL WITH THE RETURN VALUE !!!!!!
bool direction_with_reverse(bool const default_direction,
                            bool const rev_first, bool const rev_last) {
  // default_direction    reversed_first    reversed_last   real_direction
  //
  // P                    True              True            P
  // P                    True              False           AP
  // P                    False             True            AP
  // P                    False             False           P
  // AP                   True              True            AP
  // AP                   True              False           P
  // AP                   False             True            P
  // AP                   False             False           AP
  //
  // Whether the direction considering reverses (real_direction)
  // can be determined by the following condition.
  // default_direction == (reversed_first == reversed_last)
  // In current situation, continue the for loop when real_direction is AP.
  // default_direction != (reversed_first == reversed_last)
  //

  return default_direction != (rev_first == rev_last);
} // protected member function direction_with_reverse()




// ************************************************************************************
// Function get_substr()
// ************************************************************************************

SubStrand get_substr(DirectedAdjacencyList const& adj, SeqIter const iter,
                     bool const first, bool const reversed) {
  auto const serial_str_id = adj.sses.serial_strand_id[*iter];
  //
  // first  reversed  return
  // false  false     n_term
  // false  true      c_term
  // true   false     c_term
  // ture   true      n_term
  //
  return first == reversed ? adj.substrs().n_term_sub(serial_str_id) :
                             adj.substrs().c_term_sub(serial_str_id);
}

// ************************************************************************************
// Public Member Operator ()
// ************************************************************************************

bool SheetFilterBase::operator()(SeqIter const first,
                                 bool const reversed_0, bool const reversed_1) {
  SubStrand ss0, ss1;

  try {
    ss0 = get_substr(adj, first, true, reversed_0);
    ss1 = get_substr(adj, first+1, false, reversed_1);

  // if at least one of the 2 strand doesn't have sufficient hbonds and erased
  } catch (SubStrandErased const& e) {
    return false;
  }

  auto const& attr = adj.attr(ss0, ss1);
  if (not attr.reachable) { return false; }

  // Save the target of this filtering
  last_target = std::make_tuple(ss0, ss1, attr.jump);

  return actual(direction_with_reverse(attr.direction, reversed_0, reversed_1),
                attr.jump);
} // public member operator ()



// ************************************************************************************
// Private Member Funcion actual()
// ************************************************************************************
bool PCCFilter::actual(bool const direction, unsigned const jump) const {
  return ! direction and jump < min_allowed_jump;
} // private member function actual()



// ************************************************************************************
// Private Member Funcion actual()
// ************************************************************************************
bool APJumpFilter::actual(bool const direction, unsigned const jump) const {
  return direction and max_allowed_jump < jump;
} // private member function actual()

} // namespace sheet

