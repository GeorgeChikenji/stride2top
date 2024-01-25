// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_FIND_CYCLE_H_
#define SHEET_FIND_CYCLE_H_

#include <array>
#include <vector>

#include "sheet/adj_list_with_sub.h"


namespace sheet {


/// Rotate the given range so that the smallest element will come to the first position.
template<typename Iter>
decltype(auto) rotate_to_smallest(Iter first, Iter last) {
  auto const min = std::min_element(first, last);
  std::vector<typename Iter::value_type> ret(min, last);
  ret.reserve(std::distance(first, last));
  ret.insert(ret.end(), first, min);
  return ret;
}


/// Generate a reversed directional path with its smallest element palaced
/// at the first position of the returned vector.
template<typename Container>
Container invert(Container const& path) {
  return rotate_to_smallest(path.crbegin(), path.crend());
}



/// Find all the cycles in the given graph-edges
class FindCycle {
public:

  // *******************************************************************
  // Public Member Functions
  // *******************************************************************

  FindCycle(SubStrandsPairKeyVec const& key_vec):
    orig_edges{key_vec}, nodes{}, edges{init_edges()}, cycles{init_cycles(init_hidden_cycles())} {}


protected:

  // *******************************************************************
  // Protected Member Variables
  // *******************************************************************

  /// Reference to a vector, of vectors of keys in the same sheet
  SubStrandsPairKeyVec const& orig_edges;
  std::vector<SubStrand> nodes;

  /// Odd index stores the index to the first element of the edge, and even to the second
  std::vector<unsigned> edges;

public:
  std::vector<std::vector<SubStrand>> cycles;

protected:

  // *******************************************************************
  // Protected Member Functions
  // *******************************************************************

  /// Initialize \c edges .
  /// If a non-registered node is found, add it in \c nodes .
  std::vector<unsigned> init_edges();

  /// @brief  Check if the node is already registered
  /// @return The index of the found element.
  /// @retval -1 If not found.
  int has_node(IndexType const str, IndexType const substr) const;


  std::vector<std::vector<unsigned>> init_hidden_cycles() const;


  void find_new_cycles(std::vector<unsigned> const& path,
                       std::vector<std::vector<unsigned>> & cycles) const;


  /// Translatee hidden_cycles into cycles
  std::vector<std::vector<SubStrand>> init_cycles(std::vector<std::vector<unsigned>> const& hidden) const;
};


} // namespace sheet

#endif // ifndef SHEET_FIND_CYCLE_H_

