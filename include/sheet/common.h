// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_COMMON_H_
#define SHEET_COMMON_H_

#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>

#include "pdb/constants.h"

namespace sheet {

// IndexType is the same as defined in libpdb
using IndexType = pdb::IndexType;


/// Store the data for undirected_adjacency_list
struct AdjStrandData {
  AdjStrandData() : id{0}, direction{true}, count{0} {}
  AdjStrandData(IndexType const id_, bool const dir) : id{id_}, direction{dir}, count{0} {}
  AdjStrandData(IndexType const id_, bool const dir, unsigned const c):
    id{id_}, direction{dir}, count{c} {}
  IndexType id;
  bool direction;
  unsigned count;


  /// @brief ensure that the pairing and direction are the same.
  /// @return true if symmetric, false otherwise.
  bool symmetry(IndexType const other_id, AdjStrandData const& other) const {
    return id == other_id and direction == other.direction;
  }
};


/// AdjList[std::make_pair(i, j)] returns the AdjStrandData for Strand j which is adjacent to i
using AdjList = std::unordered_map<std::pair<IndexType, IndexType>,
                                   AdjStrandData,
                                   boost::hash<std::pair<IndexType, IndexType>>>;



} // namespace sheet

#endif // ifndef SHEET_COMMON_H_

