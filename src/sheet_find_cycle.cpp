// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cassert>
#include <algorithm>

#include "sheet/find_cycle.h"

namespace sheet {


std::vector<unsigned> FindCycle::init_edges() {
  std::vector<unsigned> ret;

  for (auto const& key : orig_edges) {
    int const id0 = has_node(key.str0, key.substr0);
    if (id0 == -1) {
      nodes.push_back(key.sub0());
      ret.push_back(nodes.size()-1);
    } else {
      ret.push_back(id0);
    }

    int const id1 = has_node(key.str1, key.substr1);
    if (id1 == -1) {
      nodes.push_back(key.sub1());
      ret.push_back(nodes.size()-1);
    } else {
      ret.push_back(id1);
    }
  }

  assert(ret.size() == orig_edges.size() * 2);
  return ret;
}



int FindCycle::has_node(IndexType const str, IndexType const substr) const {
  auto const itr = std::find_if(nodes.cbegin(), nodes.cend(),
                                [str, substr] (auto const& a) {
                                  return str == a.str and substr == a.substr;
                                });
  if (itr == nodes.cend()) {
    return -1;
  }
  return std::distance(nodes.cbegin(), itr);
}




std::vector<std::vector<unsigned>> FindCycle::init_hidden_cycles() const {
  std::vector<std::vector<unsigned>> ret;
  for (auto const edge : edges) {
    find_new_cycles(std::vector<unsigned>{edge}, ret);
  }
  return ret;
}




void FindCycle::find_new_cycles(std::vector<unsigned> const& path,
                                std::vector<std::vector<unsigned>> & cycles_) const {

  auto const size = edges.size();
  for (unsigned i = 0; i < size; i += 2) {
    if (path[0] == edges[i]) {
      auto const next_node = edges[i+1];
      if (std::find(path.cbegin(), path.cend(), next_node) == path.cend()) {
        std::vector<unsigned> sub{next_node};
        sub.reserve(path.size() + 1);
        sub.insert(sub.end(), path.cbegin(), path.cend());
        find_new_cycles(sub, cycles_);
      } else if (path.size() > 2 and next_node == path.back()) {
        auto const p = rotate_to_smallest(path.cbegin(), path.cend());
        auto const inv = invert(p);
        if (std::find(cycles_.cbegin(), cycles_.cend(), p) == cycles_.cend() and
            std::find(cycles_.cbegin(), cycles_.cend(), inv) == cycles_.cend()) {
          cycles_.push_back(inv);
        }
      }
    }
  }
}



std::vector<std::vector<SubStrand>> FindCycle::init_cycles(std::vector<std::vector<unsigned>> const& hidden) const {
  std::vector<std::vector<SubStrand>> ret;
  ret.reserve(hidden.size());

  auto const& node_vec = nodes;

  // tranform a vector of the index of nodes (SubStrands) into a vector of SubStrands
  for (auto const& one_cycle : hidden) {
    std::vector<SubStrand> tmp(one_cycle.size());
    std::transform(one_cycle.cbegin(), one_cycle.cend(), tmp.begin(),
                   [&node_vec](auto const a) {return node_vec[a];});
    ret.push_back(rotate_to_smallest(tmp.cbegin(), tmp.cend()));
  }
  return ret;
}


} // namespace sheet

