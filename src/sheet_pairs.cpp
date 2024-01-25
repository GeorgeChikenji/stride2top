// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cassert>
#include <cstdio>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "pdb/sses.h"
#include "pdb/tools.h"
#include "sheet/common.h"
#include "sheet/pairs.h"
#include "sheet/exceptions.h"

namespace sheet {

// **************************************************************************************
// Public Member function Pairs::resort_involved_pairs()
// **************************************************************************************


Pairs::PairsVec Pairs::resort_involved_pairs(IndexType const serial_str_id,
                                             AdjList const& adj_list, pdb::SSES const& sses) const {
  PairsVec ret = involved_pairs[serial_str_id];

  auto const sse_id = sses.gen_index_vec('E', true)[serial_str_id];

  std::sort(ret.begin(), ret.end(), [sse_id, &adj_list, &sses](auto const& a, auto const& b) {
    if (a[0] < b[0]) { return true; }
    if (a[0] == b[0]) {
      try {
        auto const sse_id_a = sses.sse_ind_of(a[1], 'E', 1, true, a[2] == 0 ? 'C' : 'N');
        auto const sse_id_b = sses.sse_ind_of(b[1], 'E', 1, true, b[2] == 0 ? 'C' : 'N');

        if (sse_id_a == sse_id_b) {
          if (adj_list.count({sse_id, sse_id_a}) == 0) {
            return false;
          }
          bool const dir = adj_list.at({sse_id, sse_id_a}).direction;
          // if parallel
          if (dir) {
            if (a[1] < b[1]) { return true; }

          // if anti-parallel
          } else {
            if (b[1] < a[1]) { return true; }
          }
          if (a[1] == b[1]) {
            return a[2] < b[2];
          }
        }
      } catch (pdb::non_sse_resnum) {
        return false;
      }
    }
    return false;
  });

  return ret;
} // public member function resort_involved_pairs()




// **************************************************************************************
// Protected Member function Pairs::read_stride_stream()
// **************************************************************************************
Pairs::PairsVec Pairs::read_stride_stream(std::istream & is) const {
  is.clear();
  is.seekg(0);

  PairsVec ret;

  for (std::string buff; std::getline(is, buff); ) {
    if (buff.substr(0, 3) == "DNR") {
      ret.push_back({{std::stoi(buff.substr(11, 4)),
                      std::stoi(buff.substr(31, 4))}});
    }
  }
  return ret;
} // private member function read_stride_stream()



// **************************************************************************************
// Protected Member function Pairs::init_involved_pairs()
// **************************************************************************************
std::vector<Pairs::PairsVec> Pairs::init_involved_pairs(PairsVec const& dnr_, pdb::SSES const& sses) const {

  std::vector<Pairs::PairsVec> involved_pairs_vec;

  for (auto const& sse_index : sses.gen_index_vec('E', true)) {
    involved_pairs_vec.push_back(involve_with(dnr_, sses[sse_index]));
  }
  return involved_pairs_vec;
} // protected member function init_involved_pairs()



// **************************************************************************************
// Protected Member function Pairs::involve_with()
// **************************************************************************************
Pairs::PairsVec Pairs::involve_with(PairsVec const& dnr_, pdb::SSE const& e, int const offset) const {
  PairsVec ret;
  for (auto const& pair : dnr_) {
    if (e.in_range(pair[0], offset, 'N')) {
      ret.push_back({{pair[0], pair[1], 0}});
    } else if (e.in_range(pair[1], offset, 'C')) {
      ret.push_back({{pair[1], pair[0], 1}});
    }
  }
  std::sort(ret.begin(), ret.end(), [](auto const& a, auto const& b)
      {
        if (a[0] < b[0]) {return true;}
        if (a[0] == b[0]) {
          if (a[2] < b[2]) { return true; }
          if (a[2] == b[2]) { return a[1] < b[1]; }
        }
        return false;
      });
  return ret;
} // protected member function involve_with()


} // namespace rperm

