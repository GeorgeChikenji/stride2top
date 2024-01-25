// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <numeric>

#include "sheet/sub_strands_range.h"

namespace sheet {


constexpr std::array<int, 2> SubStrandsRange::default_range;
constexpr IndexType SubStrandsRange::min_sub_str_len;


// ********************************************************************
// Public Member Function extend_substrand()
// ********************************************************************
void SubStrandsRange::extend_substrand(ZoneResidue const& res) {
  auto & sub_strand = data[res.serial_str_id].back();
  if (sub_strand == default_range) {
    sub_strand = std::array<int, 2>({{res.resnum, res.resnum}});
  } else {
    if (res.resnum < sub_strand[0]) {
      sub_strand[0] = res.resnum;
    } else if (sub_strand[1] < res.resnum) {
      sub_strand[1] = res.resnum;
    }
  }
} // public member function extend_substrand()




// ********************************************************************
// Public Member Function last_substr_id()
// ********************************************************************
SubStrandSet SubStrandsRange::cleanup_sheet() {
  SubStrandSet too_short_subs;

  IndexType serial_str_id = 0;
  for (auto & ranges : data) {

    // only for chenged substrands
    if (ranges.back() != default_range) {
      // if this sub-strand is too short (less than min_sub_str_len), remove it.
      if (ranges.back()[1] - ranges.back()[0] + 1 < min_sub_str_len) {
        too_short_subs.insert(SubStrand{serial_str_id,
                                        static_cast<IndexType>(ranges.size() - 1)});
        ranges.pop_back();
      }

      // for the next iteration
      ranges.push_back(default_range);
    }
    ++serial_str_id;
  }
  return too_short_subs;
} // public member function extend_substrand()



// ********************************************************************
// Public Member Function finish()
// ********************************************************************

std::unordered_map<SubStrand, SubStrand, SubStrandHasher> SubStrandsRange::finish() {
  std::vector<std::vector<IndexType>> changed_index_map;
  for (auto & ranges : data) {
    ranges.pop_back();

    auto const sorted_ids = sorted_indices(ranges);
    auto const ranges_cpy = ranges;

    std::transform(sorted_ids.cbegin(), sorted_ids.cend(), ranges.begin(),
                   [&ranges_cpy](auto const& a) { return ranges_cpy[a]; });
    changed_index_map.push_back(sorted_ids);
  }

  // Initialize index vector.
  init_index_vec();

  init_sub_strands_iters_vec();

  return convert_to_substrand(changed_index_map);
} // public member function finish()



// ********************************************************************
// Public Member Function last_substr_id()
// ********************************************************************

IndexType SubStrandsRange::last_substr_id(IndexType const str_id) const {
  if (data.at(str_id).size()) {
    return  data.at(str_id).size() - 1;
  }
  throw std::out_of_range{"No substrand: Serial Strand ID = " + std::to_string(str_id)};
} // public member function last_substr_id()




// ********************************************************************
// Public Member Function n_term_res()
// ********************************************************************

int SubStrandsRange::n_term_res(SubStrand const& s) const {
  return data.at(s.str).at(s.substr)[0];
} // public member function n_term_res()



// ********************************************************************
// Public Member Function c_term_res()
// ********************************************************************

int SubStrandsRange::c_term_res(SubStrand const& s) const {
  return data.at(s.str).at(s.substr)[1];
} // public member function term_res()




// ********************************************************************
// Public Member Function n_term_sub()
// ********************************************************************
SubStrand SubStrandsRange::n_term_sub(IndexType const str) const {
  auto const& iters = sub_strands_iters_vec.at(str);
  return iters.erased ? throw SubStrandErased{str, "In n_term_sub()"} : *iters.begin();
} // public member function n_term_sub()



// ********************************************************************
// Public Member Function c_term_sub()
// ********************************************************************
SubStrand SubStrandsRange::c_term_sub(IndexType const str) const {
  auto const& iters = sub_strands_iters_vec.at(str);
  return iters.erased ? throw SubStrandErased{str, "In c_term_sub()"} : *(iters.end() - 1);
} // public member function c_term_sub()




// ********************************************************************
// Public Member Function term_res()
// ********************************************************************

int SubStrandsRange::term_res(SubStrand const& s, IndexType const i) const {
  return data.at(s.str).at(s.substr).at(i);
} // public member function c_term_res()



#if defined(LOGGING) && defined(LOG_SHEET)
// ********************************************************************
// Public Member Function PrintSubStrands()
// ********************************************************************
void SubStrandsRange::print(std::vector<IndexType> const& strand_indices) const {
  std::clog << "\n";
  pdb::log("LIST OF SUB-STRANDS");

  auto const size = data.size();
  for (IndexType i = 0; i < size; ++i) {
    std::ostringstream iss;
    iss << "SSE[" << +strand_indices[i] << "] | ";
    for (auto const& sub_strand : data[i]) {
      iss << (boost::format("%4d ~ %4d | ") % sub_strand[0] % sub_strand[1]).str();
    }
    pdb::log(iss.str());
  }
  std::clog << std::endl;

} // public member function PrintSubStrands()

#endif // ifdef LOGGING




// ********************************************************************
// Protected Member Function sorted_indices()
// ********************************************************************
std::vector<IndexType> SubStrandsRange::sorted_indices(std::vector<std::array<int, 2>> const& substrs) {
  std::vector<IndexType> i(substrs.size());
  std::iota(i.begin(), i.end(), 0);
  std::sort(i.begin(), i.end(), [&substrs](auto const& a, auto const& b) {
                                    return substrs[a][0] < substrs[b][0]; });
  return i;
} // protected member function sorted_indices()


// ********************************************************************
// Protected Member Function convert_to_substrand()
// ********************************************************************
std::unordered_map<SubStrand, SubStrand, SubStrandHasher> SubStrandsRange::convert_to_substrand(std::vector<std::vector<IndexType>> const& changed_index_map) const {
  std::unordered_map<SubStrand, SubStrand, SubStrandHasher> ret;

  auto const n_str = changed_index_map.size();
  for (IndexType str = 0; str < n_str; ++str) {
    auto const n_substr = changed_index_map[str].size();
    for (IndexType substr = 0; substr < n_substr; ++substr) {
      if (substr != changed_index_map[str][substr]) {
        ret.insert({SubStrand{str, changed_index_map[str][substr]},
                    SubStrand{str, substr}});
      }
    }
  }

  return ret;
} // protected member function convert_to_substrand()





// ********************************************************************
// Protected Member Function init_index_vec()
// ********************************************************************

void SubStrandsRange::init_index_vec() {
  auto const n_str = data.size();
  index_vec.reserve(n_str);

  for (IndexType i = 0; i < n_str; ++i) {
    auto const n_substr = data[i].size();
    for (IndexType j = 0; j < n_substr; ++j) {
      index_vec.push_back(SubStrand{i, j});
    }
  }
  index_vec.shrink_to_fit();
} // protected member function init_index_vec()




// ********************************************************************
// Protected Member Function init_sub_strands_iters_vec()
// ********************************************************************

void SubStrandsRange::init_sub_strands_iters_vec() {
  auto const n_str = data.size();
  sub_strands_iters_vec.reserve(n_str);

  for (IndexType i = 0; i < n_str; ++i) {
    sub_strands_iters_vec.push_back(find_sub_strands(i));
  }
} // protected member function init_sub_strands_iters_vec()



// *******************************************************************************
// Protected Member Function find_sub_strands()
// *******************************************************************************

SubStrandsRange::SubStrandsIters SubStrandsRange::find_sub_strands(IndexType const serial_str_id) const {
  struct Comp {
    bool operator() (SubStrand const& s, IndexType const i) { return s.str < i; }
    bool operator() (IndexType const i, SubStrand const& s) { return i < s.str; }
  };
  auto const range = std::equal_range(index_vec.cbegin(), index_vec.cend(),
                                      serial_str_id, Comp());
  return SubStrandsIters{range.first, range.second};
} // protected memberr function find_sub_strands()



} // namespace sheet

