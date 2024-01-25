// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cstdlib>

#include "sheet/filter.h"
#include "bab/filter.h"
#include "bab/side.h"

namespace bab {


// ************************************************************************************
// Public Member Function operator()()
// ************************************************************************************

bool BabFilter::operator()(SeqIter const first, SeqIter const last, unsigned const reverse) {
  // if not strands
  if (sses[*first].type != 'E' or sses[*(last-1)].type != 'E') {

    #ifdef WITH_STAT
    last_result = BabFilterResult{};
    last_result.non_bab_reason = 1;
    #endif // WITH_STAT

    return false;
  }

  for (auto const& sub_first: adj.substrs().vec(sses.serial_strand_id[*first])) {
    for (auto const& sub_last: adj.substrs().vec(sses.serial_strand_id[*(last-1)])) {
      // If filtering failed.
      if (operator()(sub_first, sub_last, first, last, reverse)) {
        return true;
      }
    }
  }
  return false;
} // public member function operator()()




// ************************************************************************************
// Public Member Function operator()()
// ************************************************************************************

bool BabFilter::operator()(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1,
                       SeqIter const first, SeqIter const last,
                       unsigned const reverse) {
  #ifdef WITH_STAT
  last_result = BabFilterResult{ss0, ss1};
  #endif // WITH_STAT

  auto const cond = non_bab_condition(first, last);
  if (cond != 0) {

    #ifdef WITH_STAT
    last_result.non_bab_reason = cond;
    #endif // WITH_STAT

    return false;
  }

  // adj.attr will also check the reachability in inversed direction (ss1, ss0).
  auto const& attr = adj.attr(ss0, ss1);
  bool const reversed_first = reverse & (1u << *first);
  bool const reversed_last = reverse & (1u << *(last - 1));
  if (not attr.reachable
      or sheet::direction_with_reverse(attr.direction, reversed_first, reversed_last)) {

    #ifdef WITH_STAT
    last_result.non_bab_reason = 3;
    #endif // WITH_STAT

    return false;
  }
  auto const result = filter_one_unit(ss0, ss1, reversed_first, reversed_last, first, last);

  #ifdef WITH_STAT
  // If the last_result has any value (such as non_bab_reason), control never reaches here.
  last_result = result;
  last_result.jump = attr.jump;
  #endif // WITH_STAT

  // If success is false, comp need not be evaluated.
  return result.success and comp(result.left_score, cut_off_left_score);

} // public member function operator()()




// ************************************************************************************
// Protected Member Function init_sides_map()
// ************************************************************************************

SidesMap BabFilter::init_sides_map() const {
  SidesMap map;

  auto const n_substr = adj.substrs().vec().size();
  for (pdb::IndexType i = 0; i < n_substr; ++i) {
    auto const& substr0 = adj.substrs().vec()[i];
    for (pdb::IndexType j = i + 1; j < n_substr; ++j) {
      auto const& substr1 = adj.substrs().vec()[j];
      sheet::SubStrandsPairKey const key{substr0, substr1};

      if (not adj.adj_attr.at(key).reachable and
          not adj.adj_attr.at(key.reverse()).reachable) {
        continue;
      }

      auto sides_vec_0 = gen_sides_vec(substr0, substr1);
      auto sides_vec_1 = gen_sides_vec(substr1, substr0);

      auto const size_0 = sides_vec_0.size();
      auto const size_1 = sides_vec_1.size();

      // if the length of strands differs, shorten the longer sides_vec.
      if (size_0 != size_1) {
        auto & v_small = size_0 < size_1 ? sides_vec_0 : sides_vec_1;
        auto & v_big =  size_0 < size_1 ? sides_vec_1 : sides_vec_0;

        std::sort(v_big.begin(), v_big.end(), [](auto const& a, auto const& b) {
                                                return a.get_angle() > b.get_angle();
                                              });
        v_big.erase(v_big.begin() + v_small.size(), v_big.end());
      }

      map.insert({key, sides_vec_0});
      map.insert({key.reverse(), sides_vec_1});
    }
  }


  return map;
} // protected member function init_sides_map()




// *****************************************************************************
// Protected Member Function gen_sides_vec()
// *****************************************************************************

std::vector<Side> BabFilter::gen_sides_vec(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1) const {
  std::vector<Side> ret;

  auto const first = adj.atom_cbegin(ss0);
  // No triangle can be created from the last element.
  auto const last = adj.atom_cend(ss0) - 1;

  auto const opp_first = adj.atom_cbegin(ss1);
  auto const opp_last = adj.atom_cend(ss1);

  for (auto iter = first; iter != last; ++iter) {
    ret.push_back(Side{iter, iter + 1, opp_first, opp_last});
  }

  return ret;
} // protected member function gen_sides_vec()







// ************************************************************************************
// Protected Member Function each_atom()
// ************************************************************************************

std::tuple<unsigned, unsigned> BabFilter::count_left_tri(sheet::SubStrand const& b0,
                                                         sheet::SubStrand const& b1,
                                                         bool const b0_reverse,
                                                         bool const b1_reverse,
                                                         ATOM_vec_iter a_begin,
                                                         ATOM_vec_iter const a_end) const {
  unsigned total_counter = 0;
  unsigned left_counter = 0;
  auto const& sides_0 = sides_map.at(sheet::SubStrandsPairKey{b0, b1});
  auto const& sides_1 = sides_map.at(sheet::SubStrandsPairKey{b1, b0});

  auto const end_sides_0 = sides_0.cend();
  auto const end_sides_1 = sides_1.cend();

  for (auto& mid_a_itr = a_begin; mid_a_itr != a_end; ++mid_a_itr) {
    // if this is a padding atom.
    if (not mid_a_itr->pdb) {
      continue;
    }

    for (auto itr = sides_0.cbegin(); itr != end_sides_0; ++itr) {
      bool on_left, is_distant;
      std::tie(on_left, is_distant) = itr->on_left_side(mid_a_itr->xyz, b0_reverse, true,
                                                        cut_off_side_min_dist);
      if (is_distant) {
        ++total_counter;
        if (on_left) {
          ++left_counter;
        }
      }
    }

    for (auto itr = sides_1.cbegin(); itr != end_sides_1; ++itr) {
      bool on_left, is_distant;
      std::tie(on_left, is_distant) = itr->on_left_side(mid_a_itr->xyz, b1_reverse, false,
                                                        cut_off_side_min_dist);
      if (is_distant) {
        ++total_counter;
        if (on_left) {
          ++left_counter;
        }
      }
    }
  }

  return std::make_tuple(left_counter, total_counter);
}



// ************************************************************************************
// Protected Member Function filter_one_unit()
// ************************************************************************************

BabFilterResult BabFilter::filter_one_unit(sheet::SubStrand const& b0,
                                           sheet::SubStrand const& b1,
                                           bool const b0_reverse, bool const b1_reverse,
                                           SeqIter const first, SeqIter const last) {
  BabFilterResult result{b0, b1};

  #ifdef WITH_LOOP
  // the first loop right after the first strand.
  auto const& first_loop = sses.loop(*first);
  auto const counts_l = count_left_tri(b0, b1, b0_reverse, b1_reverse,
                                     first_loop.atoms.cbegin(), first_loop.atoms.cend());
  result.left_score += std::get<0>(counts_l);
  result.tri_atom_count += std::get<1>(counts_l);
  result.n_pdb_atoms += first_loop.n_pdb;

  result.mid_res_len += first_loop.atoms.size();
  // ON: Loop Bit, if loop is not empty
  if (first_loop.n_pdb) {
    result.connection_type |= 2u;
  }
  #endif // WITH_LOOP

  // past the end iterator to the sse_id of mid-SSEs (pointing to the last Strand)
  auto const mid_sse_end = last - 1;
  for (auto sse_id_itr = first + 1; sse_id_itr != mid_sse_end; ++sse_id_itr) {

    result.mid_res_len += sses[*sse_id_itr].atoms.size();
    // If the total residue length of mid-part exceeds the limit
    if (cut_off_res_len < result.mid_res_len) {

      #ifdef WITH_STAT
      result.non_bab_reason = 4;
      #endif // WITH_STAT

      return result;
    }

    // Filter SSEs
    auto const& target_sse = sses[*sse_id_itr];
    // For Helices
    if (target_sse.type == 'H') {
      auto const counts_h = count_left_tri(b0, b1, b0_reverse, b1_reverse,
                                           target_sse.atoms.cbegin(),
                                           target_sse.atoms.cend());
      result.left_score += std::get<0>(counts_h);
      result.tri_atom_count += std::get<1>(counts_h);
      result.n_pdb_atoms += target_sse.n_pdb;

      // ON: Helix Bit, if this Helix has any tri_atom.
      if (std::get<1>(counts_h)) {
        result.connection_type |= 1u;
      }

    // For Strands
    } else if (target_sse.type == 'E') {
      for (auto const& ss : adj.substrs().vec(sses.serial_strand_id[*sse_id_itr])) {
        // skip the strands in the same sheet
        if (adj.attr(b0, ss).reachable) {
          ++result.n_mid_str;

          if (cut_off_mid_str < result.n_mid_str) {

            #ifdef WITH_STAT
            result.non_bab_reason = 5;
            #endif // WITH_STAT

            return result;
          }
          continue;
        }

        // the last arguemnt to count_left_tri() must be a past the end iterator
        auto const counts_s = count_left_tri(b0, b1, b0_reverse, b1_reverse,
                               target_sse.atom_vec_iter(adj.substrs().n_term_res(ss)),
                               target_sse.atom_vec_iter(adj.substrs().c_term_res(ss))+1);
        result.left_score += std::get<0>(counts_s);
        result.tri_atom_count += std::get<1>(counts_s);
        result.n_pdb_atoms += target_sse.n_pdb;

        // ON: Strand Bit, if this Strand has any tri_atom.
        if (std::get<1>(counts_s)) {
          result.connection_type |= 4u;
        }
      }
    }


    #ifdef WITH_LOOP
    // Filter LOOPs
    auto const& target_loop = sses.loop(*sse_id_itr);
    auto const counts = count_left_tri(b0, b1, b0_reverse, b1_reverse,
                                       target_loop.atoms.cbegin(),
                                       target_loop.atoms.cend());
    result.left_score += std::get<0>(counts);
    result.tri_atom_count += std::get<1>(counts);
    result.n_pdb_atoms += target_loop.n_pdb;
    result.mid_res_len += target_loop.atoms.size();
    #endif // WITH_LOOP
  }

  // Mid-res-len check
  if (cut_off_res_len < result.mid_res_len) {

    #ifdef WITH_STAT
    result.non_bab_reason = 4;
    #endif // WITH_STAT

    return result;
  }

  // Divide the sum of the number of atoms on the left-handed side
  // by the total number of atoms times the total number of triangles.
  // If all the CA atoms are at the middle of right and left region, score will be 0;
  if (result.tri_atom_count != 0) {
    result.left_score /= result.tri_atom_count;
  } else {
    result.left_score = 0.0;
  }
  result.success = true;

  return result;
} // protected member function filter_one_unit()




// ************************************************************************************
// Protected Member Function bab_condition()
// ************************************************************************************

unsigned BabFilter::non_bab_condition(SeqIter first, SeqIter last) const {
  auto const dist = std::distance(first, last);
  switch (dist) {
    // if the first == last or only first sse.
    case 0:
    case 1: return 2;

    #ifdef WITH_LOOP
    // loop between *first and *(last-1)
    case 2: return 0;

    #else
    case 2: return 2;
    #endif // WITH_LOOP

    default: return 0;
  }
} // protected member function bab_condition()

} // namespace bab

