// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "pdb/constants.h"
#include "pdb/atom.h"
#include "pdb/sse.h"

namespace pdb {

// *********************************************************************
// Public Member Function atom_vec_iter()
// *********************************************************************

ATOM_vec_iter SSE_Base::atom_vec_iter(int const resnum) const {
  auto const idx = resnum - init;
  if (idx < 0 or atoms.size() <= static_cast<unsigned>(idx)) {
    throw resnum_out_of_range{resnum};
  }
  return atoms.cbegin() + idx;
} // public member function atom_vec_iter().



// *********************************************************************
// Protected Member Function entirety_check()
// *********************************************************************

bool SSE_Base::entirety_check(std::unordered_map<int, std::string> const& ca_lines) const {
  for (int i = init; i <= end; ++i) {
    if (ca_lines.count(i) == 0) {
      log("RESNUM '" + std::to_string(i) + "' IS MISSING IN SSE or LOOP: index = " +
          std::to_string(index));
      return false;
    }
  }
  return true;
} // protected member function entirety_check()



// *********************************************************************
// Protected Member Function read_pdb_atom()
// *********************************************************************

std::vector<ATOM> SSE_Base::read_pdb_atom(std::unordered_map<int, std::string> const& ca_lines) const {
  std::vector<ATOM> tmp_atoms;
  tmp_atoms.reserve(static_cast<unsigned>(end - init + 1));

  for (int i = init; i <= end; ++i) {
    if (ca_lines.count(i) == 0) {
      tmp_atoms.push_back(ATOM{});
    } else {
      tmp_atoms.push_back(ATOM{ca_lines.at(i)});
    }
  }

  return tmp_atoms;
} // protected member function SSE_Base::read_pdb_atom()



// *********************************************************************
// Protected Member Function count_real_atoms()
// *********************************************************************

unsigned SSE_Base::count_real_atoms() const {
  unsigned count = 0;
  for (auto const& atom : atoms) {
    if (atom.pdb) {
      ++count;
    }
  }
  return count;
} // protected member function count_real_atoms()




// *********************************************************************
// Public Member Function rep functions
// *********************************************************************


Eigen::Vector3d const& SSE::rep_outer_head(unsigned const reverse) const {
  return reverse & index_bit ? rep_atoms[1][1].xyz : rep_atoms[0][0].xyz;
}
Eigen::Vector3d const& SSE::rep_inner_head(unsigned const reverse) const {
  return reverse & index_bit ? rep_atoms[1][0].xyz : rep_atoms[0][1].xyz;
}
Eigen::Vector3d const& SSE::rep_inner_tail(unsigned const reverse) const {
  return reverse & index_bit ? rep_atoms[0][1].xyz : rep_atoms[1][0].xyz;
}
Eigen::Vector3d const& SSE::rep_outer_tail(unsigned const reverse) const {
  return reverse & index_bit ? rep_atoms[0][0].xyz : rep_atoms[1][1].xyz;
}



// *********************************************************************
// Public Member Function stamp()
// *********************************************************************

std::string SSE::stamp(SSE const& other, unsigned const reverse, bool const alt) const {
  std::string ret(2, '\0');
  if (alt) {
    ret[1] = !(reverse & index_bit) ? (index | 128u) : index;
    ret[0] = !(reverse & other.index_bit) ? (other.index | 128u) : other.index;
  } else {
    ret[0] = (reverse & index_bit) ? (index | 128u) : index;
    ret[1] = (reverse & other.index_bit) ? (other.index | 128u) : other.index;
  }
  return ret;
} // public member function SSE::stamp()



// *********************************************************************
// Public Member Function in_range()
// *********************************************************************

bool SSE::in_range(int const resnum, int const offset, char const HAtom) const noexcept {
  int const first = HAtom != 'N' ? init - offset : init;
  int const last = HAtom != 'C' ? end + offset : end;
  return first <= resnum and resnum <= last;
} // public member function SSE::in_range()




// *********************************************************************
// Public Member Function distance()
// *********************************************************************

double SSE::distance(SSE const& other) const {
  double dist = std::numeric_limits<double>::max();

  for (auto const& a_this : atoms) {
    for (auto const& a_other : other.atoms) {
      auto const tmp_dist = (a_this.xyz - a_other.xyz).norm();
      if (tmp_dist < dist) {
        dist = tmp_dist;
      }
    }
  }

  return dist;
} // public member function SSE::distance()


// *********************************************************************
// Protected Member Function too_short_check()
// *********************************************************************

bool SSE::too_short_check() const noexcept {
  return type == 'H' ?
           end - init + 1 < HELIX_MIN_LEN :
           end - init + 1 < STRAND_MIN_LEN;
} // protected member function SSE::too_short_check()


// *********************************************************************
// Protected Member Function gen_representative()
// *********************************************************************

std::array<ATOM, 2> SSE::gen_representative(int const first) const {
  if (rep.len + rep.intrvl <= static_cast<unsigned>(end - init + 1)) {
    try {
      return std::array<ATOM, 2>({{gen_representative_atom(first),
                                   gen_representative_atom(first + rep.intrvl)}});
    } catch (padding_atom_found const& p) {
      log(p.what());
      return std::array<ATOM, 2>();
    }
  } else {
    return std::array<ATOM, 2>();
  }
} // protected member function SSE::gen_representative()



// *********************************************************************
// Protected Member Function gen_representative_atom()
// *********************************************************************

ATOM SSE::gen_representative_atom(int const first) const {
  Eigen::Vector3d rep_atom(0.0, 0.0, 0.0);

  for (unsigned i = 0; i < rep.len; ++i) {
    if (atoms[first+i].pdb ==false) {
      throw padding_atom_found{};
    }

    rep_atom += atoms[first+i].xyz * rep.coeff[i];
  }
  return ATOM(rep_atom/rep.div);
} // protected member function SSE::gen_representative_atom()



// *********************************************************************
// Protected Member Function rep_check()
// *********************************************************************

bool SSE::rep_check() const {
  return rep_atoms[0] != std::array<ATOM, 2>{} and rep_atoms[1] != std::array<ATOM, 2>{};
} // protected member function SSE::valid_check()





} // namespace pdb
