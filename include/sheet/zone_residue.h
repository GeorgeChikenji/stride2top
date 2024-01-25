// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_ZONE_RESIDUE_H_
#define SHEET_ZONE_RESIDUE_H_

#include <cassert>
#include <boost/functional/hash.hpp>

#include "pdb/sses.h"
#include "sheet/common.h"
#include "sheet/exceptions.h"


namespace sheet {

/// Structure that holds the SSE id and the residue number
struct ZoneResidue {
  ZoneResidue() = default;

  ZoneResidue(IndexType const sse_id_, int const resnum_, pdb::SSES const& sses):
    sse_id{sse_id_ == sses.sse_ind_of(resnum_, 'E', 0, true, 'A') ?
            sse_id_ : throw zone_info_failure{sse_id_, resnum_}},
    resnum{resnum_}, serial_str_id{sses.serial_strand_id[sse_id]},
    serial_res_id{0 <= resnum - sses.data[sse_id].init ?
                    static_cast<unsigned>(resnum - sses.data[sse_id].init) :
                    throw zone_info_failure{sse_id, resnum}
    }, has_value{true} {}

  bool operator==(ZoneResidue const& other) const {
    return sse_id == other.sse_id and resnum == other.resnum;
  }

  /// Increment the residue numbers while on the same sse.
  /// @return whether the increment was successful.
  bool increment(pdb::SSES const& sses) {
    if (sses.data[sse_id].end == resnum) {
      return false;
    }
    ++resnum;
    ++serial_res_id;
    return true;
  }


  bool decrement(pdb::SSES const& sses) {
    if (sses.data[sse_id].init == resnum) {
      return false;
    }
    assert(serial_res_id != 0);
    --resnum;
    --serial_res_id;
    return true;
  }


  IndexType sse_id{0};
  int resnum{0};

  IndexType serial_str_id{0};
  unsigned serial_res_id{0};

  /// true if initialized through the second version of constructor
  bool has_value{false};
};



/// KeyHasher for ZoneResidue class
struct ZoneResidueHasher {
  std::size_t operator()(ZoneResidue const& k) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, boost::hash_value(k.sse_id));
    boost::hash_combine(seed, boost::hash_value(k.resnum));
    return seed;
  }
};





} // namespace sheet

#endif // ifndef SHEET_ZONE_RESIDUE_H_
