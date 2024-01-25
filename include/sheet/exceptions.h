// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEET_EXCEPTIONS_H_
#define SHEET_EXCEPTIONS_H_

#include <boost/format.hpp>

#include "pdb/exceptions.h"
#include "sheet/common.h"

namespace sheet {



// *******************************************************************
// Exception class non_symmetric
// *******************************************************************

class non_symmetric: public pdb::fatal_error_base {
public:
  non_symmetric(IndexType const sse_id_0, IndexType const sse_id_1):
    pdb::fatal_error_base("NOKEY: undirected_adj_list is not symmetric. (SSE[" +
                          std::to_string(sse_id_0) + "], SSE[" +
                          std::to_string(sse_id_1) + "]: CONNECTION DOESN'T EXIST.)") {}
};




// *******************************************************************
// Exception class non_symmetric_dir_opposite
// *******************************************************************

class non_symmetric_dir_opposite: public pdb::fatal_error_base {
public:
  non_symmetric_dir_opposite(IndexType const sse_id_0, IndexType const sse_id_1):
    pdb::fatal_error_base("NOT EQUAL: undirected_adj_list is not symmetric. (SSE[" +
                          std::to_string(sse_id_0) + "], SSE[" +
                          std::to_string(sse_id_1) + "]: DIRECTION IS OPPOSITE.)") {}
};



// *******************************************************************
// Exception class third_pair_found
// *******************************************************************

class third_pair_found: public pdb::fatal_error_base {
public:
  third_pair_found(IndexType const sse_id, int const resnum,
                   IndexType const new_sse_id, int const new_resnum,
                   IndexType const current_sse_id_0, int const current_resnum0,
                   IndexType const current_sse_id_1, int const current_resnum1):
    pdb::fatal_error_base{(boost::format("THIRD PAIR FOUND: on SSES[%d] RESNUM = %d to [%d] %d"
                                         ", but already has PAIR to [%d] %d and [%d] %d")
                                          % +sse_id
                                          % resnum
                                          % +new_sse_id
                                          % new_resnum
                                          % +current_sse_id_0
                                          % current_resnum0
                                          % +current_sse_id_1
                                          % current_resnum1).str()} {}
};



// *******************************************************************
// Exception class zone_info_failure
// *******************************************************************
//
class zone_info_failure : public pdb::fatal_error_base {
public:
  zone_info_failure(IndexType const sse_id, int const resnum):
    pdb::fatal_error_base((boost::format("CONSTRUCTION OF ZONE INFO FAILED: "
                                         "RESNUM %d NOT IN SSE[%d]")
                                          % resnum
                                          % sse_id).str()) {}
};


// *******************************************************************
// Exception class direction_failure
// *******************************************************************
class direction_failure: public pdb::fatal_error_base {
public:
  direction_failure(IndexType const sse0, IndexType const sse1, IndexType const count0, IndexType const count1):
    pdb::fatal_error_base((boost::format("CANNOT DECIDE DIRECTION BETWEEN SSE[%2d] (COUNT = %3d) "
                                         "AND SSE[%2d] (COUNT = %3d)")
                                         % sse0 % count0 % sse1 % count1
                                         ).str()) {}
};




// *******************************************************************
// Exception class substrand_cleanup_failure
// *******************************************************************
class substrand_cleanup_failure: public pdb::fatal_error_base {
public:
  substrand_cleanup_failure(IndexType const str0, IndexType const sstr0, unsigned const count0,
                            IndexType const str1, IndexType const sstr1, unsigned const count1):
    pdb::fatal_error_base((boost::format("CANNOT SET SUB_STRAND PAIR SSE[%d]_%d and SSE[%d]_%d "
                                         "(COUNT = %d and %d)")
                                         % +str0 % +sstr0 % +str1 % +sstr1
                                         % count0 % count1
                                        ).str()) {}
};



// *******************************************************************
// Exception class TargetRelativeDirectionNotSet
// *******************************************************************
class TargetRelativeDirectionNotSet: public pdb::fatal_error_base {
public:
  TargetRelativeDirectionNotSet(IndexType const sse_id):
    pdb::fatal_error_base("RELATIVE DIRECTION NOT SET: SSES[" +
                          std::to_string(sse_id) + "]") {}
};



// *******************************************************************
// Exception class PairedResidueNotFound
// *******************************************************************
class PairedResidueNotFound: public pdb::fatal_error_base {
public:
  PairedResidueNotFound(IndexType const sse_id, int const resnum,
                        IndexType const target_sse_id, IndexType const target_substr_id):
    pdb::fatal_error_base("GAP CALC: PAIRED RESIDUE NOT FOUND: SSE[" +
                          std::to_string(sse_id) + "], RESNUM = " + std::to_string(resnum)
                          + ", TARGET = SSE[" + std::to_string(target_sse_id) +
                          "], SUBSTR = " + std::to_string(target_substr_id)) {}
};



// *******************************************************************
// Exception class AdjacentSubStrandNotFound
// *******************************************************************

class AdjacentSubStrandNotFound: pdb::warning_base {
public:
  explicit AdjacentSubStrandNotFound(std::string const& msg_):
    pdb::warning_base{"Adjacent Sub-Strand Not Found: " + msg_} {}
};

} // namespace sheet



class SubStrandErased: public pdb::warning_base {
public:
  explicit SubStrandErased(std::size_t const sse_id, std::string const& msg_=""):
    pdb::warning_base{"Ignoring access for erased SubStrand. Serial Strand ID = " +
                      std::to_string(sse_id) + " :" + msg_} {}
};


#endif // ifndef SHEET_EXCEPTIONS_H_
