// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef RPERM_PDB_SSES_H_
#define RPERM_PDB_SSES_H_

#include <string>
#include <vector>

#include "pdb/constants.h"
#include "pdb/loop.h"
#include "pdb/sse.h"
#include "pdb/stride_stream.h"

namespace pdb {

/// Compare the range of an SSE and the resnum
class ResnumSSEComp {
public:
  ResnumSSEComp(std::vector<SSE> const& data_, int const n_offset_, int const c_offset_):
    data{data_}, n_offset{n_offset_}, c_offset{c_offset_} {}

  bool operator() (IndexType const id, int const resnum) const {
    return data[id].end + c_offset < resnum;
  }

  bool operator() (int const resnum, IndexType const id) const {
    return resnum < data[id].init - n_offset;
  }

protected:
  std::vector<SSE> const& data;
  int const n_offset;
  int const c_offset;
};




class SSES {
public:

  // *********************************************************************************
  // Public Member Functions
  // *********************************************************************************
  SSES() = default;

  explicit SSES(std::string const& pdb_file_, stride_stream & stride=empty_stride_stream) :
    pdb_file{pdb_file_},
    #if defined(DRYRUN) || defined(NATIVE_DRYRUN)
    pdb_basename{basename(pdb_file)},
    answer{""},
    #endif // DRYRUN

    #ifdef WITH_LOOP
    loops_data{},
    #endif // WITH_LOOP

    data{read_pdb(stride)},
    size{data.size()},
    index_vec{init_index_vec()},
    serial_strand_id{init_serial_strand_id()}
  {
    #ifdef NATIVE_DRYRUN
    answer = "";
    for (IndexType i = 0; i < size; ++i) {
      answer += "+" + std::to_string(i);
    }
    #endif //ifdef NATIVE_DRYRUN
  }


  /// Member Access Operator
  SSE const& operator[](IndexType const i) const {
    return data[i];
  }

  /// @brief  Search an SSE which contains the ATOM whose residue number is resnum
  /// @return SSE_index an index of this SSE in vector<SSE> (including helices)
  /// @param  resnum         A residue number to look for an sse index.
  /// @param  type           Specify the target SSE type to search. \c 'A' for all SSE
  ///                        types (both \c 'H' and \c 'E'), \c 'H' for helices,
  ///                        and \c 'E' for strands.
  /// @param  offset         The number of residues to extend from the default range.
  ///                        Specifying the offset value of higher than 1 is not supported.
  /// @param  with_too_short If true, include too short SSEs in the search targets.
  /// @param  hbond_atom     \c 'N' to prefer an SSE on the Nitrogen side of the \s resnum,
  ///                        \c 'C' for Carbon.\n
  ///                        'A' for all types (Any charactor other than those above will
  ///                        be treated as 'A'.).
  /// @throw  pdb::unkown_sse_type If parameter 'type' is none of {'H', 'E', or 'A'}.
  /// @throw  pdb::non_sse_resnum   If parameter 'resnum' is not found in any SSEs
  /// @sa     pdb::SSE::in_range
  IndexType sse_ind_of(int const resnum, char const type, int const offset,
                    bool const with_too_short, char const hbond_atom) const;


  /// @brief  Generate a vector of \c IndexType s for the specified type
  /// @param type             Type of SSEs to generate the index vector. Must be one of
  ///                         <tt> {'H', 'E', 'A'}</tt>.
  /// @param with_too_short   whether to include the too_short SSEs in output
  std::vector<IndexType> const& gen_index_vec(char const type, bool const with_too_short=true) const;

  #ifdef WITH_LOOP
  /// Access loops_data with range check.
  LOOP const& loop(IndexType const n) const;
  std::vector<LOOP> const& loops() const { return loops_data; }
  #endif // WITH_LOOP


  // *********************************************************************************
  // Protected Member Variables
  // *********************************************************************************

  /// The file path of the pdb_file
  std::string const pdb_file{""};

  #if defined(DRYRUN) || defined(NATIVE_DRYRUN)
  /// The basename of the given pdb_file.
  std::string const pdb_basename{""};
  std::string answer{""};
  #endif // ifdef DRYRUN

  #ifdef WITH_LOOP
protected:
  std::vector<LOOP> loops_data{};

public:
  #endif // RWITH_LOOP


  std::vector<SSE> const data{};
  std::size_t const size{0};

  std::vector<std::vector<IndexType>> const index_vec{};
  std::vector<IndexType> const serial_strand_id{};

protected:

  // *********************************************************************************
  // Protected Member Functions
  // *********************************************************************************


  /// @brief  Read the SSE header (or stride file) and ATOM lines of pdb_file.
  /// @return sse_vec
  std::vector<SSE> read_pdb(stride_stream & stride);


  /// @brief  A helper function of read_pdb() to properly get SSE headers.
  ///         Get from PDB file if a stride file is not available.
  ///         If available, read from the stride file.
  SSE::SSE_Header_vec read_sse_header(stride_stream & stride) const;


  /// @brief  Read the SSE header
  /// @return SSE_Header_vec
  /// @param  ifs input_file_stream of the input pdb_file
  SSE::SSE_Header_vec read_sse_header_pdb(std::istream & is) const;

  /// @brief  A stride version of read_sse_header_pdb()
  /// @return SSE_Header_vec
  /// @param  ifs input_file_stream of the input stride_stream
  SSE::SSE_Header_vec read_sse_header_stride(std::istream & is) const;



  /// @brief  Rewad the CA ATOM records from the pdb_file
  /// @return ca_line_dict a dictionary (key = residue number, value = ATOM line)
  /// @param  ifs in_file_stream of the input pdb_file
  std::unordered_map<int, std::string> read_atom_lines(std::ifstream & ifs);


  /// @brief  Initialize serial_strand_id
  /// @return vector like {255, 0, 1, 2} if There is 4 SSEs, {'H', 'E', 'E', 'E'}.
  std::vector<IndexType> init_serial_strand_id() const;

  /// @brief  Help initialize the \c index_vec
  std::vector<IndexType> gen_index_vec_helper(char const type, bool const with_too_short) const;

  /// @brief  Initialize \c index_vec .
  std::vector<std::vector<IndexType>> init_index_vec() const;

};

} // namespace pdb

#endif //ifndef RPERM_PDB_SSES_H_

