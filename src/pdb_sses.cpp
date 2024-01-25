// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "pdb/exceptions.h"
#include "pdb/sses.h"
#include "pdb/tools.h"

namespace pdb {


// *****************************************************************************
// Public Member gen_index_vec()
// *****************************************************************************

std::vector<IndexType> const& SSES::gen_index_vec(char const type, bool const with_too_short) const {
  IndexType index;
  switch (type) {
    case 'H': index = 0;
              break;
    case 'E': index = 1;
              break;
    case 'A': index = 2;
              break;
    default: throw unknown_sse_type(type, "In SSES::gen_index_vec()");
  }
  index += with_too_short ? 3 : 0;
  return index_vec[index];
} // public member function gen_index_vec()


#ifdef WITH_LOOP

// *****************************************************************************
// Public Member function loop()
// *****************************************************************************

LOOP const& SSES::loop(IndexType const n) const {
  try {
    return loops_data.at(n);
  } catch (std::out_of_range const& e) {
    throw loop_access_out_of_range{n, loops_data.size()-1};
  }
} // public member function loop()

#endif // WITH_LOOP



// *****************************************************************************
// Public Member function sse_ind_of()
// *****************************************************************************

IndexType SSES::sse_ind_of(int const resnum, char const type, int const offset,
                        bool const with_too_short, char const hbond_atom) const {


  // type check
  if (type != 'H' and type != 'E' and type != 'A') {
    throw unknown_sse_type(type, "In SSES::sse_ind_of()");
  }

  // Prepare the vector of IndexType.
  auto const indices = gen_index_vec(type, with_too_short);

  // Set the actual offset for N-term and C-term.
  int const n_offset = offset == 0 ? 0 : (hbond_atom != 'N' ? offset : offset - 1);
  int const c_offset = offset == 0 ? 0 : (hbond_atom != 'C' ? offset : offset - 1);

  // Comparison function object
  ResnumSSEComp const comp{data, n_offset, c_offset};

  // Run binary search
  auto const r = std::equal_range(indices.cbegin(), indices.cend(), resnum, comp);

  // If the first IndexType that is not less than resnum equals
  // the first IndexType that is greater than resnum, resnum is in the loop region.
  if (r.first == r.second) {
    throw pdb::non_sse_resnum(resnum, "In sse_ind_of().");
  }

  return *r.first;

} // private member function sse_ind_of()




// *****************************************************************************
// Public Member function read_pdb()
// *****************************************************************************

std::vector<SSE> SSES::read_pdb(stride_stream & stride) {

  // read SSE headers
  auto const headers = read_sse_header(stride);


  // if there is no SSEs in this structure, just return.
  // (To avoid invalid memory access. Otherwise, n_loops might be UNSIGNED_MAX.)
  if (headers.size() == 0) { return std::vector<SSE>{}; }

  // a vector to be returned
  std::vector<SSE> sses;

  // prepare the memory for sses
  sses.reserve(headers.size());

  // pdb file stream
  std::ifstream ifs_pdb;
  open_input(ifs_pdb, pdb_file);

  // read atom lines
  auto const atom_lines = read_atom_lines(ifs_pdb);

  #ifdef WITH_LOOP
  // Initialize loops_data
  auto const n_loops = headers.size() - 1;
  for (IndexType i = 0; i < n_loops; ++i) {
    try {
      loops_data.push_back({headers[i].end + 1, headers[i+1].init - 1, i, atom_lines});
    } catch (invalid_sse_range const&) {
      loops_data.push_back({headers[i].end, headers[i+1].init, i});
    }
  }
  #endif // WITH_LOOP


  // generate the vector of SSEs
  auto const max = headers.size();
  for (IndexType i = 0; i < max; ++i) {
    sses.emplace_back(headers[i], i, atom_lines);
  }
  return sses;
} // private member function read_pdb()




// *****************************************************************************
// Protected Member function read_sse_header()
// *****************************************************************************

SSE::SSE_Header_vec SSES::read_sse_header(stride_stream & stride) const {
  SSE::SSE_Header_vec sse_headers;

  if (stride.empty) {
    std::ifstream ifs;
    open_input(ifs, pdb_file);
    sse_headers = read_sse_header_pdb(ifs);
  } else {
    sse_headers = read_sse_header_stride(stride.ss);
  }

  // sort sse_vec by its initial residue number
  std::sort(sse_headers.begin(), sse_headers.end(),
              [](SSE_Header const& a, SSE_Header const& b) { return a.init < b.init; });
  return sse_headers;
} // private member function read_sse_header()



// *****************************************************************************
// Protected Member function read_sse_header_pdb()
// *****************************************************************************

SSE::SSE_Header_vec SSES::read_sse_header_pdb(std::istream & is) const {

  SSE::SSE_Header_vec sse_headers;

  // put back the reading point to the head of the input stream
  is.seekg(0);

  for (std::string buff; std::getline(is, buff);) {
    if (buff.substr(0, 5) == "HELIX") {
      auto const resnum_init = std::stoi(buff.substr(21, 4));
      auto const resnum_end  = std::stoi(buff.substr(33, 4));
      sse_headers.emplace_back('H', resnum_init, resnum_end);

    } else if (buff.substr(0, 5) == "SHEET") {
      auto const resnum_init = std::stoi(buff.substr(22, 4));
      auto const resnum_end  = std::stoi(buff.substr(33, 4));
      sse_headers.emplace_back('E', resnum_init, resnum_end);
    } else if (buff.substr(0, 4) == "ATOM") {
      break;
    }
  }
  return sse_headers;
} // private member function read_sse_header_pdb()



// *****************************************************************************
// Protected Member function read_sse_header_stride()
// *****************************************************************************

SSE::SSE_Header_vec SSES::read_sse_header_stride(std::istream & is) const {
  SSE::SSE_Header_vec sse_headers;

  // put back the reading point to the head of the input stream
  is.seekg(0);

  for (std::string buff; std::getline(is, buff);) {
    // if LOC finished
    if (buff.substr(0, 3) == "ASG") {
      break;
    } else if (buff.substr(0, 3) != "LOC") {
      continue;
    }

    if (buff.substr(5, 10) == "AlphaHelix") {
      auto const resnum_init = std::stoi(buff.substr(22, 27));
      auto const resnum_end = std::stoi(buff.substr(40, 45));
      sse_headers.push_back(SSE_Header{'H', resnum_init, resnum_end});

    } else if (buff.substr(5, 6) == "Strand") {
      auto const resnum_init = std::stoi(buff.substr(22, 27));
      auto const resnum_end = std::stoi(buff.substr(40, 45));
      sse_headers.push_back(SSE_Header{'E', resnum_init, resnum_end});
    }
  }
  return sse_headers;
} // private member function read_sse_header_stride()



// *****************************************************************************
// Protected Member function read_atom_lines()
// *****************************************************************************

std::unordered_map<int, std::string> SSES::read_atom_lines(std::ifstream & ifs) {
  ifs.seekg(0);
  std::unordered_map<int, std::string> atom_lines;

  for (std::string buff; std::getline(ifs, buff);) {
    if (buff.substr(0, 4) == "ATOM" and buff.substr(12, 4) == " CA ") {
      atom_lines.insert({std::stoi(buff.substr(22, 4)), buff});

    #ifdef DRYRUN // This part prevents this function to be const.
    } else if (buff.substr(0, 6) == "ANSWER") {
      answer = split(buff, "\\s+")[1];
    #endif // ifdef DRYRUN

    }
  }
  return atom_lines;
} // protected member function read_atom_lines()




// *****************************************************************************
// Protected Member function init_serial_strand_id()
// *****************************************************************************

std::vector<IndexType> SSES::init_serial_strand_id() const {
  std::vector<IndexType> ret(size, 255u);

  IndexType serial_counter = 0;
  for (auto const i : gen_index_vec('E')) {
    ret[i] = serial_counter;
    ++serial_counter;
  }

  return ret;
} // protected member function init_serial_strand_id()




// *****************************************************************************
// Protected Member function gen_index_vec_helper()
// *****************************************************************************
std::vector<IndexType> SSES::gen_index_vec_helper(char const type, bool const with_too_short) const {
  std::vector<IndexType> ret;

  for (IndexType i = 0; i < size; ++i) {
    if (type == 'A' or data[i].type == type) {
      if (with_too_short == false and data[i].too_short) {
        continue;
      }
      ret.push_back(i);
    }
  }
  return ret;
}




// *****************************************************************************
// Protected Member function init_index_vec()
// *****************************************************************************
std::vector<std::vector<IndexType>> SSES::init_index_vec() const {
  return std::vector<std::vector<IndexType>>{
    gen_index_vec_helper('H', false), // 0
    gen_index_vec_helper('E', false), // 1
    gen_index_vec_helper('A', false), // 2
    gen_index_vec_helper('H', true),  // 3
    gen_index_vec_helper('E', true),  // 4
    gen_index_vec_helper('A', true),  // 5
  };
} // protected member function init_index_vec()

} // namespace pdb

