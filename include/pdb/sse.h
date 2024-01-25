// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PDB_SSE_H_
#define PDB_SSE_H_

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "pdb/constants.h"
#include "pdb/exceptions.h"
#include "pdb/atom.h"


namespace pdb {


using ATOM_vec_iter = std::vector<ATOM>::const_iterator;


// ************************************************************************
// Struct SSERepAtom
// ************************************************************************

/// Representative Information for an SSE class object
struct SSEReprInfo {
public:
  /// @brief Constructor store the representative data based on the SSE type
  /// @exception pdb::unknown_sse_type
  SSEReprInfo(char const type, unsigned const intrvl_) :
    len{type == 'H' ? 4u : 2u},
    coeff{type == 'H' ?
            std::vector<double>{{0.74, 1.0, 1.0, 0.74}} :
            std::vector<double>{{1.0, 1.0}}},
    div{type == 'H' ? 3.48 : 2.0},
    intrvl{intrvl_}
  {
    if (type != 'H' and type != 'E') {
      throw unknown_sse_type(type, "In Construction of SSEReprInfo.\n");
    }
  }

  unsigned const len{0};
  std::vector<double> const coeff{0};
  double const div{0.0};
  unsigned const intrvl{0};
}; // struct SSEReprInfo


// ************************************************************************
// Struct SSE_header
// ************************************************************************
/// SSE header information
struct SSE_Header {
  SSE_Header(char const t, int const i, int const e) : type{t}, init{i}, end{e} {}
  // Member variables cannot be const
  // because the vector of SSE_header need to be sorted after reading the PDB file.
  char type{'H'};
  int init{0};
  int end{1};

  bool operator==(SSE_Header const& other) const {
    return type == other.type and init == other.init and end == other.end;
  }
};




// ************************************************************************
// Class SSE_Base
// ************************************************************************
class SSE_Base {
public:
  SSE_Base(int const i, int const e, IndexType const ind,
           std::unordered_map<int, std::string> const& ca_lines):
    init{i},
    end{e < init ?
      throw invalid_sse_range(init, e, "In constructor of class SSE"): e},
    index{ind}, entire{entirety_check(ca_lines)},
    atoms{read_pdb_atom(ca_lines)},
    n_pdb{count_real_atoms()}
  {}

  /// This version is used when the size of the LOOP is 0.
  SSE_Base(int const l, int const n, IndexType const ind):
    init{l}, end{n}, index{ind}, entire{false} {}


  /// Given a residue number, return a const_iterator to the ATOM with that residue
  /// number in the vector atoms.
  ATOM_vec_iter atom_vec_iter(int const resnum) const;


  int const init{0};
  int const end{1};
  IndexType const index{0};

  /// - true   if this SSE has all the ATOMs in range [init, end]
  /// - false  if there is at least one residue which is missing in the input PDB file.
  /// if an SSE is not entire, loop cross filtering cannot performed
  /// @sa check_entirety()
  bool const entire{true};


  /// Stores the ATOM class objects of the input PDB file.
  std::vector<ATOM> const atoms{};

  /// The total number of real ATOMs (not padding atoms, ATOM::pdb == true)
  unsigned const n_pdb{0};

protected:
  /// Check whether this SSE has all the ATOMs or not
  bool entirety_check(std::unordered_map<int, std::string> const& ca_lines) const;


  /// Read ca_lines and generate a vector of ATOMs
  std::vector<ATOM> read_pdb_atom(std::unordered_map<int, std::string> const& ca_lines) const;

  /// Count the number of pdb ATOMs to initialize n_pdb.
  unsigned count_real_atoms() const;

};


// ************************************************************************
// Class SSE
// ************************************************************************
/// Stores the information about an SSE
class SSE: public SSE_Base{
public:

  // ************************************************************************
  // Public Member Types
  // ************************************************************************

  /// A vector of SSE_Header (std::tuple<char, int, int>)
  using SSE_Header_vec = std::vector<SSE_Header>;

  using vec_str_itr = std::vector<std::string>::const_iterator;


  // ************************************************************************
  // Public Member Functions
  // ************************************************************************

  /// @brief constructor
  SSE(SSE_Header const& header, IndexType const index_,
      std::unordered_map<int, std::string> const& ca_lines, unsigned const intrvl=1) :
    SSE_Base{header.init, header.end, index_, ca_lines},

    type{header.type != 'H' and header.type != 'E' ?
      throw unknown_sse_type(header.type, "In constructor of class SSE") : header.type},
    index_bit{1u << index_}, rep{type, intrvl},
    too_short(too_short_check()),
    rep_atoms({{gen_representative(0),
                gen_representative(atoms.size() - static_cast<int>(rep.len)
                                   - static_cast<int>(rep.intrvl))}}),
    with_rep{rep_check()}

    {}


  /// return the coordinates of the representative ATOM
  /// reverse : true for reversed, false for non-reversed
  /// head : n-term
  /// tail : c-term

  Eigen::Vector3d const& rep_outer_head(unsigned const reverse) const;
  Eigen::Vector3d const& rep_inner_head(unsigned const reverse) const;
  Eigen::Vector3d const& rep_inner_tail(unsigned const reverse) const;
  Eigen::Vector3d const& rep_outer_tail(unsigned const reverse) const;



  /// @brief Generate a stamp specific to the loop which consists of this and other SSE.
  /// @return stamp_string which is specific to this loop
  /// @retval "+2-5" if the index for this SSE (non-reversed) is 2 and 5 for other (reversed)
  /// @retval "-5+2" if the same input as above example and alt = true
  /// @param other an SSE after this SSE in the loop
  /// @param reverse reverse flags for this and other SSE
  /// @param alt if true, the order of SSE indices returned will be swapped. default is false.
  std::string stamp(SSE const& other, unsigned const reverse, bool const alt=false) const;

  /// @brief          Check if the given \c resnum is in the range of this SSE or not.
  /// @retval true    If \c resnum is in this SSE
  /// @retval false   If \c resnum is not in this SSE
  /// @param  resnum  A residue number to decide if in this SSE.
  /// @param  offset  Number of residues to virtually extend from the default range.\n
  ///                 For example, when the range of an SSE is [10, 15]; \n
  ///                 \li By default, \c resnum such as {10, 11, ... 15} are in range.\n
  ///                 \li But with offset=1, \c resnum of {9, 10, ... 16} will be in range.
  /// @param  HAtom   Specify which side of given residue has a hbond. By default, both.
  ///                 Charactors other than \c 'N' or \c 'C' will be treated as \c 'A' .\n
  ///                 If the \c resnum is at the edge of offset range, decide whether to
  ///                 include that \c resnum in this SSE. \n
  ///                 For example, when the range of an SSE is [10, 15] and offset=1; \n
  ///                 \li \c HAtom='N' and \c resnum=16 : returns \c true .
  ///                 \li \c HAtom='C' and \c resnum=16 : returns \c false .
  ///                 \li \c HAtom='C' and \c resnum=9  : returns \c true .
  ///                 \li \c HAtom='N' and \c resnum=9  : returns \c false .
  bool in_range(int const resnum, int const offset=0, char const HAtom='A') const noexcept;



  /// @brief Calculate the distance between the closest CA atoms in this and other.
  /// @return the smallest distance
  /// @param other a SSE object to calculate the distance to
  double distance(SSE const& other) const;



  // ************************************************************************
  // Public Member Variables
  // ************************************************************************

  char const type{'H'};
  unsigned const index_bit{1};

  SSEReprInfo const rep{'H', 1};


  /// @brief  store whether this SSE is too short of not
  /// - true  This SSE is too short (Indicates that the length is less than a specific value 
  ///         defined in constants.h. ).
  /// - false This SSE has enough length to perform something.
  ///
  /// Most features which handle SSEs SHOULD NOT include the SSE if its too_short is true.
  /// Sheet topology analysis will use the SSE object even if this value is true.
  bool const too_short{false};


  /// @brief rep_atoms : representative atoms for calculating dihedral angles or disntance ...
  /// {
  ///   {ATOM_1, ATOM_2},  ( n-term side representative atoms {n-term, c-term})
  ///   {ATOM_1, ATOM_2}   ( c-term side representative atoms {n-term, c-term})
  /// }
  /// if not entirety, all ATOMs are zero
  /// @sa gen_representative()
  std::array<std::array<ATOM, 2>, 2> const rep_atoms;


  /// - true  This SSE has representative atoms.
  /// - false Representative atoms are zeros.
  bool const with_rep{false};

   
protected:


  // ************************************************************************
  // Protected Member Functions
  // ************************************************************************

  // **************************************
  // Valid check functions
  // **************************************


  /// Check if this SSE has an enough length.
  /// This method will be used to initialize too_short.
  bool too_short_check() const noexcept;

  // **************************************
  // Helper functions for initialization
  // **************************************

  /// @brief Generate an array of the representative ATOM objects using gen_representative_atom().
  /// If entire == false, return ATOMs whose all coordinates are zero
  std::array<ATOM, 2> gen_representative(int const first) const;


  /// Generate one representative ATOM object begins with resnum first.
  ATOM gen_representative_atom(int const first) const;


  /// @brief  Check if the representative atoms were generated properly.
  /// @retval true  If the representative atoms are generated.
  /// @retval false If false, representative atoms MUST NOT be accessed.
  bool rep_check() const;

};


} // namespace pdb

#endif //ifndef PDB_SSE_H_

