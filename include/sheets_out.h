// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SHEETS_OUT_H_
#define SHEETS_OUT_H_

#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sheet/adj_list_with_sub.h"
#include "sheet/sheets.h"
#include "sheet/directed_adjacency_list.h"

#include "data_store.h"
#include "table.h"

namespace sheets_out {

// **************************************************************************
// Print Sheet Information
// **************************************************************************

/// @brief  Output information of each sheet in a mmcif-like format.
void print_sheet(table::Table<table::Sheet> & tbl,
                 sheet::DirectedAdjacencyList const& adj);

// Helper Functions for print_sheet()

/// @brief  Given a sheet object, return a sorted vector of Sub-Strands that belong
///         to the sheet.
std::vector<sheet::SubStrand> sort_sheet_members(sheet::Sheet const& sheet);


/// @brief  Given a sheet object, return whether the sheet is composed of only
///         sequentially consecutive beta-strands.
bool is_all_consec(sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj);


/// @brief  Check if the given sheet is all-Parallel or all-Anti-Parallel, or neither.
std::tuple<bool, bool> check_all_pap(sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj);


// **************************************************************************
// Print Extracted Sheets and Topology String
// **************************************************************************

void extracted_adjacent_substr_out(table::TBLExtractedSheet & tbl,
                                   unsigned const n,
                                   sheet::DirectedAdjacencyList const& adj);


/// @brief  Check if the given vector of Sub-Strands has any cycles in it.
/// @return True if ss_vec has any cycle, false otherwise.
bool cycle_checker(std::vector<sheet::SubStrand> const& ss_vec,
                   sheet::DirectedAdjacencyList const& adj);

sheet::SubStrandsPairKeyVec adj_sub_vec2pair_key_vec(sheet::AdjSubVec const& adj_sub_vec);


/// @brief  Extract adjacent \c n SubStrands from \c sheet .
std::vector<std::vector<sheet::SubStrand>> extract_adjacent_substr(unsigned const n, sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj);


/// @brief  A helper function for extract_adjacent_substrands()
std::vector<std::vector<sheet::SubStrand>> extract_from_one_substr(sheet::SubStrand const& start_ss, unsigned const n, sheet::Sheet const& sheet, sheet::DirectedAdjacencyList const& adj);

/// @brief  A helper function for extract_from_one_substr()
void recursive_extract(std::vector<sheet::SubStrand> const& current_path,
                       std::vector<std::vector<sheet::SubStrand>> & found_path,
                       sheet::Sheet const& sheet,
                       unsigned const n, sheet::DirectedAdjacencyList const& adj);

/// @brief  Return a sorted copy of given vector
std::vector<sheet::SubStrand> sort_substr_vec(std::vector<sheet::SubStrand> const& ss_vec);

// **************************************************************************
// Topology String Generator
// **************************************************************************

/// @brief  For sheets that are directed and without branch, return "+1x+2x-1...".
class topology_string {
public:
  // *******************************************************************
  // Public Member Types
  // *******************************************************************

  /// One object for a consecutive pair of Sub-Strands.
  struct ss_pair_arrangement {
    ss_pair_arrangement(int const n, bool const d):
      to_next{n != 0 ? n : throw std::logic_error{"Internal Logical Error. In topology_string class."}},
      direction{d} {}

    /// The relative position of the second Sub-Strand relative to the first.
    int const to_next{0};

    /// True for Parallel, false for Anti-Parallel.
    bool const direction{true};
  };

  struct ss_position {
    ss_position(int const p, bool const d): pos{p}, direction{d} {}
    ss_position(unsigned const s, int const p, bool const d):
      seq_id{s}, pos{p}, direction{d} {}

    /// Sequential position of this Sub-Strand.
    std::size_t seq_id{0};

    /// The position of this Sub-Strand. 1 if the most left side.
    int pos{0};

    /// True for relatively Parallel to the sequentially first Sub-Strand in the sheet.
    /// False for relatively Anti-Parallel.
    bool direction{true};
  };


  // *******************************************************************
  // Public Member Functions
  // *******************************************************************

  topology_string(std::vector<sheet::SubStrand> const& ss_vec, bool const with_cycle,
                  sheet::DirectedAdjacencyList const& a):
    adj{a},
    pair_style{init_pair_style(ss_vec,
                               a.adj_sub().substr_vec2adj_sub_vec(ss_vec.cbegin(),
                                                                  ss_vec.cend()),
                               with_cycle)},
    ss_position_style{init_position_style()}
  {}


  topology_string(sheet::Sheet const& s, sheet::DirectedAdjacencyList const& a):
    adj{a},
    pair_style{init_pair_style(s)}, ss_position_style{init_position_style()}
  {
    assert((s.undirected() or s.size() != s.member().size()) ==
           (pair_style.size() == 0 and ss_position_style.size() == 0));

    assert((pair_style.size() != 0) ==
           (pair_style.size() + 1 == ss_position_style.size()));
  }

  /// Generate a string from this object.
  std::string str(unsigned const style=0) const;



private:
  // *******************************************************************
  // Private Member Functions
  // *******************************************************************

  /// Initialize pair_style from a sequentially sorted vector of sheet::SubStrands.
  std::vector<ss_pair_arrangement> init_pair_style(std::vector<sheet::SubStrand> const& ss_vec, sheet::AdjSubVec const& adj_sub_vec, bool const with_cycle) const;


  /// Wrapper function to initialize from a sheet::Sheet object.
  std::vector<ss_pair_arrangement> init_pair_style(sheet::Sheet const& sheet) const;


  /// @brief  If the first sign ('+' or '-') is '-', invert the all signs.
  std::vector<ss_pair_arrangement> modify_to_plus(std::vector<ss_pair_arrangement> const& orig) const;

  std::vector<ss_position> init_position_style() const;

  /// For pair style topology string output. Helper function for str().
  std::string pair_style_str() const;

  /// For position style topology string output. Helper function for str().
  std::string position_style_str(std::vector<ss_position> const& pos_vec,
                                 bool const use_pos=false) const;

  /// New position style topology string
  std::string position_style_2_str(std::vector<ss_position> const& pos_vec) const;

  // *******************************************************************
  // Private Member Variables
  // *******************************************************************

  sheet::DirectedAdjacencyList const& adj;

  /// Pair style data (empty vector if NA)
  std::vector<ss_pair_arrangement> const pair_style{};
  /// Each Sub-Strand style data (empty vector if NA)
  std::vector<ss_position> const ss_position_style{};
};





} // namespace sheets_out

#endif // ifndef SHEETS_OUT_H_
