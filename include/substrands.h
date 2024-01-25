// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SUBSTRANDS_H_
#define SUBSTRANDS_H_

#include <ostream>
#include <unordered_map>

#include "functions.h"
#include "table.h"
#include "pdb/exceptions.h"
#include "sheet/directed_adjacency_list.h"
#include "bab/filter.h"

namespace substrands {
/// A conversion map
/// from a std::string that corresponds to a Sub-Strand
/// to an index of the sheet that the Sub-Strand belongs to
using SubStrandStr2SheetIdxMap = std::unordered_map<std::string, std::size_t>;


/// @brief  Output the Sub-Strands and the resnum range of them.
void substrands_out(table::Table<table::SubStrand> & tbl,
                    sheet::DirectedAdjacencyList const& adj,
                    SubStrandStr2SheetIdxMap const& sheet_id_map);


/// @brief  Generate a map from SubStrandIDs to SheetIDs.
SubStrandStr2SheetIdxMap gen_sheet_id_map(sheet::DirectedAdjacencyList const& adj);

/// @brief  Output helices
void helices_out(table::Table<table::Helix> & tbl,
                 sheet::DirectedAdjacencyList const& adj);

/// @brief  Output Sub-Strand Pairs
void substrands_pair_out(table::TBLSubStrandsPair & tbl,
                         sheet::DirectedAdjacencyList const& adj,
                         SubStrandStr2SheetIdxMap const& sheet_id_map,
                         bab::BabFilter & bab);

/// @brief  Return true if both ss0 and ss1 belong to the same cycle
bool in_cycle(sheet::SubStrand const& ss0, sheet::SubStrand const& ss1,
              sheet::Sheet const& sheet);


/// @brief  Check the connection type from ss0 to ss1. Internnaly call gen_sses_lbts.
std::string check_connection_type(sheet::SubStrand const& ss0,
                                  sheet::SubStrand const& ss1,
                                  sheet::DirectedAdjacencyList const& adj,
                                  SubStrandStr2SheetIdxMap const& sheet_id_map,
                                  out::substr2str const& ss_writer);


/// @brief  Check if the substrands between itr0 and itr1 go to other sheets.
/// @retval true if cmp returns true, false otherwise.
bool check_middle_ss_sheet(std::vector<sheet::SubStrand>::const_iterator itr0,
                           std::vector<sheet::SubStrand>::const_iterator itr1,
                           std::function<bool(std::size_t, std::size_t)> cmp,
                           SubStrandStr2SheetIdxMap const& sheet_id_map,
                           out::substr2str const& ss_writer);



/// @brief  Given bab::BabFilterResult::connection_type, return connection type string.
std::string gen_sses_lbts(std::size_t const c_type, bool const on_same_sheet);


class one_directional_cycle_exception : public pdb::fatal_error_base {
public:
  one_directional_cycle_exception(std::string const& ss0, std::string const& ss1):
    pdb::fatal_error_base{ss0 + " and " + ss1 + "are in the same cycle, "
                          "but there is not path from " + ss0 + " to " + ss1 + " !!\n"
                          "There may be some bugs in substrands::substrands_pair_out()."} {}
};

} // namespace substrands

#endif // ifndef SUBSTRANDS_H_
