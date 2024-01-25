// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <algorithm>
#include <array>
#include <deque>
#include <limits>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cstdlib>

#include <boost/format.hpp>

#include "sheet/adj_list_with_sub.h"
#include "sheet/common.h"
#include "sheet/exceptions.h"
#include "sheet/cb_side.h"
#include "pdb/tools.h"

namespace sheet {


// ******************************************************************************
// Member Function Definition for Class ZoneInfo
// ******************************************************************************

// ********************************************************************
// Public Member Function add_pair()
// ********************************************************************
bool ZoneInfo::add_pair(ZoneResidue const& residue, bool const hbonded,
                        BridgeType const b_type) {
  IndexType const index = hbonded ? 1 : 0;

  // if already registered
  if (adj_residues[index] == residue) {
    return true;
  }

  // if already has the pair of the same hbond status.
  if (adj_set[index]) {
    return false;
  }

  // mark this residue as colored
  colored = true;
  adj_residues[index] = residue;
  adj_set[index] = true;
  bridge_type[index] = b_type;

  return true;
} // public member function add_pair()



// ******************************************************************************
// Member Function Definition for Class StrictZone
// ******************************************************************************

// ********************************************************************
// Public Member Function On()
// ********************************************************************
void StrictZone::On(IndexType const sse_id, int const resnum,
                    IndexType const paired_sse_id, int const paired_resnum,
                    bool const hbonded, ZoneInfo::BridgeType const bridge_type) {
  #if defined(LOGGING) && defined(LOG_SHEET)
  pdb::log((boost::format("ADDING RESIDUE PAIR: SSE[%d](%d) %3d <==> SSE[%d](%d) %3d")
                          % +sse_id % +sses.serial_strand_id[sse_id] % resnum
                          % +paired_sse_id % +sses.serial_strand_id[paired_sse_id] % paired_resnum).str());
  #endif

  try {
    On_one(ZoneResidue{sse_id, resnum, sses}, ZoneResidue{paired_sse_id, paired_resnum, sses}, hbonded, bridge_type);
  } catch (pdb::non_sse_resnum const& e) {

    #if defined(LOGGING) && defined(LOG_SHEET)
    pdb::log(e.what() + std::string("\nADDING SKIPPED"));
    #endif

    return ;
  } catch (zone_info_failure const& e) {

    #if defined(LOGGING) && defined(LOG_SHEET)
    pdb::log(e.what() + std::string("\nADDING SKIPPED"));
    #endif

    return ;
  }

} // public member function On()




// ********************************************************************
// Public Member Function decide_side()
// ********************************************************************
AdjListWithSub StrictZone::decide_side(AdjList const& undirected_adj_list) {
  AdjListWithSub adj_sub;

  // get the set of colored ZoneResidues
  ZoneResidueSet remainder;
  collect_colored(remainder);

  // Each iteration corresponds to each consistent sheet region.
  while (remainder.size() != 0) {
    bfs(adj_sub, remainder, undirected_adj_list);

    #if defined(LOGGING) && defined(LOG_SHEET)
    sub_strands.print(strand_indices);
    #endif // ifdef LOGGING

  }

  // * Remove default_range in the last element of each strand.
  // * Sort SubStrands inside their strand based on the residue numbers.
  // * Update the SubStrand ids in the adjacency map keys.
  adj_sub.update_key_substr(sub_strands.finish());

  // Remove small strand pairs
  adj_sub.cleanup(strand_indices);

  #if defined(LOGGING) && defined(LOG_SHEET)
  pdb::log("After Sorting Sub-Strands");
  sub_strands.print(strand_indices);
  #endif // ifdef LOGGING


  calc_deltas(adj_sub);


  #if defined(LOGGING) && defined(LOG_SHEET)
  // Log: Output the side (Upper, Lower, or Undefined) of Cb atoms for all residues.
  PrintSide();
  sub_strands.print(strand_indices);
  PrintList(adj_sub);
  #endif // ifdef LOGGING

  // Generate adjacency list vec with substrands
  adj_sub.gen_adj_sub_vec();

  return adj_sub;
} // public member function decide_side()




// ********************************************************************
// Public Member Function strict_info()
// ********************************************************************
ZoneInfo const& StrictZone::strict_c_info(ZoneResidue const& zone_res) const {
  return strict[zone_res.serial_str_id][zone_res.serial_res_id];
}




// ********************************************************************
// Public Member Function strict_info()
// ********************************************************************
ZoneInfo & StrictZone::strict_info(ZoneResidue const& zone_res) {
  return strict[zone_res.serial_str_id][zone_res.serial_res_id];
}




#if defined(LOGGING) && defined(LOG_SHEET)
// ********************************************************************
// Public Member Function PrintZone()
// ********************************************************************
void StrictZone::PrintZone() const {
  unsigned counter = 0;
  for (auto const& zone : strict) {
    auto const sse_id = strand_indices[counter];
    pdb::log("Zone for SSE[" + std::to_string(sse_id) + "]");

    unsigned resnum_counter = 0;
    for (auto const zone_info : zone) {
      pdb::log((boost::format("  RESNUM = %3d : %s") 
                                % (sses[sse_id].init + resnum_counter)
                                % (zone_info.colored ? "True" : "False")).str());
      ++resnum_counter;
    }
    ++counter;
  }
} // public member function PrintZone()




// ********************************************************************
// Public Member Function PrintSide()
// ********************************************************************

void StrictZone::PrintSide() const {
  unsigned counter = 0;
  for (auto const& one_sse : strict) {
    auto const sse_id = strand_indices[counter];

    pdb::log("Side for SSE[" + std::to_string(sse_id) + "]");

    unsigned resnum_counter = 0;
    for (auto const& resnum : one_sse) {
      pdb::log((boost::format(" SIDE RESNUM = %3d %2d") 
                                % (sses[sse_id].init + static_cast<int>(resnum_counter))
                                % (+resnum.side)).str());
      ++resnum_counter;
    }
    std::clog << "\n";
    ++counter;
  }
} // public member function PrintSide()



// ********************************************************************
// Public Member Function PrintList()
// ********************************************************************
void StrictZone::PrintList(AdjListWithSub const& adj_list) const {
  pdb::log("AdjListWithSub");

  std::set<SubStrandsPairKey> keys;
  for (auto const& pair : adj_list.map()) {
    keys.insert(pair.first);
  }
  for (auto const& pair : keys) {
    pdb::log((boost::format(" ADJACENT PAIR WITH SUB: %2d_%d  --->  %2d_%d %3d %3d %3d\t%s")
                                  % +pair.str0 % +pair.substr0
                                  % +pair.str1 % +pair.substr1
                                  % adj_list.map().at(pair).delta_1
                                  % adj_list.map().at(pair).delta_2
                                  % adj_list.map().at(pair).residue_pairs
                                  % (adj_list.map().at(pair).direction ?
                                                  "Parallel" : "Anti-Parallel")).str());
  }

} // public member function PrintList()



#endif // ifdef LOGGING



// ********************************************************************
// Protected Member Function On_one()
// ********************************************************************

void StrictZone::On_one(ZoneResidue const& res0, ZoneResidue const& res1,
                        bool const hbonded, ZoneInfo::BridgeType const bridge_type) {

  if (strict_info(res0).add_pair(res1, hbonded, bridge_type)
      == false) {
    auto const& zone_info = strict_c_info(res0);
    throw third_pair_found{res0.sse_id, res0.resnum,
                           res1.sse_id, res1.resnum,
                           zone_info.adj_residues[0].sse_id, zone_info.adj_residues[0].resnum,
                           zone_info.adj_residues[1].sse_id, zone_info.adj_residues[1].resnum};
  }
} // protected member function On_one




// ********************************************************************
// Protected Member Function init_strict()
// ********************************************************************

std::vector<std::vector<ZoneInfo>> StrictZone::init_strict() const {

  std::vector<std::vector<ZoneInfo>> ret(strand_indices.size());
  unsigned serial_counter = 0;
  for (auto const& sse_id : strand_indices) {
    ret[serial_counter].resize(static_cast<unsigned>(sses[sse_id].end - sses[sse_id].init + 1));
    ++serial_counter;
  }
  return ret;
} // protected member function init_strict()




// ********************************************************************
// Protected Member Function collect_colored()
// ********************************************************************

void StrictZone::collect_colored(ZoneResidueSet & set) const {

  IndexType serial_str_id = 0;
  for (auto const& one_sse : strict) {
    unsigned serial_res_id = 0;
    for (auto const& zone_info : one_sse) {
      if (zone_info.colored) {
        auto const sse_id = strand_indices[serial_str_id];
        set.insert(ZoneResidue{sse_id, static_cast<int>(serial_res_id) + sses[sse_id].init, sses});
      }
      ++serial_res_id;
    }

    ++serial_str_id;
  }
} // protected member function collect_colored()




// ********************************************************************
// Protected Member Function calc_deltas()
// ********************************************************************
void StrictZone::calc_deltas(AdjListWithSub & adj_list) const {
  // for all pairs of adjacent Sub-Strands
  for (auto & pair : adj_list.map()) {
    auto const& pair_key = pair.first;
    auto const& pair_node = pair.second;
    bool const dir = pair_node.direction;

    int const delta_1 = count_delta_1(pair_key.sub0(), pair_key.sub1(), dir);
    int const delta_2 = count_delta_2(pair_key.sub0(), pair_key.sub1(), dir);

    adj_list.update_delta(pair_key, delta_1, delta_2);
  }
} // protected member function calc_deltas()




// ********************************************************************
// Protected Member Function count_delta_base()
// ********************************************************************
std::tuple<int, ZoneResidue> StrictZone::count_delta_base(ZoneResidue const& zres_base_start, bool const toward_c_term, SubStrand const& ss_adj) const {

  int delta_base = 0;
  bool go_next = true;
  // for zres from [res] to [the most C-term residue of adj_ss]
  // zres : on the base Sub-Strand
  for (ZoneResidue zres = zres_base_start; go_next;) {

    // for each residue adjacent (counterpart of the bridge) to zres
    // adj_candidate_res : on Sub-Strand adjacent to the base Sub-Strand
    for (auto const& adj_candidate_res : strict_c_info(zres).adj_residues) {
      auto const resnum = adj_candidate_res.resnum;

      // if all the following conditions are true. 2. may not be necessary.
      // 1. adj_candidate_res belongs to Sub-Strand 'ss'
      // 2. the residue number of adj_candidate_res is in range of Sub-Strand 'ss'
      if (adj_candidate_res.serial_str_id == ss_adj.str and
          sub_strands.n_term_res(ss_adj) <= resnum and
          resnum <= sub_strands.c_term_res(ss_adj)) {
        return std::make_tuple(delta_base, adj_candidate_res);
      }
    }

    // Do not count the residues outside of strict hydrogen bonding zone.
    if (strict_c_info(zres).colored) {
      ++delta_base;
    }
    go_next = toward_c_term ? zres.increment(sses) : zres.decrement(sses);
  }

  throw PairedResidueNotFound{zres_base_start.sse_id, zres_base_start.resnum,
                              sses.serial_strand_id[ss_adj.str], ss_adj.substr};
} // protected member function count_delta_base()



// ********************************************************************
// Protected Member Function count_delta_adj()
// ********************************************************************
int StrictZone::count_delta_adj(ZoneResidue const& start, int const step,
                                ZoneResidue const& last) const {

  assert(start.serial_str_id == last.serial_str_id);
  assert(std::abs(step) == 1);

  auto const serial_strand_id_adj = start.serial_str_id;

  int delta = 0;

  // for ser_res_id in range [start, last]
  for (int ser_res_id = start.serial_res_id;
       static_cast<unsigned>(ser_res_id) != last.serial_res_id;
       ser_res_id += step) {

    assert(0 <= ser_res_id);

    // if the target residue is inside the strict hydrogen-bonding pattern,
    // count up the delta.
    if (strict.at(serial_strand_id_adj).at(ser_res_id).colored) {
      ++delta;
    }
  }

  return delta;
} // protected member function count_delta_adj




// ********************************************************************
// Protected Member Function count_delta_1()
// ********************************************************************

int StrictZone::count_delta_1(SubStrand const& ss_base, SubStrand const& ss_adj,
                              bool const dir) const {

  // The most N-term residue (inside the strict hydrogen bonding pattern)
  // of the base Sub-Strand
  ZoneResidue const zres_term_base{strand_indices[ss_base.str],
                                   sub_strands.n_term_res(ss_base), sses};

  // Find the most N-term side (when seen in direction of base Sub-Strand) residue
  // (on the ss_adj) that consists a bridge with a residue on the ss_base.
  int delta_base;
  ZoneResidue zres_bridge_adj;
  std::tie(delta_base, zres_bridge_adj) = count_delta_base(zres_term_base, true, ss_adj);


  // The most N-term residue if Parallel, C-term if Anti-Parallel
  auto const term_res_adj = dir ? sub_strands.n_term_res(ss_adj):
                                  sub_strands.c_term_res(ss_adj);
  ZoneResidue const zres_term_adj{strand_indices[ss_adj.str], term_res_adj, sses};

  // Count the number of residues on ss_adj from the common residue to the edge.
  int const delta_adj = count_delta_adj(zres_bridge_adj, dir ? -1 : 1, zres_term_adj);
  return delta_base - delta_adj;

} // protected member function count_delta_1()




// ********************************************************************
// Protected Member Function count_delta_2()
// ********************************************************************


int StrictZone::count_delta_2(SubStrand const& ss_base, SubStrand const& ss_adj,
                              bool const dir) const {

  // The most C-term residue (inside the strict hydrogen bonding pattern)
  // of the base Sub-Strand
  ZoneResidue const zres_term_base{strand_indices[ss_base.str],
                                   sub_strands.c_term_res(ss_base), sses};

  // Find the most C-term side (when seen in direction of base Sub-Strand) residue
  // (on the ss_adj) that consists a bridge with a residue on the ss_base.
  int delta_base;
  ZoneResidue zres_bridge_adj;
  std::tie(delta_base, zres_bridge_adj) = count_delta_base(zres_term_base, false, ss_adj);


  // residue number of the most N-term residue if Parallel, C-term if Anti-Parallel
  auto const term_res_adj = dir ? sub_strands.c_term_res(ss_adj):
                                  sub_strands.n_term_res(ss_adj);
  ZoneResidue const zres_term_adj{strand_indices[ss_adj.str], term_res_adj, sses};

  // Count the number of residues on ss_adj from the common residue to the edge.
  int const delta_adj = count_delta_adj(zres_bridge_adj, dir ? 1 : -1, zres_term_adj);
  return delta_adj - delta_base;

} // protected member function delta_2_parallel()



// ********************************************************************
// Protected Member Function bfs()
// ********************************************************************

void StrictZone::bfs(AdjListWithSub & adj_list,
                     ZoneResidueSet & remainder,
                     AdjList const& undirected_adj_list) {

  // Relative Direction Vector
  RelativeDirsVec relative_directions(strand_indices.size(), NotSet);

  // Set of serial SSE IDs which are involved with this consistent sheet region.
  std::unordered_set<IndexType> involved_serial_str_ids;

  // Initialize the queue (Push the first residue)
  std::deque<ZoneResidue> queue{*(remainder.begin())};
  ZoneResidueSet queue_contents;
  queue_contents.insert(queue.front());

  // Add new sheet in AdjListWithSub
  adj_list.add_sheet();

  // whether this is the fallback mode (generate an undirected graph instead of directed)
  bool fallback = false;

  // define the first residue as 'Upper'
  strict_info(queue.front()).side = ZoneInfo::SideStatus::Upper;

  // Initialize the direction of the base strand
  relative_directions[queue.front().serial_str_id] = Parallel;

  // Start Breadth First Search
  while (queue.size() != 0) {
    // target : the first element of the queue
    auto const& target = queue.front();

    involved_serial_str_ids.insert(target.serial_str_id);

    // extend the sub_strands of this strand
    sub_strands.extend_substrand(target);

    #if defined(DEBUG) && defined(LOG_SHEET)
    pdb::debug_log("TARGET = SSE[" + std::to_string(target.sse_id) + "], "
                   "RESNUM = " + std::to_string(target.resnum));
    #endif

    // ************************************************************************
    // Run search on paired residues (on other strands)
    // ************************************************************************
    // i == 0 : non-hbonding
    // i == 1 : hbonding
    for (IndexType i = 0; i < 2; ++i) {
      auto const& adj = strict_c_info(target).adj_residues[i];


#if defined(LOGGING) && defined(LOG_SHEET)
      if (adj.has_value) {
        pdb::log("  ADJACENT = SSE[" + std::to_string(adj.sse_id) + "], RESNUM = " +
                 std::to_string(adj.resnum) +
                 (remainder.count(adj) == 0 ? "  (ALREADY EXAMINED)" : ""));
      }
#endif

      // skip if default initialized empty residue, already examined
      if (adj.has_value == false or remainder.count(adj) == 0) {
        continue;
      }

      // set the side of adjacent residue
      strict_info(adj).side = strict_c_info(target).side;

      // extend the sub_strand region of adj
      sub_strands.extend_substrand(adj);

      // ********************************
      // Decide the relative directions
      // ********************************

      // If adj is parallel to the target strand,
      // relative direction of adj strand is the same as target strand.
      bool const direction_to_adj = undirected_adj_list.at({target.sse_id, adj.sse_id}).direction;


      // if the directed graph cannot be generated, fallback to undirected graph.
      if (not fallback) {
        fallback = update_rel_dir(target, adj, direction_to_adj, relative_directions);
      }

      if (fallback) {
        auto const id0 = target.serial_str_id;
        auto const id1 = adj.serial_str_id;
        auto const key0 = SubStrandsPairKey(id0, sub_strands.last_substr_id(id0),
                                            id1, sub_strands.last_substr_id(id1));
        auto const key1 = key0.reverse();

        add_adj_list_count(key0, adj_list, direction_to_adj);
        add_adj_list_count(key1, adj_list, direction_to_adj);

      } else {
        // i == 0 : non-hbonded and i == 1 : hbonded
        auto const key = gen_list_key(target, adj, i, relative_directions);
        // Add Count for adj_list
        add_adj_list_count(key, adj_list, direction_to_adj);
      }

      // Add to AdjacenctSubStrands
      add_adj_substrands(target, adj, i, direction_to_adj);


      // push into the search queue if not in the queue
      push_into_queue(queue, queue_contents, adj);

    } // for adjacent residues



    // *********************************************************************
    // for the residues before and after this reside (on the same strand)
    // *********************************************************************
    std::vector<ZoneResidue> before_after;
    before_after.reserve(2);

    // Add a residue sequentially before and after the target if they're not out of range.
    before_after.push_back(ba_check(queue.front(), -1, remainder));
    before_after.push_back(ba_check(queue.front(), 1, remainder));

    // for Befor and After the target residue
    for (auto const& ba : before_after) {
      // Skip if out of range or not colored
      if (ba.has_value == false or strict_c_info(ba).colored == false) {
        continue;
      }

      if (strict_c_info(target).side == ZoneInfo::SideStatus::Upper) {
        strict_info(ba).side = ZoneInfo::SideStatus::Lower;

      } else if (strict_c_info(target).side == ZoneInfo::SideStatus::Lower) {
        strict_info(ba).side = ZoneInfo::SideStatus::Upper;

      } else {
        throw std::runtime_error("In StrictZone::bfs()");
      }

      // push into the search queue if not in the queue
      push_into_queue(queue, queue_contents, ba);

    } // for Befor and After the target residue

    // delete the used target before poping from the queue.
    // (target is a reference to the first element of queue.)
    remainder.erase(target);
    queue_contents.erase(target);
    queue.pop_front();
  } // while queue.size() != 0


  // **************************
  // Clean Up: Sub-Strands
  // **************************

  // Store the sub-strands whose length is too short.
  SubStrandSet const too_short_subs = sub_strands.cleanup_sheet();

  // for all adjacent-pairs, delete the strands pair that contains at least 1
  // too short sub-strands.
  adj_list.erase_too_short(too_short_subs);
  if (fallback) {
    adj_list.ensure_undirected();
  }


  // Push the fallback status into adj_list
  adj_list.register_sheet_directed(not fallback);

} // protected member function bfs()



// ********************************************************************
// Protected Member Function add_adj_list_count()
// ********************************************************************

void StrictZone::add_adj_list_count(SubStrandsPairKey const& key,
                                    AdjListWithSub & adj_list,
                                    bool const direction_to_adj) {

#if defined(LOGGING) && defined(LOG_SHEET)
  auto const key_str = (boost::format("KEY = (%d-%d, %d-%d)")
                                      % +key.str0 % +key.substr0
                                      % +key.str1 % +key.substr1).str();
  std::string head, tail;
  if (adj_list.map().count(key) == 0) {
    head = "      INSERT: ";
    tail = direction_to_adj ? ", Parallel" : ", Anti-Parallel";
  } else if (adj_list.map(key).direction == direction_to_adj) {
    head = "      ADD PAIR: ";
    tail = ", PAIRS = " + std::to_string(adj_list.map(key).residue_pairs);
  } else {
    head = "\n\nFAILED ";
    tail = ", PAIRS = " + std::to_string(adj_list.map(key).residue_pairs);
  }
  pdb::log(head + key_str + tail);
#endif


  // if this is the first residue pair for this key
  if (adj_list.map().count(key) == 0) {
    adj_list.insert_map(key, direction_to_adj);
    adj_list.register_key_vec(key);

  // if the key already exists, AND if the direction is the same
  } else if (adj_list.map(key).direction == direction_to_adj) {
    adj_list.add_count_pairs_map(key);

  // if the direction was changed
  } else {
    throw std::runtime_error(direction_to_adj ?
                               "DIRECTION CHANGED: Anti-Parallel -> Parallel" :
                               "DIRECTION CHANGED: Parallel -> Anti-Parallel");
  }

} // protected member function add_adj_list_count()



// ********************************************************************
// Protected Member Function resolve_rel_dir()
// ********************************************************************

bool StrictZone::update_rel_dir(ZoneResidue const& target, ZoneResidue const& adj,
                                bool const direction_to_adj,
                                RelativeDirsVec & relative_directions) const {

  // ******************************************************
  // #rel_dir_target  #direction_to_adj  #rel_dir_adj
  // Parallel         Parallel           Parallel
  // Parallel         Anti-Parallel      Anti-Parallel
  // Anti-Parallel    Parallel           Anti-Parallel
  // Anti-Parallel    Anti-Parallel      Parallel
  // ******************************************************

  RelDir adj_rel_direction;
  // Parallel
  if (relative_directions[target.serial_str_id] == Parallel) {
    adj_rel_direction = direction_to_adj ? Parallel : AntiParallel;

  // Anti-Parallel
  } else if (relative_directions[target.serial_str_id] == AntiParallel) {
    adj_rel_direction = direction_to_adj ? AntiParallel : Parallel;

  // Not set
  } else {
    throw TargetRelativeDirectionNotSet(strand_indices[target.serial_str_id]);
  }


  // if rel_dir_adj is not set yet
  if (relative_directions[adj.serial_str_id] == NotSet) {
    relative_directions[adj.serial_str_id] = adj_rel_direction;

  // if rel_dir_adj is already set and differs from the adj_rel_direction value
  } else if (relative_directions[adj.serial_str_id] != adj_rel_direction) {
    pdb::warning("RELATIVE DIRECTION CHANGED: SSE[" +
                 std::to_string(strand_indices[adj.serial_str_id]) + "]; "
                 "FALLING BACK TO AN UNDIRECTED GRAPH");
    return true;
  }

#if defined(LOGGING) && defined(LOG_SHEET)
  pdb::log("      RELATIVE DIRECTION OF ADJACENT SSE[" + std::to_string(adj.sse_id) +
           "] = " + (adj_rel_direction ? "Parallel" : "Anti-Parallel"));
#endif
  return false;
} // protected member function update_rel_dir()



// ********************************************************************
// Protected Member Function gen_list_key()
// ********************************************************************
SubStrandsPairKey StrictZone::gen_list_key(ZoneResidue const& res0,
                                           ZoneResidue const& res1, bool const hbonded,
                                           RelativeDirsVec const& rel_dir) const {

  // if the relative direction is not set yet
  if (rel_dir[res0.serial_str_id] == NotSet) {
    throw TargetRelativeDirectionNotSet(strand_indices[res0.serial_str_id]);
  }

  auto const side = strict_c_info(res0).side;

  // **********************************************
  // # Whether res1 is on the right side of res0 can be determined
  // # by evaluating the following expression.
  //   (Rel_Direction == (Side == Hbonded))
  //
  // # Truth Table
  // #Rel_Direction   #Side    #Hbonded    #Left -> Right
  //
  // Parallel         Upper    True                res0 -> res1
  // Parallel         Upper    False       res1 -> res0
  // Parallel         Lower    True        res1 -> res0
  // Parallel         Lower    False               res0 -> res1
  //
  // Anti-Parallel    Upper    True        res1 -> res0
  // Anti-Parallel    Upper    False               res0 -> res1
  // Anti-Parallel    Lower    True                res0 -> res1
  // Anti-Parallel    Lower    False       res1 -> res0
  //
  // # These variables are treated as a binary value True and False respectively.
  //   * Relative Direction : Parallel        Anti-Parallel
  //   * Side               : Upper           Lower
  //   * Hbonded            : True            False
  //   * Left -> Right      : res0 -> res1    res1 -> res0
  //
  // **********************************************

  // if res1 is on the right of res0
  if (rel_dir[res0.serial_str_id] == (side == hbonded)) {
    return SubStrandsPairKey(res0.serial_str_id,
                             sub_strands.last_substr_id(res0.serial_str_id),
                             res1.serial_str_id,
                             sub_strands.last_substr_id(res1.serial_str_id));
  } else {
    return SubStrandsPairKey(res1.serial_str_id,
                             sub_strands.last_substr_id(res1.serial_str_id),
                             res0.serial_str_id,
                             sub_strands.last_substr_id(res0.serial_str_id));
  }
} // protected member function gen_list_key()



// ********************************************************************
// Protected Member Function add_adj_substrands()
// ********************************************************************

void StrictZone::add_adj_substrands(ZoneResidue const& res0, ZoneResidue const& res1,
                                    bool const hbonded, bool const rel_dir) {

  auto const str0 = res0.serial_str_id;
  auto const str1 = res1.serial_str_id;

  // true   if (Upper and HBonded) or (Lower and Non-HBonded)
  // false  if (Upper and Non-HBonded) or (Lower and HBonded)
  bool const right = strict_c_info(res0).side == hbonded;
  adj_substrands.add(SubStrand{str0, sub_strands.last_substr_id(str0)},
                     SubStrand{str1, sub_strands.last_substr_id(str1)},
                     right, rel_dir);

  #if defined(DEBUG) && defined(LOG_SHEET)
  auto const str_ind = sses.gen_index_vec('E');
  pdb::debug_log("SubStrand Add: " +
                 SubStrand{str1, sub_strands.last_substr_id(str1)}.string(str_ind) +
                 " is on " + (right ? "right" : "left") + " of " +
                 SubStrand{str0, sub_strands.last_substr_id(str0)}.string(str_ind) +
                 " : " + (rel_dir ? "Parallel" : "Anti-Parallel")
                 );
  #endif

} // protected member function add_adj_substrands()




// ********************************************************************
// Protected Member Function push_into_queue()
// ********************************************************************

void StrictZone::push_into_queue(std::deque<ZoneResidue>& queue,
                                 ZoneResidueSet & contents,
                                 ZoneResidue const& new_res) const {

  if (contents.count(new_res) == 0) {
    queue.push_back(new_res);
    contents.insert(new_res);

#if defined(LOGGING) && defined(LOG_SHEET)
    pdb::log("  ADD QUEUE: SSE[" + std::to_string(new_res.sse_id) + "], RESNUM = " +
             std::to_string(new_res.resnum));
#endif // ifdef LOGGING

  }
} // protected member function push_into_queue()




// ********************************************************************
// Protected Member Function ba_check()
// ********************************************************************

ZoneResidue StrictZone::ba_check(ZoneResidue const& target, int const diff,
                                 ZoneResidueSet const& remainder) const {

  // Empty ZoneResidue to be returned on failure.
  ZoneResidue const empty{};

  // if the adjacent residue is out of range
  int const new_serial_res_id = static_cast<int>(target.serial_res_id) + diff;
  if (new_serial_res_id < 0 or
      static_cast<int>(strict[target.serial_str_id].size()) <= new_serial_res_id) {

#if defined(LOGGING) && defined(LOG_SHEET)
      pdb::log((boost::format("OUT OF RANGE: SSE[%d] RESNUM = %d")
                              % +target.sse_id
                              % (target.resnum + diff)).str());
#endif // ifdef LOGGING

    return empty;
  }

  ZoneResidue const ba{target.sse_id, target.resnum + diff, sses};

  // if already examined
  if (remainder.count(ba) == 0) {

#if defined(LOGGING) && defined(LOG_SHEET)
    pdb::log((boost::format("NOT IN REMAINDER: SSE[%d] RESNUM = %d")
                            % +target.sse_id
                            % (target.resnum + diff)).str());
#endif // ifdef LOGGING

    return empty;
  }

  auto const& info_t = strict_c_info(target);
  auto const& info_ba = strict_c_info(ba);

  // if the bridges of same type continue
  for (unsigned char i_t_adj = 0; i_t_adj < 2; ++i_t_adj) {
    if (not info_t.adj_set[i_t_adj]) {
      continue;
    }
    auto const adj_sse_id = info_t.adj_residues[i_t_adj].sse_id;
    for (unsigned char i_ba_adj = 0; i_ba_adj < 2; ++i_ba_adj) {
      if (not info_ba.adj_set[i_ba_adj]) {
        continue;
      }

      // if the adjacent sses are identical and the type of bridge is the same.
      if (info_ba.adj_residues[i_ba_adj].sse_id == adj_sse_id and
          info_t.bridge_type[i_t_adj] == info_ba.bridge_type[i_ba_adj]) {

          pdb::warning((boost::format("IGNORING CONTINUOUS BRIDGES OF THE SAME TYPE: "
                                      "BETWEEN TARGET RESIDUE and SSE[%d], RESNUM = %d")
                                      % +ba.sse_id
                                      % ba.resnum).str());
        return empty;
      }
    }
  }

  return ba;
} // protected member function ba_check()

} // namespace sheet

