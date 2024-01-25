// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "pdb/tools.h"
#include "pdb/exceptions.h"

#include "sheet/adj_list_with_sub.h"
#include "sheet/cb_side.h"
#include "sheet/directed_adjacency_list.h"
#include "sheet/exceptions.h"


namespace sheet {

// **********************************************************************************
// Member Functions for class DirectedAdjacencyList
// **********************************************************************************



// *******************************************************************************
// Public Member Function atom_cbegin()
// *******************************************************************************

ATOM_vec_iter DirectedAdjacencyList::atom_cbegin(SubStrand const& ss) const {
  auto const& sse = sses[strand_indices[ss.str]];
  auto const pos = sub_strands_range.n_term_res(ss) - sse.init;
  assert(0 <= pos and static_cast<unsigned>(pos) < sse.atoms.size());
  return sse.atoms.cbegin() + pos;
} // public member function atom_cbegin()



// *******************************************************************************
// Public Member Function atom_cend()
// *******************************************************************************

ATOM_vec_iter DirectedAdjacencyList::atom_cend(SubStrand const& ss) const {
  auto const& sse = sses[strand_indices[ss.str]];
  auto const pos = sub_strands_range.c_term_res(ss) - sse.init;
  assert(0 <= pos and static_cast<unsigned>(pos) < sse.atoms.size());
  return sse.atoms.cbegin() + pos + 1;
} // public member function atom_cend()




// *******************************************************************************
// Public Member Function edge_info()
// *******************************************************************************


SubStrandsPairNode const& DirectedAdjacencyList::edge_info(IndexType const sse_id_0,
                                                           IndexType const sse_id_1,
                                                           IndexType const substrand_id_0,
                                                           IndexType const substrand_id_1) const {
  return adj_list_with_sub.map(SubStrandsPairKey{sses.serial_strand_id[sse_id_0],
                                                 substrand_id_0,
                                                 sses.serial_strand_id[sse_id_1],
                                                 substrand_id_1});
}


// *******************************************************************************
// Public Member Function substr()
// *******************************************************************************
SubStrand DirectedAdjacencyList::substr(IndexType const sse_id,
                                        IndexType const substr_id) const {
  return SubStrand{sses.serial_strand_id[sse_id], substr_id};
} // public member function substr()



// *******************************************************************************
// Public Member Function search()
// *******************************************************************************

StrandsPairAttribute DirectedAdjacencyList::search(SubStrand const& ss0,
                                                   SubStrand const& ss1) const {
  return search(ss0, ss1, adj_list_with_sub.adj_substr_vec());
} // public member function search()




// *******************************************************************************
// Public Member Function search()
// *******************************************************************************
StrandsPairAttribute DirectedAdjacencyList::search(SubStrand const& ss0,
                                                   SubStrand const& ss1,
                                                   AdjSubVec const& adj_sub_vec) const {
  StrandsPairAttribute ret{ss0, ss1};
  auto const parents = search_bfs(ss0, ss1, adj_sub_vec);

  // if not reachable
  if (parents.size() == 0) {
    return ret;
  }

  ret.reachable = true;
  search_backtrace(ret, parents);
  return ret;
} // public member function search()



// *******************************************************************************
// Public Member Function search()
// *******************************************************************************
StrandsPairAttribute DirectedAdjacencyList::search(IndexType const sse_id_0,
                                                   IndexType const sse_id_1) const {
  return search(substr(sse_id_0), substr(sse_id_1));
} // public member function search()




// *******************************************************************************
// Public Member Function attr()
// *******************************************************************************
StrandsPairAttribute const& DirectedAdjacencyList::attr(SubStrand const& ss0,
                                                        SubStrand const& ss1) const {
  SubStrandsPairKey const key{ss0, ss1};
  return adj_attr.at(key).reachable ? adj_attr.at(key) : adj_attr.at(key.reverse());
} // public member function attr()




// *******************************************************************************
// DEBUG Public Member Function
// *******************************************************************************

#if defined(DEBUG) && defined(LOG_SHEET)
void DirectedAdjacencyList::print_list(AdjList const& adj_list_, std::string const& title) const {
  std::clog << "\n";
  pdb::debug_log(title);
  pdb::debug_log("TARGET : PAIRED (   DIRECTION , COUNT)");

  for (auto const sse_id_0 : strand_indices) {
    std::string out = (boost::format("%6d :    ") % +sse_id_0).str();
    for (auto const& sse_id_1 : strand_indices) {
      auto const key = std::make_pair(sse_id_0, sse_id_1);

      // if doesn't exist
      if (adj_list_.count(key) == 0) {
        continue;
      }

      out += (boost::format("%3d (%13s, %5d), ") % +adj_list_.at(key).id
                                                 % (adj_list_.at(key).direction ?
                                                    "parallel" : "anti-parallel")
                                                 % adj_list_.at(key).count).str();
    }
    pdb::debug_log(out);
  }
  std::clog << std::endl;
}
#endif // ifdef DEBUG



// *************************************************************************************
// Protected Member Function init_list()
// *************************************************************************************

AdjListWithSub DirectedAdjacencyList::init_list(Pairs const& pairs) {

  // Initialize Undirected Adjacency List

  #ifdef LOGGING
  pdb::log("INITIALIZING: Undirected Adjacency List");
  #endif

  std::unordered_map<IndexType, std::unordered_set<IndexType>> undirected_adj_index_map;

  AdjList undirected_adj_list;

  IndexType sse_serial_id = 0;
  for (auto const target_sse_id : strand_indices) {
    auto const& involved_pairs = pairs.involved_pairs.at(sse_serial_id);

    // Add undirectred adjacency list
    auto const adj_index_set = add_undirected_adj_list(undirected_adj_list, involved_pairs,
                                                       sses[target_sse_id]);
    // Add the adjacent set of only ids
    undirected_adj_index_map.insert({target_sse_id, adj_index_set});

    ++sse_serial_id;
  }


  // debug print
  #if defined(DEBUG) && defined(LOG_SHEET)
  print_list(undirected_adj_list, "UNDIRECTED ADJACENCY LIST");
  #endif // ifdef DEBUG


  // valid check
  undirected_adj_symmetry_check(undirected_adj_list);

  // *********************************************************
  // Initialize Directed Adjacency List
  // *********************************************************
  #ifdef LOGGING
  pdb::log("INITIALIZING: Directed Adjacency List");
  #endif

  // indices of strands
  auto strands = strand_indices;


  // Decide the strict zone
  while (strands.size() != 0) {
    gen_strict_zone(undirected_adj_list, pairs, strands, undirected_adj_index_map);
  }



  #if defined(LOGGING) && defined(LOG_SHEET)
  strict_zone.PrintZone();
  #endif

  // Determine the side of CB atom for all defined residues
  // and generates a directed adjacency list with Sub-Strands
  return strict_zone.decide_side(undirected_adj_list);
}




// *************************************************************************************
// Protected Member Function add_undirected_adj_list()
// *************************************************************************************

/// create a undirected adjacency list
std::unordered_set<IndexType> DirectedAdjacencyList::add_undirected_adj_list(AdjList & undirected_adj_list, Pairs::PairsVec const& involved_pairs, pdb::SSE const& e) const {

  #if defined(DEBUG) && defined(LOG_SHEET)
  pdb::debug_log("TARGET: SSE[" + std::to_string(e.index) + "]");
  #endif // ifdef DEBUG

  auto adj_ids = adj_id_set(involved_pairs);



  // initialize the last pairs and add the empty AdjStrandData to adj_list
  std::unordered_map<IndexType, LastPair> last_pairs;
  std::unordered_map<IndexType, std::unordered_map<bool, unsigned>> direction_counter;

  for (auto const i : adj_ids) {
    last_pairs.insert({i, LastPair()});
    undirected_adj_list.insert({std::make_pair(e.index, i), AdjStrandData(i, true)});
    direction_counter.insert({i, std::unordered_map<bool, unsigned>({{true, 0}, {false, 0}})});
  }

  // for each hbond
  for (auto const pair : involved_pairs) {

    #if defined(DEBUG) && defined(LOG_SHEET)
    pdb::debug_log("pair = " + std::to_string(pair[0]) + ", " +
                   std::to_string(pair[1]) + ", " + std::to_string(pair[2]));
    #endif

    // get the SSE_ID of the paired residue (condition: this SSE is a strand)
    IndexType sse_id;
    try {
      sse_id = sses.sse_ind_of(pair[1], 'E', 1, true, pair[2] == 0 ? 'C' : 'N');

    // if not found (loop residue)
    } catch (pdb::non_sse_resnum) {

      #if defined(DEBUG) && defined(LOG_SHEET)
      pdb::debug_log(" : NON_SSE_RESNUM " + std::to_string(pair[1]) + "\n");
      #endif // ifdef DEBUG

      continue;
    }

    // if the target is a HELIX
    if (sses[sse_id].type == 'H') {
      continue;
    }

    #if defined(DEBUG) && defined(LOG_SHEET)
    pdb::debug_log(" SSE[" + std::to_string(sse_id) + "]\n");
    #endif // ifdef DEBUG


    // set key to access the AdjList
    auto const key = std::make_pair(e.index, sse_id);

    // if this is not the first time (last_pairs stores the data for the last hbond)
    if (!last_pairs[sse_id].first) {

      // store the temporal direction
      bool tmp_direction;

      if (last_pairs[sse_id].hbond_dir == pair[2]) {
#if defined(LOGGING) && defined(LOG_SHEET)
        pdb::log((boost::format("CONSECUTIVE DIRECTION OF HBONDS: "
                                "%2d %3s %2d <==> %2d %3s %2d")
                                % last_pairs[sse_id].res0
                                % (pair[2] == 0 ? "-->" : "<--")
                                % last_pairs[sse_id].res1
                                % pair[0]
                                % (pair[2] == 0 ? "-->" : "<--")
                                % pair[1]).str());
#endif
        continue;
      }


      // res 0 unchanged
      if (last_pairs[sse_id].res0 == pair[0]) {
        if (last_pairs[sse_id].res1 < pair[1]) {
          tmp_direction = true;
        } else if (last_pairs[sse_id].res1 == pair[1]) {
          tmp_direction = false;
        } else {
          tmp_direction = false;
        }
      // res 0 increased
      } else {
        if (last_pairs[sse_id].res1 < pair[1]) {
          tmp_direction = true;
        } else if (last_pairs[sse_id].res1 == pair[1]) {
          // direction of hbond: res1 -> res0 ==>> res0 -> res1
          if (last_pairs[sse_id].hbond_dir == 1 and pair[2] == 0) {
            tmp_direction = true;
          } else if (last_pairs[sse_id].hbond_dir == 0 and pair[2] == 1) {
            tmp_direction = false;
          } else {
            throw std::runtime_error("This pattern is not considered!!");
          }
        } else {
          tmp_direction = false;
        }
      }

      ++direction_counter[sse_id][tmp_direction];
    }

    ++undirected_adj_list[key].count;
    last_pairs[sse_id] = LastPair(pair[0], pair[1], pair[2], false);
  }


  // decide directions
  for (auto const sse_id : adj_ids) {
    auto const key = std::make_pair(e.index, sse_id);
    undirected_adj_list[key].direction = direction_by_majority(sse_id, e,  direction_counter);
  }


  // remove if the count is 0 or 1
  auto const last_adj_id = adj_ids.end();
  for (auto itr = adj_ids.begin(); itr != last_adj_id;) {

    // set key to access the AdjList
    auto const key = std::make_pair(e.index, *itr);

    // if the count is too small
    if (undirected_adj_list[key].count < 2) {

      #if defined(LOGGING) && defined(LOG_SHEET)
      pdb::log((boost::format("UNPAIRED: SSE_ID = (%2d, %2d), COUNT = %d")
                          % +e.index
                          % +*itr
                          % undirected_adj_list[key].count).str());
      #endif // ifdef LOGGING

      undirected_adj_list.erase(key);
      itr = adj_ids.erase(itr);
      continue;
    }

    // if not erased

    #if defined(DEBUG) && defined(LOG_SHEET)
    pdb::debug_log((boost::format("ADJACENT PAIR: "
                             "SSE_ID = %2d <==> %2d, DIR = %13s, COUNT = %3d")
                            % +e.index
                            % +*itr
                            % (undirected_adj_list[key].direction ? "parallel" : "anti-parallel")
                            % undirected_adj_list[key].count).str());
    #endif // ifdef DEBUG

    // increment iterator if not erased
    ++itr;
  }

  // return a set of SSE_IDs which is adjacency to e
  return adj_ids;
} // Protected Member Function add_undirected_adj_list()




// *************************************************************************************
// Protected Member Function adj_id_set()
// *************************************************************************************

/// generate a vector of sse_ids adjacent to a target SSE
std::unordered_set<IndexType> DirectedAdjacencyList::adj_id_set(Pairs::PairsVec const& involved_pairs) const {
  std::unordered_set<IndexType> ret;
  for (auto const& pair : involved_pairs) {
    try {
      // search the sse_id of pair[1] for strands
      // (if Helix, pdb::non_sse_resnum will be thrown.)
      auto const sse_id = sses.sse_ind_of(pair[1], 'E', 1, true, pair[2] == 0 ? 'C' : 'N');

      ret.insert(sse_id);
    } catch (pdb::non_sse_resnum) {
      continue;
    }
  }
  return ret;
} // Protected Member Function adj_id_set()




// *************************************************************************************
// Protected Member Function direction_by_majority()
// *************************************************************************************

/// Given the number of votes for parallel and anti-parallel, return the direction
bool DirectedAdjacencyList::direction_by_majority(IndexType const sse_id, pdb::SSE const& e,
                          std::unordered_map<IndexType, std::unordered_map<bool, unsigned>> direction_counter) const {

  bool const dir = direction_counter[sse_id][false] < direction_counter[sse_id][true];

  // check the ratio
  double const max = std::max(direction_counter[sse_id][false], direction_counter[sse_id][true]);
  double const min = std::min(direction_counter[sse_id][false], direction_counter[sse_id][true]);

  // the ratio of minority against the entire hbonds
  auto const frac = min / (max+min);
  if (0.5 <= frac) {
    throw std::runtime_error("[ERROR] CANNOT DECIDE DIRECTION : SSE PAIR = (" +
        std::to_string(e.index) + ", " + std::to_string(sse_id) + "), parallel = " +
        std::to_string(direction_counter[sse_id][true]) + ", anti-parallel = " +
        std::to_string(direction_counter[sse_id][false]));
  }
  return dir;
}




// *************************************************************************************
// Protected Member Function undirected_adj_symmetry_check()
// *************************************************************************************
void DirectedAdjacencyList::undirected_adj_symmetry_check(AdjList const& undirected_adj_list) const {

  auto const size = strand_indices.size();
  for (IndexType i = 0; i < size; ++i) {
    for (IndexType j = i + 1; j < size; ++j) {

      auto const key_0 = std::make_pair(strand_indices[i], strand_indices[j]);
      auto const key_1 = std::make_pair(strand_indices[j], strand_indices[i]);

      auto const count_0 = undirected_adj_list.count(key_0);
      auto const count_1 = undirected_adj_list.count(key_1);

      // if no such strands connection
      if (count_0 == 0 and count_1 == 0) {
        continue;
      }

      // existence check
      if (count_0 == 0) {
        throw non_symmetric(strand_indices[i], strand_indices[j]);
      }
      if (count_1 == 0) {
        throw non_symmetric(strand_indices[j], strand_indices[i]);
      }


      // direction check
      if (! undirected_adj_list.at(key_0).symmetry(strand_indices[j],
                                                   undirected_adj_list.at(key_1))) {
        throw non_symmetric_dir_opposite(strand_indices[i], strand_indices[j]);
      }
      if (! undirected_adj_list.at(key_1).symmetry(strand_indices[i],
                                                   undirected_adj_list.at(key_0))) {
        throw non_symmetric_dir_opposite(strand_indices[j], strand_indices[i]);
      }
    }
  }
} // Protected Member Function undirected_adj_symmetry_check()





// *************************************************************************************
// Protected Member Function gen_strict_zone()
// *************************************************************************************

void DirectedAdjacencyList::gen_strict_zone(AdjList const& undirected_adj_list,
                                            Pairs const& pairs, std::vector<IndexType>& strands,
                                            std::unordered_map<IndexType, std::unordered_set<IndexType>> const& adj_index_map) {

  std::unordered_set<IndexType> finished;
  std::deque<IndexType> queue;
  queue.push_back(strands.front());

  while (queue.size() != 0) {
    auto const target_str = queue.front();
    queue.pop_front();

    // If the input is a barrel, the last sse_id may be pushed into the queue twice.
    // If target_str is already finished, just continue.
    if (finished.count(target_str) == 1) {
      continue;
    }

    #if defined(LOGGING) && defined(LOG_SHEET)
    pdb::debug_log("TARGET = " + std::to_string(target_str));
    #endif // ifdef DEBUG


    // Mark target_str as finished. (Insert target_str into the finished set.)
    finished.insert(target_str);

    std::unordered_map<IndexType, LastPair> last_pairs;
    for (auto const i : adj_index_map.at(target_str)) {
      last_pairs.insert({i, LastPair()});
      queue.push_back(i);
    }

    for (auto const& hbond : pairs.resort_involved_pairs(sses.serial_strand_id[target_str],
                                                         undirected_adj_list, sses)) {

      #if defined(LOGGING) && defined(LOG_SHEET)
      pdb::log((boost::format("GEN_STRICT_ZONE: PROCESSING %3d %3d : %d") % hbond[0] % hbond[1] % hbond[2] ).str());
      #endif

      // get the SSE_ID of the paired residue
      IndexType adj_str;
      try {
        adj_str = sses.sse_ind_of(hbond[1], 'E', 1, true, hbond[2] == 0 ? 'C' : 'N');

      // if not found (loop residue)
      } catch (pdb::non_sse_resnum) {
        continue;
      }

      // if the adj_str is unpaired or not strand
      if (adj_index_map.at(target_str).count(adj_str) == 0) {
        continue;
      }

      
      // set key to access the AdjList
      auto const key = std::make_pair(target_str, adj_str);

      // Parallel(true) or Anti-Parallel(false)
      auto const strand_direction = undirected_adj_list.at(key).direction;

      // if this is not the first time (last_pairs stores the data for the last hbond)
      if (!last_pairs[adj_str].first) {
        // Parallel
        if (strand_direction) {
          if (hbond[0] == last_pairs[adj_str].res0 and hbond[1] == last_pairs[adj_str].res1 + 2
              and hbond[2] != last_pairs[adj_str].hbond_dir) {
            strict_zone.On(target_str, hbond[0], adj_str, hbond[1] - 1, true,
                           ZoneInfo::BridgeType::ParallelHbonds);
          } else if (hbond[0] == last_pairs[adj_str].res0 + 2 and hbond[1] == last_pairs[adj_str].res1
                     and hbond[2] != last_pairs[adj_str].hbond_dir) {
            strict_zone.On(target_str, hbond[0] - 1, adj_str, hbond[1], false,
                           ZoneInfo::BridgeType::ParallelNoHbonds);
          }

          
        // Anti-Parallel
        } else {
          // Small Ring
          if (hbond[0] == last_pairs[adj_str].res0 and hbond[1] == last_pairs[adj_str].res1
              and hbond[2] != last_pairs[adj_str].hbond_dir) {
            strict_zone.On(target_str, hbond[0], adj_str, hbond[1], true,
                           ZoneInfo::BridgeType::SmallRing);
          } else if (hbond[0] == last_pairs[adj_str].res0 + 2
                     and hbond[1] == last_pairs[adj_str].res1 - 2
                     and hbond[2] != last_pairs[adj_str].hbond_dir) {
            strict_zone.On(target_str, hbond[0] - 1, adj_str, hbond[1] + 1, false,
                           ZoneInfo::BridgeType::LargeRing);
          }

        }
      }
      last_pairs[adj_str] = LastPair(hbond[0], hbond[1], hbond[2], false);
    }


    // Delete the finished strand from the strands vector
    auto const to_erase = std::find(strands.begin(), strands.end(), target_str);
    if (to_erase != strands.cend()) {
      strands.erase(to_erase);
    }
  } // while (queue.size() != 0)
} // Protected Member Function gen_strict_zone()



// *************************************************************************************
// Protected Member Function translate_sub()
// *************************************************************************************

AdjList DirectedAdjacencyList::translate_sub() const {
  AdjList tmp_adj_list;

  for (auto const& adj_pair : adj_list_with_sub.map()) {
    auto const new_key = std::make_pair(strand_indices[adj_pair.first.str0],
                                        strand_indices[adj_pair.first.str1]);
    if (tmp_adj_list.count(new_key) != 0) {
      if (tmp_adj_list.at(new_key).count < adj_pair.second.residue_pairs * 2) {
        tmp_adj_list[new_key] = AdjStrandData{new_key.second, adj_pair.second.direction,
                                          adj_pair.second.residue_pairs * 2};
      }
    } else {
      tmp_adj_list.insert({new_key, AdjStrandData{new_key.second, adj_pair.second.direction,
                                                              adj_pair.second.residue_pairs * 2}});
    }
  }

  // debug print
  #if defined(DEBUG) && defined(LOG_SHEET)
  print_list(tmp_adj_list, "TRADITIONAL DIRECTED ADJACENCY LIST");
  #endif

  return tmp_adj_list;
} // protected member function translate_sub()




// *************************************************************************************
// Protected Member Function init_adj_index_list()
// *************************************************************************************

DirectedAdjacencyList::AdjIndList DirectedAdjacencyList::init_adj_index_list() const {
  AdjIndList ind_list;

  for (auto const& strand_pair : adj_list) {
    if (ind_list.count(strand_pair.first.first) == 0) {
      ind_list.insert({strand_pair.first.first, std::vector<IndexType>()});
    }
    ind_list[strand_pair.first.first].push_back(strand_pair.first.second);
  }

  return ind_list;
}



// *************************************************************************************
// Protected Member Function init_sheets()
// *************************************************************************************
Sheets DirectedAdjacencyList::init_sheets() {
  Sheets tmp_sheets{};

  for (auto const sub0 : sub_strands_range.vec()) {
    for (auto const sub1 : sub_strands_range.vec()) {
      SubStrandsPairKey const key{sub0, sub1};
      auto const attr = search(sub0, sub1);
      if (attr.reachable) {
        tmp_sheets.add(key, attr);
      }
    }
  }

  /// Initialize cycles inside sheets
  tmp_sheets.cycle_check(adj_list_with_sub.map());

  // Fix undirected paths that is out of cycles. (Such as d3vlaa_)
  adj_substrands.fix_undirected_paths(adj_list_with_sub, tmp_sheets);

  tmp_sheets.add_key_vec(adj_list_with_sub);
  tmp_sheets.sort_sheets();

  return tmp_sheets;
} // protected member function init_sheets()




// *************************************************************************************
// Protected Member Function init_adj_attr()
// *************************************************************************************

DirectedAdjacencyList::AdjAttrMap DirectedAdjacencyList::init_adj_attr() const {
  AdjAttrMap map;

  for (auto const sub0 : sub_strands_range.vec()) {
    for (auto const sub1 : sub_strands_range.vec()) {
      SubStrandsPairKey const key{sub0, sub1};
      auto const attr = search(sub0, sub1);
      map.insert({key, attr});
    }
  }

  return map;
} // protected member function init_adj_attr()



// *************************************************************************************
// Protected Member Function search_backtrace()
// *************************************************************************************

void DirectedAdjacencyList::search_backtrace(StrandsPairAttribute & attr,
                                             PathParents const& parents) const {
  auto last = attr.ss1;
  while (last != attr.ss0) {
    ++attr.jump;

    auto const next = parents.at(last);
    // Only next -> last is available because adj_list_with_sub is directed.
    auto const& data = adj_list_with_sub.map(SubStrandsPairKey{next, last});

    attr.direction = attr.direction == data.direction;

    attr.jumped_substrs.push_back(next);

    last = next;
  }
  --attr.jump;

} // protected member function search_backtrace()




// *************************************************************************************
// Protected Member Function search_bfs()
// *************************************************************************************

DirectedAdjacencyList::PathParents DirectedAdjacencyList::search_bfs(SubStrand const& first, SubStrand const& last, AdjSubVec const& adj_sub_vec) const {
  PathParents parents;
  SubStrandSet examined;
  std::deque<SubStrand> queue;
  queue.push_back(first);

  while (queue.size() != 0) {
    auto const node = queue.front();
    queue.pop_front();
    examined.insert(node);

    // if found
    if (node == last) {
      return parents;
    }

    // continue if no adjacent strand
    if (adj_sub_vec.count(node) == 0) {
      continue;
    }

    for (auto const& adj_key : adj_sub_vec.at(node)) {
      auto const adj = adj_key.sub1();

      // if already examined
      if (examined.count(adj) == 1) {
        continue;
      }
      queue.push_back(adj);
      parents.insert({adj, node});
    }
  }
  return PathParents();
} // protected member function search_bfs()


} // namespace rperm


