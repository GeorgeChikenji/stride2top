// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "functions.h"
#include "color.h"

namespace graphviz {

// *****************************************************************************************
// Function adj_list_to_dot()
// *****************************************************************************************

void adj_list_to_dot(std::ostream & ofs, sheet::DirectedAdjacencyList const& adj) {

  constexpr unsigned max_penwith = 5;

  // Write Graph Header
  ofs << "digraph G {\n";

  // Initialize Node list and write

  // index : sse_id in adj_list
  // data  : serial number of the index sse_id
  std::unordered_map<sheet::SubStrand, unsigned, sheet::SubStrandHasher> nodes;

  unsigned nodes_counter = 0;
  auto const n_strands = adj.strand_indices.size();
  for (sheet::IndexType serial_sse_id = 0; serial_sse_id < n_strands; ++serial_sse_id) {
    auto const sse_id = adj.strand_indices[serial_sse_id];

    auto const& substr_vec = adj.substrs().vec(serial_sse_id);

    // if all Sub-Strands of this SSE have been erased
    if (std::distance(substr_vec.begin(), substr_vec.end()) == 0) {
      ofs << nodes_counter << "[label=\"" << +sse_id << "\\n[Erased]"
          << "\", fillcolor=\""
          << color::color_split_blue_red(adj.sses.size, sse_id).to_RGB().hex_str()
          << "99\", style=filled];\n";

      nodes.insert({sheet::SubStrand{serial_sse_id, 0}, nodes_counter});
      ++nodes_counter;
      continue;
    }

    bool const one_substr = std::distance(substr_vec.begin(), substr_vec.end()) == 1;


    for (auto const& sub : substr_vec) {
      nodes.insert({sub, nodes_counter});
      ofs << +nodes_counter << "[label=\"" << +sse_id;

        // if this strand has more than 1 sub-strands
      if (not one_substr) {
        ofs << "-" << +sub.substr;
      }

      ofs << "\\n[" << adj.substrs().n_term_res(sub)
          << " ~ " << adj.substrs().c_term_res(sub) << "]"
          << "\", fillcolor=\""
          << color::color_split_blue_red(adj.sses.size, sse_id).to_RGB().hex_str()
          << "99\", style=filled];\n";

      ++nodes_counter;
    }
  }


  // Find the max count of pair residues
  auto const max_itr = std::max_element(adj.adj_sub().map().cbegin(),
                                        adj.adj_sub().map().cend(),
                                         [](auto const& a, auto const& b) {
                                            return a.second.residue_pairs < b.second.residue_pairs;
                                           });
  // if the list is empty
  if (max_itr == adj.adj_sub().map().cend()) {
    ofs << "}";
    return ;
  }
  double const max_count_residue_pairs = max_itr->second.residue_pairs;

  // Initialize Edge list
  std::unordered_set<sheet::SubStrandsPairKey, sheet::SubStrandsPairKeyHasher> drawn;

  auto const size = adj.sheets.size();
  for (unsigned i = 0; i < size; ++i) {
    for (auto const key: adj.sheets[i].substr_keys()) {

      // True if directed, false otherwise.
      bool const directed = adj.adj_sub().map().count(key.reverse()) == 0;

      // if undirected, and edge with reverse direction has already been drawn.
      if (not directed and drawn.count(key.reverse()) != 0) {
        continue;
      }

      auto const& data = adj.adj_sub().map().at(key);

      // if undirected, the delta value to display is from N-term side SubStrand to C-term.
      std::string delta_str{""};
      delta_str = std::to_string(key.substr0 < key.substr1 ? data.delta_2 : data.delta_1);
      delta_str += ":";
      delta_str += std::to_string(key.substr0 < key.substr1 ? data.delta_1 : data.delta_2);


      ofs << +nodes.at(key.sub0()) << "->"
          << +nodes.at(key.sub1())
          << " [label=\"" << (data.direction ? "" : "Anti-")
          << "Parallel " << delta_str
          << "\", "
          << "labeldistance=2.0, "
          << "penwidth=" << (max_penwith * data.residue_pairs / max_count_residue_pairs)
          << (directed == false ? ", dir=none" : "")
          << "];\n";

      drawn.insert(key);
    }
  } // for each sheet

  // Finalize the graph
  ofs << "}\n";
}

} // namespace mmcif
