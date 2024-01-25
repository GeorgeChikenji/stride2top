// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <functional>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "sheet/adj_list_with_sub.h"
#include "sheet/directed_adjacency_list.h"

#include "data_store.h"
#include "table.h"

namespace out {

/// @brief  Output Sub-Strand as a string.
class substr2str {
public:
  explicit substr2str(std::shared_ptr<sheet::DirectedAdjacencyList> const& p): ptr{p} {}

  /// Wrapper operator for static member function str()
  std::string operator()(sheet::SubStrand const& ss) const;

  /// A static member function to get a string that represents the given SubStrand 'ss'.
  static std::string str(sheet::SubStrand const& ss,
                         sheet::DirectedAdjacencyList const& adj);

private:
  std::shared_ptr<sheet::DirectedAdjacencyList> ptr;
};


/// @brief  Concatenate the elements in a container into one string
template <class Iter, class UnaryOperator>
std::string join(Iter first, Iter last, std::string const& delm,
                 UnaryOperator up = [](Iter a){return *a;}) {
  std::string joined = up(first);
  ++first;
  for (auto iter = first; iter != last; ++iter) {
    joined += delm + up(iter);
  }
  return joined;
}

} // namespace out




// ****************************************************************************
// Graphviz
// ****************************************************************************
namespace graphviz {

/// @brief  Output a Dot file of the given adjacency list into the ofstream.
///         The line width of the edge will be thick if there are many hbonds, otherwise thin.
/// @param  ofs       ostream to write the dot file.
///                   This also can be a file, cout, cerr or clog.
/// @param  adj_list  Target adjacency list to output.
void adj_list_to_dot(std::ostream & ofs, sheet::DirectedAdjacencyList const& adj);

} // namespace graphviz

// ****************************************************************************
// Cycles
// ****************************************************************************
namespace cycles {

/// Output All Cycles Pathes to the out stream.
void output_cycles(table::Table<table::Cycle> & tbl,
                   sheet::DirectedAdjacencyList const& adj);

// Helper Function for output_cycles()
using CycleMembers = std::vector<sheet::SubStrand>;
using CyclesVec = std::vector<std::tuple<std::size_t, CycleMembers>>;

/// Given an adjacency list, return a vector of lists of Sub-Strands.
CyclesVec gen_cycles_vec(sheet::DirectedAdjacencyList const& adj);

} // namespace cycles

// ****************************************************************************
// Common Utility
// ****************************************************************************

namespace mmcif {

class mmcif_like {
public:
  // **********************************************************************
  // Public Member Functions
  // **********************************************************************

  /// Constructor: Initializes this->class_name.
  mmcif_like(std::ostream& out_stream, std::string const& name):
    class_name{name}, key_head{"_" + class_name + "."}, os{out_stream} {}

  /// Given a key-value set, output a mmcif-like key-value string to the ostream.
  template<typename T>
  void key_value(std::string const& key, T value) const {
    os << "#\n"
       << key_head << key << "\t"
       << value << "\n";
  }

  /// Given a vector of keys, output a mmcif-like loop syntax to the ostream.
  void loop_head(std::vector<std::string> const& keys) const;

  // **********************************************************************
  // Public Member Variables
  // **********************************************************************

  /// The name of this output class.
  std::string const class_name{""};

  /// A class name string with a prefix and a suffix. ("_classname.")
  std::string const key_head{""};
  std::ostream& os;
};

} // namespace mmcif



namespace rpo {
int get_resnum(pdb::SSES const& sses, std::size_t const serial_str_id, std::size_t const serial_res_id);
void residue_pair_out(table::TBLResiduePair & tbl,
                      sheet::DirectedAdjacencyList const& adj);


/// Exception class
class invalid_zone_info_exception: public pdb::fatal_error_base {
public:
  explicit invalid_zone_info_exception(std::string const& _msg): pdb::fatal_error_base{_msg} {}
};

} // namespace rpo

#endif // ifndef FUNCTIONS_H_
