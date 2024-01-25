// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PDB_STRIDE_STREAM_H_
#define PDB_STRIDE_STREAM_H_

#include <sstream>
#include "tools.h"
#include "exceptions.h"

namespace pdb {

/// Store the std::stringstream object that contains the stride output
class stride_stream {
public:

  /// Construct Empty Object
  stride_stream() = default;

  /// Read directly from a stream.
  stride_stream(std::istream & is):
    empty{false}, ss{} {
      is2ss(is, ss);
    }

  /// Stride file mode
  explicit stride_stream(std::string const& stride_file):
    empty{false}, ss{} {
    std::ifstream ifs;
    open_input(ifs, stride_file);
    is2ss(ifs, ss);
  }

  /// True if empty.
  bool const empty{true};

  /// Data containing all stride output.
  std::stringstream ss{};
};

static stride_stream empty_stride_stream{};

/// @brief  Run stride command for pdb_file and redirect the results into /tmp directory,
///         then generates a stride_stream class object and return it.
/// @return stride_stream
/// @param  pdb_file a file path to the input pdb file
stride_stream pdb2stride_stream(std::string const& pdb_file);

} // namespace pdb

#endif // ifndef PDB_STRIDE_STREAM_H_
