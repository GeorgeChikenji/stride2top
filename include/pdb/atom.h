// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef PDB_ATOM_H_
#define PDB_ATOM_H_

#include <string>

#include <Eigen/Core>

namespace pdb {

// ************************************************************************
// Class ATOM
// ************************************************************************
class ATOM{
public:

  // *****************************************************************************
  // Public Member Functions
  // *****************************************************************************

  ATOM() = default;

  /// @brief Constructs from a PDB line string
  /// @param line ATOM line in a PDB file
  /// @exception std::out_of_range The line is too short.
  /// @exception std::invalid_argument No str-to-double convertion has occoured.
  explicit ATOM(std::string const& line):
    xyz{std::stod(line.substr(30, 8)),
        std::stod(line.substr(38, 8)),
        std::stod(line.substr(46, 8))},
    pdb{true} {}

  /// Constructs from the xyz coordinates.
  ATOM(double const x, double const y, double const z):
    xyz{Eigen::Vector3d(x, y, z)}, pdb{false} {}

  /// Constructs from a 3D Eigen vector
  explicit ATOM(Eigen::Vector3d const& xyz_):
    xyz{xyz_}, pdb{false} {}


  bool operator==(ATOM const& other) const {
    return xyz == other.xyz and pdb == other.pdb;
  }

  bool operator!=(ATOM const& other) const {
    return ! (*this == other);
  }



  // *****************************************************************************
  // Public Member Variables
  // *****************************************************************************

  /// xyz coordinates of this atom
  Eigen::Vector3d const xyz{0.0, 0.0, 0.0};

  /// whether this atom has real pdb coordinates.
  /// If true, this is constructed throught the first version of constructor.
  bool const pdb{false};

}; // class ATOM


} // namespace pdb

#endif // ifndef PDB_ATOM_H_

