// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef BAB_SIDE_H_
#define BAB_SIDE_H_

#include <array>
#include <vector>

#include <Eigen/Core>

#include "pdb/atom.h"
#include "sheet/directed_adjacency_list.h"

namespace bab {

/// Given 3 coordinates p0, p1, p2 return an angle p0p2p1 in radian.
/// Calculation of this angle can be done by following steps.
/// * Generate v0 and v1 which points from p2 to p0 and p1, respectively.
/// * v0 dot v1 (dot product) == |v0| * |v1| * cos(p0p2p1)
/// * p0p2p1 = acos((v0.dot(v1)) / (v0.norm() * v1.norm()))
double angle(Eigen::Vector3d const& p0, Eigen::Vector3d const& p1, Eigen::Vector3d const& p2);


using ATOM_vec_iter = sheet::ATOM_vec_iter;


// *****************************************************************************
// Class Definition
// *****************************************************************************

/// Store the information about the triangulation of strand CA atoms
class Side {
public:

  // *************************************************************************
  // Public Member Functions
  // *************************************************************************

  /// Constructor
  Side(ATOM_vec_iter const a0_, ATOM_vec_iter const a1_,
       ATOM_vec_iter const first, ATOM_vec_iter const last) :
    a0{a0_},
    a1{a1_},
    max_angle{init_max_angle(first, last)},
    normal_vec(init_normal_vec())
  {}

  auto const& opp() const { return std::get<0>(max_angle)->xyz; }
  double get_angle() const { return std::get<1>(max_angle); }


  /// @brief            A method to access the normal vector
  /// @param reversed   Pass true if the main strand (to which a0 and a1 belong) is reversed,
  ///                   false otherwise
  Eigen::Vector3d const& normal(bool const reversed) const;



  /// @brief            Return a reference to the coordinate depend on reversed or not
  /// @param reverse    True if myside is reversed, false otherwise.
  Eigen::Vector3d const& base_point(bool const reversed) const;



  /// @brief            Detect which side the given point v is located.
  /// @param v          A coordinate of target ATOM
  /// @param reversed   Whether myside is reversed or not. Pass true if reversed, false otherwise.
  /// @param myside     If true, return false when v is on the same side as normal_vec. @n
  ///                   If false, return false when v is on the opposite side of normal_vec.
  /// @retval true      v is on the left-handed side
  /// @retval false     v is on the right-handed side
  std::tuple<bool, bool> on_left_side(Eigen::Vector3d const& v,
                                      bool const reversed, bool const myside,
                                      double const min_dist=1.0) const;


  // for prototypes
#ifdef PROTO
  ATOM_vec_iter iter_a0() const { return a0; }
  ATOM_vec_iter iter_a1() const { return a1; }
  ATOM_vec_iter iter_opp() const { return std::get<0>(max_angle); }
#endif // ifdef PROTO


protected:

  // *************************************************************************
  // Protected Member Functions
  // *************************************************************************

  /// @brief  decide which atom on the opposite strand to use as one of the triangle
  /// @param  first  An iterator to the first element of the target atoms.
  /// @param  last   Past the end iterator of the target atoms.
  std::tuple<ATOM_vec_iter, double> init_max_angle(ATOM_vec_iter const first,
                                                   ATOM_vec_iter const last) const;

  /// Initialize the normal_vec
  std::array<Eigen::Vector3d, 2> init_normal_vec() const;



  // *************************************************************************
  // Protected Member Variables
  // *************************************************************************

  /// An iterator to the coordinate data a0.
  ATOM_vec_iter a0;

  /// An iterator to the coordinate data a1.
  ATOM_vec_iter a1;

  /// A tuple of the index and the angle which returned maximum angle
  std::tuple<ATOM_vec_iter, double> max_angle;


  /// An array of normal vectors against a plane which a triangle of a0, a1, opp is on.
  /// normal_vec[0] : an unit vector from a0 pointing to the right-handed side.
  /// normal_vec[1] : an unit vector from a1 pointing to the right-handed side.
  /// This value will be initialized in a member function init_normal_vec().
  std::array<Eigen::Vector3d, 2> normal_vec;
};

} // namespace bab

#endif //ifndef BAB_SIDE_H_

