// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "bab/side.h"

namespace bab {

// *****************************************************************************
// Function Definition angle()
// *****************************************************************************

// Unit of retern value is RADIAN
double angle(Eigen::Vector3d const& p0, Eigen::Vector3d const& p1, Eigen::Vector3d const& p2) {
  auto const v0 = p0 - p2;  // a vector from p2 to p0
  auto const v1 = p1 - p2;  // a vector from p2 to p1
  return std::acos(v0.dot(v1) / (v0.norm() * v1.norm()));
}


// *****************************************************************************
// class Side
// *****************************************************************************
//

// *************************************************************************
// Public Member Function normal()
// *************************************************************************

Eigen::Vector3d const& Side::normal(bool const reversed) const {
  return reversed ? normal_vec[1] : normal_vec[0];
} // public member function normal()



// *************************************************************************
// Public Member Function base_point()
// *************************************************************************

Eigen::Vector3d const& Side::base_point(bool const reversed) const {
  return reversed ? a1->xyz : a0->xyz;
} // public member function base_point()


// *************************************************************************
// Public Member Function on_left_side()
// *************************************************************************

std::tuple<bool, bool> Side::on_left_side(Eigen::Vector3d const& v,
                                          bool const reversed, bool const myside,
                                          double const min_dist) const {
  auto const dist = normal(reversed).dot(v - base_point(reversed));
  bool const is_distant = min_dist < std::fabs(dist);
 if (myside) {
    return std::make_tuple(dist < - min_dist, is_distant);
  } else {
    return std::make_tuple(dist > min_dist, is_distant);
  }
} // public member function on_left_side()


// *************************************************************************
// Protected Member Function init_max_angle()
// *************************************************************************
std::tuple<ATOM_vec_iter, double> Side::init_max_angle(ATOM_vec_iter const first,
                                                       ATOM_vec_iter const last) const {
  double tmp_max_angle = 0.0;
  ATOM_vec_iter tmp_max_iter = first;

  for (auto iter = first; iter != last; ++iter) {
    // Ignore the padding ATOMs
    if (iter->pdb == false) {
      continue;
    }


    auto const tmp_angle = angle(a0->xyz, a1->xyz, iter->xyz);
    if (tmp_max_angle < tmp_angle) {
      tmp_max_angle = tmp_angle;
      tmp_max_iter = iter;
    }
  }
  return std::make_tuple(tmp_max_iter, tmp_max_angle);
} // protected member function init_max_angle()



// *************************************************************************
// Protected Member Function init_normal_vec()
// *************************************************************************

std::array<Eigen::Vector3d, 2> Side::init_normal_vec() const {
  return {{
            (opp() - a0->xyz).cross(a1->xyz - a0->xyz).normalized(), // for sequential
            (opp() - a1->xyz).cross(a0->xyz - a1->xyz).normalized()  // for reversed
          }};
} // protected member function init_normal_vec()

} // namespace bab


