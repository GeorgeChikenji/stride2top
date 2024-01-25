// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef COLOR_H_
#define COLOR_H_

#include <cmath>
#include <cstdio>

#include <iostream>
#include <stdexcept>
#include <string>

namespace color {

struct RGB {

  /// Construct from a 0 ~ 255 integer
  RGB(unsigned char const r, unsigned char const g, unsigned char const b):
    R{r}, G{g}, B{b} {}

  /// Construct from 0.0 ~ 1.0 double
  RGB(double const r, double const g, double const b):
    R{0.0 <= r and r <= 1.0 ?
      static_cast<unsigned char>(255 * r) : throw std::invalid_argument{"Invalid r"}},
    G{0.0 <= r and r <= 1.0 ?
      static_cast<unsigned char>(255 * g) : throw std::invalid_argument{"Invalid g"}},
    B{0.0 <= r and r <= 1.0 ?
      static_cast<unsigned char>(255 * b) : throw std::invalid_argument{"Invalid b"}} {}

  std::string hex_str() const {
    char hexstring[] = "#000000";
    std::sprintf(hexstring, "#%02x%02x%02x", R, G, B);
    return std::string(hexstring);
  }

  unsigned char R{0};
  unsigned char G{0};
  unsigned char B{0};
};


struct HSV {

  HSV(double const h, double const s, double const v):
    H{0.0 <= h and h < 1.0 ? h : throw std::invalid_argument{"Invalid h"}},
    S{0.0 <= s and s <= 1.0 ? s : throw std::invalid_argument{"invalid s"}},
    V{0.0 <= v and s <= 1.0 ? v : throw std::invalid_argument{"invalid v"}} {}


  /// Convert this object into RGB and return the RGB class object
  RGB to_RGB() const {
    if (S == 0.0) {
      return RGB{V, V, V};
    }

    unsigned const i = H * 6;
    double const f = H * 6.0 - i;

    double const p = V * (1.0 - S);
    double const q = V * (1.0 - (S * f));
    double const t = V * (1.0 - (S * (1.0 - f)));

//    std::clog << "H = " << H << "\n"
//              << "S = " << S << "\n"
//              << "V = " << V << "\n"
//              << "i = " << i << "\n"
//              << "f = " << f << "\n"
//              << "p = " << p << "\n"
//              << "q = " << q << "\n"
//              << "t = " << t << "\n" ;

    switch(i) {
      case 0:
        return RGB{V, t, p};
      case 1:
        return RGB{q, V, p};
      case 2:
        return RGB{p, V, t};
      case 3:
        return RGB{p, q, V};
      case 4:
        return RGB{t, p, V};
      case 5:
      default:
        return RGB{V, p, q};
    };
  } // Public member function to_RGB

  double H{0.0};
  double S{0.0};
  double V{0.0};
};


/// @param N_node Number of nodes to split.
/// @param i_node the sequential number of the target node. (First is 0, last is N_node - 1)
HSV color_split_blue_red(unsigned const N_node, unsigned const i_node) {

  // if i_node exceeds N_node
  if (N_node <= i_node) {
    throw std::invalid_argument("Invalid i_node");
  }

  // Edge values of Hue
  constexpr double const H_RED = 0.0;
  constexpr double const H_BLUE = 240.0 / 360.0;

  // default values of Saturation and Value
  constexpr double const S = 1.0;
  constexpr double const V = 1.0;


  // if there is only one node
  if (N_node == 1) {
    return HSV{H_BLUE, S, V};
  }

  if (i_node == N_node - 1) {
    return HSV{H_RED, S, V};
  }

  double const step = H_BLUE / (N_node - 1);
  return HSV{H_BLUE - step * i_node, S, V};
}

} // namespace color

#endif // ifndef COLOR_H_
