// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef TABLE_H_
#define TABLE_H_

#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/format.hpp>

namespace table {

// ********************************************************************
// Data type typedefs
// ********************************************************************

using SubStrand = std::tuple<std::string, std::size_t, int, int>;

using Helix = std::tuple<std::size_t, int, int>;

using Sheet = std::tuple<std::size_t, std::size_t, std::size_t,
                         char, char, char, char, char,
                         std::string, std::string, std::string>;

using ExtractedSheet = std::tuple<std::size_t, std::size_t, char,
                                  std::string, std::string>;


using Cycle = std::tuple<std::size_t, std::size_t, std::string>;

// B1, B2      : std::string
// Sheet       : std::string; "same", "not"
// Dir         : std::string; "-->", "?"
// PorA        : std::string; "para", "anti". ? if Sheet == "not"
// SSEs_LBTS   : std::string; one of following if not Sheet == "not"
//                b-c-b
//                b-a-b
//                b-b-b
//                b-b'-b
//                b-ab-b
//                b-ab'-b
//                b-bb'-b
//                b-abb'-b
// Jump        : std::size_t; ? if Sheet == "not"
// D1, D2      : int; ? if Sheet == "not" or jump != 0
// Bridge      : std::size_t; ? if Sheet == "not" jump != 0
// Score       : double; ? if Sheet == "not" PorA == "para"
// NumRes_LBTS : std::size_t; ? if Sheet == "not"
using SubStrandsPair = std::tuple<std::string, std::string,
                                  std::string, std::string, std::string,
                                  std::size_t, int, int, std::size_t,
                                  double, std::string, std::size_t>;

// resnum 1, 2
// PorA
// type
// ForB
using ResiduePair = std::tuple<int, int, std::string, std::string, std::string>;


// ********************************************************************
// Helper Templates for Formatting Tuples
// ********************************************************************

template <class Tuple, std::size_t N, std::size_t FIRST, std::size_t LAST>
struct TupleFormatter {
  static void add_format(boost::format & fmt, Tuple const& t) {
    if (FIRST+1 == N) {
      fmt % std::get<FIRST>(t);
    } else {
      TupleFormatter<Tuple, N-1, FIRST, LAST>::add_format(fmt, t);

      if (N <= LAST) {
        fmt % std::get<N-1>(t);
      }
    }
  }
};

template <class Tuple, std::size_t LAST>
struct TupleFormatter<Tuple, 1, 0, LAST> {
  static void add_format(boost::format & fmt, Tuple const& t) {
    fmt % std::get<0>(t);
  }
};

template <class Tuple, std::size_t FIRST=0, std::size_t LAST=std::tuple_size<Tuple>::value>
std::string format_tuple(std::string const& fmt_str, Tuple const& t) {
  boost::format fmt{fmt_str};
  // Turn the too_many_args exception off
  fmt.exceptions(boost::io::all_error_bits ^ boost::io::too_many_args_bit);

  TupleFormatter<decltype(t), std::tuple_size<Tuple>::value, FIRST, LAST>::add_format(fmt, t);
  return fmt.str();
}



// ********************************************************************
// Class FormatStrings
// ********************************************************************

struct FormatStrings {
  FormatStrings(std::string const& _header, std::string const& _data,
                std::size_t const _type=0, std::vector<std::string> const& _opts={}):
    header{_header}, data{_data}, type{_type}, opts{_opts} {}
  std::string const header{""};
  std::string const data{""};
  /// 0 for PDB like, 1 for mmcif_like
  std::size_t const type{0};
  /// Optional format strings
  std::vector<std::string> const opts{};
};




// ********************************************************************
// Class Template TableBase
// ********************************************************************
/// A base class for class template Table
struct TableBase {
  virtual ~TableBase();

  virtual std::string format_header(std::string const&) const;

  virtual std::string format_data(std::string const&) const;

  virtual std::string format(std::string const&,
                             std::string const&) const;

  virtual std::string format(FormatStrings const&) const;
};



// ********************************************************************
// Class Template Table
// ********************************************************************

template<typename Tuple>
class Table : public TableBase {
public:

  Table(std::string const& _name, std::vector<std::string> const& _col_names):
    name(_name), col_names{_col_names}
  {
    assert(std::tuple_size<Tuple>::value == col_names.size());
  }

  void add(Tuple const& d) {
    data.push_back(d);
  }

  std::string format_header(std::string const& fmt_str) const override {
    return _format_header(fmt_str, 0, std::tuple_size<Tuple>::value);
  }

  std::string format_data(std::string const& fmt_str) const override {
    return _format_data(fmt_str);
  }

  /// If this table has some data, return the formateed string. Otherwise, empty string.
  std::string format(std::string const& fmt_str_header,
                     std::string const& fmt_str_data) const override {
    return data.size() ? format_header(fmt_str_header) + format_data(fmt_str_data) : "";
  }


  std::string format(FormatStrings const& fmt) const override {
    return format(fmt.header, fmt.data);
  }


  std::string const name{""};
  std::vector<std::string> const col_names{};

protected:

  void add_header_format(boost::format & fmt, std::string const& str) const {
    fmt % str;
  }


  /// Acctually generates formatted header string
  std::string _format_header(std::string const& fmt_str,
                            std::size_t const first,
                            std::size_t const end) const {
    boost::format fmt{fmt_str};
    assert(static_cast<std::size_t>(fmt.expected_args()) <= col_names.size());

    for (std::size_t i = first; i < end; ++i) {
      add_header_format(fmt, col_names[i]);
    }
    return fmt.str();
  }


  /// Acctually generates formatted data string
  template<std::size_t FIRST=0, std::size_t LAST=std::tuple_size<Tuple>::value>
  std::string _format_data(std::string const& fmt_str) const {
#ifndef NDEBUG
    // just for assertion
    boost::format fmt{fmt_str};
    assert(static_cast<std::size_t>(fmt.expected_args()) <= col_names.size());
#endif // ifndef NDEBUG

    std::string str{""};
    for (auto const& t : data) {
      str += format_tuple<Tuple, FIRST, LAST>(fmt_str, t);
    }
    return str;
  }

  /// Stores the actual data
  std::vector<Tuple> data{};
};

/// Given one tuple from Table::data, return a suitable boost::format object
std::string get_substrands_pair_format(SubStrandsPair const& t,
                                       FormatStrings const& fmt_strs);

// Class template member function specialization declarations

template<>
std::string Table<Sheet>::format(FormatStrings const& fmt) const;

template <>
std::string Table<SubStrandsPair>::format(FormatStrings const& fmt) const;




// Typedefs

using TBLSubStrand = Table<SubStrand>;
using TBLHelix = Table<Helix>;
using TBLSheet = Table<Sheet>;
using TBLExtractedSheet = Table<ExtractedSheet>;
using TBLCycle = Table<Cycle>;
using TBLSubStrandsPair = Table<SubStrandsPair>;
using TBLResiduePair = Table<ResiduePair>;

using Set = std::tuple<TBLSubStrand, TBLHelix, TBLSheet,
                       TBLExtractedSheet, TBLCycle, TBLSubStrandsPair,
                       TBLResiduePair>;

constexpr std::size_t const SET_TUPLE_SIZE = std::tuple_size<Set>::value;
} // namespace table

#endif // ifndef TABLE_H_
