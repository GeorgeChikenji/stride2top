// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <string>
#include "parse_argument.h"
#include "table.h"

namespace table {
TableBase::~TableBase() {}
std::string TableBase::format_header(std::string const&) const {return "";}
std::string TableBase::format_data(std::string const&) const {return "";}
std::string TableBase::format(std::string const&,
                              std::string const&) const {return "";}
std::string TableBase::format(FormatStrings const&) const {return "";}

// *************************************************************
// Function get_substrands_pair_format()
// *************************************************************

std::string get_substrands_pair_format(SubStrandsPair const& t,
                                       FormatStrings const& fmt_strs) {
  if (std::get<2>(t).substr(0, 4) != "same") {
    return fmt_strs.opts.at(4);
  }

  std::size_t flg = 0u;
  if (std::get<9>(t) < 0.0) {
    flg |= 1u;
  }
  if (std::get<5>(t) != 0) {
    flg |= 2u;
  }
  return fmt_strs.opts.at(flg);
}


// *************************************************************
// Specializations for Template class Table
// *************************************************************

template<>
std::string Table<Sheet>::format(FormatStrings const& fmt) const {
  if (fmt.type == 0) {
    assert(fmt.opts.size() == 1);
    return _format_header(fmt.header, 0, 8) +
           _format_data<0, 8>(fmt.data) + "\n" +
           "REMARK            Sheet  Description\n" +
           _format_data(fmt.opts[0]);
  } else if (fmt.type == 1) {
    return format_header(fmt.header) + format_data(fmt.data);
  } else {
    throw arg::argument_error{"unknown type '" + std::to_string(fmt.type) + "'"};
  }
}


template <>
std::string Table<SubStrandsPair>::format(FormatStrings const& fmt) const {
  if (fmt.type == 0) {
    assert(fmt.opts.size() == 5);

    auto str = format_header(fmt.header);
    for (auto const& t : data) {
      auto const fmt_str = get_substrands_pair_format(t, fmt);
      str += format_tuple<SubStrandsPair, 0,
                          std::tuple_size<SubStrandsPair>::value>(fmt_str, t);
    }
    return str;
  } else if (fmt.type == 1) {
    return "";
  } else {
    throw arg::argument_error{"unknown type '" + std::to_string(fmt.type) + "'"};
  }
}


} // namespace table
