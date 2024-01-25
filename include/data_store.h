// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef DATA_STORE_H_
#define DATA_STORE_H_

#include <array>
#include <iostream>
#include <string>
#include <tuple>

#include "parse_argument.h"
#include "table.h"

namespace data_store {

// *********************************
// Templates
// *********************************

template<std::size_t I = 0, class Tuple=table::Set>
typename std::enable_if<I == std::tuple_size<Tuple>::value, std::string>::type
print_tables(Tuple const&,
             std::array<table::FormatStrings, std::tuple_size<Tuple>::value> const&) {
  return "";
}

template<std::size_t I = 0, class Tuple=table::Set>
typename std::enable_if<I < std::tuple_size<Tuple>::value, std::string>::type
print_tables(Tuple const& t,
             std::array<table::FormatStrings, std::tuple_size<Tuple>::value> const& f) {
  return std::get<I>(t).format(f[I]) + "\n" + print_tables<I + 1, Tuple>(t, f);
}


// *****************************************************************************
// Class Data
// *****************************************************************************

template <class Tuple>
class Data {
public:
  // *********************************
  // Public Member Functions
  // *********************************

  Data(Tuple const& t) : tables{t} {
    static_assert(std::tuple_size<Tuple>::value == table::SET_TUPLE_SIZE,
                  "Tuple size is invalid.");
    static_assert(std::tuple_size<typename std::decay<decltype(t)>::type>::value ==
                  table::SET_TUPLE_SIZE,
                  "size of t is invalid.");
    static_assert(std::tuple_size<typename std::decay<decltype(tables)>::type>::value ==
                  table::SET_TUPLE_SIZE,
                  "size of tables is invlaid.");
  }

  void format_out(std::ostream & os, std::size_t const type) const {
    if (formats.size() <= type) {
      std::cerr << "Fatal error: Unknown format type '" << type << "'." << std::endl;
      throw arg::argument_error{"unknown type"};
    }
    static_assert(std::tuple_size<decltype(tables)>::value == table::SET_TUPLE_SIZE,
                  "size of tables is invalid");
    os << print_tables(tables, formats[type]);
  }

  template <class T>
  auto& table() {
    return std::get<table::Table<T>>(tables);
  }

protected:

  // *********************************
  // Protected Member Variables
  // *********************************

  Tuple tables{};

  std::array<std::array<table::FormatStrings, std::tuple_size<Tuple>::value>, 2> const formats{{
    {{
      table::FormatStrings{"REMARK %|18t|%s  %s   %s   %s\n",
                         "SUBSTRAND %|18t|%12d  %8d  %4d  %4d\n"},

      table::FormatStrings{"REMARK %|34t|%s   %s   %s\n",
                           "HELIX  %|32t|%8d  %4d  %4d\n"},

      table::FormatStrings{"REMARK %|18t|%s  %s  %s  %s  %s  %s  %s  %s\n",
                           "SHEET_INFO %|18t|%8d  %9d  %5d  %10d  %11c  %11c  %8c  %8c"
                           "\n", 0,
                           {"MEMBER %|18t|%1$5d  %9$s\n"
                            "NOMENCLATURE_R    %1$5d  %10$s\n"
                            "NOMENCLATURE_C    %1$5d  %11$s\n"}},

       table::FormatStrings{"REMARK    %|18t|%s  %s  %s  %s  %|70t|%s\n",
                            "EXT_SHEET %|18t|%8d  %9d  %1s                 %s"
                            "  %|70t|%s\n"},

       table::FormatStrings{"REMARK %|18t|%s  %s  %s\n",
                            "CYCLE %|18t|%8d  %9d  %s\n"},

       table::FormatStrings{"REMARK %|22t|%5s %5s %17s %3s %4s "
                            "%3s %3s %3s %3s %4s %8s %5s\n",
                            "", 0,

                            {
                            // Default
                             "STRAND_PAIR %|22t|%5s %5s %17s %3s %4s "
                             "%4d %3d %3d %6d  %4.2f %9s %11d\n",

                            // score < 0.0
                             "STRAND_PAIR %|22t|%1$5s %2$5s %3$17s %4$3s %5$4s "
                             "%6$4d %7$3d %8$3d %9$6d %|85t|? %11$9s %12$11d\n",

                            // Jump != 0
                             "STRAND_PAIR %|22t|%1$5s %2$5s %3$17s %4$3s %5$4s "
                             "%6$4d %|68t|? %|72t|? %|79t|?  %10$4.2f %11$9s %12$11d\n",

                            // score < 0.0 and Jump != 0
                             "STRAND_PAIR %|22t|%1$5s %2$5s %3$17s %4$3s %5$4s "
                             "%6$4d %|68t|? %|72t|? %|79t|? %|85t|? %11$9s %12$11d\n",

                            // Sheet != "same"
                             "STRAND_PAIR %|22t|%1$5s %2$5s %3$17s %|54t|? %|59t|?"
                             "%|64t|? %|68t|? %|72t|? %|79t|? %|85t|? %|95t|? %|107t|?\n"
                            }},

       table::FormatStrings{"REMARK %|18t|%s  %s  %s  %12s  %5s\n",
                            "RESIDUE_PAIR %|18t|%7d  %7d  %s  %12s  %5s\n"}
     }},
     {{
       table::FormatStrings{"#\nloop_\n_substrand.%s\n_substrand.%s\n_substrand.%s\n"
                            "_substrand.%s\n",
                            "%4s  %4d  %4d  %4d\n", 1},
       table::FormatStrings{"#\nloop_\n_helix.%s\n_helix.%s\n_helix.%s\n",
                            "%8d  %4d  %4d\n", 1},
       table::FormatStrings{"#\nloop_\n_sheet.%s\n_sheet.%s\n_sheet.%s\n_sheet.%s\n"
                            "_sheet.%s\n_sheet.%s\n_sheet.%s\n_sheet.%s\n_sheet.%s\n"
                            "_sheet.%s\n_sheet.%s\n",
                            "%3d  %3d  %3d %1c %1c %1c %1c %1c %s %|40t| %s %|65t| %s"
                            "\n", 1},

       table::FormatStrings{"#\nloop_\n_extracted_sheet.%s\n_extracted_sheet.%s\n"
                            "_extracted_sheet.%s\n_extracted_sheet.%s\n"
                            "_extracted_sheet.%s\n",
                            "%3d %3d %1s %s %|40t|%s\n", 1},

       table::FormatStrings{"#\nloop_\n_cycle.%s\n_cycle.%s\n_cycle.%s\n",
                            "%3d %3d %s\n", 1},

       table::FormatStrings{"", "", 1},

       table::FormatStrings{"", "", 1}
      }}
  }};

};

}

#endif // #ifndef DATA_STORE_H_

