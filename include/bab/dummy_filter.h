// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef BAB_DUMMY_FILTER_H_
#define BAB_DUMMY_FILTER_H_

#include "bab/filter.h"

namespace bab {

struct FalseBab {
  bool operator()(SeqIter const, SeqIter const, unsigned const) {
    return false;
  }
};


} // namespace bab


#endif // #ifndef BAB_DUMMY_FILTER_H_
