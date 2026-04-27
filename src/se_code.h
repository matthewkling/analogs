#pragma once

namespace analogs {

// SE variant code (R->C++ mapping for `se` parameter)
// 0 = none, 1 = ess, 2 = design
enum class SeCode : int {
      NONE   = 0,
            ESS    = 1,
            DESIGN = 2
};

} // namespace analogs
