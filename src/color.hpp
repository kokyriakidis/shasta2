#pragma once

#include "array.hpp"
#include "cstdint.hpp"
#include "string.hpp"

namespace shasta {

    array<double, 2> hsvToHsl(double SV, double V);

    // Convert HSL values in [0.,1.] to RGB values also in [0.,1.].
    // https://en.wikipedia.org/wiki/HSL_and_HSV#HSL_to_RGB
    array<double, 3> hslToRgb(double H, double S, double L);

    // Generate a RGB color string #RRGGBB corresponding to a HSL
    // color with a random hue computed by hashing the given id
    // and the given SL values.
    string randomHslColor(uint64_t id, double S, double L);
}



