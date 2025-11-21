#include "color.hpp"
using namespace shasta;

#include "MurmurHash2.hpp"
#include "SHASTA2_ASSERT.hpp"

#include <algorithm.hpp>
#include <array.hpp>
#include <cmath>
#include <format>
// #include <iostream.hpp>
#include <limits>
#include "string.hpp"
#include "utility.hpp"



// HSV to HSL conversion.
// Given S and V for an HSV color, returns the S and L
// values for the corresponding HSL color.
// See https://en.wikipedia.org/wiki/HSL_and_HSV#HSV_to_HSL
std::array<double, 2> shasta::hsvToHsl(double SV, double V)
{
    const double L = V * (1. - 0.5 * SV);
    const double denominator = min(L, 1. - L);
    const double SL = (denominator == 0.) ? 0. : (V - L)/denominator;
    return {SL, L};
}



// Convert HSL values in [0.,1.] to RGB values also in [0.,1.].
// https://en.wikipedia.org/wiki/HSL_and_HSV#HSL_to_RGB
std::array<double, 3> shasta::hslToRgb(double H, double S, double L)
{
    using std::fabs;
    using std::fmod;

    const double C = (1. - fabs(2. * L - 1.)) * S;
    const double Hprime = 6 * H;
    const double X = C * (1 - fabs(fmod(Hprime, 2.) - 1.));

    /*
    cout << "HSL " << H << " " << S << " " << L << endl;
    cout << "C " << C << endl;
    cout << "H' " << Hprime << endl;
    cout << "X " << X << endl;
    */

    array<double, 3> rgb;
    if(Hprime >= 0. and Hprime <1.) {
        rgb = {C, X, 0.};
    } else
    if(Hprime >= 1. and Hprime <2.) {
        rgb = {X, C, 0.};
    } else if(Hprime >= 2. and Hprime <3.) {
        rgb = {0., C, X};
    } else if(Hprime >= 3. and Hprime <4.) {
        rgb = {0., X, C};
    } else if(Hprime >= 4. and Hprime <5.) {
        rgb = {X, 0., C};
    } else if(Hprime >= 5. and Hprime <6.) {
        rgb = {C, 0., X};
    } else {
        SHASTA2_ASSERT(0);
    }

    // cout << "RGB1 " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;

    const double m = L - C / 2.;
    // cout << "m " << m << endl;
    for(uint64_t i=0; i<3; i++) {
        rgb[i] += m;
    }
    // cout << "RGB " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;

    return rgb;
}



string shasta::hslToRgbString(double H, double S, double L)
{
    using std::format;
    const array<double, 3> rgb = hslToRgb(H, S, L);

    // Normalize to 255 and construct the RGB string.
    string s = "#";
    for(const double value: rgb) {
        const int iValue = min(255, int(value * 255.));
        s.append(format("{:02x}", iValue));
    }

    return s;
}



// Generate a RGB color string #RRGGBB corresponding to a HSL
// color with a random hue computed by hashing the given id
// and the given SL values.
string shasta::randomHslColor(uint64_t id, double S, double L)
{
    using std::format;

    // Compute the random hue by hashing the id.
    const double H = double(MurmurHash2(&id, sizeof(id), 759)) / double(std::numeric_limits<uint32_t>::max());

    // Find the RGB values in [0.,1.].
    const array<double, 3> rgb = hslToRgb(H, S, L);


    // Normalize to 255 and construct the RGB string.
    string s = "#";
    for(const double value: rgb) {
        const int iValue = min(255, int(value * 255.));
        s.append(format("{:02x}", iValue));
    }

    return s;
}

