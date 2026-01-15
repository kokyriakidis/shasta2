#pragma once

#include "fstream.hpp"
#include "string.hpp"

// The html log is used to write an assembly log
// in html format that is cleaner and more readable
// that stdout.log.

namespace shasta2 {
    extern ofstream htmlLog;
    void openHtmlLog(const string& fileName);
}

