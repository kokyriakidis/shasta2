#pragma once

#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    namespace filesystem {

        // Return the extension of a path - that is, everything following
        // the last dot after the last slash.
        // If there is no dot after the last slash, throw an exception.
        string extension(const string&);

        // Find the absolute path.
        string getAbsolutePath(const string& path);

    }
}

