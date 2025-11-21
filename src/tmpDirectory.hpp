#pragma once

#include "string.hpp"

namespace shasta2 {
    
    // Return the path to a usable temporary directory, including the final "/".
    inline string tmpDirectory()
    {
        return "/dev/shm/";
    }
    
}

