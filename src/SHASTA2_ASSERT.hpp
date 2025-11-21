// Definition of macro SHASTA2_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#pragma once

namespace shasta2 {
    void handleFailedAssertion(
        const char* expression,
        const char* function,
        const char* file,
        int line
    ) __attribute__ ((__noreturn__));
}


#define SHASTA2_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (shasta2::handleFailedAssertion(#expression, __PRETTY_FUNCTION__,  __FILE__ , __LINE__)))



