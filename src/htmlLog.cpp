// The html log is used to write an assembly log
// in html format that is cleaner and more readable
// that stdout.log.

#include "htmlLog.hpp"

namespace shasta2 {
    ofstream htmlLog;
}



void shasta2::openHtmlLog(const string& fileName)
{
    htmlLog.open(fileName);
}
