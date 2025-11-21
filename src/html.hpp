#pragma once

#include "iosfwd.hpp"
#include "string.hpp"

// Miscellaneous html related functions.

namespace shasta2 {

    void writeHtmlBegin(ostream&, const string& title);
    void writeHtmlEnd(ostream&);
    void writeStyle(ostream&);

    void addSvgDragAndZoom(ostream& html);

    void writeInformationIcon(ostream& html, const string& message);
}

