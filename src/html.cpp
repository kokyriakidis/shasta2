#include "html.hpp"
using namespace shasta2;

#include "iostream.hpp"



void shasta2::writeHtmlBegin(ostream& html, const string& title)
{
    html <<
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<meta charset='UTF-8'>"
        "<title>" << title << "</title>";
    writeStyle(html);
    html << "</head>";
}



void shasta2::writeHtmlEnd(ostream& html)
{
    html << "</html>";
}



void shasta2::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 1px solid #b8b5c7d9;
        padding: 2px;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.right {
        text-align: right;
    }
    a {
        color: DarkSlateBlue;
    }

    /* This can be used to get vertical text in table cells. */
    span.rotated 
    {
      writing-mode: vertical-rl;
      transform: rotate(180deg);
    }
</style>
    )%";

}



void shasta2::addSvgDragAndZoom(ostream& html)
{
    html << R"zzz(  
<script>

var svg = document.querySelector('svg');
svg.scrollIntoView();
svg.addEventListener('pointerdown', onPointerDown); 
svg.addEventListener('pointerup', onPointerUp); 
svg.addEventListener('pointerleave', onPointerLeave); 
svg.addEventListener('pointermove', onPointerMove); 
svg.addEventListener('wheel', onMouseWheel); 

var pointerIsDown = false;

var xOrigin = 0;
var yOrigin = 0;

// The current viewbox.
var x = svg.viewBox.baseVal.x;
var y = svg.viewBox.baseVal.y;
var width = svg.viewBox.baseVal.width;
var height = svg.viewBox.baseVal.height;

var xNew = 0;
var yNew = 0;

var ratio = width / svg.getBoundingClientRect().width;

function onPointerDown(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    pointerIsDown = true;
    xOrigin = event.clientX;
    yOrigin = event.clientY;

    return false;
}

function onPointerMove (event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    if (!pointerIsDown) {
        return;
    }

    xNew = x - (event.clientX - xOrigin) * ratio;
    yNew = y - (event.clientY - yOrigin) * ratio;

    svg.setAttribute('viewBox', `${xNew} ${yNew} ${width} ${height}`);

    return false;
}

function onPointerUp(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    pointerIsDown = false;
    x = xNew;
    y = yNew;

    return false;
}

function onPointerLeave(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();
    if(pointerIsDown) {
        pointerIsDown = false;
        x = xNew;
        y = yNew;
    }

    return false;
}

function onMouseWheel(event) {
    // If the CTRL key is not pressed, don't do anything.
    if(!event.ctrlKey) {
        return;
    }
    event.stopPropagation();
    event.preventDefault();  
    var value = event.wheelDelta / 120.;
    var factor = Math.pow(1.1, value);
    zoomSvg(factor);
}

function zoomSvg(factor)
{
    ratio /= factor;

    // Adjust the viewbox so the center does not move.
    var xCenter = x + 0.5 * width;
    var yCenter = y + 0.5 * height;  
    width /= factor;
    height /= factor;
    x = xCenter - 0.5 * width;
    y = yCenter - 0.5 * height;

    svg.setAttribute('viewBox', `${x} ${y} ${width} ${height}`);
    svg.setAttribute('font-size', svg.getAttribute('font-size') / factor);

    return false;
}

</script>
    )zzz";
}




void shasta2::writeInformationIcon(ostream& html, const string& message)
{
    html << "<span style='color:Blue;font-weight:bold' title=\"" <<
            message << "\">&#9432;</span>";
}



void shasta2::writeMakeAllTablesCopyable(ostream& html)
{
    html << R"###(
    <script>

    // Copy to the clipboard the table that generated the event.
    function copyToClipboard(event)
    {
        // If the CTRL key is not pressed, don't do anything.
        if(!event.ctrlKey) {
            return;
        }

        // Prevent default behavior.
        // event.preventDefault();
        // event.stopPropagation();
        // event.returnValue = false;
        
        // Get the table element.
        var element = event.currentTarget;
         
        // Remove any previous selection.
        var selection = window.getSelection();
        selection.removeAllRanges();
        
        // Select the table.
        var range = document.createRange();
        range.selectNodeContents(element);
        selection.addRange(range);
        
        // Copy it to the clipboard.
        document.execCommand("copy");

        // Unselect it.
        selection.removeAllRanges();

        window.alert("The table was copied to the clipboard");
    }

    // Make a table copyable by Ctrl-click.
    function makeCopyable(element)
    {
        element.addEventListener('click', copyToClipboard);
        element.title = 'Ctrl-click anywhere on the table to copy the entire table to the clipboard';
    }

    // Make all tables copyable by Ctrl-click.
    function makeAllTablesCopyable()
    {
        var tables = document.getElementsByTagName('table');
        var i;
        for(i=0; i<tables.length; i++) {
            makeCopyable(tables[i]);
        }
    }
    </script>
    )###";


#if 0
    html << R"###(
<script>

// Make all tables selectable by double click.
// This must be called after all tables have
// already been created, so it can be called during onload.

// This function is called when the user double clicks on a table.
function selectElement(table)
{
    var selection = window.getSelection();
    selection.removeAllRanges();
    var range = document.createRange();
    range.selectNode(table);
    selection.addRange(range);
}

// Attach the above function to the double click event
// for all tables in the document.
// Also add to each table a title that displays a tooltip 
// explaining that the table can be selected via double click.
function makeAllTablesSelectableByDoubleClick()
{
    var allTables = document.getElementsByTagName("table");
    for (var i=0; i<allTables.length; i++) {
        var table = allTables[i];
        table.ondblclick = function() {selectElement(this);};
        table.setAttribute("title", 
        "Double click to select the entire table. You can then paste it into a spreadsheet.");
    }
}
</script>
    )###";
#endif
}
