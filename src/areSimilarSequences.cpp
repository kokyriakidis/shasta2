// Shasta.
#include "areSimilarSequences.hpp"
#include "abpoaWrapper.hpp"
#include "Base.hpp"
#include "isPeriodic.hpp"
using namespace shasta2;

// Standard library.
#include "algorithm.hpp"
#include "iterator.hpp"

// See areSimilarSequences.hpp for comments.



bool shasta2::areSimilarSequences(
        const vector<Base>& x,
        const vector<Base>& y,
        const vector<uint64_t>& minRepeatCount,
        ostream& html)
{
    // Prepare arguments for abpoa.
    const vector< vector<Base> > sequences = {x, y};
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;

    // Use abpoa to compute the alignment.
    abpoa(sequences, consensus, alignment, alignedConsensus, true);

    // Sanity checks.
    SHASTA2_ASSERT(alignment.size() == 2);
    const vector<AlignedBase>& X = alignment[0];
    const vector<AlignedBase>& Y = alignment[1];
    const uint64_t alignmentLength = X.size();
    SHASTA2_ASSERT(Y.size() == alignmentLength);



    // Debug output.
    if(html) {
        html <<
            "<h3>Alignment</h3><div style='font-family:monospace;white-space:nowrap;'>";

        // Position scale labels.
        for(size_t position=0; position<alignmentLength; ) {
            if((position%10)==0) {
                const string label = to_string(position);
                html << label;
                for(size_t i=0; i<10-label.size(); i++) {
                    html << "&nbsp;";
                }
                position += 10;
            } else {
                html << "&nbsp;";
                ++position;
            }
        }
        html<< "<br>";

        // Position scale
        for(size_t position=0; position<alignmentLength; position++) {
            if((position%10)==0) {
                html << "|";
            } else if((position%5)==0) {
                html << "+";
            } else {
                html << ".";
            }
        }

        // Alignment lines.
        for(uint64_t i=0; i<2; i++) {
            html << "<br>";
            for(uint64_t position=0; position<alignmentLength; position++) {
                const AlignedBase b = alignment[i][position];
                const AlignedBase bOther = alignment[1-i][position];
                if(b != bOther) {
                    html << "<span style='background-color:pink'>";
                }
                html << b;
                if(b != bOther) {
                    html << "</span>";
                }
            }
        }
        html << "</div><p>";
    }



    // If the alignment contains any mismatches, return false.
    for(uint64_t j=0; j<alignmentLength; j++) {
        const AlignedBase b0 = X[j];
        const AlignedBase b1 = Y[j];
        if((not b0.isGap()) and (not b1.isGap()) and (b0 != b1)) {
            if(html) {
                html << "<br>Mismatches are present, sequences are not similar." << endl;
            }
            return false;
        }
    }



    // If getting here, we did not find any mismatches.
    // Check for indels.
    // At j=0 we check for deletions in X.
    // At j=1 we check for deletions in Y.
    for(uint64_t j=0; j<2; j++) {
        const vector<AlignedBase>& A = alignment[j];
        const vector<AlignedBase>& B = alignment[1 - j];

        // Loop over streaks of gaps in A.
        for(uint64_t streakBegin=0; streakBegin<alignmentLength;) {
            if(not A[streakBegin].isGap()) {
                ++streakBegin;
                continue;
            }
            uint64_t streakEnd = streakBegin + 1;
            for(; streakEnd<alignmentLength; ++streakEnd) {
                if(not A[streakEnd].isGap()) {
                    break;
                }
            }
            const uint64_t streakLength = streakEnd - streakBegin;
            if(html) {
                html << "<p>Found a " << streakLength << " base deletion in the " << (j == 0 ? "first" : "second") <<
                    " sequence at alignment positions " << streakBegin << " " << streakEnd << "."
                    "<br>The deleted sequence is ";
                for(uint64_t position=streakBegin; position!=streakEnd; position++) {
                    html << B[position];
                }
                html << ".";
            }



            // Check all allowed periods.
            for(uint64_t period=1; period<minRepeatCount.size(); period++) {
                if(period > streakLength) {
                    break;
                }
                if(html) {
                    html << "<br>Checking period " << period;
                }
                if(not isPeriodic(B.begin()+streakBegin, B.begin()+streakEnd, period)) {
                    if(html) {
                        html << "<br>The deleted sequence does not have period " << period << endl;
                    }
                    continue;   // Continue the loop over periods.
                }

                const uint64_t periodicSequenceBegin = streakBegin;
                // const uint64_t periodicSequenceEnd = periodicSequenceBegin + period;


                // See how many times this repeats on the right.
                uint64_t copyNumberOnRight = 1;
                for(; ; ++copyNumberOnRight) {
                    bool copyIsIntact = true;
                    for(uint64_t i=0; i<period; i++) {
                        const uint64_t shiftedPosition = periodicSequenceBegin + i + period * copyNumberOnRight;
                        if(shiftedPosition >= alignmentLength) {
                            copyIsIntact = false;
                            break;
                        }
                        if(B[shiftedPosition] != B[periodicSequenceBegin + i]) {
                            copyIsIntact = false;
                            break;
                        }
                    }
                    if(not copyIsIntact) {
                        break;
                    }
                }
                --copyNumberOnRight;
                if(html) {
                    // html << "<br>Found " << copyNumberOnRight << " copies on the right of the deletion." << endl;
                }

                // See how many times this repeats on the left.
                uint64_t copyNumberOnLeft = 1;
                for(; ; ++copyNumberOnLeft) {
                    bool copyIsIntact = true;
                    if(period > streakBegin) {
                        copyIsIntact = false;
                    } else {
                        for(uint64_t i=0; i<period; i++) {
                            const uint64_t shiftedPosition = periodicSequenceBegin + i - period * copyNumberOnLeft;
                            if(B[shiftedPosition] != B[periodicSequenceBegin + i]) {
                                copyIsIntact = false;
                                break;
                            }
                        }
                    }
                    if(not copyIsIntact) {
                        break;
                    }
                }
                --copyNumberOnLeft;
                if(html) {
                    // html << "<br>Found " << copyNumberOnLeft << " copies on the left of the deletion." << endl;
                }

                const uint64_t totalCopyNumber = copyNumberOnRight + copyNumberOnLeft + 1;
                if(html) {
                    html << "<br>Found a total " << totalCopyNumber << " copies of this repeat of period " << period <<
                        " (" << totalCopyNumber * period << " bases).";
                }

                if(totalCopyNumber < minRepeatCount[period]) {
                    if(html) {
                        html << "<br>This repeat is short. This difference is significant. These sequences are not similar.";
                    }
                    return false;
                } else {
                    // Don't check other periods.
                    break;
                }
            }

            // Prepare for the next loop iteration.
            streakBegin = streakEnd;
        }

    }

    // If getting here, we did not encounter any significant differences,
    // as defined by the minRepeatCount vector. The input sequences are "similar".
    if(html) {
        html << "<br>No significant differences were found. These sequences are similar.";
    }
    return true;
}
