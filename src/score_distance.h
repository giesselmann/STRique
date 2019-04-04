// \HEADER\-------------------------------------------------------------------------
//
//  CONTENTS      : Class score_distance
//
//  DESCRIPTION   : Distance score of two nanopore signal values
//
//  RESTRICTIONS  : none
//
//  REQUIRES      : none
//
// ---------------------------------------------------------------------------------
// Copyright (c) 2018-2019,  Pay Giesselmann, Max Planck Institute for Molecular Genetics
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Written by Pay Giesselmann
// ---------------------------------------------------------------------------------
#ifndef SCORE_DISTANCE_H
#define SCORE_DISTANCE_H
// -- required headers -------------------------------------------------------------
#include "seqan/score.h"
#include <math.h>

// -- forward declarations ---------------------------------------------------------

// -- exported constants, types, classes -------------------------------------------
namespace seqan
{
    // Tag for distance score
    struct Distance_;
    typedef Tag<Distance_> Distance;

    /*!
    * @class Score
    * @brief Signal distance scoring scheme, opening gaps and extending gaps.
    *
    * @signature template <typename TValue>
    *            class Score<TValue, Distance>;
    *
    * @tparam TValue The score value to use.
    */
    template <typename TValue>
    class Score<TValue, Distance> {
    public:
        // The gap extension score.
        TValue data_gap_extend_horizontal = -1;
        TValue data_gap_extend_vertical = -1;

        // The gap open score.
        TValue data_gap_open_horizontal = -1;
        TValue data_gap_open_vertical = -1;

        // Offset to split match/mismatch score
        TValue data_dist_offset = 1;
        
        // Clip min score to this value
        TValue data_dist_min = -1;

        // default constructor
        Score(){}

        Score(TValue _gap)
            : data_gap_extend_horizontal(_gap), data_gap_extend_vertical(_gap),
                data_gap_open_horizontal(_gap), data_gap_open_vertical(_gap) {
        }

        Score(TValue _gap_extend, TValue _gap_open, TValue _dist_offset, TValue _dist_min)
            : data_gap_extend_horizontal(_gap_extend), data_gap_extend_vertical(_gap_extend), data_gap_open_horizontal(_gap_open), data_gap_open_vertical(_gap_open),
            data_dist_offset(_dist_offset), data_dist_min(_dist_min) {
        }
        
        Score(TValue _gap_extend_h, TValue _gap_extend_v, TValue _gap_open_h, TValue _gap_open_v, TValue _dist_offset, TValue _dist_min)
            : data_gap_extend_horizontal(_gap_extend_h), data_gap_extend_vertical(_gap_extend_v), data_gap_open_horizontal(_gap_open_h), data_gap_open_vertical(_gap_open_v),
            data_dist_offset(_dist_offset), data_dist_min(_dist_min) {
        }        
    };

    /*!
    * @typedef DistanceScoreTypedef DistanceScore
    * @brief Absolute distance scoring scheme.
    *
    * @signature typedef Score<int, Distance> DistanceScore;
    */
    typedef Score<int, Distance> DistanceScore;

    /*!
    * @fn Score#score
    * @brief Returns distance score for two sequence entries.
    *
    * @signature TValue score(score, entryH, entryV);
    *
    * @param[in] score  The Score to use for comparing the two sequence entries.
    * @param[in] entryH The entry in the first/horizontal sequence.
    * @param[in] entryV The entry in the second/vertical sequence.
    *
    * @return TValue The score similarity cost for gaps at the given position/entry.  TValue is the value type of
    *                score.
    */
    template <typename TValue, typename TSeqHVal, typename TSeqVVal>
    inline TValue
    score(Score<TValue, Distance> const & me, TSeqHVal valH, TSeqVVal valV)
    {
        //TValue score = me.data_dist_offset - static_cast<TValue>(valH > valV ? valH - valV : valV - valH);
        TValue score = me.data_dist_offset - static_cast<TValue>(pow(valH > valV ? valH - valV : valV - valH, 1.2));
        return score > me.data_dist_min ? score : me.data_dist_min;
    }
    
    /*!
     * @fn Score#scoreGapOpenHorizontal
     * @brief Returns the score for opening a gap in horizontal direction.
     *
     * @signature TValue scoreGapOpenHorizontal(score, entryH, entryV);
     *
     * @param[in] score  The Score to query.
     * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
     * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
     *
     * @return TValue The score open cost for gaps at the given position/entry.  TValue is the value type of score.
     *
     * @section Remarks
     *
     * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
     */
    template <typename TValue, typename TSeqHValue, typename TSeqVValue>
    inline TValue
    scoreGapOpenHorizontal(
        Score<TValue, Distance> const & me,
        TSeqHValue const & /*seqHVal*/,
        TSeqVValue const & /*seqHVal*/)
    {
        return me.data_gap_open_horizontal;
    }

    /*!
     * @fn Score#scoreGapOpenVertical
     * @brief Returns the score for opening a gap in vertical direction.
     *
     * @signature TValue scoreGapOpenVertical(score, entryH, entryV);
     *
     * @param[in] score  The Score to query.
     * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
     * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
     *
     * @return TValue The score open cost for gaps at the given position/entry.  TValue is the value type of score.
     *
     * @section Remarks
     *
     * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
     */
    template <typename TValue, typename TSeqHValue, typename TSeqVValue>
    inline TValue
    scoreGapOpenVertical(
        Score<TValue, Distance> const & me,
        TSeqHValue const & /*seqHVal*/,
        TSeqVValue const & /*seqHVal*/)
    {
        return me.data_gap_open_vertical;
    }

    /*!
     * @fn Score#scoreGapExtendHorizontal
     * @brief Returns the score for extending a gap in horizontal direction.
     *
     * @signature TValue scoreGapExtendHorizontal(score, entryH, entryV);
     *
     * @param[in] score  The Score to query.
     * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
     * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
     *
     * @return TValue The score extension cost for gaps at the given position/entry.  TValue is the value type of score.
     *
     * @section Remarks
     *
     * Corresponds to a deletion event in sequence two and an insertion event in sequence one, respectively.
     */
    template <typename TValue, typename TSeqHValue, typename TSeqVValue>
    inline TValue
    scoreGapExtendHorizontal(
        Score<TValue, Distance> const & me,
        TSeqHValue const & /*seqHVal*/,
        TSeqVValue const & /*seqHVal*/)
    {
        return me.data_gap_extend_horizontal;
    }

    /*!
     * @fn Score#scoreGapExtendVertical
     * @brief Returns the score for extending a gap in vertical direction.
     *
     * @signature TValue scoreGapExtendVertical(score, entryH, entryV);
     *
     * @param[in] score  The Score to query.
     * @param[in] entryH Entry in sequence one (horizontal), type from Score#SequenceEntryForScore.
     * @param[in] entryV Entry in sequence two (vertical), type from Score#SequenceEntryForScore.
     *
     * @return TValue The score extension cost for gaps at the given position/entry.  TValue is the value type of score.
     *
     * @section Remarks
     *
     * Corresponds to a deletion event in sequence one and an insertion event in sequence two, respectively.
     */
    template <typename TValue, typename TSeqHValue, typename TSeqVValue>
    inline TValue
    scoreGapExtendVertical(
        Score<TValue, Distance> const & me,
        TSeqHValue const & /*seqHVal*/,
        TSeqVValue const & /*seqHVal*/)
    {
        return me.data_gap_extend_vertical;
    }

}  // namespace seqan

// -- exported functions - declarations --------------------------------------------

// -- exported global variables - declarations (should be empty)--------------------
#endif  // SCORE_DISTANCE_H