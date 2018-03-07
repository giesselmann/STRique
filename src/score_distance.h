// \HEADER\---------------------------------------------------------------
//
//  CONTENTS      : Class score_distance
//
//  DESCRIPTION   : Distance score of two nanopore signal values
//
//  RESTRICTIONS  : none
//
//  REQUIRES      : none
//
// -----------------------------------------------------------------------
//  All rights reserved to Max Planck Institute for Molecular Genetics
//  Berlin, Germany
//  Written by Pay Giesselmann
// -----------------------------------------------------------------------
#ifndef SCORE_DISTANCE_H
#define SCORE_DISTANCE_H
// -- required headers ---------------------------------------------------
#include "seqan/score.h"

// -- forward declarations -----------------------------------------------

// -- exported constants, types, classes ---------------------------------
namespace seqan
{
    // Tag for distance score
    struct Distance_;
    typedef Tag<Distance_> Distance;

    /*!
    * @class DistanceScore
    * @extends Score
    * @headerfile <seqan/score.h>
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
        TValue data_gap_extend = -1;

        // The gap open score.
        TValue data_gap_open = -1;

        // Offset to split match/mismatch score
        TValue data_dist_offset = 1;
        
        // Clip min score to this value
        TValue data_dist_min = -1;

        // default constructor
        Score(){}

        Score(TValue _gap)
            : data_gap_extend(_gap), data_gap_open(_gap) {
        }

        Score(TValue _gap_extend, TValue _gap_open, TValue _dist_offset, TValue _dist_min)
            : data_gap_extend(_gap_extend), data_gap_open(_gap_open), data_dist_offset(_dist_offset), data_dist_min(_dist_min) {
        }
    };

    /*!
    * @typedef DistanceScoreTypedef DistanceScore
    * @headerfile <seqan/score.h>
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
        TValue score = me.data_dist_offset - static_cast<TValue>(valH > valV ? valH - valV : valV - valH);
        return score > me.data_dist_min ? score : me.data_dist_min;
    }

}  // namespace seqan

// -- exported functions - declarations ----------------------------------

// -- exported global variables - declarations (should be empty)----------
#endif  // SCORE_DISTANCE_H