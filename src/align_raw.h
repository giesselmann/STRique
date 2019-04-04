// \HEADER\-------------------------------------------------------------------------
//
//  CONTENTS      : Class align_raw
//
//  DESCRIPTION   : Align raw data to reference sequence
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
#ifndef ALIGN_RAW_H
#define ALIGN_RAW_H
// -- required headers -------------------------------------------------------------
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <iostream>
#include "score_distance.h"
#include "seqan/align.h"
#include "seqan/graph_msa.h"

using namespace seqan;

// -- forward declarations ---------------------------------------------------------

// -- exported constants, types, classes -------------------------------------------
template<typename TScore>
struct align_raw_settings
{
    TScore m_gap_open_h = static_cast<TScore>(-2.0L);
    TScore m_gap_open_v = static_cast<TScore>(-2.0L);
    TScore m_gap_extension_h = static_cast<TScore>(-8.0L);
    TScore m_gap_extension_v = static_cast<TScore>(-8.0L);
    TScore m_dist_offset = static_cast<TScore>(8.0L);
    TScore m_dist_min = static_cast<TScore>(-16.0);
};


template<typename TScore, typename TVal>
class align_raw
{
public:
    // default constructor
    align_raw() 
    {
        align_raw_settings<TScore> default_settings{};
        init(default_settings);
    };

    // construct with settings
    align_raw(align_raw_settings<TScore> const& settings)
    {
        init(settings);
    }

    // virtual destructor
    virtual ~align_raw() {};

    // getter and setter
    void set_gap_open(const TScore gap_open) { m_gap_open_h = gap_open; m_gap_open_v = gap_open; }
    const TScore get_gap_open(void) const { return m_gap_open_h; }
    void set_gap_extension(const TScore gap_extension) { m_gap_extension_h = gap_extension; m_gap_extension_v = gap_extension; }
    const TScore get_gap_extension(void) const { return m_gap_extension_h; }
    
    void set_gap_open_h(const TScore gap_open) { m_gap_open_h = gap_open; }
    const TScore get_gap_open_h(void) const { return m_gap_open_h; }
    void set_gap_extension_h(const TScore gap_extension) { m_gap_extension_h = gap_extension; }
    const TScore get_gap_extension_h(void) const { return m_gap_extension_h; }
    
    void set_gap_open_v(const TScore gap_open) { m_gap_open_v = gap_open; }
    const TScore get_gap_open_v(void) const { return m_gap_open_v; }
    void set_gap_extension_v(const TScore gap_extension) { m_gap_extension_v = gap_extension; }
    const TScore get_gap_extension_v(void) const { return m_gap_extension_v; }

    void set_dist_offset(const TScore dist_offset) { m_dist_offset = dist_offset; }
    const TScore get_dist_offset(void) const { return m_dist_offset; }

    void set_dist_min(const TScore dist_min) { m_dist_min = dist_min; }
    const TScore get_dist_min(void) const { return m_dist_min; }

    // semi global alignment, a >> b
    auto
    semiglobal
    (
        const std::vector<TVal>& a,
        const std::vector<TVal>& b
    )
    {
        return align<true, false>(a, b);
    }
    
    // sequence to reference alignment
    template<bool a_gap, bool b_gap>
    inline
    std::tuple<TScore, std::vector<uint64_t>, std::vector<uint64_t>>
    align
    (
        const std::vector<TVal>& a,
        const std::vector<TVal>& b,
        bool returnViewPosition = true
    )
    {
        typedef std::vector<TVal> TSequence;
        typedef Align<TSequence, ArrayGaps> TAlign;
        typedef typename Row<TAlign>::Type TRow;
        TAlign algn;
        resize(rows(algn), 2);
        assignSource(row(algn, 0), a);
        assignSource(row(algn, 1), b);
        auto score = globalAlignment(algn, Score<TScore, Distance>(m_gap_extension_h, m_gap_extension_v, m_gap_open_h, m_gap_open_v, m_dist_offset, m_dist_min),
            AlignConfig<a_gap, b_gap, b_gap, a_gap>(), AffineGaps());
        TRow & row1 = row(algn, 0);
        TRow & row2 = row(algn, 1);
        std::vector<uint64_t> a_idx, b_idx;
        if (returnViewPosition)
        {
            a_idx.resize(length(source(row1)));
            b_idx.resize(length(source(row2)));
            for (uint64_t i = 0; i < length(source(row1)); ++i)
                a_idx[i] = toViewPosition(row1, i);
            for (uint64_t i = 0; i < length(source(row2)); ++i)
                b_idx[i] = toViewPosition(row2, i);
        }
        else
        {
            a_idx.resize(length(row1));
            b_idx.resize(length(row2));
            for (uint64_t i = 0; i < length(row1); ++i)
                a_idx[i] = toSourcePosition(row1, i);
            for (uint64_t i = 0; i < length(row2); ++i)
                b_idx[i] = toSourcePosition(row2, i);
        }
        return std::make_tuple(std::move(score), std::move(a_idx), std::move(b_idx));
    }

protected:

private:
    // methods
    // Copy constructor must not be used
    align_raw(const align_raw& object);

    // Assignment operator must not be used
    const align_raw& operator=(const align_raw& rhs);

    // init
    void init(align_raw_settings<TScore> const& settings)
    {
        m_gap_open_h = settings.m_gap_open_h;
        m_gap_open_v = settings.m_gap_open_v;
        m_gap_extension_h = settings.m_gap_extension_h;
        m_gap_extension_v = settings.m_gap_extension_v;
        m_dist_offset = settings.m_dist_offset;
        m_dist_min = settings.m_dist_min;
    }
    
    // member
    TScore m_gap_open_h;
    TScore m_gap_open_v;
    TScore m_gap_extension_h;
    TScore m_gap_extension_v;
    TScore m_dist_offset;
    TScore m_dist_min;
};


// -- exported functions - declarations --------------------------------------------

// -- exported global variables - declarations (should be empty)--------------------

#endif  // ALIGN_RAW_H