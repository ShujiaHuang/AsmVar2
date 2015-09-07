#ifndef ALIGN_H
#define ALIGN_H

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009, Nov 2014
 It may not be distributed, made public, or used in other software without the permission of the copyright holder

 Modify By Shujia Huang 2015-04-08
******************************************************************************************************************/

// Add by Shujia Huang 2015-04-08
typedef struct tAlignTag {
    int score;
    int pos;
} AlignTag, *AlignTagPointer;
// Add End

// The score of 'fastAlignmentRoutine' is the smaller the better and exactly mapping will be 0, others will larger than 0.
AlignTagPointer fastAlignmentRoutine(const char* seq1, const char* seq2, const char* qual2, int len1, int len2, int gapextend, int nucprior, const char* localgapopen, char* aln1, char* aln2);

int calculateFlankScore(int hapLen, int hapFlank, const char* quals, const char* localGapOpen, int gapExtend, int nucprior, int firstpos, const char* aln1, const char* aln2);

#endif

