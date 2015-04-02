#include "align.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

#include <emmintrin.h>
#include <stdio.h>
#include <assert.h>

// defined here because it is used in fastAlignmentRoutine and calculateFlankScore
const short n_score = 0*4;

// define the DEBUG symbol to print the DP table and various running variables
//#define DEBUG 1

#define printxmm_0(label,var) \
  printf(" %s [0]=%4x %5i [1]=%4x %5i [2]=%4x %5i [3]=%4x %5i [4]=%4x %5i [5]=%4x %5i [6]=%4x %5i [7]=%4x %5i\n", \
    label, \
    _mm_extract_epi16(var,0),((_mm_extract_epi16(var,0))&0xffff) >> 2, \
    _mm_extract_epi16(var,1),((_mm_extract_epi16(var,1))&0xffff) >> 2, \
    _mm_extract_epi16(var,2),((_mm_extract_epi16(var,2))&0xffff) >> 2, \
    _mm_extract_epi16(var,3),((_mm_extract_epi16(var,3))&0xffff) >> 2, \
    _mm_extract_epi16(var,4),((_mm_extract_epi16(var,4))&0xffff) >> 2, \
    _mm_extract_epi16(var,5),((_mm_extract_epi16(var,5))&0xffff) >> 2, \
    _mm_extract_epi16(var,6),((_mm_extract_epi16(var,6))&0xffff) >> 2, \
    _mm_extract_epi16(var,7),((_mm_extract_epi16(var,7))&0xffff) >> 2);

#define printxmm(label,var) \
  printf(" %s [0]=%4x %5i [1]=%4x %5i [2]=%4x %5i [3]=%4x %5i [4]=%4x %5i [5]=%4x %5i [6]=%4x %5i [7]=%4x %5i\n", \
    label, \
    _mm_extract_epi16(var,0),((_mm_extract_epi16(var,0)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,1),((_mm_extract_epi16(var,1)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,2),((_mm_extract_epi16(var,2)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,3),((_mm_extract_epi16(var,3)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,4),((_mm_extract_epi16(var,4)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,5),((_mm_extract_epi16(var,5)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,6),((_mm_extract_epi16(var,6)-32768)&0xffff) >> 2, \
    _mm_extract_epi16(var,7),((_mm_extract_epi16(var,7)-32768)&0xffff) >> 2);


//_________________________________________________________________________________________________


int fastAlignmentRoutine(const char* seq1, const char* seq2, const char* qual2, int len1, int len2, int gapextend, int nucprior,
			 const char* localgapopen, char* aln1, char* aln2, int* firstpos) {

  // seq2 is the read; the shorter of the sequences
  // no checks for overflow are done

  // the bottom-left and top-right corners of the DP table are just
  // included at the extreme ends of the diagonal, which measures
  // n=8 entries diagonally across.  This fixes the length of the
  // longer (horizontal) sequence to 14 (2*8-2) more than the shorter

  assert( len1 == len2 + 2*8 - 1);

  // make sure that special cases at beginning and end (initialization!)
  // do not mix.  (Todo: Check this assertion.)
  assert( len1 > 8);

  const short gap_extend = gapextend*4;
  const short nuc_prior = nucprior*4;
  const char traceback = (aln1 != NULL);
  const short pos_inf = 0x7800;
  const int match_label = 0;
  const int insert_label = 1;
  const int delete_label = 3;

  register __m128i _m1;
  register __m128i _i1;
  register __m128i _d1;
  register __m128i _m2;
  register __m128i _i2;
  register __m128i _d2;
  __m128i _seq1win;
  __m128i _seq2win;
  __m128i _qual2win;
  __m128i _seq1nqual;   // 1 if N, pos_inf if not

  __m128i _gap_extend = _mm_set1_epi16( gap_extend );
  __m128i _nuc_prior = _mm_set1_epi16( nuc_prior );
  __m128i _three = _mm_set1_epi16( 3 );
  __m128i _initmask = _mm_set_epi16( 0,0,0,0,0,0,0,-1 );
  __m128i _initmask2 = _mm_set_epi16( 0,0,0,0,0,0,0,-0x8000 );
  __m128i _backpointers[ 2*(len1+8) ];

  // initialization
  _m1 = _mm_set1_epi16( pos_inf );
  _i1 = _m1;
  _d1 = _m1;
  _m2 = _m1;
  _i2 = _m1;
  _d2 = _m1;

  // sequence 1 is initialized with the n-long prefix, in forward direction
  // sequence 2 is initialized as empty; reverse direction
  _seq1win = _mm_set_epi16( seq1[7], seq1[6], seq1[5], seq1[4], seq1[3], seq1[2], seq1[1], seq1[0] );
  _seq2win = _m1;
  _qual2win = _mm_set1_epi16(64*4);

  // if N, make n_score; if != N, make pos_inf
  _seq1nqual = _mm_add_epi16( _mm_and_si128( _mm_cmpeq_epi16( _seq1win, 
							      _mm_set1_epi16( 'N' ) ), 
					     _mm_set1_epi16( n_score - pos_inf ) ), 
			      _mm_set1_epi16( pos_inf ) );


  __m128i _gap_open = _mm_set_epi16(4*localgapopen[7],4*localgapopen[6],4*localgapopen[5],4*localgapopen[4],
				    4*localgapopen[3],4*localgapopen[2],4*localgapopen[1],4*localgapopen[0]);

  short _score = 0;
  short minscore = pos_inf;
  short minscoreidx = -1;
  
  // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
  // of the seq2win/qual arrays, to simplify getting back pointers
  int s;

#ifdef DEBUG
  for(s=0;s<len1;s++) printf("%c",seq1[s]); printf("\n");
  for(s=0;s<len2;s++) printf("%c",seq2[s]); printf("\n");
  printxmm_0("gapextend ",_gap_extend);
#endif
  
  for (s=0; s<2*(len2+8-1 + 1); s+=2) {

#ifdef DEBUG
    printf("--(even)---- s = %u ------\n", s);
#endif
    
    // seq1 is current; seq2 needs updating
    _seq2win = _mm_slli_si128( _seq2win, 2 );
    _qual2win = _mm_slli_si128( _qual2win, 2 );
    if (s/2 < len2) {
      _seq2win = _mm_insert_epi16( _seq2win, seq2[ s/2 ], 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 4*qual2[ s/2 ], 0 );
    } else {
      _seq2win = _mm_insert_epi16( _seq2win, '0', 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 64*4, 0 );
    }

    //
    // S even
    //

#ifdef DEBUG
    printxmm("initmask  ",_initmask);
    printxmm("seq1win   ",_seq1win);
    printxmm_0("seq1nqual ",_seq1nqual);
    printxmm("seq2win u ",_seq2win);
    printxmm_0("qual2win u",_qual2win);
    printxmm_0("gapopen u ",_gap_open);
#endif
    
    // initialize to -0x8000
    _m1 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m1 ) );
    _m2 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m2 ) );
    _m1 = _mm_min_epi16( _m1, _mm_min_epi16( _i1, _d1 ) );

#ifdef DEBUG
    printxmm("min m1    ",_m1);
#endif
    
    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2

    if (s/2 >= len2) {
      switch (s/2 - len2) {
      case 0: _score = _mm_extract_epi16( _m1, 0 ); break;
      case 1: _score = _mm_extract_epi16( _m1, 1 ); break;
      case 2: _score = _mm_extract_epi16( _m1, 2 ); break;
      case 3: _score = _mm_extract_epi16( _m1, 3 ); break;
      case 4: _score = _mm_extract_epi16( _m1, 4 ); break;
      case 5: _score = _mm_extract_epi16( _m1, 5 ); break;
      case 6: _score = _mm_extract_epi16( _m1, 6 ); break;
      case 7: _score = _mm_extract_epi16( _m1, 7 ); break;
      }

      if (_score < minscore) {
          minscore = _score;
          minscoreidx = s;     // point back to the match state at this entry, so as not to
      }                        // have to store the state at s-2
    }

    _m1 = _mm_add_epi16( _m1, 
			 _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, 
									   _seq1win ), 
							  _qual2win ),
					_seq1nqual ) );

    _d1 = _mm_min_epi16( _mm_add_epi16( _d2,
					_gap_extend ),
			 _mm_add_epi16( _mm_min_epi16( _m2,
						       _i2 ),                // allow I->D
					_mm_srli_si128( _gap_open, 2 ) ) );

    _d1 = _mm_insert_epi16( _mm_slli_si128( _d1, 
					    2 ), 
			    pos_inf, 
			    0 );

    _i1 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i2,
						       _gap_extend ),
					_mm_add_epi16( _m2,
						       _gap_open ) ),
			 _nuc_prior );

#ifdef DEBUG
    printxmm("m1 &ptr   ",_m1);    
    printxmm("d1 &ptr   ",_d1);    
    printxmm("i1 &ptr   ",_i1);    
#endif
    
    // get back-pointers and store
    if (traceback) {
      _backpointers[ s ] = _mm_or_si128( _mm_or_si128( _mm_and_si128( _three, _m1 ),
						       _mm_slli_epi16( _mm_and_si128( _three, _i1 ), 2*insert_label ) ),
					 _mm_slli_epi16( _mm_and_si128( _three, _d1 ), 2*delete_label ) );
    
      // set state labels
      _m1 = _mm_andnot_si128( _three, _m1 );
      _i1 = _mm_or_si128( _mm_andnot_si128( _three, _i1 ), _mm_srli_epi16( _three, 1 ) );           
      _d1 = _mm_or_si128( _mm_andnot_si128( _three, _d1 ), _three );                                
    }

    //
    // S odd
    //

#ifdef DEBUG
    printf("--(odd)----- s = %u ------\n", s+1);
#endif
    
    // seq1 needs updating; seq2 is current
    char c = (8 + s/2 < len1) ? seq1[ 8+(s/2) ] : 'N';
    _seq1win   = _mm_insert_epi16( _mm_srli_si128( _seq1win,   2 ), c, 8-1 );
    _seq1nqual = _mm_insert_epi16( _mm_srli_si128( _seq1nqual, 2 ), (c=='N') ? n_score : pos_inf, 8-1 );
    _gap_open  = _mm_insert_epi16( _mm_srli_si128( _gap_open, 2 ),
				   4*localgapopen[ 8 + s/2 < len1 ? 8 + s/2 : len1-1 ],
				   8-1 );

#ifdef DEBUG
    printxmm("initmask  ",_initmask);
    printxmm("seq1win u ",_seq1win);
    printxmm_0("seq1nqul u",_seq1nqual);
    printxmm("seq2win   ",_seq2win);
    printxmm_0("qual2win  ",_qual2win);
    printxmm_0("gapopen u ",_gap_open);
#endif
    
    _initmask = _mm_slli_si128( _initmask, 2 );
    _initmask2 = _mm_slli_si128( _initmask2, 2 );
    _m2 = _mm_min_epi16( _m2, _mm_min_epi16( _i2, _d2 ) );

#ifdef DEBUG
    printxmm("min m2    ",_m2);
#endif
    
    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2
    if (s/2 >= len2) {
      switch (s/2 - len2) {
      case 0: _score = _mm_extract_epi16( _m2, 0 ); break;
      case 1: _score = _mm_extract_epi16( _m2, 1 ); break;
      case 2: _score = _mm_extract_epi16( _m2, 2 ); break;
      case 3: _score = _mm_extract_epi16( _m2, 3 ); break;
      case 4: _score = _mm_extract_epi16( _m2, 4 ); break;
      case 5: _score = _mm_extract_epi16( _m2, 5 ); break;
      case 6: _score = _mm_extract_epi16( _m2, 6 ); break;
      case 7: _score = _mm_extract_epi16( _m2, 7 ); break;
      }

      if (_score < minscore) {
          minscore = _score;
          minscoreidx = s+1;
      }
    }

    _m2 = _mm_add_epi16( _m2, 
			 _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, 
									   _seq1win ), 
							  _qual2win ),
					_seq1nqual ) );

    _d2 = _mm_min_epi16( _mm_add_epi16( _d1,
					_gap_extend ),
			 _mm_add_epi16( _mm_min_epi16( _m1,
						       _i1 ),  // allow I->D 
					_gap_open ) );

    _i2 = _mm_insert_epi16( _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _mm_srli_si128( _i1, 2 ),
									 _gap_extend ),
							  _mm_add_epi16( _mm_srli_si128( _m1, 2 ),
									 _gap_open ) ),
					   _nuc_prior ),
			    pos_inf,
			    8-1 );

#ifdef DEBUG
    printxmm("m2 &ptr   ",_m2);    
    printxmm("d2 &ptr   ",_d2);    
    printxmm("i2 &ptr   ",_i2);    
#endif

    // get back-pointers and store
    if (traceback) {
      _backpointers[ s+1 ] = _mm_or_si128( _mm_or_si128( _mm_and_si128( _three, _m2 ),
							 _mm_slli_epi16( _mm_and_si128( _three, _i2 ), 2*insert_label ) ),
					   _mm_slli_epi16( _mm_and_si128( _three, _d2 ), 2*delete_label ) );
      
      // set state labels
      _m2 = _mm_andnot_si128( _three, _m2 );
      _i2 = _mm_or_si128( _mm_andnot_si128( _three, _i2 ), _mm_srli_epi16( _three, 1 ) );          
      _d2 = _mm_or_si128( _mm_andnot_si128( _three, _d2 ), _three );                               
    }

  }

  // Backtrace.
  if (!traceback) {
    return (minscore + 0x8000) >> 2;
  }

  s = minscoreidx;    // point to the dummy match transition
  short i = s/2 - len2;
  short y = len2;
  short x = s - y;
  short alnidx = 0;
  short state = ((((short*)( _backpointers + s))[i]) >> (2*match_label)) & 3;
  s -= 2;

  // this is 2*y (s even) or 2*y+1 (s odd)
  while (y > 0)
  {
    short newstate = ((((short*)( _backpointers + s))[i]) >> (2*state)) & 3;

    if (state == match_label)
    {
      s -= 2;
      aln1[alnidx] = seq1[--x];
      aln2[alnidx] = seq2[--y];
    }
    else if (state == insert_label)
    {
      i += s&1;
      s -= 1;
      aln1[alnidx] = '-';
      aln2[alnidx] = seq2[--y];
    }
    else
    {
      s -= 1;
      i -= s&1;
      aln1[alnidx] = seq1[--x];
      aln2[alnidx] = '-';
    }

    state = newstate;
    alnidx++;
  }

  aln1[alnidx] = 0;
  aln2[alnidx] = 0;

  if (firstpos) *firstpos = x;

  // reverse them
  int j;

  for (i=0, j=alnidx-1; i<j; i++, j--)
  {
    x = aln1[i];
    y = aln2[i];
    aln1[i]=aln1[j];
    aln2[i]=aln2[j];
    aln1[j] = x;
    aln2[j] = y;
  }

#ifdef DEBUG
  printf("aln1 %s\n", aln1);
  printf("aln2 %s\n", aln2);
  printf("score %u\n", (minscore + 0x8000) >> 2 );
#endif
  
  return (minscore + 0x8000) >> 2;
}



//_________________________________________________________________________________________________


int calculateFlankScore(int hapLen, int hapFlank, const char* quals, const char* localgapopen, int gapextend, int nucprior,
			int firstpos, const char* aln1, const char* aln2) {

  char prevstate = 'M';
  int x = firstpos;     // index into haplotype
  int y = 0;            // index into read
  int i = 0;            // index into alignment
  int score = 0;        // alignment score (within flank)
  while (aln1[i]) {
    char newstate = 'M';
    if (aln1[i] == '-') newstate = 'I';
    if (aln2[i] == '-') newstate = 'D';  // can't be both '-'
    switch (newstate) {
    case 'M':
      if ( (aln1[i] != aln2[i]) &&
	   (x < hapFlank || x >= hapLen - hapFlank) ) {
	if (aln1[i] == 'N') {
	  score += n_score / 4;
	} else {
	  score += quals[y];
	}
      }
      ++x;
      ++y;
      break;
    case 'I':
      if (x < hapFlank || x >= hapLen - hapFlank) {
	if (prevstate == 'I') {
	  score += gapextend + nucprior;
	} else {
	  // gap open score is charged for insertions just after the corresponding base, hence the -1
	  score += localgapopen[x-1] + nucprior;
	}
      }
      ++y;
      break;
    case 'D':
      if (x < hapFlank || x >= hapLen - hapFlank) {
	if (prevstate == 'D') {
	  score += gapextend;
	} else {
	  score += localgapopen[x];
	}
      }
      ++x;
      break;
    }
    ++i;
    prevstate = newstate;
  }
  return score;
}
