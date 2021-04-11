/*
+--------------------------------------------------------------------------+
| CHStone : a suite of benchmark programs for C-based High-Level Synthesis |
| ======================================================================== |
|                                                                          |
| * Collected and Modified : Y. Hara, H. Tomiyama, S. Honda,               |
|                            H. Takada and K. Ishii                        |
|                            Nagoya University, Japan                      |
|                                                                          |
| * Remark :                                                               |
|    1. This source code is modified to unify the formats of the benchmark |
|       programs in CHStone.                                               |
|    2. Test vectors are added for CHStone.                                |
|    3. If "main_result" is 0 at the end of the program, the program is    |
|       correctly executed.                                                |
|    4. Please follow the copyright of each benchmark program.             |
+--------------------------------------------------------------------------+
*/
/*************************************************************************/
/*                                                                       */
/*   SNU-RT Benchmark Suite for Worst Case Timing Analysis               */
/*   =====================================================               */
/*                              Collected and Modified by S.-S. Lim      */
/*                                           sslim@archi.snu.ac.kr       */
/*                                         Real-Time Research Group      */
/*                                        Seoul National University      */
/*                                                                       */
/*                                                                       */
/*        < Features > - restrictions for our experimental environment   */
/*                                                                       */
/*          1. Completely structured.                                    */
/*               - There are no unconditional jumps.                     */
/*               - There are no exit from loop bodies.                   */
/*                 (There are no 'break' or 'return' in loop bodies)     */
/*          2. No 'switch' statements.                                   */
/*          3. No 'do..while' statements.                                */
/*          4. Expressions are restricted.                               */
/*               - There are no multiple expressions joined by 'or',     */
/*                'and' operations.                                      */
/*          5. No library calls.                                         */
/*               - All the functions needed are implemented in the       */
/*                 source file.                                          */
/*                                                                       */
/*                                                                       */
/*************************************************************************/
/*                                                                       */
/*  FILE: adpcm.c                                                        */
/*  SOURCE : C Algorithms for Real-Time DSP by P. M. Embree              */
/*                                                                       */
/*  DESCRIPTION :                                                        */
/*                                                                       */
/*     CCITT G.722 ADPCM (Adaptive Differential Pulse Code Modulation)   */
/*     algorithm.                                                        */
/*     16khz sample rate data is stored in the array test_data[SIZE].    */
/*     Results are stored in the array compressed[SIZE] and result[SIZE].*/
/*     Execution time is determined by the constant SIZE (default value  */
/*     is 2000).                                                         */
/*                                                                       */
/*  REMARK :                                                             */
/*                                                                       */
/*  EXECUTION TIME :                                                     */
/*                                                                       */
/*                                                                       */
/*************************************************************************/
#include <stdio.h>

int encode (int, int);
void decode (int);
int filtez (int *bpl, int *dlt);
void upzero (int dlt, int *dlti, int *bli);
int filtep (int rlt1, int al1, int rlt2, int al2);
int quantl (int el, int detl);
int logscl (int il, int nbl);
int scalel (int nbl, int shift_constant);
int uppol2 (int al1, int al2, int plt, int plt1, int plt2);
int uppol1 (int al1, int apl2, int plt, int plt1);
int logsch (int ih, int nbh);
void reset ();


int tqmf[24];


const int h[24] = {
  12, -44, -44, 212, 48, -624, 128, 1448,
  -840, -3220, 3804, 15504, 15504, 3804, -3220, -840,
  1448, 128, -624, 48, 212, -44, -44, 12
};

int xl, xh;


int accumc[11], accumd[11];


int xout1, xout2;

int xs, xd;



int il, szl, spl, sl, el;

const int qq4_code4_table[16] = {
  0, -20456, -12896, -8968, -6288, -4240, -2584, -1200,
  20456, 12896, 8968, 6288, 4240, 2584, 1200, 0
};


const int qq6_code6_table[64] = {
  -136, -136, -136, -136, -24808, -21904, -19008, -16704,
  -14984, -13512, -12280, -11192, -10232, -9360, -8576, -7856,
  -7192, -6576, -6000, -5456, -4944, -4464, -4008, -3576,
  -3168, -2776, -2400, -2032, -1688, -1360, -1040, -728,
  24808, 21904, 19008, 16704, 14984, 13512, 12280, 11192,
  10232, 9360, 8576, 7856, 7192, 6576, 6000, 5456,
  4944, 4464, 4008, 3576, 3168, 2776, 2400, 2032,
  1688, 1360, 1040, 728, 432, 136, -432, -136
};

int delay_bpl[6];

int delay_dltx[6];

const int wl_code_table[16] = {
  -60, 3042, 1198, 538, 334, 172, 58, -30,
  3042, 1198, 538, 334, 172, 58, -30, -60
};

const int ilb_table[32] = {
  2048, 2093, 2139, 2186, 2233, 2282, 2332, 2383,
  2435, 2489, 2543, 2599, 2656, 2714, 2774, 2834,
  2896, 2960, 3025, 3091, 3158, 3228, 3298, 3371,
  3444, 3520, 3597, 3676, 3756, 3838, 3922, 4008
};

int nbl;
int al1, al2;
int plt, plt1, plt2;
int dlt;
int rlt, rlt1, rlt2;


const int decis_levl[30] = {
  280, 576, 880, 1200, 1520, 1864, 2208, 2584,
  2960, 3376, 3784, 4240, 4696, 5200, 5712, 6288,
  6864, 7520, 8184, 8968, 9752, 10712, 11664, 12896,
  14120, 15840, 17560, 20456, 23352, 32767
};

int detl;


const int quant26bt_pos[31] = {
  61, 60, 59, 58, 57, 56, 55, 54,
  53, 52, 51, 50, 49, 48, 47, 46,
  45, 44, 43, 42, 41, 40, 39, 38,
  37, 36, 35, 34, 33, 32, 32
};


const int quant26bt_neg[31] = {
  63, 62, 31, 30, 29, 28, 27, 26,
  25, 24, 23, 22, 21, 20, 19, 18,
  17, 16, 15, 14, 13, 12, 11, 10,
  9, 8, 7, 6, 5, 4, 4
};


int deth;
int sh;
int eh;

const int qq2_code2_table[4] = {
  -7408, -1616, 7408, 1616
};

const int wh_code_table[4] = {
  798, -214, 798, -214
};


int dh, ih;
int nbh, szh;
int sph, ph, yh, rh;

int delay_dhx[6];

int delay_bph[6];

int ah1, ah2;
int ph1, ph2;
int rh1, rh2;


int ilr, rl;
int dec_deth, dec_detl, dec_dlt;

int dec_del_bpl[6];

int dec_del_dltx[6];

int dec_plt, dec_plt1, dec_plt2;
int dec_szl, dec_spl, dec_sl;
int dec_rlt1, dec_rlt2, dec_rlt;
int dec_al1, dec_al2;
int dl;
int dec_nbl, dec_dh, dec_nbh;


int dec_del_bph[6];

int dec_del_dhx[6];

int dec_szh;

int dec_rh1, dec_rh2;
int dec_ah1, dec_ah2;
int dec_ph, dec_sph;

int dec_sh;

int dec_ph1, dec_ph2;



int
abs (int n)
{
  int m;

  if (n >= 0)
    m = n;
  else
    m = -n;
  return m;
}

int
encode (int xin1, int xin2)
{
  #pragma HLS ALLOCATION instances=add limit=4 operation
  #pragma HLS ALLOCATION instances=sub limit=4 operation
  int i;
  const int *h_ptr;
  int *tqmf_ptr, *tqmf_ptr1;
  long int xa, xb;
  int decis;


  h_ptr = h;
  tqmf_ptr = tqmf;
  xa = (long) (*tqmf_ptr++) * (*h_ptr++);
  xb = (long) (*tqmf_ptr++) * (*h_ptr++);

  for (i = 0; i < 10; i++)
    {
      #pragma HLS pipeline
      xa += (long) (*tqmf_ptr++) * (*h_ptr++);
      xb += (long) (*tqmf_ptr++) * (*h_ptr++);
    }

  xa += (long) (*tqmf_ptr++) * (*h_ptr++);
  xb += (long) (*tqmf_ptr) * (*h_ptr++);

  tqmf_ptr = tqmf;
   for(i=23;i>=2;i--){
      #pragma HLS unroll
  	  tqmf_ptr[i] = tqmf_ptr[i-2];
   }

    tqmf_ptr[1] = xin1;
    tqmf_ptr[0] = xin2;



  xl = (xa + xb) >> 15;
  xh = (xa - xb) >> 15;


  szl = filtez (delay_bpl, delay_dltx);


  spl = filtep (rlt1, al1, rlt2, al2);


  sl = szl + spl;
  el = xl - sl;


  il = quantl (el, detl);


  dlt = ((long) detl * qq4_code4_table[il >> 2]) >> 15;


  nbl = logscl (il, nbl);
  detl = scalel (nbl, 8);
  plt = dlt + szl;
  upzero (dlt, delay_dltx, delay_bpl);


  al2 = uppol2 (al1, al2, plt, plt1, plt2);


  al1 = uppol1 (al1, al2, plt, plt1);


  rlt = sl + dlt;


  rlt2 = rlt1;
  rlt1 = rlt;
  plt2 = plt1;
  plt1 = plt;



  szh = filtez (delay_bph, delay_dhx);

  sph = filtep (rh1, ah1, rh2, ah2);


  sh = sph + szh;

  eh = xh - sh;


  if (eh >= 0)
    {
      ih = 3;
    }
  else
    {
      ih = 1;
    }
  decis = (564L * (long) deth) >> 12L;
  if (abs (eh) > decis)
    ih--;


  dh = ((long) deth * qq2_code2_table[ih]) >> 15L;


  nbh = logsch (ih, nbh);


  deth = scalel (nbh, 10);


  ph = dh + szh;


  upzero (dh, delay_dhx, delay_bph);


  ah2 = uppol2 (ah1, ah2, ph, ph1, ph2);


  ah1 = uppol1 (ah1, ah2, ph, ph1);


  yh = sh + dh;


  rh2 = rh1;
  rh1 = yh;
  ph2 = ph1;
  ph1 = ph;


  return (il | (ih << 6));
}



void
decode (int input)
{
 #pragma HLS ALLOCATION instances=add limit=3 operation
  int i;
  long int xa1, xa2;
  const int *h_ptr;
  int *ac_ptr, *ac_ptr1, *ad_ptr, *ad_ptr1;


  ilr = input & 0x3f;
  ih = input >> 6;




  dec_szl = filtez (dec_del_bpl, dec_del_dltx);


  dec_spl = filtep (dec_rlt1, dec_al1, dec_rlt2, dec_al2);

  dec_sl = dec_spl + dec_szl;


  dec_dlt = ((long) dec_detl * qq4_code4_table[ilr >> 2]) >> 15;


  dl = ((long) dec_detl * qq6_code6_table[il]) >> 15;

  rl = dl + dec_sl;


  dec_nbl = logscl (ilr, dec_nbl);


  dec_detl = scalel (dec_nbl, 8);


  dec_plt = dec_dlt + dec_szl;


  upzero (dec_dlt, dec_del_dltx, dec_del_bpl);


  dec_al2 = uppol2 (dec_al1, dec_al2, dec_plt, dec_plt1, dec_plt2);


  dec_al1 = uppol1 (dec_al1, dec_al2, dec_plt, dec_plt1);


  dec_rlt = dec_sl + dec_dlt;


  dec_rlt2 = dec_rlt1;
  dec_rlt1 = dec_rlt;
  dec_plt2 = dec_plt1;
  dec_plt1 = dec_plt;


  dec_szh = filtez (dec_del_bph, dec_del_dhx);


  dec_sph = filtep (dec_rh1, dec_ah1, dec_rh2, dec_ah2);


  dec_sh = dec_sph + dec_szh;


  dec_dh = ((long) dec_deth * qq2_code2_table[ih]) >> 15L;


  dec_nbh = logsch (ih, dec_nbh);


  dec_deth = scalel (dec_nbh, 10);


  dec_ph = dec_dh + dec_szh;


  upzero (dec_dh, dec_del_dhx, dec_del_bph);


  dec_ah2 = uppol2 (dec_ah1, dec_ah2, dec_ph, dec_ph1, dec_ph2);


  dec_ah1 = uppol1 (dec_ah1, dec_ah2, dec_ph, dec_ph1);


  rh = dec_sh + dec_dh;


  dec_rh2 = dec_rh1;
  dec_rh1 = rh;
  dec_ph2 = dec_ph1;
  dec_ph1 = dec_ph;




  xd = rl - rh;
  xs = rl + rh;


    h_ptr = h;
    ac_ptr = accumc;
    ad_ptr = accumd;
    xa1 = (long) xd *(h_ptr[0]);
    xa2 = (long) xs *(h_ptr[1]);


    int j;
    for(i=0,j=2;i<=10;i++,j=j+2){
          #pragma HLS unroll
          xa1 += (long) (ac_ptr[i])*(h_ptr[j]);
          xa2 += (long) (ad_ptr[i])*(h_ptr[j+1]);
    }


  xout1 = xa1 >> 14;
  xout2 = xa2 >> 14;



  for(i = 10;i>=1;i--){
      #pragma HLS unroll
  	  ac_ptr[i] = ac_ptr[i-1];
  	  ad_ptr[i] = ad_ptr[i-1];
    }
    ac_ptr[0] = xd;
    ad_ptr[0] = xs;



}

void
reset ()
{
  int i;

  detl = dec_detl = 32;
  deth = dec_deth = 8;
  nbl = al1 = al2 = plt1 = plt2 = rlt1 = rlt2 = 0;
  nbh = ah1 = ah2 = ph1 = ph2 = rh1 = rh2 = 0;
  dec_nbl = dec_al1 = dec_al2 = dec_plt1 = dec_plt2 = dec_rlt1 = dec_rlt2 = 0;
  dec_nbh = dec_ah1 = dec_ah2 = dec_ph1 = dec_ph2 = dec_rh1 = dec_rh2 = 0;

  for (i = 0; i < 6; i++)
    {
      #pragma HLS unroll
      delay_dltx[i] = 0;
      delay_dhx[i] = 0;
      dec_del_dltx[i] = 0;
      dec_del_dhx[i] = 0;
      delay_bpl[i] = 0;
      delay_bph[i] = 0;
      dec_del_bpl[i] = 0;
      dec_del_bph[i] = 0;
    }

    for (i = 0; i < 24; i++)
     {
       #pragma HLS unroll
       tqmf[i] = 0;// i<23
       if(i<12){
     	  accumc[i] = 0;
     	  accumd[i] = 0;
       }
     }
}


int filtez(int bpl[],int dlt[]){
  #pragma HLS array_map variable=bpl instance=A vertical
  #pragma HLS array_map variable=dlt instance=A vertical
	long int z1,a1,a2,a3;
	z1 = 0;
	for(int i=0;i<6;i++){
        #pragma HLS unroll factor=2
		a1 = (long)(bpl[i])*(dlt[i]);
        z1 += a1;
	}
	return ((int)(z1 >> 14));
}




int
filtep (int rlt1, int al1, int rlt2, int al2)
{
  long int pl, pl2;
  pl = 2 * rlt1;
  pl = (long) al1 *pl;
  pl2 = 2 * rlt2;
  pl += (long) al2 *pl2;
  return ((int) (pl >> 15));
}


int
quantl (int el, int detl)
{
	  int ril, mil;
	  long int wd, decis,decis1;


	  wd = abs (el);

	    for (mil = 0; mil < 30; mil=mil+2)
	    {
		  decis = (decis_levl[mil] * (long) detl) >> 15L;
		  decis1 = (decis_levl[mil+1] * (long) detl) >> 15L;
		   if (wd <= decis)
		   	   break;
		         else if(wd <= decis1){
		      	   mil=mil+1;
		      	   break;
		         }
	    }

	  if (el >= 0)
	    ril = quant26bt_pos[mil];
	  else
	    ril = quant26bt_neg[mil];
	  return (ril);
}



int
logscl (int il, int nbl)
{
 #pragma HLS inline
  long int wd;
  wd = ((long) nbl * 127L) >> 7L;	/* leak factor 127/128 */
  nbl = (int) wd + wl_code_table[il >> 2];
  if (nbl < 0)
    nbl = 0;
  if (nbl > 18432)
    nbl = 18432;
  return (nbl);
}



int
scalel (int nbl, int shift_constant)
{
  #pragma HLS inline
  int wd1, wd2, wd3;
  wd1 = (nbl >> 6) & 31;
  wd2 = nbl >> 11;
  wd3 = ilb_table[wd1] >> (shift_constant + 1 - wd2);
  return (wd3 << 3);
}



void
upzero (int dlt, int *dlti, int *bli)
{
#pragma HLS allocation instances=sub limit=2 operation
#pragma HLS allocation instances=add limit=2 operation
  int i, wd2, wd3;

  if (dlt == 0)
    {
      for (i = 0; i < 6; i++)
	{
      #pragma HLS pipeline
	  bli[i] = (int) ((255L * bli[i]) >> 8L);
	}
    }
  else
    {
      for (i = 0; i < 6; i++)
	{
      #pragma HLS unroll
	  if ((long) dlt * dlti[i] >= 0)
	    wd2 = 128;
	  else
	    wd2 = -128;
	  wd3 = (int) ((255L * bli[i]) >> 8L);
	  bli[i] = wd2 + wd3;
	}
    }

  dlti[5] = dlti[4];
  dlti[4] = dlti[3];
  dlti[3] = dlti[2];
  dlti[2] = dlti[1];
  dlti[1] = dlti[0];
  dlti[0] = dlt;
}


int
uppol2 (int al1, int al2, int plt, int plt1, int plt2)
{
 #pragma HLS inline
  long int wd2, wd4;
  int apl2;
  wd2 = 4L * (long) al1;
  if ((long) plt * plt1 >= 0L)
    wd2 = -wd2;
  wd2 = wd2 >> 7;
  if ((long) plt * plt2 >= 0L)
    {
      wd4 = wd2 + 128;
    }
  else
    {
      wd4 = wd2 - 128;
    }
  apl2 = wd4 + (127L * (long) al2 >> 7L);


  if (apl2 > 12288)
    apl2 = 12288;
  if (apl2 < -12288)
    apl2 = -12288;
  return (apl2);
}



int
uppol1 (int al1, int apl2, int plt, int plt1)
{
  #pragma HLS inline
  long int wd2;
  int wd3, apl1;
  wd2 = ((long) al1 * 255L) >> 8L;
  if ((long) plt * plt1 >= 0L)
    {
      apl1 = (int) wd2 + 192;
    }
  else
    {
      apl1 = (int) wd2 - 192;
    }

  wd3 = 15360 - apl2;
  if (apl1 > wd3)
    apl1 = wd3;
  if (apl1 < -wd3)
    apl1 = -wd3;
  return (apl1);
}



int
logsch (int ih, int nbh)
{
  int wd;
  wd = ((long) nbh * 127L) >> 7L;
  nbh = wd + wh_code_table[ih];
  if (nbh < 0)
    nbh = 0;
  if (nbh > 22528)
    nbh = 22528;
  return (nbh);
}


#define SIZE 100
#define IN_END 100



const int test_data[SIZE] = {
  0x44, 0x44, 0x44, 0x44, 0x44,
  0x44, 0x44, 0x44, 0x44, 0x44,
  0x44, 0x44, 0x44, 0x44, 0x44,
  0x44, 0x44, 0x43, 0x43, 0x43,
  0x43, 0x43, 0x43, 0x43, 0x42,
  0x42, 0x42, 0x42, 0x42, 0x42,
  0x41, 0x41, 0x41, 0x41, 0x41,
  0x40, 0x40, 0x40, 0x40, 0x40,
  0x40, 0x40, 0x40, 0x3f, 0x3f,
  0x3f, 0x3f, 0x3f, 0x3e, 0x3e,
  0x3e, 0x3e, 0x3e, 0x3e, 0x3d,
  0x3d, 0x3d, 0x3d, 0x3d, 0x3d,
  0x3c, 0x3c, 0x3c, 0x3c, 0x3c,
  0x3c, 0x3c, 0x3c, 0x3c, 0x3b,
  0x3b, 0x3b, 0x3b, 0x3b, 0x3b,
  0x3b, 0x3b, 0x3b, 0x3b, 0x3b,
  0x3b, 0x3b, 0x3b, 0x3b, 0x3b,
  0x3b, 0x3b, 0x3b, 0x3b, 0x3b,
  0x3b, 0x3b, 0x3c, 0x3c, 0x3c,
  0x3c, 0x3c, 0x3c, 0x3c, 0x3c
};

void
adpcm_main (int compressed[SIZE],int result[SIZE])
{

  #pragma HLS inline
  int i, j;
  reset ();

  j = 10;

    for (i = 0; i < IN_END; i += 2)
    {
      compressed[i / 2] = encode (test_data[i], test_data[i + 1]);
    }
  for (i = 0; i < IN_END; i += 2)
    {
      decode (compressed[i / 2]);
      result[i] = xout1;
      result[i + 1] = xout2;
    }
}

