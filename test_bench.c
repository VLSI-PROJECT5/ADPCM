#define SIZE 100
#define IN_END 100
#include<stdio.h>

int compressed[SIZE], result[SIZE];
const int test_compressed[SIZE] = {
  0xfd, 0xde, 0x77, 0xba, 0xf2,
  0x90, 0x20, 0xa0, 0xec, 0xed,
  0xef, 0xf1, 0xf3, 0xf4, 0xf5,
  0xf5, 0xf5, 0xf5, 0xf6, 0xf6,
  0xf6, 0xf7, 0xf8, 0xf7, 0xf8,
  0xf7, 0xf9, 0xf8, 0xf7, 0xf9,
  0xf8, 0xf8, 0xf6, 0xf8, 0xf8,
  0xf7, 0xf9, 0xf9, 0xf9, 0xf8,
  0xf7, 0xfa, 0xf8, 0xf8, 0xf7,
  0xfb, 0xfa, 0xf9, 0xf8, 0xf8
};
const int test_result[SIZE] = {
  0, 0xffffffff, 0xffffffff, 0, 0,
  0xffffffff, 0, 0, 0xffffffff, 0xffffffff,
  0, 0, 0x1, 0x1, 0,
  0xfffffffe, 0xffffffff, 0xfffffffe, 0, 0xfffffffc,
  0x1, 0x1, 0x1, 0xfffffffb, 0x2,
  0x2, 0x3, 0xb, 0x14, 0x14,
  0x16, 0x18, 0x20, 0x21, 0x26,
  0x27, 0x2e, 0x2f, 0x33, 0x32,
  0x35, 0x33, 0x36, 0x34, 0x37,
  0x34, 0x37, 0x35, 0x38, 0x36,
  0x39, 0x38, 0x3b, 0x3a, 0x3f,
  0x3f, 0x40, 0x3a, 0x3d, 0x3e,
  0x41, 0x3c, 0x3e, 0x3f, 0x42,
  0x3e, 0x3b, 0x37, 0x3b, 0x3e,
  0x41, 0x3b, 0x3b, 0x3a, 0x3b,
  0x36, 0x39, 0x3b, 0x3f, 0x3c,
  0x3b, 0x37, 0x3b, 0x3d, 0x41,
  0x3d, 0x3e, 0x3c, 0x3e, 0x3b,
  0x3a, 0x37, 0x3b, 0x3e, 0x41,
  0x3c, 0x3b, 0x39, 0x3a, 0x36
};


int main ()
{
  int i;
  int main_result;

      main_result = 0;
      adpcm_main (compressed,result);
      for (i = 0; i < IN_END / 2; i++)
	{
	  if (compressed[i] != test_compressed[i])
	    {
	      main_result += 1;
	    }
	}
      for (i = 0; i < IN_END; i++)
	{
	  if (result[i] != test_result[i])
	    {
	      main_result += 1;
	    }
	}
      printf ("main result = %d\n", main_result);
      return main_result;
    }
