#include <stdlib.h>
#include "ugfsr.h"
#include "unif01.h"
#include "bbattery.h"

unif01_Gen * unif_CreateCPS (int *argc, char ***argv, char *doarg_name,  int take, int drop);
void unif_DeleteCPS( unif01_Gen *gen);



int main (int argc, char **argv)
{
   unif01_Gen *gen;
//`   gen = ulcg_CreateLCG (2147483647, 16807, 0, 12345);
//   gen = ucarry_CreateRanlux (223, 0);
//   gen = ugfsr_CreateMT19937_02(83647, NULL, 0);
   gen = unif_CreateCPS(&argc, &argv, "do_arg.vml", 8,8);
   bbattery_SmallCrush (gen);
   bbattery_BigCrush (gen);
//   ugfsr_DeleteGen (gen);
	unif_DeleteCPS (gen);
   return 0;
}
