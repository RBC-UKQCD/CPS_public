



#include <stdio.h>
#include <stdlib.h>





// needed???
#include <sstream>






int main(int argc,char *argv[])
{


  
  /* get parameter from command line */

  /*  parameters: */
  /* x.x  filename X Y Z T S */
          

  if ( argc != 7 )
    {
      printf("x.x filename X Y Z T S \n");
      exit(-13);
    }
  
 
  char inNum[5], outNum[5];
  int xdim = atoi(argv[2]);
  int ydim = atoi(argv[3]); 
  int zdim = atoi(argv[4]);
  int tdim = atoi(argv[5]);
  int sdim = atoi(argv[6]);
  int in, out, count, vol, svol;	


  count=0;
  svol=xdim*ydim*zdim;
  vol= svol*tdim*sdim;

  for(int ss(0); ss < sdim; ++ss)  
   for(int tt(0); tt < tdim; ++tt)
    for(int rr(0); rr < svol; ++rr)  
	{

        out= count;
     	in  = rr + tt*svol*sdim + ss*svol;

        if( in < 1000 )
         {
          if( in < 100 )
           {
            if( in < 10 )
	     sprintf(inNum,"000%1i",in);
            else
             sprintf(inNum,"00%2i",in);
           }
          else
           sprintf(inNum,"0%3i",in);
         }
        else
         sprintf(inNum,"%4i",in);
      

        if( out < 1000 )
         {
          if( out < 100 )
           {
            if( out < 10 )
             sprintf(outNum,"000%1i",out);
            else
             sprintf(outNum,"00%2i",out);
           }
          else
           sprintf(outNum,"0%3i",out);
         }
        else
         sprintf(outNum,"%4i",out);



	//printf("%s\n",inNum);
 
	//printf("%s\n",outNum); 

        printf("cp  wrongorder.%s.vol%s %s.vol%s\n", argv[1], inNum, argv[1], outNum);

	++count;

        }






}


