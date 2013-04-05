// ======================================================================
/*


  
    * * * * * * *  N O T I C E  * * * * * * * *  

    On  2012-01-01  

    This code is now temporaliry modifed for the 4dim even/odd preconditioning.

    The modified parts code is distinguished by following preprocessor flag.


    
*/

#define FOUR_DIM_EVEN_ODD_PRECONDITION

// ======================================================================


#ifdef  FOUR_DIM_EVEN_ODD_PRECONDITION

#include "dwf_mdagm_4dpc.C"

#else

#include "dwf_mdagm_5dpc.C"

#endif
