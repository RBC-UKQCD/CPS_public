#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! commstest.C

  A series of tests of the functionality supplied by the MPI
  implementations of the QCDSP comms.

  $Id: commstest.C,v 1.4 2008-02-08 18:35:08 chulwoo Exp $

  A.N.Jackson: ajackson@epcc.ed.ac.uk                       */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdlib.h>
#include<comms/sysfunc_cps.h>
CPS_START_NAMESPACE

/*! Global integer to store the PE identity: */
int i_am;

/*! tokenizer_test: 
   Takes the teststring and splits it into tokens, delimited by any of
   the characters in the delim string.  Tests that the number of
   tokens found is equal to that expected (num_tokens).  Of course, if
   you change the test string, the expected number of tokens must also
   be changed. Anj. */
/*----------------------------------------------------------*/
void tokenizer_test(char* teststring, char* delim, int num_tokens) { 
    char *token, *tokpos, n_tok = 0;

    // Break the string up using the given delimiters:
    printf("Tokenising this: '%s'\n", teststring );
    tokpos = teststring;
    while( token = MPICommsStringTokenizer(teststring, delim, &tokpos ) ) {
	printf("%i: (%p %p) %s\n", n_tok, teststring, tokpos, token);
	n_tok++;
    }

    // Check if it worked:
    if( n_tok != num_tokens ) {
	printf("Error: the string tokenizer failed the test!\n");
	exit(EXIT_FAILURE);
    }

    return;
}

/*! Creates a unique identification number for an element of an array.
  The number is based on the coordinates of a processor in the grid
  combined with the array index.  The only compication is that the
  function actually takes arguments to specify the PE number in terms
  of its relative position to this PE.  ie, 0,0,0,0 specifies this
  PE. */
/*----------------------------------------------------------*/
int create_id( int Dt, int Dx, int Dy, int Dz, int num ) {
    int danum;
    danum = num;
    danum *= 1000; danum += (Dt+CoorT()+SizeT())%SizeT();
    danum *= 1000; danum += (Dx+CoorX()+SizeX())%SizeX();
    danum *= 1000; danum += (Dy+CoorY()+SizeY())%SizeY();
    danum *= 1000; danum += (Dz+CoorZ()+SizeZ())%SizeZ();
    return( danum );
}

/*! Looks up the RNG seeds using all four methods: */
/*----------------------------------------------------------*/
void rngseedtest( void ) {
    unsigned int rngseed, rngseedS, rngseedT, rngseedST;

    // Look-up RNG seeds:
    rngseed   = Seed();
    rngseedS  = SeedS();
    rngseedT  = SeedT();
    rngseedST = SeedST();

    // Synchronise the processors, and then print the seeds:
    sync();
    printf("[%i] seed=%i seedS=%i seedT=%i seedST=%i.\n", 
	   i_am, rngseed, rngseedS, rngseedT, rngseedST );
    return;
}

/*! Performs a simple test of contiguous and block strided int and
  IFloat transfers. */
/*----------------------------------------------------------*/
void basictransfertest(void) {
    int errorcount = 0;
    int iarray[100],iarray2[100];
    IFloat farray[100], farray2[100];

    // Set up some test data-types:
    printf("[%i] creating datatypes...\n", i_am);
    SCUDirArg iargS( (void*)iarray, SCU_XM, SCU_SEND, 100 );
    SCUDirArg iargRx( (void*)iarray2, SCU_XP, SCU_REC, 100 );
    SCUDirArg fargS, fargRx;
    fargS.Init((void*)farray, SCU_XP, SCU_SEND, 5, 10, 6);
    fargRx.Init((void*)farray2, SCU_XM, SCU_REC, 5, 10, 6);
    /* i.e. the strided sets look like this:
       #####-----#####-----#####-----.....
       where the #'s will be transmitted and the -'s won't */

    // Fill the outward arrays with something useful:
    printf("[%i] filling the arrays...\n", i_am);
    for( int i=0; i < 100; i++ ) {
	iarray[i] = create_id(0,0,0,0,i);
	farray[i] = (IFloat)create_id(0,0,0,0,100-i);
    }

    // Pass the data about:
    printf("[%i] launching comms...\n", i_am);
    printf("[%i] 1...\n", i_am);
    SCUTrans(&iargS); SCUTrans(&iargRx);
    printf("[%i] 2...\n", i_am);
    SCUTrans(&fargS); SCUTrans(&fargRx);
    SCUTransComplete();
    printf("[%i] comms have finished...\n", i_am);

    // Check the recieved numbers, and die if there are errors:
    printf("[%i] testing the data...\n", i_am);
    errorcount = 0;
    for( int i=0; i < 100; i++ ) {
      if( iarray2[i] != create_id(0,+1,0,0,i) ) errorcount++;
      if( farray2[i] != (IFloat)create_id(0,-1,0,0,100-i) && i%10<5 ) errorcount++;
      // Clear the send data, redy to be send back again:
      iarray[i] = 0; if( i%10<5 ) farray[i] = 0;
    }
    if( errorcount > 0 ) {
	printf("[%i] simple send+recieve failed, %i bad data chunks.\n",
	       i_am, errorcount );
	exit( EXIT_FAILURE );
    }
    printf("[%i] All data passed the test.\n", i_am);

    /* Now pass the data back, but use the SCUSetDMA and SCUTransAddr calls */
    iargS.Init( (void*)iarray2, SCU_XP, SCU_SEND, 100 );
    iargRx.Init( (void*)iarray, SCU_XM, SCU_REC, 100 );
    fargS.Init((void*)farray2, SCU_XM, SCU_SEND, 5, 10, 6);
    fargRx.Init((void*)farray, SCU_XP, SCU_REC, 5, 10, 6);

    // The comms:
    printf("[%i] launching comms...\n", i_am);
    printf("[%i] 1...\n", i_am);
    SCUSetDMA(&iargS); SCUSetDMA(&iargRx);
    printf("[%i] 2...\n", i_am);
    SCUTransAddr(&iargS); SCUTransAddr(&iargRx);
    printf("[%i] 3...\n", i_am);
    SCUTrans(&fargS); SCUTrans(&fargRx);
    SCUTransComplete();
    printf("[%i] comms have finished...\n", i_am);

    // Check the recieved numbers, and die if there are errors:
    printf("[%i] testing the data...\n", i_am);
    errorcount = 0;
    for( int i=0; i < 100; i++ ) {
      if( iarray[i] != create_id(0,0,0,0,i) ) errorcount++;
      if( farray[i] != (IFloat)create_id(0,0,0,0,100-i) ) errorcount++;
    }
    if( errorcount > 0 ) {
	printf("[%i] simple send+recieve failed, %i bad data chunks.\n",
	       i_am, errorcount );
	exit( EXIT_FAILURE );
    }
    printf("[%i] All data passed the test.\n", i_am);

    return;
}

/* Uses the appropriate MPI-layer added functionality to transfer data
   of the double kind in all four directions: */
void doubletransfertest( void ) {
    SCUDirArg **sendd, **recvd;
    double **darray, checkval;
    int errorcount, idim, idir, i, nmax = 250;

    // Create and fill the arrays on each PE:
    darray = new double*[NDIM];
    for( idim = 0; idim < NDIM; idim++ ) {
	darray[idim] = new double[nmax];
	for( i = 0; i < nmax; i++ ) {
	    darray[idim][i] = (double)create_id(0,0,0,0,i);
	}
    }

    // Set-up the double data objects:
    sendd = new SCUDirArg*[NDIM]; recvd = new SCUDirArg*[NDIM];
    for( idir = 0; idir < 2*NDIM; idir+=2 ) {
	idim = idir/2;
	sendd[idim] = new SCUDirArg();
	sendd[idim]->SetDataSize(sizeof(double));
	sendd[idim]->Init( (void*)darray[idim], (SCUDir)idir, SCU_SEND, nmax );
	recvd[idim] = new SCUDirArg();
	recvd[idim]->Init( (void*)darray[idim], (SCUDir)(idir+1), SCU_REC, nmax );
	recvd[idim]->SetDataSize(sizeof(double));
    }

    // Send the data out:
    printf("[%i] 4-way out: Sending double-data in all four directions.\n", i_am);
    SCUTrans( sendd, NDIM );
    SCUTrans( recvd, NDIM );
    SCUTransComplete();

    // Test it:
    errorcount = 0;
    for( idim = 0; idim < NDIM; idim++ ) {
	for( i = 0; i < nmax; i++ ) {
	    if( idim == 0 ) {
		checkval = (double)create_id(-1,0,0,0,i);
	    } else if( idim == 1 ) {
		checkval = (double)create_id(0,-1,0,0,i);
	    } else if( idim == 2 ) {
		checkval = (double)create_id(0,0,-1,0,i);
	    } else {
		checkval = (double)create_id(0,0,0,-1,i);
	    }
	    if( darray[idim][i] != checkval ) errorcount++;
	}
    }
    if( errorcount > 0 ) {
	printf("[%i] double-datatype send+recieve failed, %i bad data chunks.\n",
	       i_am, errorcount );
	exit( EXIT_FAILURE );
    }
    printf("[%i] 4-way out: All data passed the test.\n", i_am);

    // Send the data back:
    for( idir = 0; idir < 2*NDIM; idir+=2 ) {
	idim = idir/2;
	sendd[idim]->Init( (void*)darray[idim], (SCUDir)(idir+1), SCU_SEND, nmax );
	recvd[idim]->Init( (void*)darray[idim], (SCUDir)idir, SCU_REC, nmax );
    }
    printf("[%i] 4-way back: Sending double-data in all four directions.\n", i_am);
    SCUSetDMA( sendd, NDIM ); SCUTransAddr( sendd, NDIM );
    SCUSetDMA( recvd, NDIM ); SCUTransAddr( recvd, NDIM );
    SCUTransComplete();

    // Test it:
    errorcount = 0;
    for( idim = 0; idim < NDIM; idim++ ) {
	for( i = 0; i < nmax; i++ ) {
	    checkval = (double)create_id(0,0,0,0,i);
	    if( darray[idim][i] != checkval ) errorcount++;
	}
    }
    if( errorcount > 0 ) {
	printf("[%i] double-datatype send+recieve failed, %i bad data chunks.\n",
	       i_am, errorcount );
	exit( EXIT_FAILURE );
    }
    printf("[%i] 4-way back: All data passed the test.\n", i_am);

    return;
}

/*! main:
   Controls the sequence of calls to test subroutines: */
/*----------------------------------------------------------*/
int main( int argc, char** argv ) {
    int i_am_root = FALSE, num_of_pelem, mypos[NDIM], gridsize[NDIM];

    // Initialise the communications:
    SCUCommsInit();

    // Determine if this is the 0,0,0,0 PE:
    if( CoorT() == 0 &&  CoorX() == 0 &&  CoorY() == 0 &&  CoorZ() == 0 ) {
	i_am_root = TRUE;
    }

    // Testing the string tokenizer, this string should split into 10 tokens.
    if( i_am_root ) {
	printf("I am the root processor.\n");
	tokenizer_test( (char*)"This {is || == } 1 very simple\ntokenizer  t\test."
			, (char*)"{} =,\n\t", 9); 
    }

    // Look-up the number of PEs, and this ones id
    num_of_pelem = NumNodes();
    i_am = UniqueID();

    // Look-up grid aspects:
    mypos[SCU_T] = CoorT(); 
    mypos[SCU_X] = CoorX(); 
    mypos[SCU_Y] = CoorY();
    mypos[SCU_Z] = CoorZ();
    gridsize[SCU_T] = SizeT();
    gridsize[SCU_X] = SizeX();
    gridsize[SCU_Y] = SizeY();
    gridsize[SCU_Z] = SizeZ();
    printf("I am #%i of %i processing elements, in a grid of dimensions [%i,%i,%i,%i].\nMy position in the grid is [%i,%i,%i,%i].\n",
	   i_am, num_of_pelem, 
	   gridsize[SCU_T], gridsize[SCU_X], gridsize[SCU_Y], gridsize[SCU_Z],
	   mypos[SCU_T], mypos[SCU_X], mypos[SCU_Y], mypos[SCU_Z]
	   );

    // RNG seed-loading tests:
    rngseedtest();

    // SCURemap test?

    // Perform some simple 4-byte data (IFloat and int) tests:
    basictransfertest();

    // Perform 8-byte 4-way transfer test:
    doubletransfertest();

    // If we get this far, the all is good:
    printf("[%i] ***** All tests passed successfully! *****\n",i_am);

    exit(EXIT_SUCCESS);
}

CPS_END_NAMESPACE
