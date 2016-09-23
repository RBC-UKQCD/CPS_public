#include<config.h>
/*----------------------------------------------------------*/
/*!\file
  \brief  The QCD I/O Interface.

  $Id: qcdio.C,v 1.9 2008/02/08 18:35:08 chulwoo Exp $
*/
/*  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
   CVS keywords
 
   $Author: chulwoo $ 
   $Date: 2008/02/08 18:35:08 $
   $Header: /space/cvs/cps/cps++/src/util/qcdio/comsrc/qcdio.C,v 1.9 2008/02/08 18:35:08 chulwoo Exp $
   $Id: qcdio.C,v 1.9 2008/02/08 18:35:08 chulwoo Exp $
   $Name: v5_0_16_hantao_io_test_v7 $
   $Locker:  $
   $RCSfile: qcdio.C,v $
   $Revision: 1.9 $
   $Source: /space/cvs/cps/cps++/src/util/qcdio/comsrc/qcdio.C,v $
   $State: Exp $  */ 
/*----------------------------------------------------------*/

#include <util/qcdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <util/gjp.h>
#include <util/error.h>
CPS_START_NAMESPACE




static int qcdio_latt[4];
static  int qcdio_normalize_the_data = 1;
void qloadsave_pump_data( int* pos, unsigned char* buf, int byte_size, int inout);

  /* ------------------------------------------------------
  	Gauge configuration I/O:
    ------------------------------------------------------ */
  
  //! Routine for performing the bytes-swapping of double precision raw data files.
  /*!
    This is a utility routine used by functions qload_unpackstrip and
    qsave_packstrip.\n
    This subroutine taken fram swap.cpp by Balint Joo.
  */
  void block_swap_double( double *buffer, int length ) {
       int i;
       union swapper {
         double double_number;
         char pos[8];
       }
       a, b;
     
       for( i=0 ; i < length ; i++) {
         a.double_number = *buffer;
         b.pos[0] = a.pos[7];
         b.pos[1] = a.pos[6];
         b.pos[2] = a.pos[5];
         b.pos[3] = a.pos[4];
         b.pos[4] = a.pos[3];
         b.pos[5] = a.pos[2];
         b.pos[6] = a.pos[1];
         b.pos[7] = a.pos[0];
         *buffer  = b.double_number;
         buffer++;
       }
       return;
  }

//! Routine for performing the bytes-swapping of single precision raw data files.
/*!
  This is a utility routine used by functions qload_unpackstrip and
  qsave_packstrip.
*/

  void block_swap( int *buffer, int length ) {
    int i;
    union swapper {
      int integer;
      char pos[4];
    }
    a, b;
    
    for( i=0 ; i < length ; i++)
    {
      a.integer = *buffer;
      b.pos[0] = a.pos[3];
      b.pos[1] = a.pos[2];
      b.pos[2] = a.pos[1];
      b.pos[3] = a.pos[0];
      *buffer  = b.integer;
      buffer++;
    }
    return;
  }

  
  // -------------------------------------------------------------------------------
  //! Normalises the specified row of an SU3 complex matrix:
  /*!
    This is a utility routine used by function qload_unpackstrip in
    the SU(3) reunitarisation.\n
    This subroutine was taken from su3_matrix.cpp by Balint Joo.
  */
  void qcdio_norm_row(Complex* data_, int row) {
     double sum;
     int index;
     sum=0.0;
     int i;
    
     for(i=0; i<COLORS; i++) {
       index=row+i*COLORS;
       sum+=data_[index].real()*data_[index].real() + 
         data_[index].imag()*data_[index].imag();
     }
    
     if ( fabs(sum) <=10e-50) {
      ERR.General("","qload_guage",
        "Unable to normalise row as it has length numerically = 0!\n");
     }
      
     sum=1.0/sqrt(sum);
     for(i=0; i<COLORS; i++) {
       index=row+(i*COLORS);
       data_[index]*=sum;
     }
    
  }
  
  // -------------------------------------------------------------------------------
  //! Orthogonalises two rows within an SU3 matrix:
  /*!
    This is a utility routine used by function qload_unpackstrip in
    the SU(3) reunitarisation.\n
    This subroutine taken from su3_matrix.cpp by Balint Joo.
  */
  void qcdio_orthog_rows( Complex* data_, int row1, int row2) {
     Complex dotprod1;
     
     Complex dotprod2;
    
     Complex factor;
     // dotprod  <row 1*, row2>
    
     dotprod1 = conj(data_[row1])*data_[row2]
               +conj(data_[row1+COLORS])*data_[row2+COLORS]
               +conj(data_[row1+(2*COLORS)]) * data_[row2+(2*COLORS)];
    
    
     dotprod2 = conj(data_[row1])*data_[row1]
              + conj(data_[row1+COLORS])*data_[row1+COLORS]
              + conj(data_[row1+(2*COLORS)])*data_[row1 + (2*COLORS)];
    
     factor = (1.0 / dotprod2.real())*dotprod1;
    
     data_[row2]=data_[row2]  - factor * data_[row1];
     data_[row2+COLORS] = data_[row2 + COLORS] - factor*data_[row1 + COLORS];
     data_[row2+(2*COLORS)] = data_[row2 + (2*COLORS)] 
                             - factor*data_[row1 + (2*COLORS)];
    
     // sum
  }  
  
  // -------------------------------------------------------------------------------
  //! Computes the cross-product of the 1st and 2nd rows in an SU3 matrix
  /*!
    This is a utility routine used by function qload_unpackstrip in
    the SU(3) reunitarisation.\n
    This subroutine taken from su3_matrix.cpp by Balint Joo:
              Nicked from GHMC.
  */
  void qcdio_cross_rows( Complex* data_, int row1, int row2, int row3) {
     Complex temp[3];
    
     // m(row1, 1)*m(row2,2) - m(row2,1)*m(row1,2)
     temp[0]=data_[row1 + COLORS]*data_[row2 + (2*COLORS)] 
            -data_[row2 + COLORS]*data_[row1 +(2*COLORS)];
    
     // m(row2, 0)*m(row1,2) - m(row1, 0)*m(row2,2)
     temp[1]=data_[row2]*data_[row1 + (2*COLORS)]
            -data_[row1]*data_[row2 + (2*COLORS)];
    
     // m(row1, 0)*m(row2,1) - m(row2, 0)*m(row1,1)
     temp[2]=data_[row1]*data_[row2 + COLORS]
            -data_[row2]*data_[row1+COLORS];
    
     
     data_[row3]=conj(temp[0]);
     data_[row3+COLORS]=conj(temp[1]);
     data_[row3+(2*COLORS)]=conj(temp[2]);
  }
  
  // -------------------------------------------------------------------------------
  /*!
    \param newval Set to 1 if the gauge field is to be normalized, or any
    other value if it is not.
  */
  void qcdio_set_normalize( int newval ) {
    qcdio_normalize_the_data = newval;
  }
  
  // -------------------------------------------------------------------------------
  //! A utility routine used by qloadsave_gauge.

  void qload_unpackstrip( unsigned char* mtxfilebuff, 
                          int filebufsize,
			  Matrix *siteaddr,
                          int prec, 
			  int swap, 
			  int transp ) {
    // Buffer definitions:
    Matrix mtxbuff;
    Float *mtxaddr; 
    
    // Swap the byte-order:
    if( swap == 1 ) {
      if( prec == sizeof(float) ) {
        block_swap( (int*)mtxfilebuff, filebufsize/sizeof(float) );
      } else {
        block_swap_double( (double*)mtxfilebuff, filebufsize/sizeof(double) );
      }
    }
    
    // loop over the mu direction,
    for( int mu = 0; mu < 4; mu++ ) {
	        
      // Construct a complete SU(3) matrix from the incomplete,
      // storing this matrix in the right spot in the internal
      // structure
      mtxaddr = (Float*) &mtxbuff;
      for( int col=0; col < COLORS; col++) {
       for( int row=0; row < COLORS-1; row++) {
         // For each row and column that we have, copy over the real
         // and the imaginary component:
         if( prec == sizeof(float) ) {
           mtxaddr[col*COLORS*2+row*2] = 
         	(Float)((float*)mtxfilebuff)[col*(COLORS-1)*2+row*2];
           mtxaddr[col*COLORS*2+row*2+1] = 
         	(Float)((float*)mtxfilebuff)[col*(COLORS-1)*2+row*2+1];
         } else {
           mtxaddr[col*COLORS*2+row*2] = 
         	(Float)((double*)mtxfilebuff)[col*(COLORS-1)*2+row*2];
           mtxaddr[col*COLORS*2+row*2+1] = 
         	(Float)((double*)mtxfilebuff)[col*(COLORS-1)*2+row*2+1]; 
         }
       }
      }
      
      // Reconstruct the missing column:
      if( qcdio_normalize_the_data == 1 ) {
        //These three are not needed as this is done before the 
	// matrix is written out anyway:
	//qcdio_norm_row( (Complex*)mtxaddr, 0 );
        //qcdio_orthog_rows( (Complex*)mtxaddr, 0, 1 );
        //qcdio_norm_row( (Complex*)mtxaddr, 1 );
	// This reconstructs the third row:
        qcdio_cross_rows( (Complex*)mtxaddr, 0, 1, 2 );
      }

      // Transpose the matrix if TRANSPOSE_THE_MATRICIES is 1:
      if( transp == 1 ) {
        siteaddr[mu].Trans(mtxbuff);
      } else {
        siteaddr[mu] = mtxbuff;
      }

      // Move file buffer on to next matrix:
      mtxfilebuff += filebufsize/4;
	
	
    // EOLoop over the mu direction	
    }

    return;
  }
  
  // --------------------------------------------------------------------------
  //! A utility routine used by qloadsave_gauge.

  void qsave_packstrip( unsigned char* mtxfilebuff, 
                        int filebufsize,
 			Matrix *siteaddr,
                        int prec, 
			int swap, 
			int transp ) {
    // Buffer definitions:
    Matrix mtxbuf;
    unsigned char* bufoff = mtxfilebuff;

    // Loop over mu:
    for( int mu = 0; mu < 4; mu++ ) {
    
      // For each SU3 matrix, copy 3x2 into buffer, byteswapping and
      // transposing if necessary, using the requested precision:
      if( transp == 1 ) {
        mtxbuf.Trans(siteaddr[mu]);
      } else {
        mtxbuf = siteaddr[mu];
      }
      
      for( int col = 0; col < COLORS; col++ ) {
        for( int row = 0; row < COLORS-1; row++ ) {
          if( prec == sizeof(float) ) {
            ((float*)bufoff)[col*(COLORS-1)*2+row*2] = 
                        ((Float*)&mtxbuf)[col*COLORS*2+row*2];
            ((float*)bufoff)[col*(COLORS-1)*2+row*2 + 1] = 
                        ((Float*)&mtxbuf)[col*COLORS*2+row*2 + 1];
          } else {
            ((double*)bufoff)[col*(COLORS-1)*2+row*2] = 
                        ((Float*)&mtxbuf)[col*COLORS*2+row*2];
            ((double*)bufoff)[col*(COLORS-1)*2+row*2 + 1] = 
                        ((Float*)&mtxbuf)[col*COLORS*2+row*2 + 1];
          }
        }
      }
  
      // Move to next spot in buf:
      bufoff += filebufsize/4;
     
    }
    
    // Byte-swap the buffer:
    if( swap == 1 ) {
      if( prec == sizeof(float) ) {
        block_swap( (int*) mtxfilebuff,  filebufsize );
      } else {
        block_swap_double( (double*) mtxfilebuff,  filebufsize );
      }
    }

    return;
  }

#ifndef USE_QMP
  // --------------------------------------------------------------------------
  //! A utility routine used by qloadsave_gauge.

  // if inout = 0, pump down to t,0,0,0 (for saving).
  // if inout = 1, pump up from t,0,0,0 (for loading).
  void qloadsave_pump_data( int* pos, unsigned char* buf, int byte_size, int inout) {
#if TARGET == QCDOC 
    // Size of buffer in appropriate units:
    // It is technically possible that the fact that the minimum unit of comms
    // is either 4 or 8 bytes could mean that the size of the transmitted data
    // is wrong.  However, there are enough factors of 4 and 2 etc to ensure
    // that in any reasonable case the number of bytes to send will map onto
    // an integer number of COMMS_DATASIZE units.
    int size = byte_size/COMMS_DATASIZE;
    SCUDirArg* datain = new SCUDirArg();
    SCUDirArg* dataout = new SCUDirArg();
    int datainset, dataoutset;
    
    //qprintf_allid(" -> [%i,%i,%i,%i] %x %i %i\n",pos[0],pos[1],pos[2],pos[3],
    //						buf, byte_size, inout);
    //sync();

    // Otherwise, pump it in/out:
    SCUXR posdir_SR, negdir_SR;
    if( inout == 0 ) {
      posdir_SR = SCU_REC;
      negdir_SR = SCU_SEND;
    } else {
      posdir_SR = SCU_SEND;
      negdir_SR = SCU_REC;
    }
    
    // Initialise the 'send/recv to be done' flags to no-op:
    datainset = 0; dataoutset = 0;
    
    // Pump along X-axis:
    if( GJP.TnodeCoor() == pos[3] &&
        GJP.YnodeCoor() == 0 && 
        GJP.ZnodeCoor() == 0 &&
	pos[0] > 0 ) {
      // Only a subset of processors are involved:
      if( GJP.XnodeCoor() <= pos[0] ) {
        // If we are not the 0'th node:
        if( GJP.XnodeCoor() > 0 ) {
	  // Send/recieve data in the -ve direction:
          datain->Init( (void*)buf, SCU_XM, negdir_SR, size );
	  datainset = 1;
        }
	// If we are anything but the last node:
        if( GJP.XnodeCoor() < pos[0] ) {
	  // Send/recieve in the +ve dirn:
          dataout->Init( (void*)buf, SCU_XP, posdir_SR, size );
	  dataoutset = 1;
        } 
      }
    }
    
    // Pump along Y-axis:
    if( GJP.TnodeCoor() == pos[3] &&
        GJP.XnodeCoor() == pos[0] && 
        GJP.ZnodeCoor() == 0 &&
	pos[1] > 0 ) {
          // Only a subset of processors are involved:
      if( GJP.YnodeCoor() <= pos[1] ) {
        // If we are not the 0'th node:
        if( GJP.YnodeCoor() > 0 ) {
	  // Send/recieve data in the -ve direction:
          datain->Init( (void*)buf, SCU_YM, negdir_SR, size );
	  datainset = 1;
        }
	// If we are anything but the last node:
        if( GJP.YnodeCoor() < pos[1] ) {
	  // Send/recieve in the +ve dirn:
          dataout->Init( (void*)buf, SCU_YP, posdir_SR, size );
	  dataoutset = 1;
        } 
      }
    }

    // Pump along Z-axis:
    if( GJP.TnodeCoor() == pos[3] &&
        GJP.YnodeCoor() == pos[1] &&  
	GJP.XnodeCoor() == pos[0] &&
	pos[2] > 0 ) {
         // Only a subset of processors are involved:
      if( GJP.ZnodeCoor() <= pos[2] ) {
        // If we are not the 0'th node:
        if( GJP.ZnodeCoor() > 0 ) {
	  // Send/recieve data in the -ve direction:
          datain->Init( (void*)buf, SCU_ZM, negdir_SR, size );
	  datainset = 1;
        }
	// If we are anything but the last node:
        if( GJP.ZnodeCoor() < pos[2] ) {
	  // Send/recieve in the +ve dirn:
          dataout->Init( (void*)buf, SCU_ZP, posdir_SR, size );
	  dataoutset = 1;
	}
      }
    }

    // Wait for all comms in this transaction to complete:
    // Make sure we recieve before we send!
    if( inout == 0 ) {
      if( dataoutset == 1 ) {
        SCUTrans(dataout);
        SCUTransComplete();
      }
      if( datainset == 1 ) {
        SCUTrans(datain);
        SCUTransComplete();
      }
    } else {
      if( datainset == 1 ) {
        SCUTrans(datain);
        SCUTransComplete();
      }
      if( dataoutset == 1 ) {
        SCUTrans(dataout);
        SCUTransComplete();
      }
    }
    //sync();

#endif
    
    return;
  }
#endif

  // --------------------------------------------------------------------------
  //! A utility routine used by qloadsave_gauge.

  void qload_parameters( char* fprefix, Lattice& lat ) {
    int fprelen;
    char dummystr[500], *parfname, *ftoload;
    FILE* fp;
    
    //qprintf_allid("loadpar 1\n");
    
    /* Parse the file fprefix.par and checks that the lattice sizes match up */
    // Find out how long the filename prefix is:
    fprelen = 0; while( fprefix[fprelen] != 0 ) fprelen++;

    // Store the file prefix, construct the parameter filename:
    // This is fairly ugly because I don't want to rely on "external" libraries.   
    char *qcdio_fprefix = new char[fprelen];
    parfname = new char[fprelen+4];
    ftoload = new char[fprelen+4];
    for( int i = 0; i <= fprelen; i++ ) {
      qcdio_fprefix[i] = fprefix[i];
      parfname[i] = fprefix[i];
    }
    parfname[fprelen] = '.';
    parfname[fprelen+1] = 'p';
    parfname[fprelen+2] = 'a';
    parfname[fprelen+3] = 'r';
    parfname[fprelen+4] = 0;
    
    // Try to open the parameter file:
    fp = fopen(parfname, "r");
    if( fp == NULL ) {
      ERR.General("","qload_guage","Could not open file %s!\n", parfname);
    }

      double qcdio_beta;
    // Parse the parameters from the file:
    fscanf(fp,"%s %lf\n", dummystr, &qcdio_beta );
    fscanf(fp,"%s %i\n", dummystr, &(qcdio_latt[0]) );
    fscanf(fp,"%s %i\n", dummystr, &(qcdio_latt[1]) );
    fscanf(fp,"%s %i\n", dummystr, &(qcdio_latt[2]) );
    fscanf(fp,"%s %i\n", dummystr, &(qcdio_latt[3]) );
    fclose(fp);
 
    // Check that the lattice size matches up, exiting if they don't.
    if( ! ( GJP.XnodeSites()*GJP.Xnodes() == qcdio_latt[0] &&
            GJP.YnodeSites()*GJP.Ynodes() == qcdio_latt[1] &&
	    GJP.ZnodeSites()*GJP.Znodes() == qcdio_latt[2] &&
	    GJP.TnodeSites()*GJP.Tnodes() == qcdio_latt[3] &&
	    GJP.SnodeSites() == 0 && GJP.Snodes() == 1 
	    ) ) {
      ERR.General("","qload_guage",
      "Gauge configuration %s does not match the simulation dimensions!\n", 
      qcdio_fprefix );
    }
    
    // Verbose output hackin in for now:
#ifdef VERBOSE
    printf("beta = %.20f\n", qcdio_beta );
    printf("system dimensions = [%i,%i,%i,%i]\n",
         qcdio_latt[0], qcdio_latt[1], qcdio_latt[2], qcdio_latt[3] );
#endif
    //qprintf_allid("loadpar 2\n");
    
    return;
  }

  // --------------------------------------------------------------------------
  //! A utility routine used by qloadsave_gauge.

  void qsave_parameters( char* fprefix, Lattice& lat ) {
    // Grab parameters from CPS and put into our own variables.
    qcdio_latt[0] = GJP.XnodeSites()*GJP.Xnodes();
    qcdio_latt[1] = GJP.YnodeSites()*GJP.Ynodes();
    qcdio_latt[2] = GJP.ZnodeSites()*GJP.Znodes();
    qcdio_latt[3] = GJP.TnodeSites()*GJP.Tnodes();

    // Save the information to a file (ONE NODE ONLY):
    if( GJP.TnodeCoor() == 0 &&
        GJP.XnodeCoor() == 0 &&
	GJP.YnodeCoor() == 0 &&
	GJP.ZnodeCoor() == 0 ) {
      // FIXME:  This should be implemented, but perhaps as XML using the Schema	
    }
    return; 
  }

  // --------------------------------------------------------------------------
  //! A utility routine used by qload_gauge and qsave_gauge..

  // mode = 0 does a load
  // mode != 0 does a save
  void qloadsave_gauge( int mode, char* fprefix, Lattice& lat, 
             int prec,
	     int swap,
	     int transp ) {

    char ftoload[200];
    FILE *fp;
    
    //qprintf_allid("loadsavegauge 1\n");

#ifdef PARALLEL
    // Syncronise all processors before beginning to load.
    //sync();
#endif

    //qprintf_allid("loadsavegauge 2\n");

    // Check the parameters prec, swap and transp:
    if( prec != sizeof(float) && prec != sizeof(double) ) {
      ERR.General("", "qload_gauge",
        "Gauge config precision must be sizeof(float) or sizeof(double)!\n");
    }
    if( swap != 1 && swap != 0 ) {
      ERR.General("", "qload_gauge",
        "Byte-order swap flag must be 0 (no swap) or 1 (swap)!\n");
    }
    if( transp != 1 && transp != 0 ) {
      ERR.General("", "qload_gauge",
        "SU3 transpose flag must be 0 (no transpose) or 1 (transpose)!\n");
    }
    
    //qprintf_allid("loadsavegauge 3\n");

    // Load/save the parameter file:
    if( mode == 0 ) { 
      qload_parameters( fprefix, lat );
    } else {
      qsave_parameters( fprefix, lat );
    }

    //qprintf_allid("loadsavegauge 4\n");

    // Loop over the timeslices corresponding to this node
    int tslice;
    for( tslice = GJP.TnodeSites()*GJP.TnodeCoor(); 
         tslice < GJP.TnodeSites()*(GJP.TnodeCoor()+1); tslice++ ) {
	        
      // only the x=y=z=0 nodes open the file:		
      if( GJP.XnodeCoor() == 0 &&
	  GJP.YnodeCoor() == 0 &&
	  GJP.ZnodeCoor() == 0 ) {

    //qprintf_allid("loadsavegauge 4.1 open the file\n");
    
        // Construct the timeslice filename and open the file
        sprintf(ftoload, "%sT%02d", fprefix, tslice ); 
        if( mode == 0 ) {
	  fp = fopen( ftoload, "rb" );
	} else {
	  fp = fopen( ftoload, "wb" );
	}
        if( fp == NULL ) {
          ERR.General("","qload_guage", 
	    "Could not open data file %s!\n", ftoload );
        }
      }
      
      int z, y, x, mu, row, col;
      int localpos[4], pe_pos[4], bytes_read, 
               bytes_to_read = prec*COLORS*(COLORS-1)*2*4;
      Matrix *siteaddr;
      unsigned char *mtxfilebuff;
      mtxfilebuff = new unsigned char[prec*COLORS*COLORS*2*4];
      
    //qprintf_allid("loadsavegauge 5\n");

      // loop over the z direction
      for( z = 0; z < qcdio_latt[2]; z++ ) {
      
    //qprintf_allid("loadsavegauge 5.1 every z\n");
        // loop over the y direction
        for( y = 0; y < qcdio_latt[1]; y++ ) {

    //qprintf_allid("loadsavegauge 5.2 every y\n");
	  // loop over the x direction
          for( x = 0; x < qcdio_latt[0]; x++ ) {

    //qprintf_allid("loadsavegauge 5.3 every x\n");
	    // identify the internal site corresponding to the file-site
	    localpos[0] = x - GJP.XnodeSites()*GJP.XnodeCoor();
	    localpos[1] = y - GJP.YnodeSites()*GJP.YnodeCoor();
	    localpos[2] = z - GJP.ZnodeSites()*GJP.ZnodeCoor(); 
	    localpos[3] = tslice - GJP.TnodeSites()*GJP.TnodeCoor();
 
            // Determine where in the PE array the current spot lies:
	    pe_pos[0] = x/GJP.XnodeSites();
	    pe_pos[1] = y/GJP.YnodeSites();
	    pe_pos[2] = z/GJP.ZnodeSites();
	    pe_pos[3] = tslice/GJP.TnodeSites();
 
	    // If we are loading:
	    if( mode == 0 ) {
	       // if this is a x,y,z=0 node:
	      if( GJP.XnodeCoor() == 0 &&
		  GJP.YnodeCoor() == 0 &&
		  GJP.ZnodeCoor() == 0 ) {
	        // Read the incomplete Matrix, cols[3] then rows[2] then [real,img]
                // into a buffer
                bytes_read=fread((char*)mtxfilebuff, sizeof(unsigned char), 
		                  bytes_to_read, fp);
                    
                // Make sure we read all the data
                if( bytes_read != bytes_to_read ) {
                  ERR.General("","qload_guage",
  		    "fread: Could not read the whole of %s!\n",ftoload);
                }
		  
	      }
    
    //qprintf_allid("loadsavegauge 5.5 pump the data\n");
	       
	      // Pump the data out to the right PE:
 	      qloadsave_pump_data( pe_pos, mtxfilebuff, bytes_to_read , 1 ); 
	       
	      // If this is the destination of the data:
	      if( GJP.XnodeCoor() == pe_pos[0] &&
	          GJP.YnodeCoor() == pe_pos[1] &&
	          GJP.ZnodeCoor() == pe_pos[2] ) {
	        // determine destination in the local lattice array:
 	        siteaddr = lat.GaugeField() + lat.GsiteOffset( (int*)localpos );
	        // Get the data from the file buffer and put it into the right
                // spot in the Lattice:
                qload_unpackstrip( mtxfilebuff, bytes_to_read, siteaddr,prec,swap,transp);
	      } 
            // Else attempt to save the data:
	    } else {
	      // IF this PE owns the data:
	      if( GJP.XnodeCoor() == pe_pos[0] &&
	          GJP.YnodeCoor() == pe_pos[1] &&
		  GJP.ZnodeCoor() == pe_pos[2] ) {
	        // determine destination in the local lattice array:
 	        siteaddr = lat.GaugeField() + lat.GsiteOffset( (int*)localpos );
	        // Get the data from the Lattice object and pack it into the
	        // buffer.
	        qsave_packstrip( mtxfilebuff, bytes_to_read, siteaddr,prec,swap,transp);
	      }	  
               
	      // Pump the data down to the x,y,z = 0 nodes:
	      qloadsave_pump_data( pe_pos, mtxfilebuff, bytes_to_read, 0 );

              // if recieved, the x,y,z nodes=0 send the data to a file:
	      if( GJP.XnodeCoor() == 0 &&
	          GJP.YnodeCoor() == 0 &&
		  GJP.ZnodeCoor() == 0 ) {
		// Write the data to the file:
                bytes_read=fwrite((char*)mtxfilebuff, sizeof(unsigned char), 
		                  bytes_to_read, fp);
                    
                // Make sure we read all the data
                if( bytes_read != bytes_to_read ) {
                  ERR.General("","qload_guage",
  		    "fread: Could not write the whole of %s!\n",ftoload);
                }

	      } 
		  
            }
	      
	  // EOLoop over x
  	  }
        
	// EOLoop over y
	}
	
      // EOLoop over z
      }

    // EOLoop over timeslice files
    }
    
    //qprintf_allid("loadsavegauge 6\n");

#ifdef PARALLEL
    // Syncronise all processors before returning control.
    //sync();
#endif
    
    //qprintf_allid("loadsavegauge 7\n");

    // If all went well, return zero:
    return;

  }

  
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  /*! 
    Loads the gauge configuration from a set of UKQCD format
    files.
    \param fprefix The prefix of the gauge configuration filenames
    \param lat The Lattice object into which to load the gauge configuration
    \param prec The precision at which the gauge configuration is stored in
    the file, either \c sizeof(float) or \c sizeof(double)
    \param swap 1 if the data need to be byte-swapped when they are read in,
    0 otherwise,
    \param transp 1 if the gauge field matrices need to be transposed when they are read in,
  */
  void qload_gauge( char* fprefix, Lattice& lat, 
             int prec,
	     int swap,
	     int transp ) {
    qloadsave_gauge( 0, fprefix, lat, prec, swap, transp );
    return;
  }  

  /*! 
    Saves the gauge configuration into a set of UKQCD format
    files.
    Lattice object lat
    \param fprefix The prefix of the gauge configuration filenames
    \param lat The Lattice object from which to write the gauge configuration
    \param prec The precision at which the gauge configuration is stored in
    the file, either \c sizeof(float) or \c sizeof(double)
    \param swap 1 if the data need to be byte-swapped when they are written out,
    0 otherwise,
    \param transp 1 if the gauge field matrices need to be transposed when they are written out,
  */
  void qsave_gauge( char* fprefix, Lattice& lat, 
             int prec, 
	     int swap, 
	     int transp ) {
    qloadsave_gauge( 1, fprefix, lat, prec, swap, transp ); 
    return;
  }


// --------------------------------------------------------------------------


CPS_END_NAMESPACE
