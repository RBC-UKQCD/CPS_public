#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/hadron_prop_s.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: hadron_prop_s.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.10  2002/03/11 22:25:42  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.7.2.1  2002/03/08 16:34:59  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.7  2001/08/16 12:53:55  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.6  2001/08/16 10:49:40  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:00:45  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:09  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:11:29  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: hadron_prop_s.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/hadron_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// hadron_prop_s.C
//debug
CPS_END_NAMESPACE
#include <stdio.h>
#include <alg/hadron_prop_s.h>
#include <util/rcomplex.h>
#include <comms/glb.h>
#include <alg/myenum.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

const int MAX_LEN = 1023;     // max. number of IFloats that can be transfered

char HadronPropS::cname[] = "HadronPropS";

static int nx[4];
static int Nx[4];
static int vsize;
static BndCndType bc[4]; 

//------------------------------------------------------------------
// protected function used in all derived classes
//------------------------------------------------------------------
int HadronPropS::X_OFFSET(const int *x)
{ return lat.FsiteOffsetChkb(x) * VECT_LEN
           + ((x[0]+x[1]+x[2]+x[3]) & 1 ? vsize : 0); }


//------------------------------------------------------------------
// non-pure virtual function, override in NLocalProp class
//------------------------------------------------------------------
void HadronPropS::getNeighbors(int t) {}

//------------------------------------------------------------------
//  CTOR & DTOR
//------------------------------------------------------------------
HadronPropS::HadronPropS( Lattice &lattice, int num, 
		        int direction, int srcslice, int incr )
: lat(lattice),
  n_props(num), 
  dir(direction),
  slice(srcslice),
  stride(incr)
{
    char *fname = "HadronPropS(Lattice &, int, int)";
    nx[0] = GJP.XnodeSites();
    nx[1] = GJP.YnodeSites();
    nx[2] = GJP.ZnodeSites();
    nx[3] = GJP.TnodeSites();

    Nx[0] = GJP.Xnodes();
    Nx[1] = GJP.Ynodes();
    Nx[2] = GJP.Znodes();
    Nx[3] = GJP.Tnodes();

    vsize = GJP.VolNodeSites() * VECT_LEN / 2;

    bc[0] = GJP.Xbc();
    bc[1] = GJP.Ybc(); 
    bc[2] = GJP.Zbc(); 
    bc[3] = GJP.Tbc(); 

    prop = (Complex *)smalloc(nx[dir] * num * sizeof(Complex));
    if(prop == 0)
      ERR.Pointer(cname,fname, "prop");
    VRB.Smalloc(cname,fname, "prop", prop, nx[dir]*num*sizeof(Complex));
}

HadronPropS::~HadronPropS() 
{ 
  char *fname = "~HadronPropS()";
  VRB.Sfree(cname,fname, "prop",prop);
  sfree(prop); 
}

BndCndType HadronPropS::bcd()
{
  return bc[dir];
}

int HadronPropS::propLen() const 
{ return CR * n_props * nx[dir]; }

int HadronPropS::propLenTotal() const 
{ return CR * n_props * nx[dir] * Nx[dir]; }

int HadronPropS::propLenFold() const 
{ return n_props * (nx[dir] * Nx[dir]/2 + 1); }
//------------------------------------------------------------------
//  Compute the contribution from each site and sum over a
//  hyperplane. For nonlocal operators scu transfers are needed
//  for each hyperplane.
//------------------------------------------------------------------

void HadronPropS::getHadronPropS()
{
  int s[4];
  for (int len = 0; len < n_props * nx[dir]; len++)
    prop[len] = 0;

  Complex *currt_p = prop; 

  //---------------------------------------------------------------
  // sum on slice orthogonal to dir (not necessary time)
  // i, j, k are the rest of directions
  //---------------------------------------------------------------

  int i = (dir+1)%4;
  int j = (dir+2)%4;
  int k = (dir+3)%4;

  //---------------------------------------------------------------
  //  For each slice s[dir], sum over the hyperplane on this node 
  //---------------------------------------------------------------

  for (s[dir] = 0; s[dir] < nx[dir]; s[dir]++) { 
    getNeighbors(s[dir]);
    for (s[i] = 0; s[i] < nx[i]; s[i] += stride) {    
      for (s[j] = 0; s[j] < nx[j]; s[j] += stride) { 
        for (s[k] = 0; s[k] < nx[k]; s[k] += stride) {

	  localVal(currt_p, s);   
        }
      }
    }
    currt_p += n_props;	
  }

  //---------------------------------------------------------------
  // sum over all sites on all slices 
  //---------------------------------------------------------------

  // added if ( ) { } else {} by manke
  // this is to protect this routine against large blocks > 1023
  // which are not dealt with properly in slice sum

  if (propLen() > MAX_LEN){
    // if vector is too long  split slice sum into seperate pieces
    int ll_max = propLen()/MAX_LEN + 1;
    int shift  = 0;

    for (int ll=0; ll < ll_max; ll++){
      Float *temp_prop = ((Float *)prop)+shift;
      int length= (ll==ll_max-1) ? (propLen()-shift) : MAX_LEN;

      slice_sum((Float *)temp_prop, length, dir);

      shift+=MAX_LEN;
    }
  } else {

    slice_sum((Float *)prop, propLen(), dir);
  }

}

void HadronPropS::download_prop(HadronType type, Float *buf)
{
  const Float *tmp = propPtr();

#ifdef PARALLEL

  // added by manke
  if (propLen() > MAX_LEN){
    collect_large_prop(type, buf, tmp, n_props * CR, dir, slice);
  } else {
    collect_prop(type, buf, tmp, n_props * CR, dir, slice);
  }

#else 

  collect_prop(type, buf, tmp, n_props * CR, dir, slice);

#endif


  /*
  if (type==SMOMMESON){
    printf("download_prop: after collect_prop\n");
 
    IFloat checksum=0.;
    for (int ii = 0; ii < 400; ii++) checksum+=buf[ii];
    printf("download_prop: checksum =  %e \n",checksum);
  }
  */
}


#ifdef PARALLEL
//---------------------------------------------------------------------
// download data from all nodes in the propagator direction to
// node 0 (origin of the lattice)
//---------------------------------------------------------------------

const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };

static int isNodeOrigin();
static void storeData(IFloat *buf, const IFloat *term, int len);

void HadronPropS::collect_prop(HadronType type, Float *sum_buf, 
			      const Float *IFloat_p, int unit, 
			      int direction, int t)
{
  char *fname ="collect_prop(HadronType, Float, Float, *IFloat_p, int, int,int)";

  int blcklength = unit*nx[direction];

  IFloat *receive_buf = (IFloat *)smalloc(blcklength*Nx[direction]*sizeof(IFloat));

  // added by manke for check

  if(receive_buf == 0)
    ERR.Pointer(cname,fname, "receive_buf");
  VRB.Smalloc(cname,fname, "receive_buf", receive_buf, blcklength*Nx[direction]*
sizeof(IFloat));

  // end added

  IFloat *receive_buf_p = receive_buf;
  IFloat *transmit_buf_p;
  
  int itmp;		  
	// loop index with 1<= itmp < Nx[i] 

  //--------------------------------------------------------------
  // store data on node 0 first 
  //--------------------------------------------------------------
  if (isNodeOrigin()) {
    storeData(receive_buf_p, (IFloat *)IFloat_p, blcklength);
    receive_buf_p += blcklength;
  }
  //--------------------------------------------------------------
  // address of buffer of to be sent first (data on this node) 
  //--------------------------------------------------------------
  transmit_buf_p = (IFloat *)IFloat_p; 

  SCUDirArg send(transmit_buf_p, neg_dir[direction], SCU_SEND, blcklength);
  SCUDirArg recv(receive_buf_p, pos_dir[direction], SCU_REC, blcklength);

  //--------------------------------------------------------------
  // tranmit & receive Nx[direction] - 1 times in snd_dir[direction] direction
  //--------------------------------------------------------------
  for ( itmp = 1; itmp < Nx[direction]; itmp++) {


	 //-----------------------------------------------------------
         // do SCU transfers
	 //-----------------------------------------------------------
         SCUTrans(&send);
         SCUTrans(&recv);
	 SCUTransComplete();

	 //-----------------------------------------------------------
         // the received data will be sent out     
	 // the free buffer will be used to receive      
	 //-----------------------------------------------------------
	 send.Addr(transmit_buf_p = receive_buf_p);
	 recv.Addr(receive_buf_p += blcklength );
  }





  // everything is in receive_buf now

  if(isNodeOrigin()) {
    int i;
    for(i = unit * t; i < unit * nx[direction] * Nx[direction]; i++)
	*sum_buf++ += receive_buf[i];
    if( (type != SMESON) && (type != SMOMMESON) && (bc[dir] == BND_CND_APRD) ) {
      for(i = 0; i < unit * t; i++)
	*sum_buf++ -= receive_buf[i];
    }
    else {
      for(i = 0; i < unit * t; i++)
	*sum_buf++ += receive_buf[i];
    }
  }

  sfree(receive_buf);
}


// added by manke to work properly
// if "blcklength" is larger than MAX_LEN=1023 --> split the transfer into 
// several pieces

void HadronPropS::collect_large_prop(HadronType type, Float *sum_buf,
                              const Float *IFloat_p, int unit,
                              int direction, int t)
{
  char *fname = "collect_large_prop(HadronType, Float, Float, *IFloat_p,  int, int, int)";

  int blcklength = unit*nx[direction];

  IFloat *receive_buf = (IFloat *)smalloc(blcklength*Nx[direction]*sizeof(IFloat));

  // added by manke for check
  if(receive_buf == 0)
    ERR.Pointer(cname,fname, "receive_buf");
  VRB.Smalloc(cname,fname, "receive_buf", receive_buf, blcklength*Nx[direction]*
sizeof(IFloat));
  // end added

  IFloat *receive_buf_p = receive_buf;
  IFloat *transmit_buf_p;

  int ll_max = 1;                        // transfer everything at once
  int length = blcklength;               // default length of block
  int shift  = 0;                        // shift in receive_buf

  if (blcklength > MAX_LEN) ll_max = blcklength/MAX_LEN + 1;

  /*
  printf("collect_large_prop: IFloat_p \n");
  for (int ii=0; ii< blcklength; ii+=8)
    printf("%d %e %e  %e %e   %e %e   %e %e   \n",ii, IFloat_p[ii],IFloat_p[ii+1],
 IFloat_p[ii+2],IFloat_p[ii+3], IFloat_p[ii+4],IFloat_p[ii+5], IFloat_p[ii+6],IFloat_p
[ii+7]);
  printf("blcklenth = %d --> run transfer %d times \n",blcklength,ll_max);
  */

  for (int ll=0; ll < ll_max; ll++){
    length= (ll==ll_max-1) ? (blcklength-shift) : MAX_LEN;  
		//  block length for transfer

    // printf("%d transfer:  shift = %d with length = %d \n",ll,shift,length);

    //---------------------------------------------------------------------

    // address of buffer of to be sent/received first
    // send/receive "length" IFloats starting from IFloat_p + shift 
		// on every processor

    //---------------------------------------------------------------------

    transmit_buf_p = (IFloat *)(IFloat_p + shift);
    receive_buf_p  = receive_buf + shift;  
		// initial offset if blcklength > MAX_LEN


    //--------------------------------------------------------------
    // store data on node 0 first
    // NOTICE: if MAX_LEN = length < blcklength then
    // receive_buf[length] ... receive_buf[blcklength] is still undefined
    // for ll=0. It will be filled in later.
    //--------------------------------------------------------------
    if (isNodeOrigin()) {
      storeData(receive_buf_p, transmit_buf_p , length);
      receive_buf_p += blcklength;
    }

    SCUDirArg send(transmit_buf_p, neg_dir[direction], SCU_SEND, length);
    SCUDirArg recv(receive_buf_p, pos_dir[direction], SCU_REC, length);

    //--------------------------------------------------------------
    // tranmit & receive Nx[direction] - 1 times in snd_dir[direction] direction
    //--------------------------------------------------------------

    for ( int itmp = 1; itmp < Nx[direction]; itmp++) {


      //-----------------------------------------------------------
      // do SCU transfers
      //-----------------------------------------------------------
      SCUTrans(&send);
      SCUTrans(&recv);
      SCUTransComplete();

      //-----------------------------------------------------------
      // the received data will be transmitted again to the next processor
      // the buffer for receive will be shifted by blcklenght to receive
      // the next data
      // the data from the chain of Nx[direction] processors in the
      // direction "direction", will be accumulated as shifted
      // entries in the (global) receive buffer
      // every processor has a copy of this, but only for
      // the node at the origin the order in receive_buf corresponds
      // to the true temporal ordering

      // PROCESSORS ARE SEPARATED BY blcklength in the TEMPORAL DIRECTION
      // THE ACTUAL NUMBER OF TRANSFERED FLOATS may be 
      // smaller length < blcklength
      //-----------------------------------------------------------

      send.Addr(transmit_buf_p = receive_buf_p);
      recv.Addr(receive_buf_p += blcklength );

    }

    shift+=MAX_LEN; // increase shift in receive_buf for next set of transfers

  } // end loop over blocks of all sizes
 
  // everything is in receive_buf now
  // let the origin do the final analysis

  if(isNodeOrigin()) {
    int i;
    for(i = unit * t; i < unit * nx[direction] * Nx[direction]; i++)
        *sum_buf++ += receive_buf[i];
    if( (type != SMESON) && (type != SMOMMESON) && (bc[dir] == BND_CND_APRD) ) {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ -= receive_buf[i];
    }
    else {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ += receive_buf[i];
    }
  }

  sfree(receive_buf);
}



int isNodeOrigin() 
{
  return (GJP.XnodeCoor() == 0) && (GJP.YnodeCoor() == 0) && 
	 (GJP.ZnodeCoor() == 0) && (GJP.TnodeCoor() == 0);
}

void storeData(IFloat *buf, const IFloat *term, int len)
{
  for(int i = 0; i < len; i++) {
    *buf++ = *term++;
  }
}

#else

void HadronPropS::collect_prop(HadronType type, Float *sum_buf, 
	      const Float *IFloat_p, int unit, int direction, int t)
{
    int i;
    for(i = unit * t; i < unit * nx[direction] * Nx[direction]; i++)
        *sum_buf++ += IFloat_p[i];
    if( (type != SMESON) && (type != SMOMMESON) && (bc[dir] == BND_CND_APRD) ) {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ -= IFloat_p[i];
    }
    else {
      for(i = 0; i < unit * t; i++)
        *sum_buf++ += IFloat_p[i];
    }
}

#endif


CPS_END_NAMESPACE
