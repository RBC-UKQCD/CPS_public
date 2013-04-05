#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  LinkBuffer class methods and some Lattice class methods involving
  the link buffer.

  $Id $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_base/link_buffer.C,v 1.9 2013-04-05 17:51:14 chulwoo Exp $
//  $Id: link_buffer.C,v 1.9 2013-04-05 17:51:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/lattice_base/link_buffer.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/verbose.h>
//#include <comms/nga_reg.h>
#include <comms/scu.h>
#include <comms/cbuf.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/link_buffer.h>
#include <util/list.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE

enum {MATRIX_SIZE = 18};

//! The basic object in the link buffer.
/*!
  This is manipulated in the internal workings of the Linkbuffer class.
*/
struct LinkEntry{
  list_head hash_entry;
  //to put this LinkEntry in the hash table
  
  list_head flush_entry; 
  //to put this LinkEntry in the free_list or the flush_list
  
  int d_node_id;
  //the node that this link belong to
  
  Matrix link; 
  // the link
};

//  CRAM temp buffer
#ifdef _TARTAN
Matrix & mat1 = * ((Matrix*) CRAM_SCRATCH_ADDR);
Matrix & mat2 = * (&mat1+1);
Matrix & mat3 = * (&mat2+1);
Matrix & mat4 = * (&mat3+1);
Matrix & mat5 = * (&mat4+1);

const Matrix *new_mp = &mat1;
Matrix *result1_mp = &mat2;
Matrix *result_mp = &mat3;
Matrix *acumulate_mp =&mat4;
#else

Matrix CRAM_SCRATCH[5] ;

Matrix & mat1 = CRAM_SCRATCH[0] ;
Matrix & mat2 = CRAM_SCRATCH[1] ;
Matrix & mat3 = CRAM_SCRATCH[2] ;
Matrix & mat4 = CRAM_SCRATCH[3] ;
Matrix & mat5 = CRAM_SCRATCH[4] ;

const Matrix *new_mp = &mat1;
Matrix *result1_mp = &mat2;
Matrix *result_mp = &mat3;
Matrix *acumulate_mp =&mat4;
#endif

//! Prints the addresses of the items in a linked list 
void print_list(list_head * lp){ 
  list_head * lh_p = lp->next;
  printf("\n\nlist at %x:\n", (lp));
  printf("%x %x %x\n", (lp), (lp)->prev, (lp)->next);
  while(lh_p != lp){
    printf("%x %x %x\n", (lh_p), (lh_p)->prev, (lh_p)->next);
    lh_p = lh_p->next;
  }
}

static int cnt;
static int hit_cnt;

//! The number of items in a linked list
int list_len(list_head * l){
  list_head * l_iter = l->next; 
   int i=0;
  while(l_iter != l){ l_iter=l_iter->next; i++;}
  return i;
}
//int
//LinkBuffer::check_lists(){
//  hash_sz=0;
//  for(int i=0;i< GJP.VolNodeSites()*4;i++){
//    hash_sz += list_len(&hash_tab[i]);
//  }
//  int flush_size = list_len(&flush_list);
//  int free_size = list_len(&free_list);
//  if(hash_sz != flush_size || hash_sz + free_size!=buf_sz){
//    printf("hash_sz %d flush_size %d free_size %d\n", hash_sze,
//            flush_size, free_size);
//    return 1;
//  }
//  return 0;
//}

LinkBuffer::LinkBuffer(Lattice &_lat, int _buf_sz):lat(_lat), buf_sz(_buf_sz){
  cname = "LinkBuffer";
  char * fname = "LinkBuffer";
  VRB.Func(cname, fname);

  cnt=0; hit_cnt=0;
  
  tab_size = GJP.VolNodeSites()*4;
  hash_tab = (list_head *)smalloc(tab_size*sizeof(list_head)); 
  if(hash_tab==0) ERR.Pointer(cname, fname, "hash_tab");
  
  int i;
  for(i=0; i<tab_size; i++) {
    INIT_LIST_HEAD(&hash_tab[i]); 
    //print_list(&hash_tab[i]);
  }
  
  INIT_LIST_HEAD(&flush_list); 
  INIT_LIST_HEAD(&free_list); 
  
  link_entry_buf = (LinkEntry*)smalloc(sizeof (LinkEntry) * buf_sz); 
  if(link_entry_buf==0) ERR.Pointer(cname, fname, "link_entry_buf");
  
  LinkEntry * tmp = link_entry_buf; 
  for(i =0 ; i < buf_sz ; i++, tmp++){
    list_add(&(tmp->hash_entry), &free_list); 
    // hash_tab and free_list use the same entry in the LinkEntry
    // put everything on the freelist initially
    // INIT_LIST_HEAD(&(tmp->flush_entry)); 
  }
  site_stride[0]=4; 
  node_stride[0]=1;
  site_size[0]=GJP.Nodes(0)*GJP.NodeSites(0);
  for(i=1;i<4;i++){
   site_stride[i]=site_stride[i-1]*GJP.NodeSites(i-1);
   site_size[i]=GJP.Nodes(i)*GJP.NodeSites(i);
   node_stride[i]=node_stride[i-1]*GJP.Nodes(i-1);
  }
}

/*!
  Prints to \c stdout the number of times a lenk has been fetched from the
  buffer and the the number of times that link was on this node.
*/
LinkBuffer::~LinkBuffer(){ 
  //char * fname = "~LinkBuffer()";
  //VRB.Func(cname, fname);
  sfree(link_entry_buf); sfree(hash_tab);
  printf("GetBufferedLink called %d times with %d hits\n", cnt, hit_cnt);
}

//-----------------------------------------------------------------------
//LinkBuffer::hash() 
// given the lattice coordinates of a link and its direction.x[0-4] can be
// any integer, mu can be 0, 1, 2 ,3.
// the returned value is between 0 and VolNodeSize*4
//-----------------------------------------------------------------------
int LinkBuffer:: hash(const int *x, int mu){
  //char * fname ="hash";
  int result=0;

  for(int i=0; i<4; i++){
    int local_x = x[i];
    while(local_x<0) local_x += GJP.NodeSites(i); 
    while(local_x>=GJP.NodeSites(i)) local_x -= GJP.NodeSites(i);
    result += local_x*site_stride[i];
  } 
  return (result + mu);
}

//----------------------------------------------------------------------
// LinkBuffer::GetNodeId()
// given the lattice coordinates of a link calculate the coordinate of the 
// node that this link resides on. These coordinates are further converted
// into an integer.
//----------------------------------------------------------------------
int LinkBuffer::GetNodeId(const int *x){
  //char * fname = "GetNodeId";

  int result=0;
  for(int i = 0 ;i< 4; i++) {
    int tmp_size =site_size[i]; 
    int tmp_x = x[i];
    while(tmp_x<0) tmp_x += tmp_size;
    while(tmp_x>=tmp_size) tmp_x -= tmp_size;
    result += node_stride[i] * (tmp_x/GJP.NodeSites(i));
  }
  //VRB.Flow(cname, fname, "node_id=%d\n", result);
  return result;
}

//----------------------------------------------------------------------
/*!
  Looks for the link \a U_mu(x) in the buffer. If it is not there it is
  brought in.
  \param x The lattice coordinates of the link.
  \param mu The direction index of the link.
  \return The link \a U_mu(x).
*/
//----------------------------------------------------------------------
Matrix * LinkBuffer::GetBufferedLink(const int *x, int mu){
  char * fname = "GetBufferedLink()";
  VRB.Flow(cname, fname, "cnt =%d hit=%d \n",cnt, hit_cnt);

  //int h=hash(x,mu);
  list_head * hash_list= &hash_tab[hash(x,mu)];
  int node_id = GetNodeId(x);
  list_head * list_iter=hash_list->next;
  list_head * list_p;
  LinkEntry * link_entry;
  //printf("%d %d %d %d %d %d %d\n", x[0],x[1],x[2],x[3],mu, h, node_id );
  cnt++;
  
  while(list_iter!=hash_list) {
    link_entry = list_entry(list_iter, LinkEntry, hash_entry); 
    if (link_entry->d_node_id == node_id){
      hit_cnt++;

      //move the link from its current position in the flush_list
      //to the end of it. so it won't be easily flushed too soon.
      //a very simple policy it may need improvements
      //-------------------------------------------------------------
      //list_p =& link_entry->flush_entry;
      //list_del(list_p);
      //list_add(list_p, flush_list.prev );

      //move the link from its current position in the hash list to the 
      //head of it, so it's easier to locate it next item.
      //-------------------------------------------------------------
      list_del(list_iter);
      list_add(list_iter, hash_list); 

      return &link_entry->link;
    }
    list_iter= list_iter->next;
  } 
  
  if(list_empty(&free_list)){
      //free_list is empty, the flush_list must not be empty. 
    list_p = flush_list.next; 
    list_del(list_p);
      //remove the first element from the flush_list
    list_add(list_p, flush_list.prev);
      //append to the end of the flush_list, so it don't get flushed too soon
    link_entry = list_entry(list_p, LinkEntry, flush_entry); 
      //get its LinkEntry
    list_p = &link_entry -> hash_entry;
    list_del(list_p);
      //remove it from its current hash_list
    list_add(list_p, hash_list); 
      //add it to the head of the new hash_list
  }
  else {
      //free_list is not empty.

    list_p = free_list.next; 
    list_del(list_p);
      //remove an element from the head of free_list

    list_add(list_p, hash_list);
      //add it to the hash_list

    link_entry = list_entry(list_p, LinkEntry, hash_entry); 
    list_add(&link_entry->flush_entry, flush_list.prev);
      //get its LinkEntry and append it to the end of the flush_list
  } 
  link_entry->d_node_id = node_id;
  
  moveMem((IFloat*)&link_entry->link,(IFloat*) lat.GetLink(x, mu),
	 MATRIX_SIZE * sizeof(IFloat)); 
  return & link_entry->link;
}

//----------------------------------------------------------------------
/*!
  Remove from the buffer all the links with local lattice site \a x
  and direction \a mu.
  \param x The lattice coordinates.
  \param mu The direction index.
*/
//----------------------------------------------------------------------

void LinkBuffer::ClearBufferedLink(const int *x, int mu){
  char * fname = "ClearBufferedLink()";
  VRB.Func(cname, fname);
  
  list_head * hash_list = &hash_tab[hash(x,mu)];
  //get the hash list in the hash table
  
  list_head * list_iter=hash_list->next;
  LinkEntry * link_entry;
  int i=0; 
  while (list_iter != hash_list ) {
    i++;
    link_entry = list_entry(list_iter, LinkEntry, hash_entry); 
    //get the LinkEntry of this object
    
    list_del(&link_entry->flush_entry);
    //delete this LinkEntry from the  flush_list 

    list_iter = list_iter->next;
  }

  //printf("clearbufferedlink() removing %d links\n", i); 
  list_splice(hash_list, & free_list);
  //move all elements on the hash_list to the free_list
  INIT_LIST_HEAD(hash_list); 
  //clean up this hash_list 
  //if(check_lists())ERR.Pointer(cname, fname, "check_list");
}

void LinkBuffer::ClearAll(){
  for(int i=0;i<tab_size;i++){
    list_splice(&hash_tab[i], &free_list);
    INIT_LIST_HEAD(&hash_tab[i]);
  }
  INIT_LIST_HEAD(&flush_list);
}

void Lattice::ClearAllBufferedLink(){
  if(LinkBufferIsEnabled())
    link_buffer->ClearAll();
}

/*!
  This function does not create a buffer if there already is one.
  \param buf_sz The size of the link buffer.
  \return True if there is a link buffer, false otherwise.
*/
int Lattice::EnableLinkBuffer(int buf_sz){
  //char * fname = "EnableLinkBuffer()";
  //VRB.Func(cname, fname);
  if(link_buffer==0 && buf_sz > 0 )
    link_buffer = new LinkBuffer(*this, buf_sz);
  if(link_buffer) return 1;
  else return 0;
}

void Lattice::DisableLinkBuffer(){
  if(link_buffer != 0){
    delete link_buffer;
    link_buffer = 0;
  }
}

/*!
  Looks for the link \a U_mu(x) in the buffer. If it is not there it is
  brought in.
  \param x The lattice coordinates of the link.
  \param mu The direction index of the link.
  \return The link \a U_mu(x).
*/
const Matrix * Lattice::
GetBufferedLink(const int *x , int mu){
  if(IsOnNode(x)) return gauge_field+GsiteOffset((int*)x) + mu;    
  if(LinkBufferIsEnabled()) return link_buffer->GetBufferedLink(x, mu);
  return GetLink(x,mu);
}

/*!
  Remove from the buffer all the links with local lattice site \a x
  and direction \a mu.
  \param x The lattice coordinates.
  \param mu The direction index.
*/
void  Lattice::
ClearBufferedLink(const int *x, int mu){
  if(LinkBufferIsEnabled())
    link_buffer->ClearBufferedLink(x, mu); 
}

/*!
  \param x The lattice site coordinates
  \return 1 if x is on this node, 0 otherwise.
 */
int Lattice::
IsOnNode(const int * x){
  for(int i=0;i<4;i++)
   if(x[i]<0 || x[i]>=node_sites[i]) return 0; 
  return 1;
}


//------------------------------------------------------------------
//Lattice::BufferedStaple()
/*!
  The staple sum around the link \f$ U_\mu(x) \f$is
  
\f[
  \sum_{v \neq \mu} [                                  
              U_\nu(x+\mu) U_\mu(x+\nu) U^\dagger_\nu(x)                  
           +  U^\dagger_\nu(x+\mu-\nu) U^\dagger_u(x-\nu) U_\nu(x-\nu) ]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedStaple(Matrix& stap, const int *x, int u){
  //char * fname = "BufferedStaple()";
  //VRB.Func(cname, fname);
  int link_site[4];

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  acumulate_mp->ZeroMatrix();
  
  for(int i=0;i<4;i++)link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[3];
  stap.ZeroMatrix();
  
  int u1= u+4;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
    dir[0]=v;     //v0
    dir[1]=u1;
    dir[2]=(v+4)&7; //v1
    PathOrdProdPlus(*acumulate_mp, link_site, dir, 3);
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp,
    //		   MATRIX_SIZE*sizeof(IFloat));
  } 
  moveMem((IFloat*)&stap,  (IFloat*) acumulate_mp, MATRIX_SIZE*sizeof(IFloat)); 
}


//------------------------------------------------------------------
//Lattice::BufferedRectStaple()
/*! The 5-link rectangle staple sum around the link U_\mu(x) is:

\f[
 \sum_{\nu \neq \mu} \left[\right.                                     
 U_\mu(x+\mu)  U_\nu(x+2\mu)  U^\dagger_\mu(x+\mu+\nu)
 U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x)       \f]\f[
 + U_\mu(x+\mu)  U^\dagger_\nu(x+2\mu-\nu) U^\dagger_\mu(x+\mu-\nu)
 U^\dagger_\mu(x-\nu)  U_\nu(x-\nu)  \f]\f[
 + U_\nu(x+\mu)  U^\dagger_\mu(x+\nu)  U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)      \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  U^\dagger_\mu(x-\mu-\nu)
 U_\nu(x-\mu-\nu) U_\mu(x-\mu)  \f]\f[
 + U_\nu(x+\mu)  U_\nu(x+\mu+\nu)  U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\nu(x+\mu-2\nu) U^\dagger_\mu(x-2\nu)
 U_\nu(x-2\nu)  U_\nu(x-\nu)
 \left.\right]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedRectStaple(Matrix& stap, const int *x, int u){
  //char * fname = "BufferedRectStaple()";
  //VRB.Func(cname, fname);
  
  int link_site[4];
  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  acumulate_mp->ZeroMatrix();
  
  for(int i=0;i<4;i++)link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[5];
  stap.ZeroMatrix();
  int u1= u+4;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
      
    int v1= (v+4)&7;
      
    dir[0] = v; 
    dir[1] = v;
    dir[2] = u1;
    dir[3] = v1;
    dir[4] = v1;
      
    PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
      
    //dir[0] = v 
    dir[1] = u1;
    //dir[2] = u1;
    //dir[3] = v1; 
    dir[4] = u;
      
    PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
      
    dir[0] = u; 
    dir[1] = v;
    //dir[2] = u1;
    dir[3] = u1; 
    dir[4] = v1;
    PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
    
  }
  moveMem((IFloat*)&stap, (IFloat*)acumulate_mp, MATRIX_SIZE*sizeof(IFloat));
}


//------------------------------------------------------------------
/*! The chair shaped 5-link staple sum around the link U_\mu(x) is:
  
\f[
 \sum_{\pm \nu, |\nu|\neq \mu} \sum_{\pm \rho, |\rho|\neq \nu, |\rho|\neq \mu}(
\left[\right.
 U_\rho(x+\mu) U_\nu(x+\mu+\rho) U_{-\rho}(x+\mu+\nu+\rho)
 U_{-\mu}(x+\mu+\nu) U_{-\nu}(x+\nu)   \f]\f[
 + U_\rho(x+\mu) U_\nu(x+\mu+\rho) U_{-\mu}(x+\mu+\nu+\rho)
 U_{-\nu}(x+\nu+\rho) U_{-\rho}(x+\rho)    \f]\f[
 + U_\rho(x+\mu) U_{-\mu}(x+\mu+\rho) U_\nu(x+\rho)
 U_{-\rho}(x+\rho+\nu) U_{-\nu}(x+\nu)        
 \left.\right]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedChairStaple(Matrix& stap, const int *x, int u){
 
  int link_site[4];
  
  for(int i=0;i<4;i++) link_site[i]=x[i]; 

  link_site[u]++;

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  acumulate_mp->ZeroMatrix();
  
  int dir[5];
  Matrix m_tmp;
  
  int u1= u+4;
  
  for(int v=0; v<8; v++){
    
    if((v&3)==u) continue;
      
    int v1= (v+4)&7;
      
    for(int w=0; w<8; w++){
      if((w&3) ==u || (w&3) == (v&3) ) continue;
	  
      int w1= (w+4)&7;
	  
      dir[0] = w; 
      dir[1] = v;
      dir[2] = w1;
      dir[3] = u1;
      dir[4] = v1;

      PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //          MATRIX_SIZE*sizeof(IFloat)); 
	  
      //dir[0] = w 
      //dir[1] = v
      dir[2] = u1;
      dir[3] = v1; 
      dir[4] = w1;
      PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //            MATRIX_SIZE*sizeof(IFloat)); 
	  
      //dir[0] = w; 
      dir[1] = u1;
      dir[2] = v;
      dir[3] = w1; 
      dir[4] = v1;
      PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //             MATRIX_SIZE*sizeof(IFloat)); 

    }
  }
  moveMem((IFloat *) &stap, (IFloat*)acumulate_mp, MATRIX_SIZE*sizeof(IFloat));
}
//------------------------------------------------------------------------
/*!
  The staple sum around the link U_\mu(x) is
\f[
\sum_{\pm \nu, |\nu|\neq \mu} \sum_{\pm \rho, |\rho|\neq \nu, |\rho|\neq \mu}

     U_\nu(x+\mu) U_\rho(x+\mu+\nu) U_{-\mu}(x+\mu+\nu+\rho) U_{-\nu}(x+\nu+\rho) U_{-\rho}(x+\rho)
    
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/     
//------------------------------------------------------------------------
void Lattice:: 
BufferedCubeStaple(Matrix &stap, const int *x, int u){ 
  int link_site[4];

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  acumulate_mp->ZeroMatrix();

  for(int i=0;i<4;i++) link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[5];
  stap.ZeroMatrix();
  Matrix m_tmp;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
    int v1 = (v+4)&7;
    for(int w=0; w<8; w++){
      if((w&3) ==u || (w&3) == (v&3) ) continue;
      dir[0] = v; 
      dir[1] = w;
      dir[2] = u+4;
      dir[3] = v1;
      dir[4] = (w+4)&7;

      PathOrdProdPlus(*acumulate_mp, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //	       MATRIX_SIZE*sizeof(IFloat)); 
      
    }
  }
  moveMem((IFloat*)&stap, (IFloat*)acumulate_mp, MATRIX_SIZE*sizeof(IFloat));
}


//-------------------------------------------------------------------
/*!
  Given the starting site x, the directions of each step on the path
  and the number of steps. calculate the path ordered product of
  all the links along the path and add it to the given matrix m.
  Each direction is one of 0, 1, 2, 3, 4, 5, 6 or 7} corresponding to
  the directions X, Y, Z, T, -X, -Y, -Z and -T respectively.

  \param m The initial matrix.
  \param x The coordinates of the starting point of the path
  \param dirs The list of directions.
  \param n The number of links in the path.
  \post \a The product along the path is added to \a m.
*/
// 
// in this implementation, each link is retrieved from other sites,
// and assembled on the local node.
//
// in another implementation, one could imagine passing a matrix
// around along the path, and multiply the links to it along the way.
// this would save a lot of communication, in
// many cases, but not optimal if there is only one offnode link.
// and the link buffer would not be of much use in this case.
//-------------------------------------------------------------------
void Lattice::
PathOrdProdPlus(Matrix & mat, const int * x, const int* dirs, int n){
  //char * fname = "PathOrdProd"; 
  //VRB.Flow(cname, fname,"(,,,%d)\n",n);

  int abs_dir;
  int dir_sign;
  int link_site[4];
  
  const Matrix * p1;
  //Matrix m1, m2, m3;

  int i;
  for(i=0;i<4;i++)link_site[i]=x[i]; 

  //deal with the first link
  //------------------------------
  abs_dir = dirs[0]&3; 
    //abs_dir is {0,1,2,3}
  dir_sign = dirs[0]>>2; 
    //dir_sign is {0, 1}

  link_site[abs_dir] -= dir_sign; 
    //if dir_sign == 1, the link is at x-n_v

  p1 = GetBufferedLink(link_site, abs_dir);  
    //get the first link 

  link_site[abs_dir] += 1-dir_sign; 
    //if dir_sign == 0, march on to the next site, if dir_sign == 1, we have
    //already moved.

  if (dir_sign) 
      result1_mp->Dagger((IFloat*)p1); 
        //if dir_sign==1 the link is going backward so get its dagger
  else 
      moveMem((IFloat*) result1_mp, (IFloat*)p1, 
              MATRIX_SIZE * sizeof(IFloat));
        //simply move to the cram
  
  for(i=1;i<n;i++){
    abs_dir = dirs[i]&3; 
    dir_sign = dirs[i]>>2; 
    
    link_site[abs_dir] -= dir_sign; 

    p1 = GetBufferedLink(link_site, abs_dir);  

    link_site[abs_dir] += 1-dir_sign; 

    //put the next link on the path in mat1
    //--------------------------------------
    if(dir_sign){  
      mat1.Dagger((IFloat*)p1); 
    }
    else 
      moveMem((IFloat*)&mat1, (IFloat*)p1, 
              MATRIX_SIZE*sizeof(IFloat));

    if(i!=n-1)
      mDotMEqual((IFloat*) result_mp, (IFloat*)result1_mp, (IFloat*)&mat1);
        //if not the last link on the path, just multiply to the earlier result
    else
      mDotMPlus((IFloat*) &mat, (IFloat*) result1_mp, (IFloat*) &mat1);
        //if the last link, multiply and add to mat.

    Matrix * tmp_p = result1_mp; result1_mp=result_mp; result_mp = tmp_p;
    //swap result_mp and result1_mp;
  }
}

//-------------------------------------------------------------------
/*!
  Given the starting site x, the directions of each step on the path
  and the number of steps. calculate the path ordered product of
  all the links along the path and add it to the given matrix m.
  Each direction is one of 0, 1, 2, 3, 4, 5, 6 or 7} corresponding to
  the directions X, Y, Z, T, -X, -Y, -Z and -T respectively.

  The idea is whenever the path hits a boundary the current partial result is 
  passed to the next processer on the path, which will calculate the part of 
  the product that is on this node and pass on, so the result ends up on the 
  last processor where the path stops.
  
  \param m The product along the path.
  \param x The coordinates of the starting point of the path
  \param dirs The list of directions.
  \param n The number of links in the path.
*/
//-------------------------------------------------------------------

void Lattice::
PathOrdProd(Matrix & mat, const int * x, const int* dirs, int n){
  //char * fname = "PathOrdProd"; 
  //VRB.Flow(cname, fname,"(,,,%d)\n",n);

  int abs_dir;
  int dir_sign;
  int link_site[4];
 
  const Matrix * p1;
  Matrix m1, m2,  m4;
  Matrix * r1_mp = &mat1;
  Matrix * r_mp  = &mat2;
  Matrix * buf1_mp=&mat3;

  int i;
  for(i=0;i<4;i++) { 
    int l_x = x[i];
    while(l_x<0) l_x += node_sites[i]; 
    while(l_x>=node_sites[i]) l_x -= node_sites[i];
    link_site[i]=l_x;
  }

  //deal with the first link
  //-------------------------
  abs_dir = dirs[0]&3; 
  dir_sign = dirs[0]>>2; 

  link_site[abs_dir] -= dir_sign; 
    //go to the site where the link is stored
  if(link_site[abs_dir]<0) link_site[abs_dir] +=node_sites[abs_dir];
    //if out of boundary move to the other boundary in this direction
    //where a path is entering the territory of this node.

  p1 = gauge_field+GsiteOffset((int*)link_site) + abs_dir;
    //get the fisrt link

  link_site[abs_dir] += 1-dir_sign;
    //march on to the next site, if haven't done so 

  if(link_site[abs_dir]==node_sites[abs_dir]){ 
      //just hit another boundary we need to pass the partial result
      //on to the positive direction and get a link from the minus direction
    link_site[abs_dir]=0;
    getMinusData((IFloat*) &m1, (IFloat*) p1, MATRIX_SIZE, abs_dir); 
    moveMem((IFloat*)buf1_mp, (IFloat*)&m1, MATRIX_SIZE);
      //move it to the cram
  }
  else
    moveMem((IFloat*)buf1_mp, (IFloat*)p1, MATRIX_SIZE);
      //didn't hit the boundary, directly move the cram

#define SWAP(a, b) {Matrix* tmp_mp = (a); (a)=(b); (b)=tmp_mp;}

  if (dir_sign) 
     //link going backward so take its hermite conjugate 
    r1_mp->Dagger((IFloat*) buf1_mp); 
  else 
     //put the most recent result in r1_mp
    SWAP(r1_mp, buf1_mp) 
      

  for(i=1;i<n;i++){
    abs_dir = dirs[i]&3; 
      //abs_dir is {0,1,2,3}
    dir_sign = dirs[i]>>2; 
      //dir_sign is {0, 1}
    
    link_site[abs_dir] -= dir_sign; 
      //go to the site of the link

    if(link_site[abs_dir]<0){
        //the next link is offnode pass the current result to that node
        //and receive one from the plus direction
      link_site[abs_dir] += node_sites[abs_dir];
      moveMem((IFloat*)&m1, (IFloat*)r1_mp, MATRIX_SIZE); 
        //move from cram to dram, so scu can access.
      getPlusData((IFloat*)&m2, (IFloat*)&m1, MATRIX_SIZE, abs_dir); 
      moveMem((IFloat*)r1_mp, (IFloat*)&m2, MATRIX_SIZE); 
        //move from dram to cram to speed matrix multiplication.
    } 
    
    p1 = gauge_field+GsiteOffset((int*)link_site) + abs_dir;
      //get the next link
    
    if(dir_sign){ 
       //going backward, take the hermite conjugate
      buf1_mp->Dagger((IFloat*)p1); 
    }
    else{ 
      moveMem((IFloat*)buf1_mp,(Matrix *)p1, MATRIX_SIZE*sizeof(IFloat));
       //simply move to the cram
    }

    mDotMEqual((IFloat*) r_mp, (IFloat*)r1_mp, (IFloat*)buf1_mp);
      //every thing is in cram, just do the multiplication

    link_site[abs_dir] += 1-dir_sign;
       //march on to the next site 

    if(link_site[abs_dir]==node_sites[abs_dir]){
       //just hit the wall in the plus direction, pass on
      link_site[abs_dir]=0;
      moveMem((IFloat*)&m1, (IFloat*)r_mp, MATRIX_SIZE); 
       //cram to dram
      getMinusData((IFloat*) &m2, (IFloat*) &m1, MATRIX_SIZE, abs_dir);
      moveMem((IFloat*)r1_mp, (IFloat*)&m2, MATRIX_SIZE*sizeof(IFloat)); 
       //dram to cram
    }
    else 
       //didn't hit the wall, swap r1_mp, r_mp so r1_mp will hold the most
       //recent result of the products.
      SWAP(r1_mp, r_mp)
        
#undef SWAP
  }
  moveMem((IFloat*)&mat, (IFloat*)r1_mp, MATRIX_SIZE* sizeof(IFloat));
}


CPS_END_NAMESPACE
