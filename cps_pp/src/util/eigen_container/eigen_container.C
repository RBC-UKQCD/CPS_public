#include<util/eig_io.h>
#include<util/eigen_container.h>


// needed to declare globally

CPS_START_NAMESPACE

#if 1
//const char *EvecReader::cname = "EvecReader";
const char* EvecReader::header =
"QCD eigenvector decompressor\n"
"Authors: Christoph Lehner, Ziyuan Bai, Chulwoo Jung\n"
"Date: 2017\n";
#endif

std::vector<EigenCache*> EigenCacheList(0);
//Search contents that match to arguments, return 0 if not found
EigenCache* EigenCacheListSearch( char* fname_root_bc, int neig )
{

  const char *fname("EigenCacheListSearch()");
  EigenCache* ecache=NULL;

  VRB.Flow("",fname,"cache_name=%s neig=%d\n",fname_root_bc,neig);
  VRB.Result("",fname,"EigenCacheList.size()=%d %p \n",EigenCacheList.size(),EigenCacheList[0]);

  for(int i=0; i< EigenCacheList.size(); ++i  ){
    VRB.Result("EigenCacheList",fname,"EigenCacheListSearch(%d): %s %d\n",
		i,fname_root_bc,neig);
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) ){
      ecache = EigenCacheList[i];
  if(VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
  printf("Node %d: %s cache_name=%s neig=%d ecache=%p \n",UniqueID(),fname, fname_root_bc,neig,ecache);fflush(stdout);
      return ecache;
    }
  }
  if(VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
  printf("Node %d: %s cache_name=%s neig=%d ecache=%p \n",UniqueID(),fname, fname_root_bc,neig,ecache);fflush(stdout);
 
  return ecache;

}


// Cleanup list EigenCache, it also destroys contents pointed by the elements.
void EigenCacheListCleanup( )
{
  for(size_t i=0;i< EigenCacheList. size(); ++i){
    EigenCacheList[i]-> dealloc(); 
   delete  EigenCacheList[i];
  }
  EigenCacheList.clear() ;
}
CPS_END_NAMESPACE
