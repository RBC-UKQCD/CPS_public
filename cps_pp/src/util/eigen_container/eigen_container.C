#include<util/eig_io.h>
#include<util/eigen_container.h>


// needed to declare globally

CPS_START_NAMESPACE

#if 1
const char *EvecReader::cname = "EvecReader";
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
  EigenCache* ecache=0;

  printf("EigenCacheList.size()=%d\n",EigenCacheList.size());

  for(int i=0; i< EigenCacheList.size(); ++i  ){
    VRB.Debug("EigenCacheList",fname,"EigenCacheListSearch(%d): %s %d\n",
		i,fname_root_bc,neig);
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) )
      ecache = EigenCacheList[i];
  }

  return ecache;
}


// Cleanup list EigenCache, it also destroies contents pointed by the elements.
void EigenCacheListCleanup( )
{
  for(size_t i=0;i< EigenCacheList. size(); ++i){
    EigenCacheList[i]-> dealloc(); 
  }
  EigenCacheList.clear() ;
}
CPS_END_NAMESPACE
