#include<util/eigen_container.h>


// needed to declare globally

CPS_START_NAMESPACE
std::vector<EigenCache*> EigenCacheList(0);
//Search contents that match to arguments, return 0 if not found
EigenCache* EigenCacheListSearch( char* fname_root_bc, int neig )
{
  EigenCache* ecache=0;

  for(int i=0; i< EigenCacheList.size(); ++i  )
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) )
      ecache = EigenCacheList[i];

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
