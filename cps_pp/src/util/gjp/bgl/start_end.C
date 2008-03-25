#include <config.h>
#include <util/gjp.h>
#include <comms/bgl_net.h>

CPS_START_NAMESPACE

void Start(int * argc, char ***argv)
{Start(NULL); GJP.setArg(argc,argv);}

void Start(){
  Start(NULL);
}

void Start(const  BGLAxisMap *axis_map){
  int i;
  if (axis_map != NULL){
    bgl_machine_dir[0] = 2*axis_map->bgl_machine_dir_x;
    bgl_machine_dir[1] = 2*axis_map->bgl_machine_dir_x+1;
    bgl_machine_dir[2] = 2*axis_map->bgl_machine_dir_y;
    bgl_machine_dir[3] = 2*axis_map->bgl_machine_dir_y+1;
    bgl_machine_dir[4] = 2*axis_map->bgl_machine_dir_z;
    bgl_machine_dir[5] = 2*axis_map->bgl_machine_dir_z+1;
    bgl_machine_dir[6] = 2*axis_map->bgl_machine_dir_t;
    bgl_machine_dir[7] = 2*axis_map->bgl_machine_dir_t+1;
  } else 
    for(i = 0; i<8 ; i++) 
      bgl_machine_dir[i] = i;
  
  for(i = 0; i<8 ; i++)
    bgl_cps_dir[bgl_machine_dir[i]] = i;
  if (UniqueID()==0)
  for(i = 0; i<8 ; i++){
    printf("bgl_machine_dir[%d]=%d bgl_cps_dir[%d]=%d\n",
    i, bgl_machine_dir[i], i, bgl_cps_dir[i]);
  }
  BGLCPSNet_Init();
}

void End(){
}
CPS_END_NAMESPACE
