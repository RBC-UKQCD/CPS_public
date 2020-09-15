// #ifndef MEM_DISTRIBUTE_H
// #define MEM_DISTRIBUTE_H

// #ifdef USE_MPI
// #warning "Using MPI RMA-distributed memory for A2A"

// class MemoryDistribute{
//   std::vector<size_t> rank_offset; //current byte offset within MPI window for each node
//   std::vector<size_t> data_offsets; //each block of data given to the manager is 

//   MPI_Win win; //the MPI window for this node
//   void* win_mem; //the memory associated with the window
//   size_t size; //window size in bytes
//  public:
//   MemoryDistribute(const size_t window_size){
//     size = window_size;
//     win_mem = malloc(window_size);
//     if(win_mem == NULL) ERR.General("MemoryDistribute","MemoryDistribute","Failed to allocate window memory\n");

//     int ret = MPI_Win_create(win_mem, window_size, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
//     if(ret != MPI_SUCCESS) ERR.General("MemoryDistribute","MemoryDistribute","Failed to create window\n");
//   }
//   //Store some data. It will be given to the node with the least memory currently allocated. Returns a tag for retrieval
//   int store


// };






// #else
// //Must be single-node


// #endif


// #endif
