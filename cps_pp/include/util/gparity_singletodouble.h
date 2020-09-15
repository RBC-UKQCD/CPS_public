#ifndef GPARITY_SINGLE_TO_DOUBLE_H
#define GPARITY_SINGLE_TO_DOUBLE_H
// #include <config.h>
// #include <util/lattice.h>
// #include <util/random.h>
// #include <util/gjp.h>

//CK 2012

//Class to take a G-parity array of objects for which separate versions exist for each G-parity flavour (e.g. RNGs and the gauge links) in the two-flavour implementation and 
//communicate them onto a lattice that has been doubled/quadrupled in the one-flavour implemenation.

//Assumes the layout has already been reinitialised with the new, doubled size.

//For usage cf. tests/gparity_1f_2f_compare/main.C
CPS_START_NAMESPACE

class SingleToDouble {
private:
  SingleToDouble(){} //cannot be created without parameters

  static int pos_quadrant(int *pos){
    //LL LR UR UL
    if(pos[0]<GJP.XnodeSites()/2 && pos[1]<GJP.YnodeSites()/2) return 0;
    if(pos[0]>=GJP.XnodeSites()/2 && pos[1]<GJP.YnodeSites()/2) return 1;
    if(pos[0]>=GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2) return 2;
    if(pos[0]<GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2) return 3;
  }
public:
  bool gparity_X;
  bool gparity_Y;

  SingleToDouble(bool _gparity_X, bool _gparity_Y): gparity_X(_gparity_X), gparity_Y(_gparity_Y){}

  virtual int SiteSize() = 0; //in units of Float
  virtual int Ndata() = 0; //total amount of data. Should be 2* the value without G-parity
  virtual void Store(Float *buf, int site)=0; //fill buffer
  virtual void LoadX(Float *buf, int* pos, int flav)=0; //load from buf into new data at position pos and flavour flav (gparity_X)
  virtual void LoadXY(Float *buf, int* pos, int flav)=0; //load from buf into new data at position pos and flavour flav (gparity_X && gparity_Y)
  virtual void IncrPos(int *pos) = 0; //increment the position vector
  
  virtual void StartPos(int *pos){
    pos[0]=0; pos[1]=0; pos[2]=0; pos[3]=0; pos[4]=0;
  }

  /*
    |0 1 2 3|         |4 5 6 7|               |8 9 10 11|               |12 13 14 15|
->
    |0 1 2 3 4 5 6 7| |8 9 10 11 12 13 14 15| |0* 1* 2* 3* 4* 5* 6* 7*| |8* 9* 10* 11* 12* 13* 14* 15*|
   */

  void RunGparityX(){
    if(GJP.Xnodes()>1){
      int ndata = Ndata();
      int sitesize =  SiteSize();
	
      int buf_size = ndata * sitesize*sizeof(Float); //send both flavours
      Float *recv_buf = (Float *) pmalloc(buf_size);
      Float *send_buf = (Float *) pmalloc(buf_size);

      if(!UniqueID()){ printf("Storing data to buffer\n"); fflush(stdout); }
      for(int site =0; site < ndata; site++){
	Store(send_buf,site);
      }
      int data_nodecoor_hf1;  //what xnode coor is this nodes data for the first halves currently stored on? (second half is always on the next node)
      int data_nodecoor_hf2;
      int data_flav = 0; 
      int x_origin = GJP.XnodeCoor()*GJP.XnodeSites(); //x position of start of first half
      if(GJP.XnodeCoor()>=GJP.Xnodes()/2){  
	x_origin = (GJP.XnodeCoor()-GJP.Xnodes()/2)*GJP.XnodeSites(); 
	data_flav = 1;
      }
      data_nodecoor_hf1 = (x_origin/(GJP.XnodeSites()/2) ) % GJP.Xnodes();
      data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Xnodes();

      bool printnode = false;
      if(printnode) printf("Node %d: need 1st half from node %d and second half from node %d (flav %d)\n",GJP.XnodeCoor(),data_nodecoor_hf1,data_nodecoor_hf2,data_flav);

      Float nodes_unhappy = 1.0;
      Float *cur_data_buf = send_buf;
      Float *send_buf_p = send_buf;
      Float *recv_buf_p = recv_buf;
      int xnode_coor_of_buf_data = GJP.XnodeCoor(); 
      int got_hf1 = 0;
      int got_hf2 = 0;

      if(!UniqueID()){ printf("Starting site loop\n"); fflush(stdout); }

      while(nodes_unhappy != 0.0){
	if(xnode_coor_of_buf_data == data_nodecoor_hf1 || xnode_coor_of_buf_data == data_nodecoor_hf2 ){
	  //if(xnode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, placing on 1st half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);
	  //else if(xnode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, placing on 2nd half (flav %d)\n",GJP.XnodeCoor(),xnode_coor_of_buf_data,data_flav);

	  int pos[5];
	  StartPos(pos);
	    
	  for(int site =0 ; site < ndata; site++){
	    if(pos[0]>=GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf2){
	      LoadX(cur_data_buf,pos,data_flav);
	      got_hf2 = 1;
	    }else if(pos[0]<GJP.XnodeSites()/2 && xnode_coor_of_buf_data == data_nodecoor_hf1){
	      LoadX(cur_data_buf,pos,data_flav);
	      got_hf1 = 1;
	    }
	    IncrPos(pos);
	  }
	}
	if(got_hf1 && got_hf2) nodes_unhappy = 0.0;
	else nodes_unhappy = 1.0;

	glb_sum(&nodes_unhappy);

	if(nodes_unhappy!=0.0){
	  cur_data_buf = recv_buf_p;

	  getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), 0); //converted ints to floats for comms. Why doesn't CPS allow sending of generic data types?
	  xnode_coor_of_buf_data = (xnode_coor_of_buf_data+1) % GJP.Xnodes();

	  //swap buffers over for next send
	  Float *tmp = recv_buf_p;
	  recv_buf_p = send_buf_p;
	  send_buf_p = tmp;
	}
      }
      pfree(recv_buf);
      pfree(send_buf);      
    }else{
      //single node
      int ndata = Ndata();
      int sitesize =  SiteSize();
	
      int buf_size = ndata * sitesize*sizeof(Float); //send both flavours
      Float *buf = (Float *) pmalloc(buf_size);
      if(!UniqueID()){ printf("Storing data to buffer\n"); fflush(stdout); }
      for(int site =0; site < ndata; site++){
	Store(buf,site);
      }
      
      int pos[5];
      StartPos(pos);
      for(int site =0 ; site < ndata; site++){
	if(pos[0]>=GJP.XnodeSites()/2){
	  LoadX(buf,pos,1);
	}else if(pos[0]<GJP.XnodeSites()/2){
	  LoadX(buf,pos,0);
	}
	IncrPos(pos);
      }
      pfree(buf);
    }

  }
    
  void RunGparityXY(){
    if( (GJP.Xnodes()>1 && GJP.Ynodes()==1) || (GJP.Ynodes()>1 && GJP.Xnodes()==1) ){
  /*
    |0 1 2 3|                 |4 5 6 7|                       |8 9 10 11|               |12 13 14 15|
->


^ shortaxis
|
|----->longaxis

    |0* 1* 2* 3* 4* 5* 6* 7*| |8* 9* 10* 11* 12* 13* 14* 15*| |0  1  2  3  4  5  6  7 | |8  9  10  11  12  13  14  15 |
    |0  1  2  3  4  5  6  7 | |8  9  10  11  12  13  14  15 | |0* 1* 2* 3* 4* 5* 6* 7*| |8* 9* 10* 11* 12* 13* 14* 15*|
   */


      int longaxis = 0; int shortaxis = 1;
      if(GJP.Ynodes()>1){ longaxis = 1; shortaxis = 0; }

      int ndata = Ndata(); //original amount of data
      int sitesize =  SiteSize();
	
      int buf_size = ndata * sitesize*sizeof(Float); //send both flavours
      Float *recv_buf = (Float *) pmalloc(buf_size);
      Float *send_buf = (Float *) pmalloc(buf_size);

      if(!UniqueID()){ printf("Storing data to buffer\n"); fflush(stdout); }
      for(int site =0; site < ndata; site++){
	Store(send_buf,site);
      }

      int data_nodecoor_hf1;  //what longaxis node coor is this nodes data for the first half currently stored on? (second half is always on the next node)
      int data_nodecoor_hf2;
      int data_flav_lower = 0; 
      int data_flav_upper = 1; //upper quadrants contain U* on first half of longaxis 
      int longaxis_origin = GJP.NodeCoor(longaxis)*GJP.NodeSites(longaxis); //longaxis position of start of first half

      if(GJP.NodeCoor(longaxis)>=GJP.Nodes(longaxis)/2){  
	longaxis_origin = (GJP.NodeCoor(longaxis)-GJP.Nodes(longaxis)/2)*GJP.NodeSites(longaxis); 
	data_flav_lower = 1;
	data_flav_upper = 0;
      }
      data_nodecoor_hf1 = (longaxis_origin/(GJP.NodeSites(longaxis)/2) ) % GJP.Nodes(longaxis);
      data_nodecoor_hf2 = (data_nodecoor_hf1+1) % GJP.Nodes(longaxis);

      bool printnode = true;
      if(printnode) printf("Node %d: along axis %d need 1st half from node %d and second half from node %d (flav %d)\n",GJP.NodeCoor(longaxis),longaxis,data_nodecoor_hf1,data_nodecoor_hf2,data_flav_lower);

      Float nodes_unhappy = 1.0;
      Float *cur_data_buf = send_buf;
      Float *send_buf_p = send_buf;
      Float *recv_buf_p = recv_buf;
      int lanode_coor_of_buf_data = GJP.NodeCoor(longaxis); 
      int got_LL = 0; //lower-left quadrant
      int got_LR = 0; //lower-right quadrant
      int got_UL = 0; //upper-left quadrant
      int got_UR = 0; //upper-right quadrant
      
      if(!UniqueID()){ printf("Starting site loop\n"); fflush(stdout); }

      while(nodes_unhappy != 0.0){
	printf("Node %d: buffer contains data from node %d\n",GJP.NodeCoor(longaxis),lanode_coor_of_buf_data); fflush(stdout);

	if(lanode_coor_of_buf_data == data_nodecoor_hf1 || lanode_coor_of_buf_data == data_nodecoor_hf2 ){
	  //if(lanode_coor_of_buf_data == data_nodecoor_hf1 && printnode) printf("Node %d: Buffer contains data from node %d, placing on 1st half (flav %d)\n",GJP.NodeCoor(longaxis),lanode_coor_of_buf_data,data_flav_lower);
	  //else if(lanode_coor_of_buf_data == data_nodecoor_hf2 && printnode) printf("Node %d: Buffer contains data from node %d, placing on 2nd half (flav %d)\n",GJP.NodeCoor(longaxis),lanode_coor_of_buf_data,data_flav_lower);

	  int pos[5];
	  StartPos(pos);
	  for(int site =0 ; site < 2*ndata; site++){
	    if(pos[longaxis]>=GJP.NodeSites(longaxis)/2 && pos[shortaxis]<GJP.NodeSites(shortaxis)/2 && lanode_coor_of_buf_data == data_nodecoor_hf2){
	      //printf("Node %d: pos (%d %d %d %d %d), LR quad\n",GJP.NodeCoor(longaxis),pos[0],pos[1],pos[2],pos[3],pos[4]); fflush(stdout);
	      LoadXY(cur_data_buf,pos,data_flav_lower);
	      got_LR = 1;
	    }else if(pos[longaxis]<GJP.NodeSites(longaxis)/2 && pos[shortaxis]<GJP.NodeSites(shortaxis)/2 && lanode_coor_of_buf_data == data_nodecoor_hf1){
	      //printf("Node %d: pos (%d %d %d %d %d), LL quad\n",GJP.NodeCoor(longaxis),pos[0],pos[1],pos[2],pos[3],pos[4]); fflush(stdout);
	      LoadXY(cur_data_buf,pos,data_flav_lower);
	      got_LL = 1;
	    }else if(pos[longaxis]>=GJP.NodeSites(longaxis)/2 && pos[shortaxis]>=GJP.NodeSites(shortaxis)/2 && lanode_coor_of_buf_data == data_nodecoor_hf2){
	      //printf("Node %d: pos (%d %d %d %d %d), UR quad\n",GJP.NodeCoor(longaxis),pos[0],pos[1],pos[2],pos[3],pos[4]); fflush(stdout);
	      LoadXY(cur_data_buf,pos,data_flav_upper);
	      got_UR = 1;
	    }else if(pos[longaxis]<GJP.NodeSites(longaxis)/2 && pos[shortaxis]>=GJP.NodeSites(shortaxis)/2 && lanode_coor_of_buf_data == data_nodecoor_hf1){
	      //printf("Node %d: pos (%d %d %d %d %d), UL quad\n",GJP.NodeCoor(longaxis),pos[0],pos[1],pos[2],pos[3],pos[4]); fflush(stdout);
	      LoadXY(cur_data_buf,pos,data_flav_upper);
	      got_UL = 1;
	    }

	    IncrPos(pos);
	  }
	}
	printf("Node %d: quad status LL %d LR %d UR %d UL %d\n",GJP.NodeCoor(longaxis),got_LL,got_LR,got_UR,got_UL); fflush(stdout);

	if(got_LL && got_LR && got_UL && got_UR) nodes_unhappy = 0.0;
	else nodes_unhappy = 1.0;

	glb_sum(&nodes_unhappy);

	if(nodes_unhappy!=0.0){
	  cur_data_buf = recv_buf_p;
	  if(!UniqueID()){ printf("Doing comms\n"); fflush(stdout); }

	  getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), longaxis); //converted ints to floats for comms. Why doesn't CPS allow sending of generic data types?
	  lanode_coor_of_buf_data = (lanode_coor_of_buf_data+1) % GJP.Nodes(longaxis);

	  //swap buffers over for next send
	  Float *tmp = recv_buf_p;
	  recv_buf_p = send_buf_p;
	  send_buf_p = tmp;
	}
      }
      pfree(recv_buf);
      pfree(send_buf);      
    }else if(GJP.Xnodes()>1 && GJP.Ynodes()>1){
  /*
    |16 17 18 19|                 |20 21 22 23|                       |24 25 26 27|               |28 29 30 31|
-------------------------------------------------------------------------------------------------------------
    |0  1  2  3 |                 |4  5  6  7 |                       |8  9  10 11|               |12 13 14 15|
->

    |16*  17*  18*  19*  20*  21*  22*  23* | |24*  25*  26*  27*  28*  29*  30*  31* |   |16   17   18   19   20   21   22   23  | |24   25   26   27   28   29   30   31 |
    |0*   1*    2*   3*   4*   5*   6*  7*  | |8*   9*   10*  11*  12*  13*  14*  15* |   |0    1    2    3    4    5    6    7   | |8    9    10   11   12   13   14   15 |
----------------------------------------------------------------------------------------------------------------------------------------------------------------
    |16   17   18   19   20   21   22   23  | |24   25   26   27   28   29   30   31  |   |16*  17*  18*  19*  20*  21*  22*  23* | |24*  25*  26*  27*  28*  29*  30*  31* |
    |0    1    2    3    4    5    6    7   | |8    9    10   11   12   13   14   15  |   |0*   1*   2*   3*   4*   5*   6*   7*  | |8*   9*   10*  11*  12*  13*  14*  15* |
   */
      
      int ndata = Ndata(); //original amount of data
      int sitesize =  SiteSize();
	
      int buf_size = ndata * sitesize*sizeof(Float); //send both flavours
      Float *recv_buf = (Float *) pmalloc(buf_size);
      Float *send_buf = (Float *) pmalloc(buf_size);

      if(!UniqueID()){ printf("Storing data to buffer\n"); fflush(stdout); }
      for(int site =0; site < ndata; site++){
	Store(send_buf,site);
      }

      int node_quadrant; //LL LR UR UL
      if(GJP.XnodeCoor()<GJP.Xnodes()/2 && GJP.YnodeCoor()<GJP.Ynodes()/2) node_quadrant = 0;
      else if(GJP.XnodeCoor()>=GJP.Xnodes()/2 && GJP.YnodeCoor()<GJP.Ynodes()/2) node_quadrant = 1;
      else if(GJP.XnodeCoor()>=GJP.Xnodes()/2 && GJP.YnodeCoor()>=GJP.Ynodes()/2) node_quadrant = 2;
      else if(GJP.XnodeCoor()<GJP.Xnodes()/2 && GJP.YnodeCoor()>=GJP.Ynodes()/2) node_quadrant = 3;

      int node_flav = 0; //U field
      if(node_quadrant == 1 || node_quadrant == 3) node_flav = 1; //U* field

      int x_global_node_ll = GJP.XnodeCoor()*GJP.XnodeSites(); //global x coordinate of lower-left corner site of node
      int y_global_node_ll = GJP.YnodeCoor()*GJP.YnodeSites();
      
      //determine original lattice coordinates of lower-left corner site
      int x_orig_global_node_ll = x_global_node_ll;
      if(node_quadrant == 1 || node_quadrant == 2) x_orig_global_node_ll -= GJP.XnodeSites()*GJP.Xnodes()/2;
      
      int y_orig_global_node_ll = y_global_node_ll;
      if(node_quadrant == 2 || node_quadrant == 3) y_orig_global_node_ll -= GJP.YnodeSites()*GJP.Ynodes()/2;

      //the four quadrants of data on this node come from four different nodes.
      //determine those nodes.
      
      int quad_nodes[4][2]; //[LL LR UR UL][x,y node coor]
      quad_nodes[0][0] = (x_global_node_ll/(GJP.XnodeSites()/2)) % GJP.Xnodes();
      quad_nodes[0][1] = (y_global_node_ll/(GJP.YnodeSites()/2)) % GJP.Ynodes();

      quad_nodes[1][0] = (quad_nodes[0][0]+1) % GJP.Xnodes();
      quad_nodes[1][1] = quad_nodes[0][1];
      
      quad_nodes[2][0] = quad_nodes[1][0];
      quad_nodes[2][1] = (quad_nodes[1][1]+1) % GJP.Ynodes();   
	
      quad_nodes[3][0] = quad_nodes[0][0];
      quad_nodes[3][1] = quad_nodes[2][1];

      printf("Node (%d,%d) has LL corner coord (%d,%d) corresponding to point (%d,%d) on original lattice. For the four quadrants on this node, we need data from LL = (%d,%d) LR = (%d,%d) UR = (%d,%d) UL = (%d,%d)\n",
	     GJP.XnodeCoor(),GJP.YnodeCoor(),x_global_node_ll,y_global_node_ll,x_orig_global_node_ll,y_orig_global_node_ll,
	     quad_nodes[0][0], quad_nodes[0][1], quad_nodes[1][0], quad_nodes[1][1], quad_nodes[2][0], quad_nodes[2][1], quad_nodes[3][0], quad_nodes[3][1]);
      fflush(stdout);

      
      Float nodes_unhappy = 1.0;
      Float *cur_data_buf = send_buf;
      Float *send_buf_p = send_buf;
      Float *recv_buf_p = recv_buf;

      int got_quadrants[4] = {0,0,0,0};

      int buf_data_from_node[2] = {GJP.XnodeCoor(),GJP.YnodeCoor()}; //x,y node coordinates of origin for data currently in buffer

      bool xrow_complete = false;
      int xcomms_count = 0;

      while(nodes_unhappy != 0.0){
	printf("Node (%d,%d): Data in buffer originated on node (%d,%d)\n",GJP.XnodeCoor(),GJP.YnodeCoor(),buf_data_from_node[0],buf_data_from_node[1]);
	fflush(stdout);

	if(buf_data_from_node[0] == quad_nodes[0][0] || buf_data_from_node[0] == quad_nodes[1][0] ||
	   buf_data_from_node[1] == quad_nodes[0][1] || buf_data_from_node[1] == quad_nodes[2][1] ){
	  printf("Node (%d,%d): Using data from buffer\n",GJP.XnodeCoor(),GJP.YnodeCoor()); fflush(stdout);

	  int pos[5];
	  StartPos(pos);
	  for(int site =0 ; site < 2*ndata; site++){
	    int quad = pos_quadrant(pos);
	    if(buf_data_from_node[0] == quad_nodes[quad][0] && buf_data_from_node[1] == quad_nodes[quad][1]){
	      LoadXY(cur_data_buf,pos,node_flav);
	      got_quadrants[quad] = 1;
	    }
	    IncrPos(pos);
	  }
	}

	if(got_quadrants[0] + got_quadrants[1] + got_quadrants[2] + got_quadrants[3] == 4) nodes_unhappy = 0.0;
	printf("Node (%d,%d): happiness level %d\n",GJP.XnodeCoor(),GJP.YnodeCoor());

	glb_sum(&nodes_unhappy);
	if(!UniqueID()){ printf("Global happiness level %f\n",nodes_unhappy); fflush(stdout); }
	

	if(nodes_unhappy!=0.0){
	  int send_dir = 0; //first send in x direction
	  if(xcomms_count == GJP.Xnodes()-1){ send_dir = 1; xcomms_count = 0; }
	  
	  cur_data_buf = recv_buf_p;

	  getPlusData((IFloat *)recv_buf_p, (IFloat *)send_buf_p, buf_size/sizeof(Float), send_dir); //converted ints to floats for comms. Why doesn't CPS allow sending of generic data types?
	  buf_data_from_node[send_dir] = (buf_data_from_node[send_dir]+1) % GJP.Nodes(send_dir);
	  if(send_dir == 0) xcomms_count++;

	  //swap buffers over for next send
	  Float *tmp = recv_buf_p;
	  recv_buf_p = send_buf_p;
	  send_buf_p = tmp;
	}
      }
      if(!UniqueID()){ printf("Finished quadrupling\n"); fflush(stdout); }

      pfree(recv_buf);
      pfree(send_buf);  
    }else{
      //single node (in directions that matter)
      int ndata = Ndata();
      int sitesize =  SiteSize();
	
      int buf_size = ndata * sitesize*sizeof(Float); //send both flavours
      Float *buf = (Float *) pmalloc(buf_size);
      if(!UniqueID()){ printf("Storing data to buffer\n"); fflush(stdout); }
      for(int site =0; site < ndata; site++){
	Store(buf,site);
      }
      
      int pos[5];
      StartPos(pos);
      for(int site =0 ; site < 2*ndata; site++){
	if(pos[0]>=GJP.XnodeSites()/2 && pos[1]<GJP.YnodeSites()/2){ //LR
	  LoadXY(buf,pos,1);
	}else if(pos[0]<GJP.XnodeSites()/2 && pos[1]<GJP.YnodeSites()/2){ //LL
	  LoadXY(buf,pos,0);
	}else if(pos[0]>=GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2){ //UR
	  LoadXY(buf,pos,0);
	}else if(pos[0]<GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2){ //UL
	  LoadXY(buf,pos,1);
	}
	IncrPos(pos);
      }
      pfree(buf);
    }

  }

  void Run(){
    if(gparity_X && !gparity_Y){
      RunGparityX();
    }else if(gparity_X && gparity_Y){
      RunGparityXY();
    }
  }
};

class RNGComms{
public:
  static void rngstore(UGrandomGenerator &what, Float *to){
    GaussianRandomGenerator * wg = (GaussianRandomGenerator *)(&what);
    int sz = wg->RNGints();
    int *tbuf = new int[sz];
    wg->store(tbuf);
    //convert to floats
    for(int i=0;i<sz;i++){
      to[i] = tbuf[i];
    }
    delete[] tbuf;
  }
  static void rngload(UGrandomGenerator &into, Float *from){
    GaussianRandomGenerator * wg = (GaussianRandomGenerator *)(&into);
    int sz = wg->RNGints();
    int *tbuf = new int[sz];
    //convert from floats
    for(int i=0;i<sz;i++){
      tbuf[i] = from[i];
    }
    wg->load(tbuf);
    delete[] tbuf;
  }

  static int rngnint(UGrandomGenerator &augrand){
    GaussianRandomGenerator * wg = (GaussianRandomGenerator *)(&augrand);
    return wg->RNGints();
  }
};

class SingleToDouble4dRNG : public SingleToDouble{
public:
  UGrandomGenerator *ugran_4d_orig;
  int rngsz;
  int n_rgen_4d_orig;
  
  SingleToDouble4dRNG(bool _gparity_X, bool _gparity_Y): SingleToDouble(_gparity_X,_gparity_Y){
    n_rgen_4d_orig = GJP.VolNodeSites()/16; 
    //factor of 2 in X direction included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) n_rgen_4d_orig/=2; //duplication in Y direction

    ugran_4d_orig = new UGrandomGenerator[Ndata()];
    for(int i=0;i<n_rgen_4d_orig;i++) ugran_4d_orig[i] = LRG.UGrandGen4D(i);
    rngsz = RNGComms::rngnint(ugran_4d_orig[0]);
  } 
  ~SingleToDouble4dRNG(){ delete[] ugran_4d_orig; }

  int SiteSize(){ return rngsz; }//in units of Float
  int Ndata(){ return n_rgen_4d_orig; }
  void Store(Float *buf, int site){
    Float* to = buf + rngsz*site;
    RNGComms::rngstore(ugran_4d_orig[site],to);
  }

  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    if(pos[0] >= GJP.XnodeSites()/2)
      orig_idx = (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + flav * n_rgen_4d_orig/2;
    else 
      orig_idx = pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + flav * n_rgen_4d_orig/2;

    Float* orig_data = buf + rngsz * orig_idx;
    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
    RNGComms::rngload(LRG.UGrandGen4D(new_idx),orig_data);
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    int posorig[4] = {pos[0],pos[1],pos[2],pos[3]}; //equivalent site on one-quarter lattice
    if(posorig[0]>=GJP.XnodeSites()/2) posorig[0]-=GJP.XnodeSites()/2;
    if(posorig[1]>=GJP.YnodeSites()/2) posorig[1]-=GJP.YnodeSites()/2;

    orig_idx = posorig[0]/2 + GJP.XnodeSites()/4*(posorig[1]/2 + GJP.YnodeSites()/4*(posorig[2]/2 + GJP.ZnodeSites()/2*posorig[3]/2)) + flav * n_rgen_4d_orig/2;

    Float* orig_data = buf + rngsz * orig_idx;
    int new_idx = pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
    RNGComms::rngload(LRG.UGrandGen4D(new_idx),orig_data);
  }

  void IncrPos(int *pos){
    //increment the position vector
    pos[0] +=2; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] +=2;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] +=2;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] +=2;
	}
      }
    }
  }

};

class SingleToDouble5dRNG : public SingleToDouble{
public:
  UGrandomGenerator *ugran_orig;
  int rngsz;
  int n_rgen_orig;
  int blocks_per_s_layer_orig;
  int blocks_per_s_layer_new;

  int stk_index_5d_off;

  SingleToDouble5dRNG(bool _gparity_X, bool _gparity_Y): SingleToDouble(_gparity_X,_gparity_Y){
    n_rgen_orig = GJP.VolNodeSites()/16;
    if (GJP.SnodeSites()>=2)
      n_rgen_orig = GJP.VolNodeSites()*GJP.SnodeSites() / 32;

    blocks_per_s_layer_new = n_rgen_orig /( GJP.SnodeSites() / 2 );

    if(gparity_X && gparity_Y) n_rgen_orig/=2; //duplication in Y direction

    blocks_per_s_layer_orig = n_rgen_orig /( GJP.SnodeSites() / 2 );
    stk_index_5d_off = blocks_per_s_layer_orig/2; //offset for R' on 5D orig latt (only double X : for XY we have to halve this)

    ugran_orig = new UGrandomGenerator[n_rgen_orig];
    for(int i=0;i<n_rgen_orig;i++) ugran_orig[i] = LRG.UGrandGen(i);
    rngsz = RNGComms::rngnint(ugran_orig[0]);
  } 
  ~SingleToDouble5dRNG(){ delete[] ugran_orig; }

  int SiteSize(){ return rngsz; }//in units of Float
  int Ndata(){ return n_rgen_orig; }
  void Store(Float *buf, int site){
    Float* to = buf + rngsz*site;
    if(!UniqueID()){ printf("Store to site %d, offset is %d, bufsize is %d\n",site,rngsz*site,n_rgen_orig*rngsz); fflush(stdout); }
    RNGComms::rngstore(ugran_orig[site],to);
  }

  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    if(pos[0] >= GJP.XnodeSites()/2)
      orig_idx = pos[4]/2*blocks_per_s_layer_orig + (pos[0]-GJP.XnodeSites()/2)/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + flav * stk_index_5d_off;
    else 
      orig_idx = pos[4]/2*blocks_per_s_layer_orig + pos[0]/2 + GJP.XnodeSites()/4*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2)) + flav * stk_index_5d_off;

    Float* orig_data = buf + rngsz * orig_idx;
    int new_idx = pos[4]/2*blocks_per_s_layer_new + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
    RNGComms::rngload(LRG.UGrandGen(new_idx),orig_data);
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    int posorig[5] = {pos[0],pos[1],pos[2],pos[3],pos[4]}; //equivalent site on one-quarter lattice
    if(posorig[0]>=GJP.XnodeSites()/2) posorig[0]-=GJP.XnodeSites()/2;
    if(posorig[1]>=GJP.YnodeSites()/2) posorig[1]-=GJP.YnodeSites()/2;

    orig_idx = posorig[4]/2*blocks_per_s_layer_orig + posorig[0]/2 + GJP.XnodeSites()/4*(posorig[1]/2 + GJP.YnodeSites()/4*(posorig[2]/2 + GJP.ZnodeSites()/2*posorig[3]/2)) + flav * stk_index_5d_off;
    Float* orig_data = buf + rngsz * orig_idx;

    int new_idx = pos[4]/2*blocks_per_s_layer_new + pos[0]/2 + GJP.XnodeSites()/2*(pos[1]/2 + GJP.YnodeSites()/2*(pos[2]/2 + GJP.ZnodeSites()/2*pos[3]/2));
    RNGComms::rngload(LRG.UGrandGen(new_idx),orig_data);
  }


  void IncrPos(int *pos){
    //increment the position vector
    pos[0] +=2; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] +=2;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] +=2;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] +=2;
	  if(pos[3]>=GJP.TnodeSites()){
	    pos[3] = 0;
	    pos[4] +=2;
	  }
	}
      }
    }
  }

};


class SingleToDoubleMatrixField : public SingleToDouble{
public:
  Matrix *orig_matrix; //as lattice is singleton user must take a copy of the lattice data before reinitialising for doubled lattice
  Matrix *dbl_matrix;
  int nsites_orig;
  int sitesize;

  SingleToDoubleMatrixField(bool _gparity_X, bool _gparity_Y,int nmat_per_site,Matrix *_orig_matrix, Matrix *_dbl_matrix): SingleToDouble(_gparity_X,_gparity_Y), orig_matrix(_orig_matrix),dbl_matrix(_dbl_matrix){
    nsites_orig = GJP.VolNodeSites(); //factor of 2 in X direction included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) nsites_orig/=2; //duplication in Y direction
    sitesize = nmat_per_site * 9 * 2; //9 matrix elements, 2 re/im
  } 

  int SiteSize(){ return sitesize; }//in units of Float
  int Ndata(){ return nsites_orig; }
  void Store(Float *buf, int site){
    Float *from = (Float*)(orig_matrix) + sitesize * site;
    Float* to = buf + sitesize*site;
    memcpy((void*)to,(void*)from,sitesize*sizeof(Float));
  }

  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    if(pos[0] >= GJP.XnodeSites()/2)
      orig_idx = (pos[0]-GJP.XnodeSites()/2) + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;
    else 
      orig_idx = pos[0] + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = (Float*)(dbl_matrix) + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    int posorig[4] = {pos[0],pos[1],pos[2],pos[3]}; //equivalent site on one-quarter lattice
    if(posorig[0]>=GJP.XnodeSites()/2) posorig[0]-=GJP.XnodeSites()/2;
    if(posorig[1]>=GJP.YnodeSites()/2) posorig[1]-=GJP.YnodeSites()/2;

    orig_idx = posorig[0] + GJP.XnodeSites()/2*(posorig[1] + GJP.YnodeSites()/2*(posorig[2] + GJP.ZnodeSites()*posorig[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = (Float*)(dbl_matrix) + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void IncrPos(int *pos){
    //increment the position vector
    pos[0] ++; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] ++;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] ++;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] ++;
	}
      }
    }
  }

};


class SingleToDoubleLattice : public SingleToDouble{
public:
  Matrix *orig_lattice; //as lattice is singleton user must take a copy of the lattice data before reinitialising for doubled lattice
  Matrix *dbl_lattice;
  int nsites_orig;
  int sitesize;

  SingleToDoubleLattice(bool _gparity_X, bool _gparity_Y, Matrix *_orig_lattice, Lattice &doubled_lattice): SingleToDouble(_gparity_X,_gparity_Y), orig_lattice(_orig_lattice){
    dbl_lattice = doubled_lattice.GaugeField();
    nsites_orig = GJP.VolNodeSites(); //factor of 2 in X direction (and 2 in Y-direction if XY) included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) nsites_orig/=2; //duplication in Y direction
    sitesize = 4 * 9 * 2; //4 directions, 9 matrix elements, 2 re/im
  } 

  int SiteSize(){ return sitesize; }//in units of Float
  int Ndata(){ return nsites_orig; }
  void Store(Float *buf, int site){
    Float *from = (Float*)(orig_lattice) + sitesize * site;
    Float* to = buf + sitesize*site;
    memcpy((void*)to,(void*)from,sitesize*sizeof(Float));
  }

  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    if(pos[0] >= GJP.XnodeSites()/2)
      orig_idx = (pos[0]-GJP.XnodeSites()/2) + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;
    else 
      orig_idx = pos[0] + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = (Float*)(dbl_lattice) + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    int posorig[4] = {pos[0],pos[1],pos[2],pos[3]}; //equivalent site on one-quarter lattice
    if(posorig[0]>=GJP.XnodeSites()/2) posorig[0]-=GJP.XnodeSites()/2;
    if(posorig[1]>=GJP.YnodeSites()/2) posorig[1]-=GJP.YnodeSites()/2;

    orig_idx = posorig[0] + GJP.XnodeSites()/2*(posorig[1] + GJP.YnodeSites()/2*(posorig[2] + GJP.ZnodeSites()*posorig[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = (Float*)(dbl_lattice) + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void IncrPos(int *pos){
    //increment the position vector
    pos[0] ++; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] ++;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] ++;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] ++;
	}
      }
    }
  }

};


class SingleToDoubleField : public SingleToDouble{
public:
  Float *orig_field;
  Float *dbl_field;
  int nsites_orig;
  int sitesize;

 SingleToDoubleField(bool _gparity_X, bool _gparity_Y,int _sitesize, Float *_orig_field, Float *_dbl_field): 
  SingleToDouble(_gparity_X,_gparity_Y), 
    orig_field(_orig_field),
    dbl_field(_dbl_field),
    sitesize(_sitesize){

    nsites_orig = GJP.VolNodeSites(); //factor of 2 in X direction included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) nsites_orig/=2; //duplication in Y direction
  } 

  int SiteSize(){ return sitesize; }//in units of Float
  int Ndata(){ return nsites_orig; }
  void Store(Float *buf, int site){
    Float *from = orig_field + sitesize * site;
    Float* to = buf + sitesize*site;
    memcpy((void*)to,(void*)from,sitesize*sizeof(Float));
  }

  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    if(pos[0] >= GJP.XnodeSites()/2)
      orig_idx = (pos[0]-GJP.XnodeSites()/2) + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;
    else 
      orig_idx = pos[0] + GJP.XnodeSites()/2*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = dbl_field + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int orig_idx;
    int posorig[4] = {pos[0],pos[1],pos[2],pos[3]}; //equivalent site on one-quarter lattice
    if(posorig[0]>=GJP.XnodeSites()/2) posorig[0]-=GJP.XnodeSites()/2;
    if(posorig[1]>=GJP.YnodeSites()/2) posorig[1]-=GJP.YnodeSites()/2;

    orig_idx = posorig[0] + GJP.XnodeSites()/2*(posorig[1] + GJP.YnodeSites()/2*(posorig[2] + GJP.ZnodeSites()*posorig[3])) + flav * nsites_orig/2;

    Float* orig_data = buf + sitesize * orig_idx;
    int new_idx = pos[0] + GJP.XnodeSites()*(pos[1] + GJP.YnodeSites()*(pos[2] + GJP.ZnodeSites()*pos[3]));

    Float* to = dbl_field + sitesize * new_idx;
    memcpy((void*)to,(void*)orig_data,sitesize*sizeof(Float));
  }

  void IncrPos(int *pos){
    //increment the position vector
    pos[0] ++; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] ++;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] ++;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] ++;
	}
      }
    }
  }

};

class SingleToDouble5dVectorField : public SingleToDouble{
public:
  Vector *orig_field;
  Vector *dbl_field;
  int nsites_orig;
  StrOrdType ord;
  int Ncb; //number of checkerboards in vector

  SingleToDouble5dVectorField(bool _gparity_X, bool _gparity_Y, Vector *_orig_field, Vector *_dbl_field, StrOrdType _ord, const int &_Ncb=2): 
    SingleToDouble(_gparity_X,_gparity_Y), 
    orig_field(_orig_field),
    dbl_field(_dbl_field),
    ord(_ord),
    Ncb(_Ncb){
      
      if(ord == WILSON){
	ERR.General("SingleToDouble5dVectorField","SingleToDouble5dVectorField(...)","WILSON ord conversion has not been tested and probably doesn't work!");
      }

    nsites_orig = GJP.VolNodeSites()*GJP.SnodeSites(); //factor of 2 in X direction included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) nsites_orig/=2; //duplication in Y direction
    if(ord == WILSON && Ncb == 1) nsites_orig/=2;
  } 

  int SiteSize(){ return 24; }//in units of Float
  int Ndata(){ return nsites_orig; }
  void Store(Float *buf, int site){
    Float *from = (Float*)orig_field + 24 * site;
    Float* to = buf + 24 * site;
    memcpy((void*)to,(void*)from,24*sizeof(Float));
  }

  int SiteOffset(int* pos, int flav, bool gparity, int *sz){
    if(ord == CANONICAL){
      if(gparity) return pos[0] + sz[0]*(pos[1]+sz[1]*(pos[2]+sz[2]*(pos[3]+sz[3]*(flav + 2*pos[4]))));
      else return pos[0] + sz[0]*(pos[1]+sz[1]*(pos[2]+sz[2]*(pos[3]+sz[3]*pos[4])));
    }else if(ord == WILSON){
      int index = pos[4]; if(gparity) index*=2;
      int vol = sz[4];
      int parity = (pos[4]+pos[3]+pos[2]+pos[1]+pos[0]+1)%2; //Odd first

      if(Ncb == 1 && parity == 1) ERR.General("SingleToDouble5dVectorField","SiteOffset","Invalid position for single checkerboard: %d,%d,%d,%d,%d\n",pos[0],pos[1],pos[2],pos[3],pos[4]);

      for(int i = 3; i>=0;i--){
	index = index*sz[i]+pos[i];
	vol *= sz[i];
      }
      int cboff = vol;
      if(gparity) cboff*=2;
      else flav =0;
      return (index + flav*vol/sz[4] + cboff*parity)/2;
    }else ERR.General("SingleToDouble5dVectorField","SiteOffset","Invalid fermion ordering\n");
  }


  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int pos_orig[5]; 
    memcpy((void*)pos_orig,(void*)pos,5*sizeof(int));
    if(pos[0] >= GJP.XnodeSites()/2) pos_orig[0] -= GJP.XnodeSites()/2;

    int sz[5] = {GJP.XnodeSites()/2, GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites(), GJP.SnodeSites()};
    int orig_idx = SiteOffset(pos_orig,flav,true,sz);
      
    Float* orig_data = buf + 24 * orig_idx;
    sz[0]*=2;
    int new_idx = SiteOffset(pos,0,false,sz);

    Float* to = (Float*)dbl_field + 24*new_idx;
    memcpy((void*)to,(void*)orig_data,24*sizeof(Float));

    //printf("LoadX from pos %d,%d,%d,%d,%d (%d) maps to %d,%d,%d,%d,%d flav %d: %f %f\n",pos[0],pos[1],pos[2],pos[3],pos[4],new_idx,pos_orig[0],pos_orig[1],pos_orig[2],pos_orig[3],pos_orig[4],flav,*orig_data,*to);
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int pos_orig[5]; 
    memcpy((void*)pos_orig,(void*)pos,5*sizeof(int));
    if(pos_orig[0]>=GJP.XnodeSites()/2) pos_orig[0]-=GJP.XnodeSites()/2;
    if(pos_orig[1]>=GJP.YnodeSites()/2) pos_orig[1]-=GJP.YnodeSites()/2;

    int sz[5] = {GJP.XnodeSites()/2, GJP.YnodeSites()/2, GJP.ZnodeSites(), GJP.TnodeSites(), GJP.SnodeSites()};
    int orig_idx = SiteOffset(pos_orig,flav,true,sz);
      
    Float* orig_data = buf + 24 * orig_idx;
    sz[0]*=2; sz[1]*=2;
    int new_idx = SiteOffset(pos,0,false,sz);

    Float* to = (Float*)dbl_field + 24*new_idx;
    memcpy((void*)to,(void*)orig_data,24*sizeof(Float));
    
    //upper-right quadrant needs (-) sign
    if( (GJP.Xnodes()>1 && GJP.Ynodes()>1 && GJP.XnodeCoor()>=GJP.Xnodes()/2 && GJP.YnodeCoor()>=GJP.Ynodes()/2) ||
	(GJP.Xnodes()>1 && GJP.Ynodes() == 1 && GJP.XnodeCoor()>=GJP.Xnodes()/2 && pos[1]>=GJP.YnodeSites()/2) ||
	(GJP.Xnodes()==1 && GJP.Ynodes()>1 && GJP.YnodeCoor()>=GJP.Ynodes()/2 && pos[0]>=GJP.XnodeSites()/2) ||
	(GJP.Xnodes()==1 && GJP.Ynodes()==1 && pos[0]>=GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2) ){
      for(int i=0;i<24;i++) to[i]*=-1;
    }

  }

  void IncrPos(int *pos){
    //increment the position vector
    pos[0] ++; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] ++;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] ++;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] ++;
	  if(pos[3]>=GJP.TnodeSites()){
	    pos[3] = 0;
	    pos[4] ++;
	  }
	}
      }
    }
    if(Ncb == 1 && (pos[4]+pos[3]+pos[2]+pos[1]+pos[0])%2 == 0) return IncrPos(pos); //skip even parity
  }
  void StartPos(int *pos){
    if( (ord == WILSON && Ncb==2) || ord == CANONICAL){
      pos[0]=0; pos[1]=0; pos[2]=0; pos[3]=0; pos[4]=0;
    }else{
      //odd sites only
      pos[0]=1; pos[1]=0; pos[2]=0; pos[3]=0; pos[4]=0;
    }
  }
};


class SingleToDouble4dVectorField : public SingleToDouble{
public:
  Vector *orig_field;
  Vector *dbl_field;
  int nsites_orig;
  StrOrdType ord;
  int Ncb; //number of checkerboards in vector

  SingleToDouble4dVectorField(bool _gparity_X, bool _gparity_Y, Vector *_orig_field, Vector *_dbl_field, StrOrdType _ord, const int &_Ncb=2): 
    SingleToDouble(_gparity_X,_gparity_Y), 
    orig_field(_orig_field),
    dbl_field(_dbl_field),
    ord(_ord),
    Ncb(_Ncb){
      
    nsites_orig = GJP.VolNodeSites(); //factor of 2 in X direction included already as VolNodeSites calculated *after* lattice doubling
    if(gparity_X && gparity_Y) nsites_orig/=2; //duplication in Y direction
    if(ord == WILSON && Ncb == 1) nsites_orig/=2;
  } 

  int SiteSize(){ return 24; }//in units of Float
  int Ndata(){ return nsites_orig; }
  void Store(Float *buf, int site){
    Float *from = (Float*)orig_field + 24 * site;
    Float* to = buf + 24 * site;
    memcpy((void*)to,(void*)from,24*sizeof(Float));
  }

  int SiteOffset(int* pos, int flav, bool gparity, int *sz){
    if(ord == CANONICAL){
      if(gparity) return pos[0] + sz[0]*(pos[1]+sz[1]*(pos[2]+sz[2]*(pos[3]+sz[3]*flav)));
      else return pos[0] + sz[0]*(pos[1]+sz[1]*(pos[2]+sz[2]*pos[3]));
    }else if(ord == WILSON){
      int index = 0;
      int vol = 1;
      int parity = (pos[3]+pos[2]+pos[1]+pos[0]+1)%2; //Odd first

      if(Ncb == 1 && parity == 1) ERR.General("SingleToDouble4dVectorField","SiteOffset","Invalid position for single checkerboard: %d,%d,%d,%d\n",pos[0],pos[1],pos[2],pos[3]);

      for(int i = 3; i>=0;i--){
	index = index*sz[i]+pos[i];
	vol *= sz[i];
      }
      int cboff = vol;
      if(gparity) cboff*=2;
      else flav =0;

      int offset = (index + flav*vol + cboff*parity)/2;
      //if(gparity) printf("2f siteOffset of site (%d %d %d %d) flav %d, parity is %d, offset is %d (float offset %d)\n",pos[0],pos[1],pos[2],pos[3],flav,parity,offset,offset*24);
      //else printf("1f siteOffset of site (%d %d %d %d), parity is %d, offset is %d (float offset %d)\n",pos[0],pos[1],pos[2],pos[3],parity,offset,offset*24);

      return offset;
    }else ERR.General("SingleToDouble4dVectorField","SiteOffset","Invalid fermion ordering\n");
  }


  void LoadX(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int pos_orig[4]; 
    memcpy((void*)pos_orig,(void*)pos,4*sizeof(int));
    if(pos[0] >= GJP.XnodeSites()/2) pos_orig[0] -= GJP.XnodeSites()/2;

    //printf("LoadX called with pos (%d %d %d %d) [flav %d]. Pos on original latt (%d %d %d %d)\n",pos[0],pos[1],pos[2],pos[3],flav,pos_orig[0],pos_orig[1],pos_orig[2],pos_orig[3]);
    

    int sz[4] = {GJP.XnodeSites()/2, GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites()};
    int orig_idx = SiteOffset(pos_orig,flav,true,sz);
      
    Float* orig_data = buf + 24 * orig_idx;
    sz[0]*=2;
    int new_idx = SiteOffset(pos,0,false,sz);

    Float* to = (Float*)dbl_field + 24*new_idx;
    memcpy((void*)to,(void*)orig_data,24*sizeof(Float));

    //printf("LoadX from pos %d,%d,%d,%d,%d (%d) maps to %d,%d,%d,%d,%d flav %d: %f %f\n",pos[0],pos[1],pos[2],pos[3],pos[4],new_idx,pos_orig[0],pos_orig[1],pos_orig[2],pos_orig[3],pos_orig[4],flav,*orig_data,*to);
  }

  void LoadXY(Float *buf, int* pos, int flav){
    //load from buf into new data at position pos and flavour flav
    int pos_orig[4]; 
    memcpy((void*)pos_orig,(void*)pos,4*sizeof(int));
    if(pos_orig[0]>=GJP.XnodeSites()/2) pos_orig[0]-=GJP.XnodeSites()/2;
    if(pos_orig[1]>=GJP.YnodeSites()/2) pos_orig[1]-=GJP.YnodeSites()/2;

    int sz[4] = {GJP.XnodeSites()/2, GJP.YnodeSites()/2, GJP.ZnodeSites(), GJP.TnodeSites()};
    int orig_idx = SiteOffset(pos_orig,flav,true,sz);
      
    Float* orig_data = buf + 24 * orig_idx;
    sz[0]*=2; sz[1]*=2;
    int new_idx = SiteOffset(pos,0,false,sz);

    Float* to = (Float*)dbl_field + 24*new_idx;
    memcpy((void*)to,(void*)orig_data,24*sizeof(Float));
    
    //upper-right quadrant needs (-) sign
    if( (GJP.Xnodes()>1 && GJP.Ynodes()>1 && GJP.XnodeCoor()>=GJP.Xnodes()/2 && GJP.YnodeCoor()>=GJP.Ynodes()/2) ||
	(GJP.Xnodes()>1 && GJP.Ynodes() == 1 && GJP.XnodeCoor()>=GJP.Xnodes()/2 && pos[1]>=GJP.YnodeSites()/2) ||
	(GJP.Xnodes()==1 && GJP.Ynodes()>1 && GJP.YnodeCoor()>=GJP.Ynodes()/2 && pos[0]>=GJP.XnodeSites()/2) ||
	(GJP.Xnodes()==1 && GJP.Ynodes()==1 && pos[0]>=GJP.XnodeSites()/2 && pos[1]>=GJP.YnodeSites()/2) ){
      for(int i=0;i<24;i++) to[i]*=-1;
    }

  }

  void IncrPos(int *pos){
    //printf("IncrPos incrementing position (%d %d %d %d)\n",pos[0],pos[1],pos[2],pos[3]);

    //increment the position vector
    pos[0] ++; 
    if(pos[0]>=GJP.XnodeSites()){
      pos[0] = 0;
      pos[1] ++;
      if(pos[1]>=GJP.YnodeSites()){
	pos[1] = 0;
	pos[2] ++;
	if(pos[2]>=GJP.ZnodeSites()){
	  pos[2] = 0;
	  pos[3] ++;
	}
      }
    }
    if(Ncb == 1 && (pos[3]+pos[2]+pos[1]+pos[0])%2 == 0){
      //printf("Odd checkerboard only and new position (%d %d %d %d) has even parity, skipping\n",pos[0],pos[1],pos[2],pos[3]);
      return IncrPos(pos); //skip even parity
    }
  }
  void StartPos(int *pos){
    if( (ord == WILSON && Ncb==2) || ord == CANONICAL){
      pos[0]=0; pos[1]=0; pos[2]=0; pos[3]=0;
    }else{
      //odd sites only
      pos[0]=1; pos[1]=0; pos[2]=0; pos[3]=0;
    }
  }
};


CPS_END_NAMESPACE
#endif
