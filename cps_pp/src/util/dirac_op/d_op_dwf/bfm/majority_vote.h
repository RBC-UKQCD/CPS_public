
//typedef bfm_qdp<double> bfm_t;
//
// Majority vote incarnation overrides the global sum  
// to instead broadcast result on node zero and do bit identical compare
//
// Any nodes that dribble a bit will print "Oops, I did it again!".  
//
//
class majorityVote : public bfm_qdp<double> {
public:
  int called;
  majorityVote() {called=0;}
  void IdentifyNode(void) { 
// intentionally crash
    double *temp = NULL;
    *temp = 1.0;
  };
  void comm_gsum(double &val) 
  {
    int me = thread_barrier();
    double myval=val;
    thread_sum(myval,me);
    val = myval;
    if ( me == 0 ) { 
//      Printf("comm_gsum\n");
//    QMP_sum_double_array(&val,1);
//    myval = val;
      QMP_broadcast(&val,sizeof(val));
      if ( myval != val ){
	 printf("%d: Oops, I did it again! %0.14g %0.14g\n",called,myval,val);
      	IdentifyNode();
      }
      called++;
    }
    thread_barrier();
  }
};
