CPS_START_NAMESPACE
#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H
CPS_END_NAMESPACE

#include<config.h>
#include<string>

CPS_START_NAMESPACE


class CorrelationFunction{
public:
  enum ThreadType { THREADED, UNTHREADED };
private:
  std::string label;
  int ncontract;
  int time_size;
  Rcomplex** wick; //unthreaded wick contractions [contraction][t]
  Rcomplex*** wick_threaded; //threaded wick contractions [thread][contraction][t]

  ThreadType threadtype;
  int max_threads; //the max amount of threads when the array of Rcomplex was created
  bool global_sum_on_write; //do a global summation (lattice sum + thread sum) before writing. Defaults to true.
public:
  inline void setLabel(const char* _label){ label = _label; }
  const std::string & getLabel() const{ return label; }

  void setNcontractions(const int &n);
  const int & nContractions() const{ return ncontract; }
  void write(const char *file);
  void write(FILE *fp);

  //Get the value at *global* time coordinate t
  Rcomplex & operator()(const int &contraction_idx, const int &t);
  Rcomplex & operator()(const int &thread_idx, const int &contraction_idx, const int &t);

  void sumThreads(); //Form the sum of wick[contraction][t] over all threads on the local node. Result is stored in wick and after sum accessible via operator() (i.e. the non-threaded version)
  void sumLattice(); //Form the sum of wick[contraction][t] over all nodes (and threads if threaded). Result is stored in wick and after sum accessible via operator() (i.e. the non-threaded version)
  void clear();

  void setGlobalSumOnWrite(const bool &b){ global_sum_on_write = b; } //Global sum on write is disabled automatically if sumLattice() is called manually

  //Set the thread type. Cannot be changed once memory has been allocated
  void setThreadType(const ThreadType &thread);
  const ThreadType & threadType() const{ return threadtype; }

  CorrelationFunction(const char *_label, const ThreadType &thread = UNTHREADED);
  CorrelationFunction(const char *_label, const int &n_contractions, const ThreadType &thread = UNTHREADED);
  CorrelationFunction();

  ~CorrelationFunction();
};

#endif
CPS_END_NAMESPACE
