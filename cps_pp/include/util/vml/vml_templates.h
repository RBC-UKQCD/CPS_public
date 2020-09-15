#ifndef VML_TEMPLATES_H
#define VML_TEMPLATES_H

CPS_END_NAMESPACE
#include <string.h>
CPS_START_NAMESPACE
/*C.Kelly 2012 - templates for rpccommands*/


/*Deep copy of vml class types. Please ensure that the rpccommand GENERATE_DEEPCOPY_METHOD is run contained
  in the .x definition of any member requiring a *non-trivial copy constructor* (e.g. classes containing pointers, arrays or unions);
  this creates a template specialization for that class type. This cannot be checked at compile time without forcing the command also 
  to be applied to all enums used within. This can be fixed if we start using boost, which has an is_enum check. Also type_traits in TR1
  can do the same thing (and presumably c++11), but I don't want to force new standards when we might continue to support older compilers
*/

template<typename T>
struct rpc_deepcopy{
  //copy routine for single object. Specialize this for stucts and classes that require a non-trivial copy constructor to perform a deep copy
  static void doit(T&into, T const&from){ into = from; }
};
template<typename T>
struct rpc_deepcopy<T*>{
  //could be single object or array, need to provide a size
  static void doit(T* &into, T const* const &from, const int &size){ 
    into = new T[size]; for(int i=0;i<size;i++) rpc_deepcopy<T>::doit(into[i],from[i]);
  }
};
template<typename T>
struct rpc_deepcopy<T const *>{
  //pointer to const object (or array of). Must allocate new memory location then reassign pointer
  //T const* const &from   is a constant pointer to a constant object!
  static void doit(T const* &into, T const* const &from, const int &size){ 
    T* tmp; rpc_deepcopy<T *>::doit(tmp,from,size); into = tmp;
  }
};

#ifndef _USE_STDLIB
#define _USE_STDLIB
#endif

#ifdef _USE_STDLIB

CPS_END_NAMESPACE
#include<iostream>
#include<string>
#include<sstream>
CPS_START_NAMESPACE

template<typename T>
struct rpc_print{
  static void doit(T const&what, const std::string &prefix="" ){ std::cout << prefix << what << std::endl; }
};
template<typename T>
struct rpc_print<T*>{
  static void doit(T const * const &what, const int &size, const std::string &prefix="" ){
    if(size == 0) return;

    if(size>1){
      std::cout << prefix << " array of size " << size << " at " << static_cast<void const* const &>(what) << std::endl;
      for(int i=0;i<size;i++){
	std::ostringstream os; os << prefix << "[" << i << "] = ";
	rpc_print<T>::doit(what[i],os.str());
      }
    }else{
      std::cout << prefix << " object at " << static_cast<void const* const &>(what) << std::endl;
      rpc_print<T>::doit(*what,prefix);
    }

  }
};
template<typename T>
struct rpc_print<T const *>{
  static void doit(T const * const &what, const int &size, const std::string &prefix="" ){ 
    rpc_print<T*>::doit(what,size,prefix);
  }
};
template<>
struct rpc_print<char *>{
  static void doit(char const * const &what, const int &size, const std::string &prefix="" );
};

#endif



#endif
