#ifndef _TEMPLATE_WIZARDRY_H
#define _TEMPLATE_WIZARDRY_H

template <int LorR, typename T, typename U>
struct _selectLR{};
template <typename T, typename U>
struct _selectLR<0,T,U>{
  typedef T Type;
};
template <typename T, typename U>
struct _selectLR<1,T,U>{
  typedef U Type;
};

template<typename T,typename U>
struct _equal{
  enum { value = 0 };
};
template<typename T>
struct _equal<T,T>{
  enum { value = 1 };
};


template<int i,int j>
struct intEq{ static const bool val = false; };
template<int i>
struct intEq<i,i>{ static const bool val = true; };

template<bool v, typename T>
struct my_enable_if{};

template<typename T>
struct my_enable_if<true,T>{ typedef T type; };


template<typename T>
struct is_double_or_float{ enum {value = 0}; };

template<>
struct is_double_or_float<double>{ enum {value = 1}; };

template<>
struct is_double_or_float<float>{ enum {value = 1}; };


//A method of asking whether a type is an std::complex<double> or std::complex<float>
template<typename T>
struct is_complex_double_or_float{ enum {value = 0}; };

template<>
struct is_complex_double_or_float<std::complex<double> >{ enum {value = 1}; };

template<>
struct is_complex_double_or_float<std::complex<float> >{ enum {value = 1}; };

#ifdef USE_GRID
//A method of asking a Grid vector type if it's scalar_type is complex, predicated on whether the type is indeed a Grid vector type
template<bool is_grid_vector, typename T>
struct is_scalar_type_complex{};

template<typename T>
struct is_scalar_type_complex<true,T>{
  //safe to call scalar type
  enum {value = Grid::is_complex<typename T::scalar_type>::value };
};
template<typename T>
struct is_scalar_type_complex<false,T>{
  enum {value = 0};
};
  
//A method of asking whether a type is a Grid *complex* vector type
template<typename T>
struct is_grid_vector_complex{
  enum {value = is_scalar_type_complex<Grid::is_simd<T>::value, T>::value };
};
#else

template<typename T>
struct is_grid_vector_complex{
  enum {value = 0};
};

#endif


//A method of providing a list of conditions and associated types for classification
struct no_mark{};

template<bool Condition, typename IfTrue, typename NextTest>
struct TestElem{};

template<typename IfTrue, typename NextTest>
struct TestElem<true,IfTrue,NextTest>{
  typedef IfTrue type;
};
template<typename IfTrue, typename NextTest>
struct TestElem<false,IfTrue,NextTest>{
  typedef typename NextTest::type type;
};

struct LastElem{
  typedef no_mark type;
};


  
//An expandable method of classifying a complex type
struct complex_double_or_float_mark{};
struct grid_vector_complex_mark{};

template<typename T>
struct ComplexClassify{
  typedef typename TestElem< is_complex_double_or_float<T>::value, complex_double_or_float_mark,
			     TestElem< is_grid_vector_complex<T>::value, grid_vector_complex_mark,				       
				       LastElem>
			     >::type type;
};



//A version of is_base_of (copied from the internet!  https://stackoverflow.com/questions/2910979/how-does-is-base-of-work)

template <typename B, typename D>
struct _is_base_of_host
{
  operator B*() const;
  operator D*();
};

template <typename B, typename D>
struct my_is_base_of
{
  template <typename T> 
  static char check(D*, T);
  static std::pair<char,char> check(B*, int);

  static const bool value = sizeof(check(_is_base_of_host<B,D>(), int())) == sizeof(char);
};

//Generate a test for if a class has an enum with a particular name
#define define_test_has_enum(ENUM)			\
template<class T> \
struct has_enum_##ENUM{ \
  typedef char yes; \
  typedef yes (&no)[2]; \
  \
  template<int> \
  struct test2; \
  \
  template<class U> \
  static yes test(test2<U:: ENUM>*); \
  template<class U> \
  static no  test(...); \
  \
  static bool const value = sizeof(test<T>(0)) == sizeof(yes); \
}

//Example usage:
//<HEAD>
// define_test_has_enum(HELLO);

// struct Astruct{
//   enum { HELLO=0 };
// };
// struct Bstruct{
// };
//
//<BODY>
  // assert( has_enum_HELLO<Astruct>::value == true );
  // assert( has_enum_HELLO<Bstruct>::value == false );


//A method of using compound templated classes and recursive functions to iterate over types in a static list
struct ListEnd{};

template<typename V, typename N>
struct Elem{
  typedef V ValueType;
  typedef N NextType;
};

//Get a type by index
template<typename TypeList, int i>
struct getTypeFromList{
  typedef typename getTypeFromList<typename TypeList::NextType, i-1>::type type;
};
template<typename TypeList>
struct getTypeFromList<TypeList,0>{
  typedef typename TypeList::ValueType type;
};

//A struct containing instances of the list of types.
template<typename TypeList>
struct ListStruct{
  typedef typename TypeList::ValueType ValueType;
  typedef ListStruct<typename TypeList::NextType> NextType;
  
  ValueType v;
  NextType n;

  ListStruct(){}
  ListStruct(const ListStruct &r): v(r.v),n(r.n){}
};
template<>
struct ListStruct<ListEnd>{};

//Access elements by compile time index
template<typename aListStruct, int i>
struct getConstElemFromListStruct{
  static inline const typename getTypeFromList<aListStruct,i>::type & get(const aListStruct &from){
    return getConstElemFromListStruct<typename aListStruct::NextType,i-1>::get(from.n);
  }
};
template<typename aListStruct>
struct getConstElemFromListStruct<aListStruct,0>{
  static inline const typename aListStruct::ValueType & get(const aListStruct &from){
    return from.v;
  }
};
    
template<typename aListStruct, int i>
struct getElemFromListStruct{
  static inline typename getTypeFromList<aListStruct,i>::type & get(aListStruct &from){
    return getElemFromListStruct<typename aListStruct::NextType,i-1>::get(from.n);
  }
};
template<typename aListStruct>
struct getElemFromListStruct<aListStruct,0>{
  static inline typename aListStruct::ValueType & get(aListStruct &from){
    return from.v;
  }
};

//Number of elements in ListStruct
template<typename aListStruct, int count = 0>
struct getSizeOfListStruct{
  enum{ value = getSizeOfListStruct<typename aListStruct::NextType,count+1>::value };
};
template<int count>
struct getSizeOfListStruct<ListStruct<ListEnd>, count>{
  enum{ value = count };
};




//Example print all elements
template<typename aListStruct>
struct _printAll{
  static void doit(std::ostream &into, const aListStruct &src){
    into << src.v << '\n';
    typedef typename aListStruct::NextType Next;
    _printAll<Next>::doit(into,src.n);
  }
};
template<>
struct _printAll< ListStruct<ListEnd> >{
  static void doit(std::ostream &into, const ListStruct<ListEnd> &src){
  }
};

//Example print single element
template<typename aListStruct, int idx>
struct _printElement{
  static void doit(std::ostream &into, const aListStruct &src){
    typedef typename aListStruct::NextType Next;
    _printElement<Next,idx-1>::doit(into,src.n);
  }
};
template<typename aListStruct>
struct _printElement<aListStruct,0>{
  static void doit(std::ostream &into, const aListStruct &src){
    into << src.v << '\n';
  }
};

//Perform an operation with templated parameter and static member 'doit' on all elements of list
template< template<typename> class Operation, typename aListStruct>
struct _operationAll{
  static void doit(const aListStruct &src){
    typedef typename aListStruct::ValueType ValueType;
    typedef typename aListStruct::NextType NextType;
    Operation<ValueType>::doit(src.v);
    _operationAll<Operation, NextType>::doit(src.n);
  }
};
template< template<typename> class Operation>
struct _operationAll<Operation, ListStruct<ListEnd> >{
  static void doit(const ListStruct<ListEnd> &src){
  }
};

//An example operation
template<typename T>
struct _printcout{
  static void doit(const T &src){
    std::cout << src << '\n';
  }
};




#endif
