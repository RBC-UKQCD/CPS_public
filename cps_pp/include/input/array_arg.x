class FloatArray{
   memfun void resize(int n_floats);
   memfun FloatArray();
   Float Floats<>;
};

class ParamArg{
  string name<>;
  Float  val;
};

class ParamArray{
  memfun void resize(u_int num);
  memfun ParamArray();
  ParamArg params<>;
};

class IntArray {
  int v<>;
};

/*
typedef char CharArray[100];
class Filenames{
   CharArray Names<>;
};
*/
