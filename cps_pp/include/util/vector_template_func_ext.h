CPS_START_NAMESPACE
template < typename AFloat, typename BFloat >
inline  void compDotProduct (IFloat * c_r, IFloat * c_i,
                       const AFloat * a, const BFloat * b, int len)
{
    *c_r = *c_i = 0.0;
    Float re=0.,im=0.;
#pragma omp parallel for reduction(+:re,im)
    for(int i = 0; i < len; i += 2)
    {
       re+= a[i] * b[i]     + a[i+1] * b[i+1];   // real part
       im+= a[i] * b[i+1] - a[i+1] * b[i];       // imag part
    }
    *c_r =re;
    *c_i =im;
}
CPS_END_NAMESPACE
