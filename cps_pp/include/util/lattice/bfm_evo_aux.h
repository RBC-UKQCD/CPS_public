/* -*- mode:c++; c-basic-offset:4 -*- */
#ifndef INCLUDED_BFM_EVO_AUX_HT_H
#define INCLUDED_BFM_EVO_AUX_HT_H

namespace bfm_evo_aux {

    template<class Float>
    static inline void trless_am(Float *p, Float coef)
    {
        p[0] = p[8] = p[16] = 0.;

        Float tmp = 0.5*(p[2] - p[6]) * coef;
        p[2]=tmp; p[6] = -tmp;

        tmp = 0.5*(p[3] + p[7]) * coef;
        p[3]=tmp; p[7] = tmp;

        tmp = 0.5*(p[4] - p[12]) * coef;
        p[4]=tmp; p[12] = -tmp;

        tmp = 0.5*(p[5] + p[13]) * coef;
        p[5]=tmp; p[13] = tmp;

        tmp = 0.5*(p[10] - p[14]) * coef;
        p[10]=tmp; p[14] = -tmp;

        tmp = 0.5*(p[11] + p[15]) * coef;
        p[11]=tmp; p[15] = tmp;

        Float c = 1./3. * (p[1] + p[9] + p[17]);

        p[1] = (p[1] - c) * coef;
        p[9] = (p[9] - c) * coef;
        p[17] = (p[17] - c) * coef;
    }

    template<class Float>
    static inline void mDotMEqual(Float* c, const Float* a, const Float* b)
    {
        *c      = *a      * *b      - *(a+1)  * *(b+1)    +
            *(a+2)  * *(b+6)  - *(a+3)  * *(b+7)    +
            *(a+4)  * *(b+12) - *(a+5)  * *(b+13);
        *(c+1)  = *a      * *(b+1)  + *(a+1)  * *b        +
            *(a+2)  * *(b+7)  + *(a+3)  * *(b+6)    +
            *(a+4)  * *(b+13) + *(a+5)  * *(b+12);

        *(c+2)  = *a      * *(b+2)  - *(a+1)  * *(b+3)    +
            *(a+2)  * *(b+8)  - *(a+3)  * *(b+9)    +
            *(a+4)  * *(b+14) - *(a+5)  * *(b+15);
        *(c+3)  = *a      * *(b+3)  + *(a+1)  * *(b+2)    +
            *(a+2)  * *(b+9)  + *(a+3)  * *(b+8)    +
            *(a+4)  * *(b+15) + *(a+5)  * *(b+14);

        *(c+4)  = *a      * *(b+4)  - *(a+1)  * *(b+5)    +
            *(a+2)  * *(b+10) - *(a+3)  * *(b+11)   +
            *(a+4)  * *(b+16) - *(a+5)  * *(b+17);
        *(c+5)  = *a      * *(b+5)  + *(a+1)  * *(b+4)    +
            *(a+2)  * *(b+11) + *(a+3)  * *(b+10)   +
            *(a+4)  * *(b+17) + *(a+5)  * *(b+16);

        *(c+6)  = *(a+6)  * *b      - *(a+7)  * *(b+1)    +
            *(a+8)  * *(b+6)  - *(a+9)  * *(b+7)    +
            *(a+10) * *(b+12) - *(a+11) * *(b+13);
        *(c+7)  = *(a+6)  * *(b+1)  + *(a+7)  * *b        +
            *(a+8)  * *(b+7)  + *(a+9)  * *(b+6)    +
            *(a+10) * *(b+13) + *(a+11) * *(b+12);

        *(c+8)  = *(a+6)  * *(b+2)  - *(a+7)  * *(b+3)    +
            *(a+8)  * *(b+8)  - *(a+9)  * *(b+9)    +
            *(a+10) * *(b+14) - *(a+11) * *(b+15);
        *(c+9)  = *(a+6)  * *(b+3)  + *(a+7)  * *(b+2)    +
            *(a+8)  * *(b+9)  + *(a+9)  * *(b+8)    +
            *(a+10) * *(b+15) + *(a+11) * *(b+14);

        *(c+10) = *(a+6)  * *(b+4)  - *(a+7)  * *(b+5)    +
            *(a+8)  * *(b+10) - *(a+9)  * *(b+11)   +
            *(a+10) * *(b+16) - *(a+11) * *(b+17);
        *(c+11) = *(a+6)  * *(b+5)  + *(a+7)  * *(b+4)    +
            *(a+8)  * *(b+11) + *(a+9)  * *(b+10)   +
            *(a+10) * *(b+17) + *(a+11) * *(b+16);

        *(c+12) = *(a+12) * *b      - *(a+13) * *(b+1)    +
            *(a+14) * *(b+6)  - *(a+15) * *(b+7)    +
            *(a+16) * *(b+12) - *(a+17) * *(b+13);
        *(c+13) = *(a+12) * *(b+1)  + *(a+13) * *b        +
            *(a+14) * *(b+7)  + *(a+15) * *(b+6)    +
            *(a+16) * *(b+13) + *(a+17) * *(b+12);

        *(c+14) = *(a+12) * *(b+2)  - *(a+13) * *(b+3)    +
            *(a+14) * *(b+8)  - *(a+15) * *(b+9)    +
            *(a+16) * *(b+14) - *(a+17) * *(b+15);
        *(c+15) = *(a+12) * *(b+3)  + *(a+13) * *(b+2)    +
            *(a+14) * *(b+9)  + *(a+15) * *(b+8)    +
            *(a+16) * *(b+15) + *(a+17) * *(b+14);

        *(c+16) = *(a+12) * *(b+4)  - *(a+13) * *(b+5)    +
            *(a+14) * *(b+10) - *(a+15) * *(b+11)   +
            *(a+16) * *(b+16) - *(a+17) * *(b+17);
        *(c+17) = *(a+12) * *(b+5)  + *(a+13) * *(b+4)    +
            *(a+14) * *(b+11) + *(a+15) * *(b+10)   +
            *(a+16) * *(b+17) + *(a+17) * *(b+16);
    }

    // a += b
    template<class Float>
    static inline void su3_add(Float *a, Float *b)
    {
        for(int i = 0; i < 18; ++i) {
            a[i] += b[i];
        }
    }

    //------------------------------------------------------------------
    // sproj with (1 + gamma_0)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrXp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++) 
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) - *(tv+19)) * *(tw   )
                        + (*(tv+ 1) + *(tv+18)) * *(tw+ 1)
                        + (*(tv+ 6) - *(tv+13)) * *(tw+ 6)
                        + (*(tv+ 7) + *(tv+12)) * *(tw+ 7)
                        + (*(tv+12) + *(tv+ 7)) * *(tw+12)
                        + (*(tv+13) - *(tv+ 6)) * *(tw+13)
                        + (*(tv+18) + *(tv+ 1)) * *(tw+18)
                        + (*(tv+19) - *(tv   )) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) + *(tv+18)) * *(tw   )
                        - (*(tv   ) - *(tv+19)) * *(tw+ 1)
                        + (*(tv+ 7) + *(tv+12)) * *(tw+ 6)
                        - (*(tv+ 6) - *(tv+13)) * *(tw+ 7)
                        + (*(tv+13) - *(tv+ 6)) * *(tw+12)
                        - (*(tv+12) + *(tv+ 7)) * *(tw+13)
                        + (*(tv+19) - *(tv   )) * *(tw+18)
                        - (*(tv+18) + *(tv+ 1)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 - gamma_0)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrXm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++)
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) + *(tv+19)) * *(tw   )
                        + (*(tv+ 1) - *(tv+18)) * *(tw+ 1)
                        + (*(tv+ 6) + *(tv+13)) * *(tw+ 6)
                        + (*(tv+ 7) - *(tv+12)) * *(tw+ 7)
                        + (*(tv+12) - *(tv+ 7)) * *(tw+12)
                        + (*(tv+13) + *(tv+ 6)) * *(tw+13)
                        + (*(tv+18) - *(tv+ 1)) * *(tw+18)
                        + (*(tv+19) + *(tv   )) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) - *(tv+18)) * *(tw   )
                        - (*(tv   ) + *(tv+19)) * *(tw+ 1)
                        + (*(tv+ 7) - *(tv+12)) * *(tw+ 6)
                        - (*(tv+ 6) + *(tv+13)) * *(tw+ 7)
                        + (*(tv+13) + *(tv+ 6)) * *(tw+12)
                        - (*(tv+12) - *(tv+ 7)) * *(tw+13)
                        + (*(tv+19) + *(tv   )) * *(tw+18)
                        - (*(tv+18) - *(tv+ 1)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 + gamma_1)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrYp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    { 
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++)
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) - *(tv+18)) * *(tw   )
                        + (*(tv+ 1) - *(tv+19)) * *(tw+ 1)
                        + (*(tv+ 6) + *(tv+12)) * *(tw+ 6)
                        + (*(tv+ 7) + *(tv+13)) * *(tw+ 7)
                        + (*(tv+12) + *(tv+ 6)) * *(tw+12)
                        + (*(tv+13) + *(tv+ 7)) * *(tw+13)
                        + (*(tv+18) - *(tv   )) * *(tw+18)
                        + (*(tv+19) - *(tv+ 1)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) - *(tv+19)) * *(tw   )
                        - (*(tv   ) - *(tv+18)) * *(tw+ 1)
                        + (*(tv+ 7) + *(tv+13)) * *(tw+ 6)
                        - (*(tv+ 6) + *(tv+12)) * *(tw+ 7)
                        + (*(tv+13) + *(tv+ 7)) * *(tw+12)
                        - (*(tv+12) + *(tv+ 6)) * *(tw+13)
                        + (*(tv+19) - *(tv+ 1)) * *(tw+18)
                        - (*(tv+18) - *(tv   )) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 - gamma_1)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrYm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++)
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) + *(tv+18)) * *(tw   )
                        + (*(tv+ 1) + *(tv+19)) * *(tw+ 1)
                        + (*(tv+ 6) - *(tv+12)) * *(tw+ 6)
                        + (*(tv+ 7) - *(tv+13)) * *(tw+ 7)
                        + (*(tv+12) - *(tv+ 6)) * *(tw+12)
                        + (*(tv+13) - *(tv+ 7)) * *(tw+13)
                        + (*(tv+18) + *(tv   )) * *(tw+18)
                        + (*(tv+19) + *(tv+ 1)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) + *(tv+19)) * *(tw   )
                        - (*(tv   ) + *(tv+18)) * *(tw+ 1)
                        + (*(tv+ 7) - *(tv+13)) * *(tw+ 6)
                        - (*(tv+ 6) - *(tv+12)) * *(tw+ 7)
                        + (*(tv+13) - *(tv+ 7)) * *(tw+12)
                        - (*(tv+12) - *(tv+ 6)) * *(tw+13)
                        + (*(tv+19) + *(tv+ 1)) * *(tw+18)
                        - (*(tv+18) + *(tv   )) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 + gamma_2)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrZp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++) 
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) - *(tv+13)) * *(tw   )
                        + (*(tv+ 1) + *(tv+12)) * *(tw+ 1)
                        + (*(tv+ 6) + *(tv+19)) * *(tw+ 6)
                        + (*(tv+ 7) - *(tv+18)) * *(tw+ 7)
                        + (*(tv+12) + *(tv+ 1)) * *(tw+12)
                        + (*(tv+13) - *(tv   )) * *(tw+13)
                        + (*(tv+18) - *(tv+ 7)) * *(tw+18)
                        + (*(tv+19) + *(tv+ 6)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) + *(tv+12)) * *(tw   )
                        - (*(tv   ) - *(tv+13)) * *(tw+ 1)
                        + (*(tv+ 7) - *(tv+18)) * *(tw+ 6)
                        - (*(tv+ 6) + *(tv+19)) * *(tw+ 7)
                        + (*(tv+13) - *(tv   )) * *(tw+12)
                        - (*(tv+12) + *(tv+ 1)) * *(tw+13)
                        + (*(tv+19) + *(tv+ 6)) * *(tw+18)
                        - (*(tv+18) - *(tv+ 7)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 - gamma_2)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrZm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++) 
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) + *(tv+13)) * *(tw   )
                        + (*(tv+ 1) - *(tv+12)) * *(tw+ 1)
                        + (*(tv+ 6) - *(tv+19)) * *(tw+ 6)
                        + (*(tv+ 7) + *(tv+18)) * *(tw+ 7)
                        + (*(tv+12) - *(tv+ 1)) * *(tw+12)
                        + (*(tv+13) + *(tv   )) * *(tw+13)
                        + (*(tv+18) + *(tv+ 7)) * *(tw+18)
                        + (*(tv+19) - *(tv+ 6)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) - *(tv+12)) * *(tw   )
                        - (*(tv   ) + *(tv+13)) * *(tw+ 1)
                        + (*(tv+ 7) + *(tv+18)) * *(tw+ 6)
                        - (*(tv+ 6) - *(tv+19)) * *(tw+ 7)
                        + (*(tv+13) + *(tv   )) * *(tw+12)
                        - (*(tv+12) - *(tv+ 1)) * *(tw+13)
                        + (*(tv+19) - *(tv+ 6)) * *(tw+18)
                        - (*(tv+18) + *(tv+ 7)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 + gamma_3)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrTp(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++) 
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) + *(tv+12)) * *(tw   )
                        + (*(tv+ 1) + *(tv+13)) * *(tw+ 1)
                        + (*(tv+ 6) + *(tv+18)) * *(tw+ 6)
                        + (*(tv+ 7) + *(tv+19)) * *(tw+ 7)
                        + (*(tv+12) + *(tv   )) * *(tw+12)
                        + (*(tv+13) + *(tv+ 1)) * *(tw+13)
                        + (*(tv+18) + *(tv+ 6)) * *(tw+18)
                        + (*(tv+19) + *(tv+ 7)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) + *(tv+13)) * *(tw   )
                        - (*(tv   ) + *(tv+12)) * *(tw+ 1)
                        + (*(tv+ 7) + *(tv+19)) * *(tw+ 6)
                        - (*(tv+ 6) + *(tv+18)) * *(tw+ 7)
                        + (*(tv+13) + *(tv+ 1)) * *(tw+12)
                        - (*(tv+12) + *(tv   )) * *(tw+13)
                        + (*(tv+19) + *(tv+ 7)) * *(tw+18)
                        - (*(tv+18) + *(tv+ 6)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }


    //------------------------------------------------------------------
    // sproj with (1 - gamma_3)
    //------------------------------------------------------------------
    template<class IFloat>
    void sprojTrTm(IFloat *f, IFloat *v, IFloat *w, int num_blk,
                   int v_stride, int w_stride)
    {
        int row, col, blk ;

        IFloat *tf, *tv, *tw ;

        tf = f ;

        for (row=0; row<18; row++) 
            *tf++ = 0.0 ;

        for (blk=0; blk<num_blk; blk++) {
            tf = f ;
            tv = v ;
            tw = w ;
            for (row=0; row<3; row++) {
                for (col=0; col<3; col++) {
                    *tf++ +=   (*(tv   ) - *(tv+12)) * *(tw   )
                        + (*(tv+ 1) - *(tv+13)) * *(tw+ 1)
                        + (*(tv+ 6) - *(tv+18)) * *(tw+ 6)
                        + (*(tv+ 7) - *(tv+19)) * *(tw+ 7)
                        + (*(tv+12) - *(tv   )) * *(tw+12)
                        + (*(tv+13) - *(tv+ 1)) * *(tw+13)
                        + (*(tv+18) - *(tv+ 6)) * *(tw+18)
                        + (*(tv+19) - *(tv+ 7)) * *(tw+19) ;

                    *tf++ +=   (*(tv+ 1) - *(tv+13)) * *(tw   )
                        - (*(tv   ) - *(tv+12)) * *(tw+ 1)
                        + (*(tv+ 7) - *(tv+19)) * *(tw+ 6)
                        - (*(tv+ 6) - *(tv+18)) * *(tw+ 7)
                        + (*(tv+13) - *(tv+ 1)) * *(tw+12)
                        - (*(tv+12) - *(tv   )) * *(tw+13)
                        + (*(tv+19) - *(tv+ 7)) * *(tw+18)
                        - (*(tv+18) - *(tv+ 6)) * *(tw+19) ;

                    tw += 2 ;
                }
                tw = w ;
                tv += 2 ;
            }
            v += 24 + v_stride ;
            w += 24 + w_stride ;
        }
    }

}

#endif
