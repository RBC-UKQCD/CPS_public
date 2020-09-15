#ifndef GP_MEAS_ENUMS_H
#define GP_MEAS_ENUMS_H
CPS_START_NAMESPACE

enum PropPrecision { Sloppy, Exact };
enum TbcCombination { Single, CombinationF, CombinationB }; //F=(P+A)/2  B=(P-A)/2
enum QuarkType { Light, Heavy };
enum MomentumOf { SrcPsiBar, SrcPsi, DaggeredProp, UndaggeredProp, Total };

CPS_END_NAMESPACE
#endif
