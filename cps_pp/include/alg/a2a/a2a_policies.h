#ifndef _A2A_POLICIES_H
#define _A2A_POLICIES_H

#include<alg/a2a/a2a_allocpolicies.h>

CPS_START_NAMESPACE

//Type policies needed for sources

#ifdef USE_GRID
struct GridSIMDSourcePolicies{
  typedef Grid::vComplexD ComplexType;
  typedef ThreeDSIMDPolicy DimensionPolicy;
  typedef Aligned128AllocPolicy AllocPolicy;
};
struct GridSIMDSourcePoliciesSingle{
  typedef Grid::vComplexF ComplexType;
  typedef ThreeDSIMDPolicy DimensionPolicy;
  typedef Aligned128AllocPolicy AllocPolicy;
};
#endif

struct StandardSourcePolicies{
  typedef cps::ComplexD ComplexType;
  typedef SpatialPolicy DimensionPolicy;
  typedef StandardAllocPolicy AllocPolicy;
};

struct A2ApoliciesBase{
  typedef cps::ComplexD ComplexTypeD;
  typedef cps::ComplexF ComplexTypeF;
  typedef GwilsonFdwf LatticeType;
};

#define INHERIT_A2A_POLICIES_BASE_TYPEDEFS \
  typedef typename A2ApoliciesBase::ComplexTypeD ComplexTypeD; \
  typedef typename A2ApoliciesBase::ComplexTypeF ComplexTypeF; \
  typedef typename A2ApoliciesBase::LatticeType LatticeType;


#ifdef USE_GRID

CPS_END_NAMESPACE
#include<util/lattice/fgrid.h>
CPS_START_NAMESPACE

//These typedefs are needed if Grid is being used at all even if the main program is not using SIMD vectorized data types
struct BaseGridPolicies{
# ifdef USE_GRID_GPARITY
  typedef FgridGparityMobius FgridFclass;
  typedef GnoneFgridGparityMobius FgridGFclass;
  typedef Grid::QCD::GparityMobiusFermionD GridDirac;
  typedef Grid::QCD::GparityMobiusFermionF GridDiracF; //single prec
  enum { FGRID_CLASS_NAME=F_CLASS_GRID_GPARITY_MOBIUS };
# else
  typedef FgridMobius FgridFclass;
  typedef GnoneFgridMobius FgridGFclass;
  typedef Grid::QCD::MobiusFermionD GridDirac;
  typedef Grid::QCD::MobiusFermionF GridDiracF;
  enum { FGRID_CLASS_NAME=F_CLASS_GRID_MOBIUS };
# endif  
  
  typedef typename GridDirac::FermionField GridFermionField;
  typedef typename GridDiracF::FermionField GridFermionFieldF;
};

#define INHERIT_BASE_GRID_TYPEDEFS \
  typedef typename BaseGridPolicies::FgridFclass FgridFclass; \
  typedef typename BaseGridPolicies::FgridGFclass FgridGFclass; \
  typedef typename BaseGridPolicies::GridDirac GridDirac; \
  typedef typename BaseGridPolicies::GridFermionField GridFermionField; \
  typedef typename BaseGridPolicies::GridDiracF GridDiracF; \
  typedef typename BaseGridPolicies::GridFermionFieldF GridFermionFieldF; \
  enum { FGRID_CLASS_NAME=BaseGridPolicies::FGRID_CLASS_NAME }

#endif


//Policy choices
struct A2ApoliciesDoubleAutoAlloc{
  INHERIT_A2A_POLICIES_BASE_TYPEDEFS;
#ifdef USE_GRID
  INHERIT_BASE_GRID_TYPEDEFS;
#endif
  
  typedef cps::ComplexD ComplexType;
  typedef StandardAllocPolicy AllocPolicy;
  typedef cps::ComplexD ScalarComplexType;
  typedef CPSfermion4D<ComplexType, FourDpolicy, DynamicFlavorPolicy, AllocPolicy> FermionFieldType;
  typedef CPScomplex4D<ComplexType, FourDpolicy, DynamicFlavorPolicy, AllocPolicy> ComplexFieldType;
  typedef StandardSourcePolicies SourcePolicies;

  SET_A2AVECTOR_AUTOMATIC_ALLOC(A2ApoliciesDoubleAutoAlloc);
};

struct A2ApoliciesDoubleManualAlloc{
  INHERIT_A2A_POLICIES_BASE_TYPEDEFS;
#ifdef USE_GRID
  INHERIT_BASE_GRID_TYPEDEFS;
#endif
  
  typedef cps::ComplexD ComplexType;
  typedef StandardAllocPolicy AllocPolicy;
  typedef cps::ComplexD ScalarComplexType;
  typedef CPSfermion4D<ComplexType, FourDpolicy, DynamicFlavorPolicy, AllocPolicy> FermionFieldType;
  typedef CPScomplex4D<ComplexType, FourDpolicy, DynamicFlavorPolicy, AllocPolicy> ComplexFieldType;
  typedef StandardSourcePolicies SourcePolicies;

  SET_A2AVECTOR_MANUAL_ALLOC(A2ApoliciesDoubleManualAlloc);
};


#ifdef USE_GRID

//Base stuff for SIMD policies
struct A2ApoliciesGridBase{
  typedef Grid::vComplexD ComplexTypeD;
  typedef Grid::vComplexF ComplexTypeF;
};

#define INHERIT_A2A_POLICIES_GRID_BASE_TYPEDEFS \
  typedef typename A2ApoliciesGridBase::ComplexTypeD ComplexTypeD; \
  typedef typename A2ApoliciesGridBase::ComplexTypeF ComplexTypeF;


//Policy choices
struct A2ApoliciesSIMDdoubleAutoAlloc{
  INHERIT_BASE_GRID_TYPEDEFS;
  INHERIT_A2A_POLICIES_GRID_BASE_TYPEDEFS;

  typedef FgridGFclass LatticeType;
  typedef Grid::vComplexD ComplexType;
  typedef Aligned128AllocPolicy AllocPolicy;
  typedef cps::ComplexD ScalarComplexType;
  typedef CPSfermion4D<ComplexType, FourDSIMDPolicy, DynamicFlavorPolicy, AllocPolicy> FermionFieldType;
  typedef CPScomplex4D<ComplexType, FourDSIMDPolicy, DynamicFlavorPolicy, AllocPolicy> ComplexFieldType;
  typedef GridSIMDSourcePolicies SourcePolicies;

  SET_A2AVECTOR_AUTOMATIC_ALLOC(A2ApoliciesSIMDdoubleAutoAlloc);
};

struct A2ApoliciesSIMDdoubleManualAlloc{
  INHERIT_BASE_GRID_TYPEDEFS;
  INHERIT_A2A_POLICIES_GRID_BASE_TYPEDEFS;

  typedef FgridGFclass LatticeType;
  typedef Grid::vComplexD ComplexType;
  typedef Aligned128AllocPolicy AllocPolicy;
  typedef cps::ComplexD ScalarComplexType;
  typedef CPSfermion4D<ComplexType, FourDSIMDPolicy, DynamicFlavorPolicy, AllocPolicy> FermionFieldType;
  typedef CPScomplex4D<ComplexType, FourDSIMDPolicy, DynamicFlavorPolicy, AllocPolicy> ComplexFieldType;
  typedef GridSIMDSourcePolicies SourcePolicies;

  SET_A2AVECTOR_MANUAL_ALLOC(A2ApoliciesSIMDdoubleManualAlloc);
};

#endif

CPS_END_NAMESPACE

#endif
