#ifndef INCLUDED_LATTICE_TYPES_H
#define INCLUDED_LATTICE_TYPES_H

#include <config.h>

CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// The following classes have double inheritance. The virtual base 
// class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The classes below inherit from one 
// gauge class and from one fermion class. All combinations are 
// present.
//
/*! \defgroup latactions  Lattice Actions */
//------------------------------------------------------------------

//------------------------------------------------------------------
//! Trivial gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFnone  
    : public virtual Lattice, 
    public Gnone, 
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFnone();
    virtual ~GnoneFnone();
};


//------------------------------------------------------------------
//! Trivial gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFasqtad 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gnone, 
    public Fasqtad
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFasqtad();
    virtual ~GnoneFasqtad();
};

//------------------------------------------------------------------
//! Trivial gauge action with improved staggered fermion action (P4)
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFp4
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gnone, 
    public Fp4
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFp4();
    virtual ~GnoneFp4();
};

//------------------------------------------------------------------
//! Trivial gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gnone, 
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFstag();
    virtual ~GnoneFstag();
};


//------------------------------------------------------------------
//! Trivial gauge action with wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFwilson 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFwilson();
    virtual ~GnoneFwilson();
};


//------------------------------------------------------------------
//! Trivial gauge action with clover Wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFclover 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFclover();
    virtual ~GnoneFclover();
};


//------------------------------------------------------------------
//! Trivial gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFdwf 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFdwf();
    virtual ~GnoneFdwf();
};

//------------------------------------------------------------------
//! Trivial gauge action with Mobius domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFmdwf 
    : public virtual Lattice, 
    public Gnone, 
    public Fmdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFmdwf();
    virtual ~GnoneFmdwf();
};


//------------------------------------------------------------------
//! Wilson gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFnone 
//    : Public virtual Lattice, 
    : public Gwilson, 
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFnone();
    virtual ~GwilsonFnone();
};


//------------------------------------------------------------------
//! Wilson gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gwilson, 
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFstag();
    virtual ~GwilsonFstag();
};


//------------------------------------------------------------------
//! Wilson gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFasqtad 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gwilson, 
    public Fasqtad
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFasqtad();
    virtual ~GwilsonFasqtad();
};

//------------------------------------------------------------------
//! Wilson gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFp4 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gwilson, 
    public Fp4
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFp4();
    virtual ~GwilsonFp4();
};

//------------------------------------------------------------------
//! Wilson gauge action with wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFwilson 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFwilson();
    virtual ~GwilsonFwilson();
};


//------------------------------------------------------------------
//! Wilson gauge action with clover Wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFclover 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFclover();
    virtual ~GwilsonFclover();
};


//------------------------------------------------------------------
//! Wilson gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFdwf 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFdwf();
    virtual ~GwilsonFdwf();
};

//------------------------------------------------------------------
//! Wilson gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFmdwf 
    : public virtual Lattice, 
    public Gwilson, 
    public Fmdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFmdwf();
    virtual ~GwilsonFmdwf();
};


//------------------------------------------------------------------
//! Power plaquette gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFnone 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFnone();
    virtual ~GpowerPlaqFnone();
};


//------------------------------------------------------------------
//! Power plaquette gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public GpowerPlaq, 
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFstag();
    virtual ~GpowerPlaqFstag();
};


//------------------------------------------------------------------
//! Power plaquette gauge action with Wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFwilson 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFwilson();
    virtual ~GpowerPlaqFwilson();
};


//------------------------------------------------------------------
//! Power plaquette gauge action with clover PowerPlaq fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFclover 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFclover();
    virtual ~GpowerPlaqFclover();
};


//------------------------------------------------------------------
//! Power plaquette gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFdwf 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFdwf();
    virtual ~GpowerPlaqFdwf();
};

//------------------------------------------------------------------
//! Power plaquette gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFmdwf 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fmdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerPlaqFmdwf();
    virtual ~GpowerPlaqFmdwf();
};

//------------------------------------------------------------------
//! Improved rectangle gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFnone
    : public virtual Lattice,
    public GimprRect,
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFnone();
    virtual ~GimprRectFnone();
};


//------------------------------------------------------------------
//! Improved rectangle gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFstag
    : public virtual Lattice,
    public virtual FstagTypes,
    public GimprRect,
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFstag();
    virtual ~GimprRectFstag();
};


//------------------------------------------------------------------
//! Improved rectangle gauge action with wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFwilson
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFwilson();
    virtual ~GimprRectFwilson();
};


//------------------------------------------------------------------
//! Improved rectangle gauge action with clover Wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFclover
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFclover();
    virtual ~GimprRectFclover();
};


//------------------------------------------------------------------
//! Improved rectangle gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFdwf
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFdwf();
    virtual ~GimprRectFdwf();
};

//------------------------------------------------------------------
//! Improved rectangle gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFmdwf
    : public virtual Lattice,
    public GimprRect,
    public Fmdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFmdwf();
    virtual ~GimprRectFmdwf();
};

//------------------------------------------------------------------
//! Improved rectangle gauge action with P4 staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFp4 : public GimprRect, public Fp4 {
 private:
  const char *cname;    // Class name.
 
 public:
  GimprRectFp4();     
  ~GimprRectFp4();
};

//------------------------------------------------------------------
//! Power rectangle gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFnone 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerRectFnone();
    virtual ~GpowerRectFnone();
};


//------------------------------------------------------------------
//! Power rectangle gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public GpowerRect, 
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerRectFstag();
    virtual ~GpowerRectFstag();
};


//------------------------------------------------------------------
//! Power rectangle gauge action with powerRect fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFwilson 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerRectFwilson();
    virtual ~GpowerRectFwilson();
};


//------------------------------------------------------------------
//! Power rectangle gauge action with clover PowerRect fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFclover 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerRectFclover();
    virtual ~GpowerRectFclover();
};


//------------------------------------------------------------------
//! Power rectangle gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFdwf 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GpowerRectFdwf();
    virtual ~GpowerRectFdwf();
};

//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with no fermions
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFnone
    : public virtual Lattice,
    public GimprOLSym,
    public Fnone
{
 private:
    const char *cname;    // Class name.

 public:
    GimprOLSymFnone();
    virtual ~GimprOLSymFnone();
};


//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFstag
    : public virtual Lattice,
    public virtual FstagTypes,
    public GimprOLSym,
    public Fstag
{
 private:
    const char *cname;    // Class name.

 public:
    GimprOLSymFstag();
    virtual ~GimprOLSymFstag();
};


//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFwilson
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fwilson
{
 private:
    const char *cname;    // Class name.

 public:
    GimprOLSymFwilson();
    virtual ~GimprOLSymFwilson();
};


//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with clover Wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFclover
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fclover
{
 private:
    const char *cname;    // Class name.

 public:
    GimprOLSymFclover();
    virtual ~GimprOLSymFclover();
};


//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with domain wall fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFdwf
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fdwf
{
 private:
    const char *cname;    // Class name.

 public:
    GimprOLSymFdwf();
    virtual ~GimprOLSymFdwf();
};

//------------------------------------------------------------------
//! Power plaquette gauge action with Asqtad staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerPlaqFasqtad : public GpowerPlaq, public Fasqtad {
  private:
    const char *cname;    // Class name.

  public:
    GpowerPlaqFasqtad();
    ~GpowerPlaqFasqtad();
};

//------------------------------------------------------------------
//! Improved rectangle gauge action with Asqtad staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFasqtad : public GimprRect, public Fasqtad {
  private:
    const char *cname;    // Class name.

  public:
    GimprRectFasqtad();
    ~GimprRectFasqtad();
};


//------------------------------------------------------------------
//! Power rectangle gauge action with Asqtad staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GpowerRectFasqtad : public GpowerRect, public Fasqtad {
  private:
    const char *cname;    // Class name.

  public:
    GpowerRectFasqtad();
    ~GpowerRectFasqtad();
};


//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with asqtad fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFasqtad : public GimprOLSym, public Fasqtad{
    
  private:
    const char *cname;    // Class name.

  public:
    GimprOLSymFasqtad();
    ~GimprOLSymFasqtad();
};

//------------------------------------------------------------------
//! One Loop Symanzik improved gauge action with Asqtad staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprOLSymFp4 : public GimprOLSym, public Fp4{
    
  private:
    const char *cname;    // Class name.

  public:
    GimprOLSymFp4();
    ~GimprOLSymFp4();
};

//------------------------------------------------------------------
//! Tadpole-improved rectangle gauge action with p4 staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GtadpoleRectFp4 : public GtadpoleRect, public Fp4{
    
  private:
    const char *cname;    // Class name.

  public:
    GtadpoleRectFp4();
    ~GtadpoleRectFp4();
};

//------------------------------------------------------------------
//! Tadpole-improved rectangle gauge action with asqtad staggered fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GtadpoleRectFasqtad : public GtadpoleRect, public Fasqtad{
    
  private:
    const char *cname;    // Class name.

  public:
    GtadpoleRectFasqtad();
    ~GtadpoleRectFasqtad();
};

//------------------------------------------------------------------
//! Tadpole-improved rectangle gauge action with no fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GtadpoleRectFnone : public GtadpoleRect, public Fnone{
    
  private:
    const char *cname;    // Class name.

  public:
    GtadpoleRectFnone();
    ~GtadpoleRectFnone();
};

//------------------------------------------------------------------
//! Trivial gauge action with twisted-mass wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GnoneFwilsonTm
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public Gnone,
    public virtual Fwilson,
    public FwilsonTm
{
 private:
    const char *cname;    // Class name.

 public:
    GnoneFwilsonTm() {
      cname = "GnoneFwilsonTm";
      const char *fname = "GnoneFwilsonTm()";
      VRB.Func(cname,fname);
    }

    ~GnoneFwilsonTm() {
      const char *fname = "~GnoneFwilsonTm()";
      VRB.Func(cname,fname);
    }

};
//------------------------------------------------------------------
//! Trivial gauge action with twisted-mass wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GwilsonFwilsonTm
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public Gwilson,
    public virtual Fwilson,
    public FwilsonTm
{
 private:
    const char *cname;    // Class name.

 public:
    GwilsonFwilsonTm() {
      cname = "GwilsonFwilsonTm";
      const char *fname = "GwilsonFwilsonTm()";
      VRB.Func(cname,fname);
    }

    ~GwilsonFwilsonTm() {
      const char *fname = "~GwilsonFwilsonTm()";
      VRB.Func(cname,fname);
    }

};
//------------------------------------------------------------------
//! Trivial gauge action with twisted-mass wilson fermion action
/*! \ingroup latactions */
//------------------------------------------------------------------
class GimprRectFwilsonTm
: public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public virtual Fwilson,
    public FwilsonTm
{
 private:
    const char *cname;    // Class name.

 public:
    GimprRectFwilsonTm() {
      cname = "GimprRectFwilsonTm";
      const char *fname = "GimprRectFwilsonTm()";
      VRB.Func(cname,fname);
    }

    ~GimprRectFwilsonTm() {
      const char *fname = "~GimprRectFwilsonTm()";
      VRB.Func(cname,fname);
    }
    
};

CPS_END_NAMESPACE
#endif
