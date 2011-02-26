
enum MeasLimits { MAX_MEAS_TASK = 24 };

enum MeasType { 

	MeasAlgPlaq, 
	MeasAlgPbp, 
	MeasAlgWspect, 
	MeasAlgEig,
	MeasAlgPot,
	MeasAlgFixGauge,
	MeasAlgFixGaugeFree,
	MeasAlgQPropW,
	MeasAlgNuc3pt,
	MeasAlgRandomGauge
};	

enum MeasIOTask {
 	MeasIOLoad,
	MeasIOSave,
	MeasIONone
};
class MeasTask {

  MeasType Measurement;
  string ArgFilename<>;
  string OutputFilestem<>;

};

class MeasArg { 

  FclassType Fermion;
  GclassType Gluon;

  string WorkDirectory<>;
  string GaugeStem<>;
  string RNGStem<>;

  MeasIOTask GaugeIO;
  MeasIOTask RNGIO;

  int TrajStart;
  int TrajIncrement;
  int TrajLessThanLimit;
  int TrajCur;

  int HdwXCsum;
  int HdwRCsum;
  int IOconcurrency;

  MeasTask TaskList<>;
  
};
