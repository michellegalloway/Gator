#if defined(__ROOTCLING__)

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ namespace Gator;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ all_function Gator;

#pragma link C++ class+protected Gator::GatorCalibClass-;

#pragma link C++ all_function Gator::GatorCalibClass;

#pragma link off class SPEdata;


#endif
