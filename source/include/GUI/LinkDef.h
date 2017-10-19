#if defined(__ROOTCLING__)

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link off all methods;

#pragma link C++ namespace Gator;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class+protected Gator::GatorCalibGUI+;
#pragma link C++ class+protected Gator::SpectraLoaderDialog+;
#pragma link C++ class+protected Gator::SaveLinesDialog+;

#pragma link C++ all_function Gator::GatorCalibGUI;
#pragma link C++ all_function Gator::SpectraLoaderDialog;
#pragma link C++ all_function Gator::SaveLinesDialog;


#endif
