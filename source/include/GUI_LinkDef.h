#if defined(__ROOTCLING__)

#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

#pragma link C++ namespace Gator;
#pragma link C++ defined_in Gator;
#pragma link C++ all_function Gator;
#pragma link C++ datamember Gator;

#pragma link C++ class+protected Gator::GatorCalibClass-;
#pragma link C++ class+protected SpectraLoaderDialog-;
#pragma link C++ class+protected SaveLinesDialog-;
#pragma link C++ class+protected Gator::GatorCalibGUI-;

#endif
