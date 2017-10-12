#include "GatorCalibClass.hh"

#include "BAT/BCAux.h"
#include "BAT/BCLog.h"
#include "BAT/BCHistogramFitter.h"
#include "BAT/BCParameter.h"
#include "BAT/BCH2D.h"

#include "TQObject.h"

#include "TGFrame.h"
#include "TGClient.h"
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGTab.h"
#include "TGComboBox.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGListTree.h"
#include "TGCanvas.h"
#include "TSystem.h"
#include "TGNumberEntry.h"

#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>


TString FullFileName(TGListTreeItem* item);
TString PathName(TString filefullname);
//TString DirName(TGListTreeItem* item){return FullFileName(item);};//Must be removed at some point
void InitDirsTree(TGListTree* pContents);
void InitFilesTree(TGListTree* pContents);
void BrowseNewDir(TGListTreeItem *item, TGListTree *pContents, bool rootfiles);

namespace Gator{
	
	class GatorCalibGUI: public TGMainFrame, public GatorCalib
	{
	private:
		map<string, BCHistogramFitter*> fFittersMap;//This is to check the outcome of the fits when they are accepted by the user.
		
		//Stuff of the spectra loader tab
		TGComboBox *fSpectSel; //This used in the Spectra Manager Tab
		TRootEmbeddedCanvas *fSpectrumCanvas;
		
		
		//Stuff of the Line manager tab
		TGTextEntry *fIsotMassInput;
		TGTextEntry *fElementInput;
		TGNumberEntry *fLineEnergyInput;
		TGNumberEntry *fLineEnergyErrInput;
		TGNumberEntry *fHistoLowEdgeInput;
		TGNumberEntry *fHistoUpEdgeInput;
		TGComboBox *fSpectraList; //This used in the Line Manager Tab
		
		TGListBox *fLinesSel;//This is for the selection of the lines to be saved
		
		TGListTree *fDirsContents;//This is to browse the directories and choose the root file from which the lines are loaded
		TGTextEntry *fInputLinesFile;//This is the full path name of the root file used to load the lines
		
		//Stuff of the Line fit tab
		TGComboBox *fFitLineSelBox; //This is for the selection of the line to be fitted
		TRootEmbeddedCanvas *fFitCanvas;
		TGNumberEntry *fMeanEntry, *fMeanErrEntry, *fSigmaEntry, *fSigmaErrEntry, *fBetaEntry, *fBetaErrEntry, *fAmplEntry, *fAmplErrEntry, *fTailEntry, *fTailErrEntry, *fConstEntry, *fConstErrEntry, *fStepEntry, *fStepErrEntry, *fRedChi2Entry, *fPvalEntry, *fFitRangeMinEntry, *fFitRangeMaxEntry;
		
		
		
		void MakeSpectraManagerTab(TGTab *pTabs);
		void MakeLinesManagerTab(TGTab *pTabs);
		void MakeLineFitTab(TGTab *pTabs);
		
		//void InitFilesTree(TGListTree* pContents);
		
	public:
		GatorCalibGUI();
		virtual ~GatorCalibGUI();
		
		
		// slots
		void OpenSpectrumLoadWin();
		void FillSpectraList();
		void SelectAndDrawSpect(const char* sourcename);
		void OnSingleClickFilesTree(TGListTreeItem* item, Int_t btn);
		void OnDoubleClickFilesTree(TGListTreeItem* item, Int_t btn);
		void AddLineFromGui();
		void FillLinesLists();
		void LoadLinesFromFile();
		void OpenSavingDialog();
		void SaveSelectedLines(const char* outfilename);
		void SelectAndDrawLine(const char* linename);
		void GuiLineFit();
		
#if defined(__CLING__)
		ClassDef(GatorCalibGUI,0)
#endif
	};
	
}//End of 'Gator' namespace



using namespace Gator;

class SpectraLoaderDialog: public TGTransientFrame
{
private:
	GatorCalibGUI *fGatorCalib;
	
	TGTextEntry *fSpectName;
	TGListTree *fContents;
	TGTextEntry *fSelectedDirText;
	
	//TString DirName(TGListTreeItem* item);
	
	
public:
	SpectraLoaderDialog(GatorCalibGUI *_GatorCalib);
	virtual ~SpectraLoaderDialog(){Cleanup();};

	//slots
	void OnSingleClick(TGListTreeItem* item, Int_t btn);
	void OnDoubleClick(TGListTreeItem* item, Int_t btn);
	void LoadSpectrum(); // *SIGNAL*
	
#if defined(__CLING__)
	ClassDef(SpectraLoaderDialog,0)
#endif
};



class SaveLinesDialog: public TGTransientFrame
{
private:
	GatorCalibGUI *fGatorCalib;
	
	TGListTree *fContents;
	TGTextEntry *fSelectedOutFile;
	
public:
	SaveLinesDialog(GatorCalibGUI *_GatorCalib);
	virtual ~SaveLinesDialog(){Cleanup();};
	
	
	//Slots
	void OnSingleClick(TGListTreeItem* item, Int_t btn);
	void OnDoubleClick(TGListTreeItem* item, Int_t btn);
	void SaveSelectedLines(); 
	
	void SaveSelectedLines(const char* outfilename); // *SIGNAL*
	
#if defined(__CLING__)
	ClassDef(SaveLinesDialog,0)
#endif
};


/*
class SaveLinesDialog: public TGTransientFrame
{
private:
	GatorCalibGUI *fGatorCalib;
	
	TGListTree *fContents;
	TGTextEntry *fSelectedOutFile;
	
public:
	SaveLinesDialog(GatorCalibGUI *_GatorCalib); // *SIGNAL*
	virtual ~SaveLinesDialog(){Cleanup();};
	
	//Slots
	void OnSingleClick(TGListTreeItem* item, Int_t btn);
	void OnDoubleClick(TGListTreeItem* item, Int_t btn);
	
	ClassDef(SaveLinesDialog,0)
};
*/