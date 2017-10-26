#include "GatorCalibClass.hh"
#include "GUI/GatorCalibGUIClass.hh"

#include "TQObject.h"
#include "TGPicture.h"

#include "TApplication.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TCollection.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TF1.h"



#include <unistd.h>


using namespace std;


const TGPicture *gIcon=NULL;

GatorCalibGUI::GatorCalibGUI(): TGMainFrame(gClient->GetRoot(), 800, 600)
{
	gIcon = gClient->GetPicture("rootdb_t.xpm");
	
	fSpectSel = NULL;
	fSpectraList = NULL;
	
	fLinesSel = NULL;
	fFitLineSelBox = NULL;
	
	
	//Create a Tab widget
	TGCompositeFrame *pTabsTopFrame = new TGCompositeFrame(this, 800, 600);
	TGTextButton *pExitButton = new TGTextButton(this, "Exit");
	
	TGTab *pTabs = new TGTab(pTabsTopFrame, 800, 600);
	
	
	AddFrame(pTabsTopFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY, 5, 1, 5, 1));
	AddFrame(pExitButton, new TGLayoutHints(kLHintsBottom, 5, 1, 5, 1));
	pTabsTopFrame->AddFrame(pTabs, new TGLayoutHints(kLHintsBottom | kLHintsExpandX | kLHintsExpandY, 5, 5, 2, 1) );
	
	
	MakeSpectraManagerTab(pTabs);
	MakeLinesManagerTab(pTabs);
	MakeLineFitTab(pTabs);
	
	
	//TGCompositeFrame *fCorrTab = pTabs->AddTab("Correlations"); //This panel is only to plot between two parameters
	
	
	pExitButton->Connect("Clicked()", "GatorCalibGUI", this, "CloseWindow()");
	
#if !defined(__CLING__)
	pExitButton->Connect("Clicked()", "TApplication", gApplication, "Terminate()");
#endif
	
	
	// Set a name to the main frame
	SetWindowName("Gator calibration");
	
	// Map all subwindows of main frame
	MapSubwindows();
	
	// Initialize the layout algorithm
	Resize(GetDefaultSize());
	
	// Map main frame
	MapWindow();
}


/*
GatorCalibGUI::~GatorCalibGUI() {
	//Clean up used widgets: frames, buttons, layout hints
	Cleanup();
}
*/

void GatorCalibGUI::MakeSpectraManagerTab(TGTab *pTabs)
{
	//Add the tab relative to the calibration spectra
	TGCompositeFrame *fSpectTab = pTabs->AddTab("Spectra"); //This panel is to load, select and draw the calibration spectra
	
	//Fill this tab with a button on top and a permanent canvas at the bottom
	TGCompositeFrame *fSpecButtFrame = new TGCompositeFrame(fSpectTab, 900, 50, kHorizontalFrame);
	
	TGTextButton *fLoadSpectButt = new TGTextButton(fSpecButtFrame, "&Load spectrum...");
	
	TGLabel *pLabel1 = new TGLabel(fSpecButtFrame,"Selected spectrum:");
	
	fSpectSel = new TGComboBox(fSpecButtFrame);
	fSpectSel->Resize(100, 20);
	FillSpectraList();
	
	// Create canvas widget where the spectra are shown
	fSpectrumCanvas = new TRootEmbeddedCanvas("Source Spectrum", fSpectTab, 900, 600);
	
	
	//"AddFrame" section
	fSpectTab->AddFrame(fSpecButtFrame, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5, 5, 2, 5) );
	fSpecButtFrame->AddFrame(fLoadSpectButt, new TGLayoutHints(kLHintsCenterY, 1, 1, 1, 1) );
	fSpecButtFrame->AddFrame(pLabel1, new TGLayoutHints(kLHintsCenterY | kLHintsLeft,15,5,5,5));
	fSpecButtFrame->AddFrame(fSpectSel, new TGLayoutHints(kLHintsCenterY | kLHintsLeft,1,5,5,5));
	fSpectTab->AddFrame(fSpectrumCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5,1,5,10) );
	
	//"Connect" section
	fLoadSpectButt->Connect("Clicked()","GatorCalibGUI",this,"OpenSpectrumLoadWin()");
	fSpectSel->Connect("Selected(const char*)","GatorCalibGUI",this,"SelectAndDrawSpect(const char*)");
}


void GatorCalibGUI::MakeLinesManagerTab(TGTab *pTabs)
{
	//Build the layout of the line manager tab
	TGCompositeFrame *fLinesManTab = pTabs->AddTab("Lines"); //This panel is to manage the lines (add, remove, etc...)
	
	
	//Add line frame
	TGGroupFrame *pSubFrame1 = new TGGroupFrame(fLinesManTab,"Add calibration line",kVerticalFrame);
	
	TGCompositeFrame *pSubFrame1_1 = new TGCompositeFrame(pSubFrame1,(pSubFrame1->GetWidth())-2,50,kHorizontalFrame);
	TGCompositeFrame *pSubFrame1_2 = new TGCompositeFrame(pSubFrame1,(pSubFrame1->GetWidth())-2,50,kHorizontalFrame);
	
	TGLabel *pIsotMassLabel = new TGLabel(pSubFrame1_1, "Isotope mass: ");
	fIsotMassInput = new TGTextEntry(pSubFrame1_1, "");
	fIsotMassInput->SetWidth(150);
	
	TGLabel *pElementLabel = new TGLabel(pSubFrame1_1, "Element: ");
	fElementInput = new TGTextEntry(pSubFrame1_1, "");
	fElementInput->SetWidth(150);
	
	TGLabel *pSpectraListLabel = new TGLabel(pSubFrame1_1, "Spectrum: ");
	fSpectraList = new TGComboBox(pSubFrame1_1);
	fSpectraList->SetWidth(150);
	fSpectraList->SetHeight(fElementInput->GetHeight());
	
	TGLabel *pLineEnergyLabel = new TGLabel(pSubFrame1_2, "Energy (keV): ");
	fLineEnergyInput = new TGNumberEntry(pSubFrame1_2, 0);
	
	TGLabel *pLineEnergyErrLabel = new TGLabel(pSubFrame1_2, "Energy error (keV): ");
	fLineEnergyErrInput = new TGNumberEntry(pSubFrame1_2, 0);
	
	TGLabel *pHistoLowEdgeLabel = new TGLabel(pSubFrame1_2, "MCA min ch: ");
	fHistoLowEdgeInput = new TGNumberEntry(pSubFrame1_2, 0);
	
	TGLabel *pHistoUpEdgeLabel = new TGLabel(pSubFrame1_2, "MCA max ch: ");
	fHistoUpEdgeInput = new TGNumberEntry(pSubFrame1_2, 0);
	
	TGTextButton *pAddLineButt = new TGTextButton(pSubFrame1_2, "&Add");
	
	
	//Load lines frame
	TGGroupFrame *pSubFrame2 = new TGGroupFrame(fLinesManTab,"Load lines from tree",kVerticalFrame);
	
	TGLabel *pSelectRootFileLabel = new TGLabel(pSubFrame2, "Select root file");
	TGCanvas* pSelectRootFileCanvas = new TGCanvas(pSubFrame2, 2*(fElementInput->GetWidth()), 8*(fElementInput->GetHeight()) );
	fDirsContents = new TGListTree(pSelectRootFileCanvas, kHorizontalFrame);
	TGListTreeItem* item = fDirsContents->AddItem(0,"/");
	BrowseNewDir(item, fDirsContents, true);
	item->SetOpen(true);
	//InitFilesTree(fDirsContents); //This doesn't work properly
	
	TGCompositeFrame *pSubFrame2_1 = new TGCompositeFrame(pSubFrame2, (pSubFrame2->GetWidth())-2, (fElementInput->GetHeight())+2, kHorizontalFrame);
	
	fInputLinesFile = new TGTextEntry(pSubFrame2_1, "");
	fInputLinesFile->SetWidth(200);
	TGTextButton *pLoadLinesButt = new TGTextButton(pSubFrame2_1, "Load lines");
	
	
	//Save lines frame
	TGGroupFrame *pSubFrame3 = new TGGroupFrame(fLinesManTab,"Save fitted lines",kHorizontalFrame);
	
	TGLabel *pSelectLineLabel = new TGLabel(pSubFrame3, "Select line: ");
	fLinesSel = new TGListBox(pSubFrame3);
	fLinesSel->SetMultipleSelections();
	fLinesSel->SetWidth(150);
	fLinesSel->SetHeight(5*fElementInput->GetHeight());
	TGTextButton *pSaveLinesButt = new TGTextButton(pSubFrame3, "Save lines...");
	
	
	//"AddFrame" section
	fLinesManTab->AddFrame(pSubFrame1, new TGLayoutHints(kLHintsExpandX, 5,10,5,1));
	
	pSubFrame1->AddFrame(pSubFrame1_1,new TGLayoutHints(kLHintsExpandX, 1,5,1,1));
	pSubFrame1->AddFrame(pSubFrame1_2,new TGLayoutHints(kLHintsExpandX, 1,1,1,1));
	
	pSubFrame1_1->AddFrame(pIsotMassLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 1, 1));
	pSubFrame1_1->AddFrame(fIsotMassInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	pSubFrame1_1->AddFrame(pElementLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 1, 1));
	pSubFrame1_1->AddFrame(fElementInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	pSubFrame1_1->AddFrame(pSpectraListLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 2, 1));
	pSubFrame1_1->AddFrame(fSpectraList, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	
	pSubFrame1_2->AddFrame(pLineEnergyLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 1, 1));
	pSubFrame1_2->AddFrame(fLineEnergyInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	pSubFrame1_2->AddFrame(pLineEnergyErrLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 1, 1));
	pSubFrame1_2->AddFrame(fLineEnergyErrInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	pSubFrame1_2->AddFrame(pHistoLowEdgeLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 2, 1));
	pSubFrame1_2->AddFrame(fHistoLowEdgeInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	pSubFrame1_2->AddFrame(pHistoUpEdgeLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 15, 1, 2, 1));
	pSubFrame1_2->AddFrame(fHistoUpEdgeInput, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 5, 1));
	
	pSubFrame1_2->AddFrame(pAddLineButt, new TGLayoutHints(kLHintsRight|kLHintsCenterY, 10, 1, 10, 1));
	
	
	fLinesManTab->AddFrame(pSubFrame2, new TGLayoutHints(kLHintsCenterX|kLHintsTop|kLHintsExpandX, 1,1,1,1));
	
	pSubFrame2->AddFrame(pSelectRootFileLabel, new TGLayoutHints(kLHintsCenterX|kLHintsTop, 1,1,1,1));
	pSubFrame2->AddFrame(pSelectRootFileCanvas, new TGLayoutHints(kLHintsCenterX|kLHintsTop, 50,1,50,1));
	pSubFrame2->AddFrame(pSubFrame2_1, new TGLayoutHints(kLHintsCenterX|kLHintsTop|kLHintsExpandX, 1,1,1,1));
	
	pSubFrame2_1->AddFrame(fInputLinesFile, new TGLayoutHints(kLHintsLeft|kLHintsCenterY|kLHintsExpandX, 1,1,1,1));;
	pSubFrame2_1->AddFrame(pLoadLinesButt, new TGLayoutHints(kLHintsRight|kLHintsCenterY, 5,1,5,1));
	
	
	fLinesManTab->AddFrame(pSubFrame3, new TGLayoutHints(kLHintsExpandX, 5, 1, 5, 1));
	pSubFrame3->AddFrame(pSelectLineLabel, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
	pSubFrame3->AddFrame(fLinesSel, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
	pSubFrame3->AddFrame(pSaveLinesButt, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
	
	
	
	//"Connect" section
	pAddLineButt->Connect("Clicked()", "GatorCalibGUI", this, "AddLineFromGui()");
	pLoadLinesButt->Connect("Clicked()", "GatorCalibGUI", this, "LoadLinesFromFile()");
	pSaveLinesButt->Connect("Clicked()", "GatorCalibGUI", this, "OpenSavingDialog()");
	
	fDirsContents->Connect("Clicked(TGListTreeItem*,Int_t)","GatorCalibGUI",this, "OnSingleClickFilesTree(TGListTreeItem*,Int_t)");
	fDirsContents->Connect("DoubleClicked(TGListTreeItem*,Int_t)","GatorCalibGUI",this, "OnDoubleClickFilesTree(TGListTreeItem*,Int_t)");
	
}


void GatorCalibGUI::MakeLineFitTab(TGTab *pTabs)
{
	TGCompositeFrame *fLineFitTab = pTabs->AddTab("Fits"); //This panel is to fit the lines
	
	//This actually goes on top of everything but I need to get the dimensions from above
	TGCompositeFrame *pSubFrame0 = new TGCompositeFrame(fLineFitTab, (fLineFitTab->GetWidth())-2, 50, kHorizontalFrame);
	TGLabel *pSelLineLabel = new TGLabel(pSubFrame0, "Select a line: ");
	fFitLineSelBox = new TGComboBox(pSubFrame0);
	fFitLineSelBox->SetWidth(150);
	fFitLineSelBox->SetHeight(pSelLineLabel->GetHeight());
	
	TGGroupFrame *pSubFrame1 = new TGGroupFrame(fLineFitTab,"Line fit setup: ",kVerticalFrame);
	TGCompositeFrame *pSubFrame1_1 = new TGCompositeFrame(pSubFrame1,(pSubFrame1->GetWidth())-2,150,kHorizontalFrame);
	TGCompositeFrame *pSubFrame1_2 = new TGCompositeFrame(pSubFrame1,(pSubFrame1->GetWidth())-2,150,kHorizontalFrame);
	TGCompositeFrame *pFitRangeFrame = new TGCompositeFrame(pSubFrame1, (pSubFrame1->GetWidth())-2, 150, kHorizontalFrame);
	TGCompositeFrame *pFitRangeSubFrame1 = new TGCompositeFrame(pFitRangeFrame, (pFitRangeFrame->GetWidth())-2, 150, kVerticalFrame);
	TGCompositeFrame *pFitRangeSubFrame1_1 = new TGCompositeFrame(pFitRangeSubFrame1, (pFitRangeSubFrame1->GetWidth())-2, 150, kHorizontalFrame);
	
	TGCompositeFrame *pMeanSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pMeanErrSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pSigmaSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pSigmaErrSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pBetaSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pBetaErrSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pAmplSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pAmplErrSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pTailSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	TGCompositeFrame *pTailErrSubFrame = new TGCompositeFrame(pSubFrame1_1,kVerticalFrame);
	
	
	TGCompositeFrame *pConstSubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	TGCompositeFrame *pConstErrSubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	TGCompositeFrame *pStepSubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	TGCompositeFrame *pStepErrSubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	TGCompositeFrame *pRedChi2SubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	TGCompositeFrame *pPvalSubFrame = new TGCompositeFrame(pSubFrame1_2,kVerticalFrame);
	
	
	TGTextButton *pFitLineButton = new TGTextButton(pFitRangeFrame, "Fit!");
	
	TGLabel *pFitRangeLabel = new TGLabel(pFitRangeSubFrame1, "Fit range");
	
	TGLabel *pFitRangeMinLabel = new TGLabel(pFitRangeSubFrame1_1, "Min: ");
	fFitRangeMinEntry = new TGNumberEntry(pFitRangeSubFrame1_1,0);
	TGLabel *pFitRangeMaxLabel = new TGLabel(pFitRangeSubFrame1_1, "Max: ");
	fFitRangeMaxEntry = new TGNumberEntry(pFitRangeSubFrame1_1,0);
	
	
	
	TGLabel *pMeanLabel = new TGLabel(pMeanSubFrame,"Mean");
	fMeanEntry = new TGNumberEntry(pMeanSubFrame,0);
	TGLabel *pMeanErrLabel = new TGLabel(pMeanErrSubFrame,"Mean error");
	fMeanErrEntry = new TGNumberEntry(pMeanErrSubFrame,0);
	fMeanErrEntry->SetState(false);
	
	TGLabel *pSigmaLabel = new TGLabel(pSigmaSubFrame,"Sigma");
	fSigmaEntry = new TGNumberEntry(pSigmaSubFrame,0);
	TGLabel *pSigmaErrLabel = new TGLabel(pSigmaErrSubFrame,"Sigma error");
	fSigmaErrEntry = new TGNumberEntry(pSigmaErrSubFrame,0);
	fSigmaErrEntry->SetState(false);
	
	TGLabel *pBetaLabel = new TGLabel(pBetaSubFrame,"Beta");
	fBetaEntry = new TGNumberEntry(pBetaSubFrame,0);
	TGLabel *pBetaErrLabel = new TGLabel(pBetaErrSubFrame,"Beta error");
	fBetaErrEntry = new TGNumberEntry(pBetaErrSubFrame,0);
	fBetaErrEntry->SetState(false);
	
	TGLabel *pAmplLabel = new TGLabel(pAmplSubFrame,"Ampl");
	fAmplEntry = new TGNumberEntry(pAmplSubFrame,0);
	TGLabel *pAmplErrLabel = new TGLabel(pAmplErrSubFrame,"Ampl error");
	fAmplErrEntry = new TGNumberEntry(pAmplErrSubFrame,0);
	fAmplErrEntry->SetState(false);
	
	TGLabel *pTailLabel = new TGLabel(pTailSubFrame,"Tail");
	fTailEntry = new TGNumberEntry(pTailSubFrame,0);
	TGLabel *pTailErrLabel = new TGLabel(pTailErrSubFrame,"Tail error");
	fTailErrEntry = new TGNumberEntry(pTailErrSubFrame,0);
	fTailErrEntry->SetState(false);
	
	TGLabel *pConstLabel = new TGLabel(pConstSubFrame,"Const");
	fConstEntry = new TGNumberEntry(pConstSubFrame,0);
	fConstEntry->SetState(false);
	TGLabel *pConstErrLabel = new TGLabel(pConstErrSubFrame,"Const error");
	fConstErrEntry = new TGNumberEntry(pConstErrSubFrame,0);
	fConstErrEntry->SetState(false);
	
	TGLabel *pStepLabel = new TGLabel(pStepSubFrame,"Step");
	fStepEntry = new TGNumberEntry(pStepSubFrame,0);
	fStepEntry->SetState(false);
	TGLabel *pStepErrLabel = new TGLabel(pStepErrSubFrame,"Step error");
	fStepErrEntry = new TGNumberEntry(pStepErrSubFrame,0);
	fStepErrEntry->SetState(false);
	
	TGLabel *pRedChi2Label = new TGLabel(pRedChi2SubFrame,"Red Chi2");
	fRedChi2Entry = new TGNumberEntry(pRedChi2SubFrame,0);
	fRedChi2Entry->SetState(false);
	
	TGLabel *pPvalLabel = new TGLabel(pPvalSubFrame,"P-value");
	fPvalEntry = new TGNumberEntry(pPvalSubFrame,0);
	fPvalEntry->SetState(false);
	
	
	
	//"AddFrame" section
	TGLayoutHints* pSimpleFrameLoh = new TGLayoutHints(kLHintsExpandX, 1, 10, 1, 5);
	TGLayoutHints* pMiddleSubFrameLoh = new TGLayoutHints(kLHintsExpandX, 1,2,1,2);
	TGLayoutHints* pLeftSubFrameLoh = new TGLayoutHints(kLHintsExpandX | kLHintsLeft, 1,1,1,2);
	TGLayoutHints* pTopWidgetsLoh = new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 1, 5);
	TGLayoutHints* pBoxFrameLoh = new TGLayoutHints(kLHintsTop|kLHintsCenterX|kLHintsExpandX, 5, 10, 5, 5);
	
	TGLayoutHints* pValueBoxLoh = new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 1, 1, 1);
	TGLayoutHints* pWidgetLoh = new TGLayoutHints(kLHintsCenterX|kLHintsCenterY, 1, 1, 1, 1);
	TGLayoutHints* pBotWidgetLoh = new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 1, 1, 1, 1);
	
	
	fLineFitTab->AddFrame(pSubFrame0, pSimpleFrameLoh);
	pSubFrame0->AddFrame(pSelLineLabel, pTopWidgetsLoh);
	pSubFrame0->AddFrame(fFitLineSelBox, pTopWidgetsLoh);
	
	
	fLineFitTab->AddFrame(pSubFrame1, pBoxFrameLoh);
	
	pSubFrame1->AddFrame(pSubFrame1_1, new TGLayoutHints(kLHintsExpandX, 1,5,1,2));
	pSubFrame1->AddFrame(pSubFrame1_2, pMiddleSubFrameLoh);
	pSubFrame1->AddFrame(pFitRangeFrame, pMiddleSubFrameLoh);
	
	
	pSubFrame1_1->AddFrame(pMeanSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pMeanErrSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pSigmaSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pSigmaErrSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pBetaSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pBetaErrSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pAmplSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pAmplErrSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pTailSubFrame, pValueBoxLoh);
	pSubFrame1_1->AddFrame(pTailErrSubFrame, pValueBoxLoh);
	
	
	pSubFrame1_2->AddFrame(pConstSubFrame, pValueBoxLoh);
	pSubFrame1_2->AddFrame(pConstErrSubFrame, pValueBoxLoh);
	pSubFrame1_2->AddFrame(pStepSubFrame, pValueBoxLoh);
	pSubFrame1_2->AddFrame(pStepErrSubFrame, pValueBoxLoh);
	pSubFrame1_2->AddFrame(pRedChi2SubFrame, pValueBoxLoh);
	pSubFrame1_2->AddFrame(pPvalSubFrame, pValueBoxLoh);
	
	
	pFitRangeFrame->AddFrame(pFitRangeSubFrame1, pLeftSubFrameLoh);
	pFitRangeFrame->AddFrame(pFitLineButton, new TGLayoutHints(kLHintsRight|kLHintsCenterY, 1, 1, 1, 1));
	
	pFitRangeSubFrame1->AddFrame(pFitRangeLabel, pWidgetLoh);
	pFitRangeSubFrame1->AddFrame(pFitRangeSubFrame1_1, new TGLayoutHints(kLHintsExpandX, 1,1,1,2));
	
	pFitRangeSubFrame1_1->AddFrame(pFitRangeMinLabel, pBotWidgetLoh);
	pFitRangeSubFrame1_1->AddFrame(fFitRangeMinEntry, pBotWidgetLoh);
	pFitRangeSubFrame1_1->AddFrame(pFitRangeMaxLabel, pBotWidgetLoh);
	pFitRangeSubFrame1_1->AddFrame(fFitRangeMaxEntry, pBotWidgetLoh);
	
	
	pMeanSubFrame->AddFrame(pMeanLabel, pWidgetLoh);
	pMeanSubFrame->AddFrame(fMeanEntry, pWidgetLoh);
	
	pMeanErrSubFrame->AddFrame(pMeanErrLabel, pWidgetLoh);
	pMeanErrSubFrame->AddFrame(fMeanErrEntry, pWidgetLoh);
	
	pSigmaSubFrame->AddFrame(pSigmaLabel, pWidgetLoh);
	pSigmaSubFrame->AddFrame(fSigmaEntry, pWidgetLoh);
	
	pSigmaErrSubFrame->AddFrame(pSigmaErrLabel, pWidgetLoh);
	pSigmaErrSubFrame->AddFrame(fSigmaErrEntry, pWidgetLoh);
	
	pBetaSubFrame->AddFrame(pBetaLabel, pWidgetLoh);
	pBetaSubFrame->AddFrame(fBetaEntry, pWidgetLoh);
	
	pBetaErrSubFrame->AddFrame(pBetaErrLabel, pWidgetLoh);
	pBetaErrSubFrame->AddFrame(fBetaErrEntry, pWidgetLoh);
	
	pAmplSubFrame->AddFrame(pAmplLabel, pWidgetLoh);
	pAmplSubFrame->AddFrame(fAmplEntry, pWidgetLoh);
	
	pAmplErrSubFrame->AddFrame(pAmplErrLabel, pWidgetLoh);
	pAmplErrSubFrame->AddFrame(fAmplErrEntry, pWidgetLoh);
	
	pTailSubFrame->AddFrame(pTailLabel, pWidgetLoh);
	pTailSubFrame->AddFrame(fTailEntry, pWidgetLoh);
	
	pTailErrSubFrame->AddFrame(pTailErrLabel, pWidgetLoh);
	pTailErrSubFrame->AddFrame(fTailErrEntry, pWidgetLoh);
	
	pConstSubFrame->AddFrame(pConstLabel, pWidgetLoh);
	pConstSubFrame->AddFrame(fConstEntry, pWidgetLoh);
	
	pConstErrSubFrame->AddFrame(pConstErrLabel, pWidgetLoh);
	pConstErrSubFrame->AddFrame(fConstErrEntry, pWidgetLoh);
	
	pStepSubFrame->AddFrame(pStepLabel, pWidgetLoh);
	pStepSubFrame->AddFrame(fStepEntry, pWidgetLoh);
	
	pStepErrSubFrame->AddFrame(pStepErrLabel, pWidgetLoh);
	pStepErrSubFrame->AddFrame(fStepErrEntry, pWidgetLoh);
	
	pRedChi2SubFrame->AddFrame(pRedChi2Label, pWidgetLoh);
	pRedChi2SubFrame->AddFrame(fRedChi2Entry, pWidgetLoh);
	
	pPvalSubFrame->AddFrame(pPvalLabel, pWidgetLoh);
	pPvalSubFrame->AddFrame(fPvalEntry, pWidgetLoh);
	
	
	
	//This is constructed here after all the other frames have been placed
	fFitCanvas = new TRootEmbeddedCanvas("Source Spectrum", fLineFitTab, (fLineFitTab->GetWidth())-10, (fLineFitTab->GetHeight())-(pSubFrame0->GetHeight()+pSubFrame1->GetHeight()) );
	
	fLineFitTab->AddFrame(fFitCanvas, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY, 5,5,5,5));
	
	
	//"Connect" section
	fFitLineSelBox->Connect("Selected(const char*)","GatorCalibGUI",this,"SelectAndDrawLine(const char*)");
	pFitLineButton->Connect("Clicked()","GatorCalibGUI",this,"GuiLineFit()");
	
}


void GatorCalibGUI::OpenSpectrumLoadWin()
{
	if( IsDebug() ){
		cout << "Debug---> GatorCalibMainFrame::OpenSpectrumLoadWin(): activated." << endl;
	}
	
	SpectraLoaderDialog *dialog = new SpectraLoaderDialog(this);
}


void GatorCalibGUI::FillSpectraList()
{
	if(fSpectSel) fSpectSel->RemoveAll();
	if(fSpectraList) fSpectraList->RemoveAll();
	
	map<string, TH1D*>::iterator It;
	for(It=fSpectramap.begin(); It!=fSpectramap.end(); It++)
	{
		if(fSpectSel) fSpectSel->AddEntry((It->first).c_str(), -1);
		if(fSpectraList) fSpectraList->AddEntry((It->first).c_str(), -1);
	}
}


void GatorCalibGUI::SelectAndDrawSpect(const char* sourcename)
{
	if(fDebug){
		cout << "Debug---> GatorCalibGUI::SelectAndDrawSpect(...): Executed with sourcename=" << sourcename << endl;
	}
	
	fSpectrumCanvas->GetCanvas()->cd();
	
	if(SelectSpectrum(sourcename))
	{
		DrawSpectrum();
	}
	fSpectrumCanvas->GetCanvas()->Update();
}


void GatorCalibGUI::AddLineFromGui()
{
	//Make line name
	stringstream lName; lName.str("");
	lName << fIsotMassInput->GetText() << fElementInput->GetText() << "_" << ((int)(fLineEnergyInput->GetNumber()+0.5)) << "keV";
	
	if(fLinesmap.find(lName.str())!=fLinesmap.end())
	{
		if(fLinesmap[lName.str()]) delete fLinesmap[lName.str()];
	}
	
	if(!fSpectraList->GetSelectedEntry()) return;
	
	string spectrumname = ((TGTextLBEntry*)fSpectraList->GetSelectedEntry())->GetText()->GetString();
	
	if(fSpectramap.find(spectrumname)==fSpectramap.end()) return;
	
	TH1D* spectrum = fSpectramap[spectrumname];
	
	int firstbin = spectrum->FindBin(fHistoLowEdgeInput->GetNumber());
	int lastbin = spectrum->FindBin(fHistoUpEdgeInput->GetNumber());
	int nbins = 1 + lastbin - firstbin;
	double xmin = spectrum->GetBinLowEdge(firstbin);
	double xmax = spectrum->GetBinLowEdge(lastbin+1);
	
	fLinesmap[lName.str()] = new CalibLine(lName.str(), fIsotMassInput->GetText(), fElementInput->GetText(), fLineEnergyInput->GetNumber(), fLineEnergyErrInput->GetNumber());
	
	stringstream ss_histoname; ss_histoname.str(""); ss_histoname << ((int)(fLineEnergyInput->GetNumber()+0.5)) << "keV" ;
	fLinesmap[lName.str()]->histo = new TH1D(ss_histoname.str().c_str(),";Channel;Counts",nbins,xmin,xmax);
	fLinesmap[lName.str()]->histo -> SetDirectory(0);
	
	for(int iBin=firstbin; iBin<=lastbin; iBin++)
	{
		fLinesmap[lName.str()]->histo->SetBinContent( iBin-firstbin+1, spectrum->GetBinContent(iBin) );
	}
	
	fLinesmap[lName.str()]->fit = new TF1("ff_MCA", &GatorCalib::peakFitFuncB, fHistoLowEdgeInput->GetNumber(), fHistoUpEdgeInput->GetNumber(), 7);
	
	//fLinesSel->AddEntry(lName.str().c_str(),-1);
	//fFitLineSelBox->AddEntry(lName.str().c_str(), -1);
	
	if(fFittersMap.find(lName.str())!=fFittersMap.end())
	{
		if(fFittersMap[lName.str()])
		{
			delete fFittersMap[lName.str()];
		}
	}
	fFittersMap[lName.str()] = NULL;
	
	FillLinesLists();
}


void GatorCalibGUI::FillLinesLists()
{
	//Clear the lists and remake it
	if(fLinesSel) fLinesSel->RemoveAll();
	if(fFitLineSelBox) fFitLineSelBox->RemoveAll();
	
	map<string, CalibLine*>::iterator It;
	int id=0;
	for(It=fLinesmap.begin(); It!=fLinesmap.end(); It++)
	{
		if(fLinesSel) fLinesSel->AddEntry((It->first).c_str(), id);
		if(fFitLineSelBox) fFitLineSelBox->AddEntry((It->first).c_str(), id);
		id++;
	}
	Layout();
	//fClient->NeedRedraw(fFitLineSelBox);
}


void GatorCalibGUI::OnSingleClickFilesTree(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item) return;
	
	TString filefullname = FullFileName(item);
	TString pathname = PathName(filefullname);
	TSystemFile selectedfile(filefullname, pathname);
	
	
	if(selectedfile.IsDirectory())
	{
		if(IsDebug())
		{
			TSystemDirectory dir( filefullname.Data(), filefullname.Data() );
			TList *files = dir.GetListOfFiles();
			cout << "Debug---> GatorCalibGUI::OnSingleClickFilesTree(...): <" << filefullname << "> is a directory containing " << files->GetEntries() << " files." << endl;
		}
		
		if( !(bool)item->GetUserData() )
		{
			if(IsDebug())
			{
				cout << "Debug---> GatorCalibGUI::OnSingleClickFilesTree(...): The directory <" << filefullname << "> was not browsed before." << endl;
			}
			BrowseNewDir(item, fDirsContents, true);
			item->SetOpen(true);
		}
	}
	
	item->SetOpen(true);
	
	fInputLinesFile->SetText(filefullname.Data(), true);
	Layout();
}


void GatorCalibGUI::OnDoubleClickFilesTree(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item) return;
	
	// use UserData to indicate that item was already browsed
	
	TString filefullname = FullFileName(item);
	TString pathname = PathName(filefullname);
	TSystemFile selectedfile(filefullname, pathname);
	
	if(IsDebug()){
		cout << "\nDebug---> GatorCalibGUI::OnDoubleClickFilesTree(...): item->GetText() = " << item->GetText() << endl;
		cout << "Debug---> GatorCalibGUI::OnDoubleClickFilesTree(...): FullFileName(item) = " << FullFileName(item) << endl;
		cout << "Debug---> GatorCalibGUI::OnDoubleClickFilesTree(...): selectedfile.GetName() = " << selectedfile.GetName() << endl;
		cout << "Debug---> GatorCalibGUI::OnDoubleClickFilesTree(...): selectedfile.GetTitle() = " << selectedfile.GetTitle() << endl;
	}
	
	if(filefullname.EndsWith(".root"))
	{
		fInputLinesFile->SetText(filefullname.Data(), true);
		LoadLinesFromFile();
	}
}


void GatorCalibGUI::LoadLinesFromFile()
{
	TString filename = fInputLinesFile->GetText();
	
	if( (filename=="") || (!filename.EndsWith(".root")) ) return;
	
	if( LoadLinesFromTree( string(filename.Data()) ) )
	{
		FillLinesLists();
	}
}


void GatorCalibGUI::SelectAndDrawLine(const char* linename)
{
	fFitCanvas->GetCanvas()->cd();
	
	if(fLinesmap.find(string(linename))==fLinesmap.end()) return;
	
	if( (!fLinesmap[string(linename)]) || (!fLinesmap[string(linename)]->histo) ) return;
	
	fLinesmap[string(linename)]->histo->Draw();
	fSpectrumCanvas->GetCanvas()->Modified();
	fSpectrumCanvas->GetCanvas()->Update();
	Layout();
}


void GatorCalibGUI::GuiLineFit()
{
	if(!fFitLineSelBox->GetSelectedEntry()) return;
	
	string selLine = ((TGTextLBEntry*)fFitLineSelBox->GetSelectedEntry())->GetText()->GetString();
	
	if(fLinesmap.find(selLine)==fLinesmap.end()) return;
	
	CalibLine *line = fLinesmap[selLine];
	
	if(!line) return;
	
	if(IsDebug())
	{
		cout << "\nDebug ---> GatorCalibGUI::GuiLineFit(): line name = " << line->linename << endl;
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): line histo address = " << line->histo << endl;
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): line fit address = " << line->fit << endl;
	}
	
	if(!( line->histo && line->fit )) return;
	
	if(IsDebug())
	{
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): initialising the fit parameters.... ";
	}
	
	line->FitInit( fMeanEntry->GetNumber(), fSigmaEntry->GetNumber(), fBetaEntry->GetNumber(), fAmplEntry->GetNumber(), fTailEntry->GetNumber());
	
	if(IsDebug())
	{
		cout << " done" << endl;
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): setting the fit range. " << endl;
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): value of fFitRangeMinEntry = " << fFitRangeMinEntry->GetNumber() << endl;
		cout << "Debug ---> GatorCalibGUI::GuiLineFit(): value of fFitRangeMaxEntry = " << fFitRangeMaxEntry->GetNumber() << endl;
	}
	
	line->SetFitRange( fFitRangeMinEntry->GetNumber(), fFitRangeMaxEntry->GetNumber() );
	
	
	fFitCanvas->GetCanvas()->cd();
	
	BCHistogramFitter *histofitter = GatorCalib::FitLine(*line);
	
	if(!histofitter)
	{
		cerr << "\nWARNING ---> GatorCalibGUI::GuiLineFit(): The fit for line <" << line->linename << "> failed!" << endl;
		return;
	}
	
	if(fFittersMap.find(selLine)!=fFittersMap.end())
	{
		if(fFittersMap[selLine]) delete fFittersMap[selLine];
	}
	fFittersMap[selLine] = histofitter;
	
	
	//Copy the parameters into the numerical input widgets
	
	fMeanEntry->SetNumber(line->mean);
	fMeanErrEntry->SetNumber(line->mean_err);
	
	fSigmaEntry->SetNumber(line->sigma);
	fSigmaErrEntry->SetNumber(line->sigma_err);
	
	fBetaEntry->SetNumber(line->beta);
	fBetaErrEntry->SetNumber(line->beta_err);
	
	fAmplEntry->SetNumber(line->ampl);
	fAmplErrEntry->SetNumber(line->ampl_err);
	
	fTailEntry->SetNumber(line->tail);
	fTailErrEntry->SetNumber(line->tail_err);
	
	fConstEntry->SetNumber(line->cost);
	fConstErrEntry->SetNumber(line->cost_err);
	
	fStepEntry->SetNumber(line->step);
	fStepErrEntry->SetNumber(line->step_err);
	
	fRedChi2Entry->SetNumber(line->chi2ndof);
	fPvalEntry->SetNumber(line->p_value_ndof);
}


void GatorCalibGUI::OpenSavingDialog()
{
	SaveLinesDialog *dialog = new SaveLinesDialog(this);
}


void GatorCalibGUI::SaveSelectedLines(const char* outfilename)
{
	
	map<string, CalibLine*> oldLinesMap = fLinesmap;
	
	set<string> fSaveLinesList;
	
	TList *pSelected = new TList;
	
	fLinesSel->GetSelectedEntries(pSelected);
	
	TIter next(pSelected);
	TGTextLBEntry *pEntry;
	
	while( (pEntry=(TGTextLBEntry*)next()) )
	{
		fSaveLinesList.insert( string(pEntry->GetText()->GetString()) );
	}
	
	if( IsDebug() )
	{
		cout << "Debug---> GatorCalibGUI::SaveSelectedLines(...): Number of lines to be saved: " << fSaveLinesList.size() << endl;
	}
	
	map<string, CalibLine*>::iterator It;
	for(It=fLinesmap.begin(); It!=fLinesmap.end(); It++)
	{
		if( fSaveLinesList.find(It->first)==fSaveLinesList.end() )
		{
			fLinesmap.erase(It->first);
		}
	}
	
	
	if( IsDebug() )
	{
		cout << "Debug---> GatorCalibGUI::SaveSelectedLines(...): fLinesmap.size() = " << fLinesmap.size() << endl;
		cout << "Debug---> GatorCalibGUI::SaveSelectedLines(...): Saving the lines in the file <" << outfilename << ">" << endl;
	}
	
	
	SaveLines(string(outfilename));
	
	fLinesmap = oldLinesMap;
}




//------------------------------------------------//
//  Definitions of the SpectraLoaderDialog class  //
//------------------------------------------------//

//using namespace Gator;

SpectraLoaderDialog::SpectraLoaderDialog(GatorCalibGUI *_GatorCalib):
	TGTransientFrame(gClient->GetRoot(), _GatorCalib, 300, 200, kVerticalFrame)
{
	
	fGatorCalib = _GatorCalib;
	
	SetCleanup(kDeepCleanup);
	Connect("CloseWindow()", "SpectraLoaderDialog", this, "CloseWindow()");
	DontCallClose();
	
	if(!_GatorCalib)
	{
		MakeZombie();
		return;
	}
	
	
	
	TGCompositeFrame *pTopFrame = new TGCompositeFrame(this, 600, 100, kHorizontalFrame);
	AddFrame(pTopFrame, new TGLayoutHints(kLHintsExpandX, 5, 1, 5, 5));
	
	TGLabel *pLabel = new TGLabel(pTopFrame, "Source name");
	fSpectName = new TGTextEntry(pTopFrame, "");
	
	pTopFrame->AddFrame(pLabel, new TGLayoutHints(kLHintsCenterY, 5, 1, 5, 1) );
	pTopFrame->AddFrame(fSpectName, new TGLayoutHints(kLHintsCenterY, 5, 1, 5, 1));
	
	
	//Make the list of directories to be selected
	TGCanvas* tgcanvas = new TGCanvas(this, 600, 200);
	fContents = new TGListTree(tgcanvas, kHorizontalFrame);
	TGListTreeItem* item = fContents->AddItem(0,"/");
	BrowseNewDir(item, fContents, true);
	item->SetOpen(true);
	//InitDirsTree(fContents);
	
	
	
	TGCompositeFrame *pBottomFrame = new TGCompositeFrame(this, 600, 100, kHorizontalFrame);
	
	
	TGLabel *pLabel2 = new TGLabel(pBottomFrame, "Selected directory: ");
	
	
	fSelectedDirText = new TGTextEntry(pBottomFrame, "");
	fSelectedDirText->SetMinWidth(150);
	
	//fSelectedDirLabel->ChangeOptions( fSelectedDirLabel->GetOptions() | kSunkenFrame);
	
	
	
	
	//Bottom frame for OK and Cancel bottons
	TGCompositeFrame *pCommandsFrame = new TGCompositeFrame(this, 600, 100, kHorizontalFrame);
	
	
	TGButton *pOkButt = new TGTextButton(pCommandsFrame, "O&k");
	TGButton *pCancelButt = new TGTextButton(pCommandsFrame, "&Cancel");
	
	
	
	
	AddFrame(tgcanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY| kLHintsCenterX, 20, 1, 20, 10));
	AddFrame(pBottomFrame, new TGLayoutHints(kLHintsExpandX, 5, 1, 5, 5));
	AddFrame(pCommandsFrame, new TGLayoutHints(kLHintsLeft, 5, 1, 5, 5));
	
	pBottomFrame->AddFrame(pLabel2, new TGLayoutHints(kLHintsNormal, 5, 1, 5, 1) );
	pBottomFrame->AddFrame(fSelectedDirText, new TGLayoutHints(kLHintsCenterY|kLHintsExpandX, 5, 1, 5, 1) );
	
	pCommandsFrame->AddFrame(pOkButt, new TGLayoutHints(kLHintsLeft, 5, 1, 5, 1));
	pCommandsFrame->AddFrame(pCancelButt, new TGLayoutHints(kLHintsLeft, 5, 1, 5, 1));
	
	
	
	
	
	fContents->Connect("Clicked(TGListTreeItem*,Int_t)","SpectraLoaderDialog",this, "OnSingleClick(TGListTreeItem*,Int_t)");
	fContents->Connect("DoubleClicked(TGListTreeItem*,Int_t)","SpectraLoaderDialog",this, "OnDoubleClick(TGListTreeItem*,Int_t)");
	
	pOkButt->Connect("Clicked()", "SpectraLoaderDialog", this, "LoadSpectrum()");
	pOkButt->Connect("Clicked()", "SpectraLoaderDialog", this, "CloseWindow()");
	//Connect("LoadSpectrum()", "GatorCalibGUI", fGatorCalib, "FillSpectraList()");
	pCancelButt->Connect("Clicked()", "SpectraLoaderDialog", this, "CloseWindow()");
	
	
	
	SetWindowName("Load spectrum dialog.");
	MapSubwindows();
	Resize(GetDefaultSize());
	CenterOnParent(kTRUE, TGTransientFrame::kRight);
	
	MapWindow();
}


void SpectraLoaderDialog::OnSingleClick(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item ) return;
	
	TString dirfullpath = FullFileName(item);
	
	TSystemDirectory dir(dirfullpath.Data(), dirfullpath.Data());
	
	
	if(fGatorCalib->IsDebug()){
		cout << "Debug---> SpectraLoaderDialog::OnSingleClick: item->GetText() = " << item->GetText() << endl;
		cout << "Debug---> SpectraLoaderDialog::OnSingleClick: FullFileName(item) = " << FullFileName(item) << endl;
	}
	
	if( !(bool)item->GetUserData() ){
		BrowseNewDir(item, fContents, false);
		item->SetOpen(true);
	}
	
	
	fSelectedDirText->SetText(dirfullpath.Data(), true);
	Layout();
}


void SpectraLoaderDialog::OnDoubleClick(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item) return;
	
	TString dirfullpath = FullFileName(item);
	
	//fSelectedDirText->SetText(dirfullpath.Data());
	Layout();
}


void SpectraLoaderDialog::LoadSpectrum()
{
	if(fSpectName->GetDisplayText()==TString("")) return;
	if( access(fSelectedDirText->GetDisplayText().Data(), R_OK|X_OK)!=0 ) return;
	
	if(fGatorCalib->IsDebug())
	{
		cout << "Debug---> SpectraLoaderDialog::LoadSpectrum: Loading spectrum \"" << fSpectName->GetDisplayText() << "\" from directory <" << fSelectedDirText->GetDisplayText() << ">" << endl;
	}
	
	string dirname = (fSelectedDirText->GetDisplayText()).Data();
	
	fGatorCalib->LoadCalibFiles( string((fSpectName->GetDisplayText()).Data()), dirname+string("/") );
	fGatorCalib->FillSpectraList();
	
	Emit("LoadSpectrum()");
}




//--------------------------------------------//
//  Definitions of the SaveLinesDialog class  //
//--------------------------------------------//

SaveLinesDialog::SaveLinesDialog(GatorCalibGUI *_GatorCalib):TGTransientFrame(gClient->GetRoot(), _GatorCalib, 300, 200, kVerticalFrame), fGatorCalib(_GatorCalib)
{
	SetCleanup(kDeepCleanup);
	Connect("CloseWindow()", "SaveLinesDialog", this, "CloseWindow()");
	DontCallClose();
	
	if(!_GatorCalib)
	{
		MakeZombie();
		return;
	}
	
	
	TGCompositeFrame *pTopFrame = new TGCompositeFrame(this, GetWidth()-2, GetHeight()-2, kHorizontalFrame);
	
	TGLabel *pLabel = new TGLabel(pTopFrame, "File name: ");
	fSelectedOutFile = new TGTextEntry(pTopFrame, "");
	
	//Make the list of directories to be selected
	TGCanvas* tgcanvas = new TGCanvas(this, 150, 100);
	fContents = new TGListTree(tgcanvas, kHorizontalFrame);
	TGListTreeItem* item = fContents->AddItem(0,"/");
	BrowseNewDir(item, fContents, true);
	item->SetOpen(true);
	//InitFilesTree(fContents); //Doesn't work properly
	
	
	//Bottom frame for OK and Cancel bottons
	TGCompositeFrame *pCommandsFrame = new TGCompositeFrame(this, 600, 100, kHorizontalFrame);
	
	
	TGButton *pOkButt = new TGTextButton(pCommandsFrame, "O&k");
	TGButton *pCancelButt = new TGTextButton(pCommandsFrame, "&Cancel");
	
	
	
	AddFrame(pTopFrame, new TGLayoutHints(kLHintsExpandX, 1, 5, 1, 5));
	AddFrame(tgcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY| kLHintsCenterX, 20, 1, 20, 10));
	AddFrame(pCommandsFrame, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 5, 1, 5, 5));
	
	pTopFrame->AddFrame(pLabel, new TGLayoutHints(kLHintsCenterY, 1, 1, 1, 5));
	pTopFrame->AddFrame(fSelectedOutFile, new TGLayoutHints(kLHintsCenterY|kLHintsExpandX, 1, 1, 1, 5));
	
	pCommandsFrame->AddFrame(pOkButt, new TGLayoutHints(kLHintsLeft, 5, 1, 5, 1));
	pCommandsFrame->AddFrame(pCancelButt, new TGLayoutHints(kLHintsLeft, 5, 1, 5, 1));
	
	
	//Connect section
	
	fContents->Connect("Clicked(TGListTreeItem*,Int_t)","SaveLinesDialog",this, "OnSingleClick(TGListTreeItem*,Int_t)");
	fContents->Connect("DoubleClicked(TGListTreeItem*,Int_t)","SaveLinesDialog",this, "OnDoubleClick(TGListTreeItem*,Int_t)");
	
	Connect("SaveSelectedLines(const char*)","GatorCalibGUI", fGatorCalib, "SaveSelectedLines(const char*)");
	pOkButt->Connect("Clicked()", "SaveLinesDialog", this, "SaveSelectedLines()");
	//pOkButt->Connect("Clicked()", "SaveLinesDialog", this, "CloseWindow()");
	pCancelButt->Connect("Clicked()", "SaveLinesDialog", this, "CloseWindow()");
	
	
	SetWindowName("Save selected lines.");
	MapSubwindows();
	Resize(GetDefaultSize());
	CenterOnParent(kTRUE, TGTransientFrame::kRight);
	
	MapWindow();
	
}


void SaveLinesDialog::SaveSelectedLines()
{
	TString outfilename = fSelectedOutFile->GetText();
	if( (outfilename=="") || (!outfilename.EndsWith(".root")) ) return;
	
	//fGatorCalib->SetOutputFile( string(outfilename.Data()) );
	
	SaveSelectedLines(outfilename.Data()); //This just emits a signal
	Emit("CloseWindow()");
}


void SaveLinesDialog::SaveSelectedLines(const char* outfilename)
{
	Emit("SaveSelectedLines(const char*)", outfilename);
}


void SaveLinesDialog::OnSingleClick(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item) return;
	
	// use UserData to indicate that item was already browsed
	
	TString filefullname = FullFileName(item);
	TString pathname = PathName(filefullname);
	
	
	TSystemFile selectedfile(filefullname, pathname);
	
	
	if(selectedfile.IsDirectory())
	{
		if( !((bool)item->GetUserData()) )
		{
			BrowseNewDir(item, fContents, true);
			item->SetOpen(true);
		}
	}
	
	fSelectedOutFile->SetText(filefullname.Data(), true);
	Layout();
}
	

void SaveLinesDialog::OnDoubleClick(TGListTreeItem* item, Int_t btn)
{
	if ((btn!=kButton1) || !item ) return;
	
	TString filefullname = FullFileName(item);
	
	fSelectedOutFile->SetText(filefullname.Data(), true);
	Layout();
	
	//if(filefullname.EndsWith(".root")) SaveSelectedLines();
}


//Global functions
TString FullFileName(TGListTreeItem* item)
{
	// Returns an absolute path.
	TGListTreeItem* parent;
	TString filename = item->GetText();
	while ((parent=item->GetParent()))
	{
		filename = gSystem->ConcatFileName(parent->GetText(), filename);
		item = parent;
	}
	return filename;
}


TString PathName(TString filefullname)
{
	string str = filefullname.Data();
	
	if(filefullname.EndsWith("/"))
	{//This should never be the case
		str = str.substr(0,str.find_last_of("/"));
	}
	
	TString pathname = (str.substr(0,str.find_last_of("/"))).c_str();
	return pathname;
}


void InitFilesTree(TGListTree* pContents)
{
	if(!pContents) return;
	
	
	
	string dirfullpath = gSystem->GetWorkingDirectory();
	
	TGListTreeItem *item=pContents->AddItem(0,"/"), *newitem=NULL;
	
	item->SetOpen(true);
	
	/*
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): Before browsing the directories." << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item address = " << item << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item->GetText() = " << item->GetText() << endl;
	*/
	
	vector<TString> parents;
	if( dirfullpath.find_last_of("/")!=string::npos )
	{
		do{
			size_t pos = dirfullpath.find_last_of("/");
			//cout << "Debug---> GatorCalibGUI::InitFilesTree(...): dirfullpath = " << dirfullpath << endl;
			TString dirname = dirfullpath.substr(pos+1);
			//cout << "Debug---> GatorCalibGUI::InitFilesTree(...): dirname = " << dirname << endl;
			if(dirname!="")
			{
				parents.push_back(dirname);
			}
			dirfullpath = dirfullpath.substr(0, pos);
		}while( dirfullpath.find_last_of("/")!=string::npos );
	}
	
	if(!(bool)item->GetUserData()) BrowseNewDir(item, pContents, true);
	
	for(int iEn=parents.size()-1; iEn>=0; iEn--)
	{
		newitem = pContents->AddItem(item, parents.at(iEn));
		if(!(bool)newitem->GetUserData()) BrowseNewDir(newitem, pContents, true);
		newitem->SetOpen(true);
		pContents->SetSelected(newitem);
		pContents->OpenItem(pContents->GetSelected());
		item = pContents->GetSelected();
	}
	
	/*
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): After browsing the directories." << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item address = " << item << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item->GetText() = " << item->GetText() << endl;
	*/
	
	/*
	//Populate the tree with the elements of the current directory
	TSystemDirectory dir(item->GetText(), (gSystem->GetWorkingDirectory()).c_str());
	TList *files = dir.GetListOfFiles();
	
	pContents->SetSelected(item);
	
	
	if(files)
	{
		TIter next(files);
		TSystemFile *file;
		TString fname;
		while((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if(file->IsDirectory())
			{
				if((fname!=".") && (fname!=".."))
				{ // skip it
					pContents->AddItem(fname);
				}
			}
			else if(fname.EndsWith(".root"))
			{
				pContents->AddItem(item,fname,fIcon,fIcon);
			}
		}
		delete files;
	}
	*/
}


void InitDirsTree(TGListTree* pContents)
{
	if(!pContents) return;
	
	string dirfullpath = gSystem->GetWorkingDirectory();
	
	TGListTreeItem *item=pContents->AddItem(0,"/"), *newitem=NULL;
	
	item->SetOpen(true);
	
	/*
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): Before browsing the directories." << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item address = " << item << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item->GetText() = " << item->GetText() << endl;
	*/
	
	vector<TString> parents;
	if( dirfullpath.find_last_of("/")!=string::npos )
	{
		do{
			size_t pos = dirfullpath.find_last_of("/");
			//cout << "Debug---> GatorCalibGUI::InitFilesTree(...): dirfullpath = " << dirfullpath << endl;
			TString dirname = dirfullpath.substr(pos+1);
			//cout << "Debug---> GatorCalibGUI::InitFilesTree(...): dirname = " << dirname << endl;
			if(dirname!="")
			{
				parents.push_back(dirname);
			}
			dirfullpath = dirfullpath.substr(0, pos);
		}while( dirfullpath.find_last_of("/")!=string::npos );
	}
	
	if(!(bool)item->GetUserData())
	{
		BrowseNewDir(item, pContents, false);
	}
	
	for(int iEn=parents.size()-1; iEn>=0; iEn--)
	{
		newitem = pContents->AddItem(item, parents.at(iEn));
		if(!(bool)newitem->GetUserData()) BrowseNewDir(newitem, pContents, false);
		newitem->SetOpen(true);
		pContents->SetSelected(newitem);
		pContents->OpenItem(pContents->GetSelected());
		item = pContents->GetSelected();
	}
	
	/*
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): After browsing the directories." << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item address = " << item << endl;
	cout << "Debug---> GatorCalibGUI::InitFilesTree(...): item->GetText() = " << item->GetText() << endl;
	*/
	
	//Populate the tree with the elements of the current directory
	/*
	TSystemDirectory dir(item->GetText(), (gSystem->GetWorkingDirectory()).c_str());
	TList *files = dir.GetListOfFiles();
	
	pContents->SetSelected(item);
	
	
	if(files)
	{
		TIter next(files);
		TSystemFile *file;
		TString fname;
		while((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if(file->IsDirectory())
			{
				if((fname!=".") && (fname!=".."))
				{ // skip it
					pContents->AddItem(fname);
				}
			}
			
		}
		delete files;
	}
	*/
}


void BrowseNewDir(TGListTreeItem *item, TGListTree *pContents, bool rootfiles)
{
	if( (!item) || (!pContents) || ((bool)item->GetUserData()) ) return;
	
	
	TString filefullname = FullFileName(item);
	TString pathname = PathName(filefullname);
	
	//TSystemDirectory dir(filefullname.Data(), pathname.Data());
	TSystemDirectory dir(filefullname.Data(), filefullname.Data());
	
	if(!dir.IsDirectory()) return;
	
	TList *files = dir.GetListOfFiles();
	
	pContents->SetSelected(item);
	
	if(files)
	{
		item->SetUserData((void*)1);
		
		TIter next(files);
		TSystemFile *file;
		TString fname;
		while((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if(file->IsDirectory())
			{
				if((fname!=".") && (fname!=".."))
				{ // skip it
					if(!fname.BeginsWith("."))
					{
						pContents->AddItem(item, fname);
					}
				}
			}
			else if( rootfiles && fname.EndsWith(".root") )
			{
				//pContents->AddItem(item, fname, gIcon, gIcon);
				pContents->AddItem(item, fname);
			}
		}
		delete files;
	}
}


