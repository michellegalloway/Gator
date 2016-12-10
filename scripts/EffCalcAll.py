import ROOT



#Uncomment the line below if you want to run in batch
#ROOT.gROOT.SetBatch()

ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetFuncColor(2)
  
#Here are the setup and initialization stuff to be changed manually
datadir = "/home/atp/fpiastra/g4sims/Screening/xebottle/"
sample="SteelBottle"
events = "1e8"

fileprefix = "merged_10PMTs-ceramic"

DataPath = datadir + sample
	
workdir = datadir + sample + "/effic"
	
outfilename = "lines.list"
	
NumInEvt =100000000
	
#Switches
U238 = True #If the 238U chain is simulated
Ra226 = True #If the 226Ra chain is simulated
Ra228 = True #If the "232Th" chain is simulated (only 228Ac lines are detectable)
Th228 = True #If the 228Th chain is simulated (Th232 late chain)
U235 = False #If want or not the 235U info
U186 = False #To calculate the 186 keV BRxEff from direct gamma simulation and not from the 235U decay
K40 = True
Co60 = False
Cs137 = False
Sc46 = False #If 46Sc is simulated
Pa1001 = False #To calculate the 1001 keV BRxEff from direct gamma simulation and not from the 238U decay chain
Pa = True #If want or not the 234mPa info
Mn54 = False #If I want Mn54 info
	
writetxt = True
	
compbins = 10 #number of bins left and right of lines in which estimate the compton BG

#Int_t leftcompbins, rightcompbins #Used in non simmetric Compton estimation


###    Define the branching ratio (BR) for the main lines  
BR186 = 0.570 # 186 keV of 235U
BR239 = 0.436  # 239 keV of 212Pb
BR339 = 0.114  # 339 keV of 228Ac
BR583 = 0.3054 # 583 keV of 208Tl (from nucleide.org)
BR911 = 0.262  # 911 keV of 228Ac
BR969 = 0.159  # 969 keV of 228Ac
BR2615 = 0.3584 # 2615 keV of 208Tl
BR295 = 0.184  # 295 keV of 214Pb
BR352 = 0.356  # 352 keV of 214Pb
BR609 = 0.4549  # 609 keV of 214Bi
BR1001 = 0.00846 # 1001 keV of 234mPa
BR1120 = 0.1491 # 1120 keV of 214Bi
BR1765 = 0.1531 # 1765 keV of 214Bi
BR1173 = 0.9985 # 1173 keV of 60Co
BR1332 = 0.9998 # 1332 keV of 60Co
BR1460 = 0.1055 # 1460 keV of 40K
BR835 = 0.999746 # 835 keV of 54Mn
BR889 = 0.99984  # 889 keV of 46Sc
BR1120s = 0.99987  # 1120 keV of 46Sc
BR662 = 0.8499  # 662 keV of 137Cs


point = 0 #Point number of the efficiency TGraphErrors

g_effic = ROOT.TGraphErrors()

##########################################################
############ Starting with the 228Ra chain ###############
##########################################################
if Ra228:
	
	c1 = ROOT.TCanvas("c1","Simulated Spectra",870,500)
	c1.SetLogy()
	c1.cd()
  
	leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
	
	FullFN1 = DataPath+"/"+fileprefix+"_228Ra_"+events+".root"
	
	description1 = fileprefix+" ^{228}Ra Simulation"
	
	print "\nOpening <"+FullFN1+">"
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h1 = ROOT.TH1F("h1",description1,27000,0,2700)
  
	h1.SetXTitle("Energy [keV]")
	h1.SetYTitle("counts")
	h1.SetStats(0)
	t1.Draw("GeEtot>>h1")
	
	


################ Calculate the efficiencies ################

	#----------  338.32 keV from 228Ac (11 % BR)  -----------#
	binx = h1.FindBin(338.320)
	
	BRxEff339 = (h1.GetBinContent(binx) - (h1.Integral(binx-(compbins+1),binx-2)/compbins + h1.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff339 = (h1.GetBinContent(339) - (h1.Integral(333,337)/5 + h1.Integral(341,355)/5)/2)/NumInEvt
 
 	RealEff339 = BRxEff339/BR339

 	print "Efficiency x BR (338): "+str(BRxEff339)
 	print "Real Efficiency (338): "+str(RealEff339)
	
	g_effic.SetPoint(point,338.320,RealEff339)
	g_effic.SetPointError(point,0.0,0.1*RealEff339)
	point=point+1
	
	#----------  911.196 keV from 228Ac (27 % BR)  -----------#
	binx = h1.FindBin(911.2) #Geant4 for doesn't know exactly the energy of this line..... this should adjust the things
	
	BRxEff911 = (h1.GetBinContent(binx) - (h1.Integral(binx-(compbins+1),binx-2)/compbins + h1.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff911 = (h1.GetBinContent(912) - (h1.Integral(906,910)/5 + h1.Integral(914,918)/5)/2)/NumInEvt

	RealEff911 = BRxEff911/BR911

	print "Efficiency x BR (911): "+str(BRxEff911)
	print "Real Efficiency (911): "+str(RealEff911)
	
	g_effic.SetPoint(point,911.196,RealEff911)
	g_effic.SetPointError(point,0.0,0.1*RealEff911)
	point=point+1

	#----------  968.7 keV from 228Ac (16 % BR)  -----------#
	binx = h1.FindBin(968.96)
	BRxEff969 = (h1.GetBinContent(binx) - (h1.Integral(binx-(compbins+1),binx-2)/compbins + h1.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff969 = (h1.GetBinContent(969) - (h1.Integral(963,967)/5 + h1.Integral(976,980)/5)/2)/NumInEvt

	RealEff969 = BRxEff969/BR969
 
	print "Efficiency x BR (969): "+str(BRxEff969)
	print "Real Efficiency (969): "+str(RealEff969)
	
	g_effic.SetPoint(point,968.96,RealEff969)
	g_effic.SetPointError(point,0.0,0.1*RealEff969)
	point=point+1

	#*****  When I simulate the entire chain and not separately	*****#
	if not Th228:
		
		#----------  238.6 keV from 212Pb (44 % BR)  -----------#
		binx = h1.FindBin(238.632)
		BRxEff239 = (h1.GetBinContent(binx)-(h1.Integral(binx-(compbins+1),binx-2)/compbins + h1.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

		RealEff239 = BRxEff239/BR239
  
		print "Efficiency x BR (239): "+str(BRxEff239)
		print "Real Efficiency (239): "+str(RealEff239)
 
		g_effic.SetPoint(point,238.632,RealEff239)
		g_effic.SetPointError(point,0.0,0.1*RealEff239)
		point=point+1
		
		#----------  583.2 keV from 208Tl (30.4 % BR)  -----------#
#NOTE: At 583.4 there the peaks of 228Ac 228Pa (visible) ==> Bias on the left the Compton counting
		rightcompbins = 1
		leftcompbins = 5
		binx = h1.FindBin(583.187)
	
		BRxEff583 = (h1.Integral(binx,binx+1) - 2*(h1.Integral(binx-(leftcompbins+1),binx-2)/leftcompbins + h1.Integral(binx+2,binx+(rightcompbins+1))/rightcompbins)/2)/NumInEvt
	
		#BRxEff583 = (h1.GetBinContent(584) - (h1.Integral(578,582)/5 + h1.Integral(586,590)/5)/2)/NumInEvt
 
		RealEff583 = BRxEff583/BR583

		print "Efficiency x BR (583): "+str(BRxEff583)
		print "Real Efficiency (583): "+str(RealEff583)
 	
		g_effic.SetPoint(point,583.187,RealEff583)
		g_effic.SetPointError(point,0.0,0.1*RealEff583)
		point=point+1
		
		#----------  2614.5 keV from 208Tl (35.6 % BR)  -----------#
		binx = h1.FindBin(2614.511)
		BRxEff2615 = (h1.GetBinContent(binx) - (h1.Integral(binx-(compbins+1),binx-2)/compbins + h1.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

		#BRxEff2615 = (h1.GetBinContent(2615) - (h1.Integral(2610,2614)/5 + h1.Integral(2617,2621)/5))/NumInEvt

		RealEff2615 = BRxEff2615/BR2615
		
		print "Efficiency x BR (2615): "+str(BRxEff2615)
		print "Real Efficiency (2615): "+str(RealEff2615)
	
		g_effic.SetPoint(point,2614.511,RealEff2615)
		g_effic.SetPointError(point,0.0,0.1*RealEff2615)
		point=point+1

	
	c1.SaveAs(workdir+"/"+fileprefix+"_228Ra.C")

	t1.Delete()
###################################################
############ Finished with 228Ra chain ############
###################################################



##########################################################
############ Starting with the 228Th chain ###############
##########################################################
if Th228:
	
	c2 = ROOT.TCanvas("c2","Simulated Spectra",870,500)
	c2.SetLogy()
	c2.cd()
  
	leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
	
	FullFN1 = DataPath+"/"+fileprefix+"_228Th_"+events+".root"
	
	description1 = fileprefix+" ^{228}Th Simulation"
	
	print "\nOpening <"+FullFN1+">"
    
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h2 = ROOT.TH1F("h2",description1,27000,0,2700)
  
	h2.SetXTitle("Energy [keV]")
	h2.SetYTitle("counts")
	h2.SetStats(0)
	t1.Draw("GeEtot>>h2")
	
	
################ Calculate the efficiencies ################
	
	#----------  238.6 keV from 212Pb (44 % BR)  -----------#
	binx = h2.FindBin(238.632)
	BRxEff239 = (h2.GetBinContent(binx)-(h2.Integral(binx-(compbins+1),binx-2)/compbins + h2.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

	RealEff239 = BRxEff239/BR239
  
	print "Efficiency x BR (239): "+str(BRxEff239)
	print "Real Efficiency (239): "+str(RealEff239)
 
	g_effic.SetPoint(point,238.632,RealEff239)
	g_effic.SetPointError(point,0.0,0.1*RealEff239)
	point=point+1

	#----------  583.2 keV from 208Tl (30.4 % BR)  -----------#
#NOTE: At 583.4 there the peaks of 228Ac 228Pa (visible) ==> Bias on the left the Compton counting
	rightcompbins = 1
	leftcompbins = 5
	binx = h2.FindBin(583.187)
	
	BRxEff583 = (h2.Integral(binx,binx+1) - 2*(h2.Integral(binx-(leftcompbins+1),binx-2)/leftcompbins + h2.Integral(binx+2,binx+(rightcompbins+1))/rightcompbins)/2)/NumInEvt
	
	#BRxEff583 = (h2.GetBinContent(584) - (h2.Integral(578,582)/5 + h2.Integral(586,590)/5)/2)/NumInEvt
 
 	RealEff583 = BRxEff583/BR583

 	print "Efficiency x BR (583): "+str(BRxEff583)
 	print "Real Efficiency (583): "+str(RealEff583)
 	
 	g_effic.SetPoint(point,583.187,RealEff583)
	g_effic.SetPointError(point,0.0,0.1*RealEff583)
	point=point+1

	#----------  2614.5 keV from 208Tl (35.6 % BR)  -----------#
	binx = h2.FindBin(2614.511)
	BRxEff2615 = (h2.GetBinContent(binx) - (h2.Integral(binx-(compbins+1),binx-2)/compbins + h2.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

	#BRxEff2615 = (h2.GetBinContent(2615) - (h2.Integral(2610,2614)/5 + h2.Integral(2617,2621)/5))/NumInEvt

	RealEff2615 = BRxEff2615/BR2615
 
	print "Efficiency x BR (2615): "+str(BRxEff2615)
	print "Real Efficiency (2615): "+str(RealEff2615)
	
	g_effic.SetPoint(point,2614.511,RealEff2615)
	g_effic.SetPointError(point,0.0,0.1*RealEff2615)
	point=point+1

	c2.SaveAs(workdir+"/"+fileprefix+"_228Th.C")

	t1.Delete()
###################################################
############ Finished with 228Th chain ############
###################################################
	
	
	
#########################################################
############ Starting with the 238U chain ###############
#########################################################
if U238:

	c3 = ROOT.TCanvas("c3","Simulated Spectra",870,500)
	c3.SetLogy()
	c3.cd()

	FullFN1 = DataPath+"/"+fileprefix+"_238U_"+events+".root"
	
	description1 = fileprefix+" ^{238}U Simulation"
	
	print "\nOpening <"+FullFN1+">"
    
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h3 = ROOT.TH1F("h3",description1,27000,0,2700)
  
	h3.SetXTitle("Energy [keV]")
	h3.SetYTitle("counts")
	h3.SetStats(0)
	t1.Draw("GeEtot>>h3")
	
	

	################ Calculate the efficiencies ################

	#----------  92.38keV and 92.80keV from 234Th (2.18% BR and 2.15% BR)  -----------#
	# Since in real data I don't resolve the two peaks, here I have to define a mean BRxEff for 92.6keV (weighted mean)
	# All the variables used are defined in this sub-section.

	bin92a = h3.FindBin(92.38)
	bin92b = h3.FindBin(92.75) # the correct position should be 92.75 but sometimes from the simulations it is shifted a bit
	
	BR92a = 0.0218
	BR92b = 0.0215
	
	BR92 = BR92a + BR92b
	
	#Here the stuff is a bit involved
	leftcompbins = 7
	rightcompbins = 1
	BRxEff92a = (h3.GetBinContent(bin92a) - (h3.Integral(bin92a-(leftcompbins+1),bin92a-2)/leftcompbins + h3.Integral(bin92a+2,bin92a+(rightcompbins+1))/rightcompbins)/2)/NumInEvt
	
	leftcompbins = 1
	rightcompbins = 7
	BRxEff92b = (h3.Integral(bin92b,bin92b+1) - 2*(h3.Integral(bin92b-(leftcompbins+1),bin92b-2)/leftcompbins + h3.Integral(bin92b+3,bin92b+(rightcompbins+2))/rightcompbins)/2)/NumInEvt
	
	BRxEff92 = BRxEff92a + BRxEff92b
	
	RealEff92 = BRxEff92 / BR92
	
	print "Efficiency x BR (92): "+str(BRxEff92)
	print "Real Efficiency (92): "+str(RealEff92)
	
	g_effic.SetPoint(point,92.6,RealEff92)
	g_effic.SetPointError(point,0.0,0.1*RealEff92)
	point=point+1
	
	if (not Pa1001) and Pa:
		#----------  1001 keV from 234mPa (0.59 % BR)  -----------#
		# This fainty line is simulated alone if here cannot be determined (visual check of the spectrum)
		binx = h3.FindBin(1001.026)
		BRxEff1001 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
		
		#BRxEff1001 = (h3.GetBinContent(1002) - (h3.Integral(996,1000)/5 + h3.Integral(1004,1008)/5)/2)/NumInEvt
		
		RealEff1001 = BRxEff1001 / BR1001
		
		print "Efficiency x BR (1001): "+str(BRxEff1001)
		print "Real Efficiency (1001): "+str(RealEff1001)
		
		g_effic.SetPoint(point,1001.026,RealEff1001)
		g_effic.SetPointError(point,0.0,0.1*RealEff1001)
		point=point+1
	
	#*****  When I simulate the entire 238U chain and not separately	*****#
	if not Ra226:
		
		#----------  295.2 keV from 214Pb (19 % BR)  -----------#
		binx = h3.FindBin(295.224)
		BRxEff295 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
		
		#BRxEff295 = (h3.GetBinContent(296) - (h3.Integral(290,294)/5 + h3.Integral(298,302)/5)/2)/NumInEvt
		
		RealEff295 = BRxEff295 / BR295
		
		print "Efficiency x BR (295): "+str(BRxEff295)
		print "Real Efficiency (295): "+str(RealEff295)
	
		g_effic.SetPoint(point,295.224,RealEff295)
		g_effic.SetPointError(point,0.0,0.1*RealEff295)
		point=point+1
	
		#----------  351.9 keV from 214Pb (37 % BR)  -----------#
		binx = h3.FindBin(351.932)
		BRxEff352 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
		#BRxEff352 = (h3.GetBinContent(352) - (h3.Integral(346,350)/5 + h3.Integral(354,358)/5)/2)/NumInEvt
		
		RealEff352 = BRxEff352 / BR352
		
		print "Efficiency x BR (352): "+str(BRxEff352)
		print "Real Efficiency (352): "+str(RealEff352)
	
		g_effic.SetPoint(point,351.932,RealEff352)
		g_effic.SetPointError(point,0.0,0.1*RealEff352)
		point=point+1
		
		#----------  609.3 keV from 214Bi (46 % BR)  -----------#
		binx = h3.FindBin(609.312)
		BRxEff609 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
		
		#BRxEff609 = (h3.GetBinContent(610) - (h3.Integral(604,608)/5 + h3.Integral(612,616)/5)/2)/NumInEvt
		
		RealEff609 = BRxEff609 / BR609
		
		print "Efficiency x BR (609): "+str(BRxEff609)
		print "Real Efficiency (609): "+str(RealEff609)
		
		g_effic.SetPoint(point,609.312,RealEff609)
		g_effic.SetPointError(point,0.0,0.1*RealEff609)
		point=point+1
		
		#----------  1120.3 keV from 214Bi (15 % BR)  -----------#
		binx = h3.FindBin(1120.287)
		BRxEff1120 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
		#BRxEff1120 = (h3.GetBinContent(1121) - (h3.Integral(1115,1119)/5 + h3.Integral(1123,1127)/5)/2)/NumInEvt
		
		RealEff1120 = BRxEff1120 / BR1120

		print "Efficiency x BR (1120): "+str(BRxEff1120)
		print "Real Efficiency (1120): "+str(RealEff1120)
		
		g_effic.SetPoint(point,1120.287,RealEff1120)
		g_effic.SetPointError(point,0.0,0.1*RealEff1120)
		point=point+1
	
		#----------  1764.5 keV from 214Bi (16 % BR)  -----------#
		binx = h3.FindBin(1764.494)
		BRxEff1765 = (h3.GetBinContent(binx) - (h3.Integral(binx-(compbins+1),binx-2)/compbins + h3.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
		#BRxEff1765 = (h3.GetBinContent(1765) - (h3.Integral(1759,1763)/5 + h3.Integral(1767,1771)/5)/2)/NumInEvt
		
		RealEff1765 = BRxEff1765 / BR1765

		print "Efficiency x BR (1765): "+str(BRxEff1765)
		print "Real Efficiency (1765): "+str(RealEff1765)
	
		g_effic.SetPoint(point,1764.494,RealEff1765)
		g_effic.SetPointError(point,0.0,0.1*RealEff1765)
		point=point+1
	
	c3.SaveAs(workdir+"/"+fileprefix+"_238U.C")
	
	t1.Delete()
##################################################
############ Finished with 238U chain ############
##################################################



##########################################################
############ Starting with the 226Ra chain ###############
##########################################################
if Ra228:
	
	c4 = ROOT.TCanvas("c4","Simulated Spectra",870,500)
	c4.SetLogy()
	c4.cd()

	FullFN1 = DataPath+"/"+fileprefix+"_226Ra_"+events+".root"
	
	description1 = fileprefix+" ^{226}Ra Simulation"
	
	print "\nOpening <"+FullFN1+">"
    
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h4 = ROOT.TH1F("h4",description1,27000,0,2700)
  
	h4.SetXTitle("Energy [keV]")
	h4.SetYTitle("counts")
	h4.SetStats(0)
	t1.Draw("GeEtot>>h4")
	
	
	################ Calculate the efficiencies ################
	
	#----------  295.2 keV from 214Pb (19 % BR)  -----------#
	binx = h4.FindBin(295.224)
	BRxEff295 = (h4.GetBinContent(binx) - (h4.Integral(binx-(compbins+1),binx-2)/compbins + h4.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff295 = (h4.GetBinContent(296) - (h4.Integral(290,294)/5 + h4.Integral(298,302)/5)/2)/NumInEvt
	
	RealEff295 = BRxEff295 / BR295
	
	print "Efficiency x BR (295): "+str(BRxEff295)
	print "Real Efficiency (295): "+str(RealEff295)

	g_effic.SetPoint(point,295.224,RealEff295)
	g_effic.SetPointError(point,0.0,0.1*RealEff295)
	point=point+1

	#----------  351.9 keV from 214Pb (37 % BR)  -----------#
	binx = h4.FindBin(351.932)
	BRxEff352 = (h4.GetBinContent(binx) - (h4.Integral(binx-(compbins+1),binx-2)/compbins + h4.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

	#BRxEff352 = (h4.GetBinContent(352) - (h4.Integral(346,350)/5 + h4.Integral(354,358)/5)/2)/NumInEvt
	
	RealEff352 = BRxEff352 / BR352
		
	print "Efficiency x BR (352): "+str(BRxEff352)
	print "Real Efficiency (352): "+str(RealEff352)

	g_effic.SetPoint(point,351.932,RealEff352)
	g_effic.SetPointError(point,0.0,0.1*RealEff352)
	point=point+1
	
	#----------  609.3 keV from 214Bi (46 % BR)  -----------#
	binx = h4.FindBin(609.312)
	BRxEff609 = (h4.GetBinContent(binx) - (h4.Integral(binx-(compbins+1),binx-2)/compbins + h4.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff609 = (h4.GetBinContent(610) - (h4.Integral(604,608)/5 + h4.Integral(612,616)/5)/2)/NumInEvt
	
	RealEff609 = BRxEff609 / BR609
	
	print "Efficiency x BR (609): "+str(BRxEff609)
	print "Real Efficiency (609): "+str(RealEff609)
	
	g_effic.SetPoint(point,609.312,RealEff609)
	g_effic.SetPointError(point,0.0,0.1*RealEff609)
	point=point+1
	
	#----------  1120.3 keV from 214Bi (15 % BR)  -----------#
	binx = h4.FindBin(1120.287)
	BRxEff1120 = (h4.GetBinContent(binx) - (h4.Integral(binx-(compbins+1),binx-2)/compbins + h4.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

	#BRxEff1120 = (h4.GetBinContent(1121) - (h4.Integral(1115,1119)/5 + h4.Integral(1123,1127)/5)/2)/NumInEvt
	
	RealEff1120 = BRxEff1120 / BR1120
	
	print "Efficiency x BR (1120): "+str(BRxEff1120)
	print "Real Efficiency (1120): "+str(RealEff1120)
	
	g_effic.SetPoint(point,1120.287,RealEff1120)
	g_effic.SetPointError(point,0.0,0.1*RealEff1120)
	point=point+1

	#----------  1764.5 keV from 214Bi (16 % BR)  -----------#
	binx = h4.FindBin(1764.494)
	BRxEff1765 = (h4.GetBinContent(binx) - (h4.Integral(binx-(compbins+1),binx-2)/compbins + h4.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt

	#BRxEff1765 = (h4.GetBinContent(1765) - (h4.Integral(1759,1763)/5 + h4.Integral(1767,1771)/5)/2)/NumInEvt
	
	RealEff1765 = BRxEff1765 / BR1765
	
	print "Efficiency x BR (1765): "+str(BRxEff1765)
	print "Real Efficiency (1765): "+str(RealEff1765)

	g_effic.SetPoint(point,1764.494,RealEff1765)
	g_effic.SetPointError(point,0.0,0.1*RealEff1765)
	point=point+1
	
	
	c4.SaveAs(workdir+"/"+fileprefix+"_226Ra.C")
	
	t1.Delete()
###################################################
############ Finished with 226Ra chain ############
###################################################



########################################################
############ Starting with the 40K chain ###############
########################################################
if K40:

	c5 = ROOT.TCanvas("c5","Simulated Spectra",870,500)
	c5.SetLogy()
	c5.cd()

	FullFN1 = DataPath+"/"+fileprefix+"_40K_"+events+".root"
	
	description1 = fileprefix+" ^{40}K Simulation"
	
	print "\nOpening <"+FullFN1+">"
    
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h5 = ROOT.TH1F("h5",description1,27000,0,2700)
  
	h5.SetXTitle("Energy [keV]")
	h5.SetYTitle("counts")
	h5.SetStats(0)
	t1.Draw("GeEtot>>h5")
	
	
	################ Calculate the efficiencies ################
	#NOTE: This is the highest energy line.... no Compton affects it.
	
	#----------  1460.882 keV from 40K (11 % BR)  -----------#
	binx = h5.FindBin(1460.882)
	BRxEff1460 = h5.GetBinContent(binx)/NumInEvt
	
	#BRxEff1460 = h5.GetBinContent(1461)/NumInEvt

	RealEff1460 = BRxEff1460 / BR1460

	print "Efficiency x BR (1460): "+str(BRxEff1460)
	print "Real Efficiency (1460): "+str(RealEff1460)
	
	g_effic.SetPoint(point,1460.882,RealEff1460)
	g_effic.SetPointError(point,0.0,0.1*RealEff1460)
	point=point+1
	
	
	c5.SaveAs(workdir+"/"+fileprefix+"_40K.C")
	
	t1.Delete()
#################################################
############ Finished with 40K chain ############
#################################################	



#########################################################
############ Starting with the 60Co chain ###############
#########################################################
if Co60:

	c6 = ROOT.TCanvas("c6","Simulated Spectra",870,500)
	c6.SetLogy()
	c6.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_60Co_"+events+".root"
	
	description1 = fileprefix+" ^{60}Co Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h6 = ROOT.TH1F("h6",description1,27000,0,2700)
  
	h6.SetXTitle("Energy [keV]")
	h6.SetYTitle("counts")
	h6.SetStats(0)
	t1.Draw("GeEtot>>h6")

	
	################ Calculate the efficiencies ################
	
	#----------  1173.2 keV from 60Co (99 % BR)  -----------#
	binx = h6.FindBin(1173.228)
	BRxEff1173 = (h6.GetBinContent(binx) - (h6.Integral(binx-(compbins+1),binx-2)/compbins + h6.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
	
	#BRxEff1173 = (h6.GetBinContent(1174) - (h6.Integral(1168,1172)/5 + h6.Integral(1176,1180)/5)/2)/NumInEvt

	RealEff1173 = BRxEff1173 / BR1173
 
	print "Efficiency x BR (1173): "+str(BRxEff1173)
	print "Real Efficiency (1173): "+str(RealEff1173)
	
	g_effic.SetPoint(point,1173.228,RealEff1173)
	g_effic.SetPointError(point,0.0,0.1*RealEff1173)
	point=point+1
	

	#----------  1332 keV from 60Co (99 % BR)  -----------#
#NOTE: Is clearly visible his own Compton in the neibur bins.... shift the compton interval more to the right
	binx = h6.FindBin(1332.495) #It should be 1332.495 but G4.9.3 sometimes gives it at 1332.52
	
	rightcompbins = 5
	leftcompbins = 3
	
	BRxEff1332 = (h6.Integral(binx,binx+1) - 2*(h6.Integral(binx-(leftcompbins+1),binx-2)/leftcompbins + h6.Integral(binx+3,binx+(rightcompbins+2))/rightcompbins)/2)/NumInEvt
	
	#BRxEff1332 = (h6.GetBinContent(1333) - (h6.Integral(1327,1331)/5 + h6.Integral(1335,1339)/5)/2)/NumInEvt

	RealEff1332 = BRxEff1332 / BR1332

	print "Efficiency x BR (1332): "+str(BRxEff1332)
	print "Real Efficiency (1332): "+str(RealEff1332)
	
	g_effic.SetPoint(point,1332.495,RealEff1332)
	g_effic.SetPointError(point,0.0,0.1*RealEff1332)
	point=point+1
	
	
	c6.SaveAs(workdir+"/"+fileprefix+"_40K.C")
	
	t1.Delete()
##################################################
############ Finished with 60Co chain ############
##################################################	


	
#########################################################
############ Starting with the 46Sc chain ###############
#########################################################
if Sc46:
	
	c7 = ROOT.TCanvas("c7","Simulated Spectra",870,500)
	c7.SetLogy()
	c7.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_46Sc_"+events+".root"
	
	description1 = fileprefix+" ^{46}Sc Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h7 = ROOT.TH1F("h7",description1,27000,0,2700)
  
	h7.SetXTitle("Energy [keV]")
	h7.SetYTitle("counts")
	h7.SetStats(0)
	t1.Draw("GeEtot>>h7")
	
	
	################ Calculate the efficiencies ################

	#----------  889.3 keV from 46Sc (85 % BR)  -----------#
	binx = h7.FindBin(889.271)
		
	rightcompbins = 5
	leftcompbins = 3
		
	BRxEff889 = (h7.GetBinContent(binx) - (h7.Integral(binx-(leftcompbins+1),binx-2)/leftcompbins + h7.Integral(binx+2,binx+(rightcompbins+1))/rightcompbins)/2)/NumInEvt
		
		
	#BRxEff889 = (h7.GetBinContent(890) - (h7.Integral(884,888)/5 + h7.Integral(892,896)/5)/2)/NumInEvt

	RealEff889 = BRxEff889 / BR889
 
	print "Efficiency x BR (889): "+str(BRxEff889)
	print "Real Efficiency (889): "+str(RealEff889)
	
	g_effic.SetPoint(point,889.271,RealEff889)
	g_effic.SetPointError(point,0.0,0.1*RealEff889)
	point=point+1


	#----------  1120.5 keV from 46Sc (85 % BR)  -----------#
#NOTE: Is clearly visible his own Compton in the neighbor bins.... shift the compton interval more to the right
	binx = h7.FindBin(1120.537)
	BRxEff1120s = (h7.GetBinContent(binx) - (h7.Integral(binx-(compbins+1),binx-2)/compbins + h7.Integral(binx+2,binx+(compbins+1))/compbins)/2)/NumInEvt
 		
	#BRxEff1120s = h7.GetBinContent(1121)/NumInEvt

	RealEff1120s = BRxEff1120s / BR1120s
 
	print "Efficiency x BR (1120): "+str(BRxEff1120s)
	print "Real Efficiency (1120): "+str(RealEff1120s)
	
	g_effic.SetPoint(point,1120.537,RealEff1120s)
	g_effic.SetPointError(point,0.0,0.1*RealEff1120s)
	point=point+1
	
	
	c7.SaveAs(workdir+"/"+fileprefix+"_46Sc.C")
  
	t1.Delete()
##################################################
############ Finished with 46Sc chain ############
##################################################	
	
	
	
##########################################################
############ Starting with the 137Cs chain ###############
##########################################################
if Cs137:

	c8 = ROOT.TCanvas("c8","Simulated Spectra",870,500)
	c8.SetLogy()
	c8.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_137Cs_"+events+".root"
	
	description1 = fileprefix+" ^{137}Cs Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h8 = ROOT.TH1F("h8",description1,27000,0,2700)
  
	h8.SetXTitle("Energy [keV]")
	h8.SetYTitle("counts")
	h8.SetStats(0)
	t1.Draw("GeEtot>>h8")
	
	
	################ Calculate the efficiencies ################
	
	#----------  661.7 keV from 137Cs (84.99 % BR)  -----------#
	#NOTE: This is the highest energy line.... no Compton affects it.
	binx = h8.FindBin(661.657)
	BRxEff662 = (h8.GetBinContent(binx) - (h8.Integral(binx-(compbins+1),binx-2)/compbins))/NumInEvt
	
	
	#BRxEff662 = (h8.GetBinContent(662) - (h8.Integral(656,660)/5 + h8.Integral(664,668)/5)/2)/NumInEvt
	
	RealEff662 = BRxEff662 / BR662
	
	print "Efficiency x BR (662): "+str(BRxEff662)
	print "Real Efficiency (662): "+str(RealEff662)
	
	g_effic.SetPoint(point,661.657,RealEff662)
	g_effic.SetPointError(point,0.0,0.1*RealEff662)
	point=point+1
	
	
	c8.SaveAs(workdir+"/"+fileprefix+"_137Cs.C")
  
	t1.Delete()
###################################################
############ Finished with 137Cs chain ############
###################################################	
	
	
	
#########################################################
############ Starting with the 235U chain ###############
#########################################################
if (not U186) and U235:
	
	c9 = ROOT.TCanvas("c9","Simulated Spectra",870,500)
	c9.SetLogy()
	c9.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_137Cs_"+events+".root"
	
	description1 = fileprefix+" ^{137}Cs Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h9 = ROOT.TH1F("h9",description1,27000,0,2700)
  
	h9.SetXTitle("Energy [keV]")
	h9.SetYTitle("counts")
	h9.SetStats(0)
	t1.Draw("GeEtot>>h9")
	
	
		
	################ Calculate the efficiencies ################
	
	#---------- 185.712 keV from 235U (57.0 % BR)  -----------#
	binx = h9.FindBin(185.712) #Sometimes Geant4 knows it as 185.6 keV ==> must be checked each time 
	BRxEff186 = (h9.GetBinContent(binx) - (h9.Integral(binx-(compbins+1),binx-2)/compbins))/NumInEvt
	
	#BRxEff186 = (h9.GetBinContent(186) - (h9.Integral(180,184)/5 + h9.Integral(188,192)/5)/2)/NumInEvt
	
	RealEff186 = BRxEff186 / BR186
	
	print "Efficiency x BR (186): "+str(BRxEff186)
	print "Real Efficiency (186): "+str(RealEff186)
	
	g_effic.SetPoint(point,185.712,RealEff186)
	g_effic.SetPointError(point,0.0,0.1*RealEff186)
	point=point+1
	
	
	c9.SaveAs(workdir+"/"+fileprefix+"_235U.C")
	
	t1.Delete()
##################################################
############ Finished with 235U chain ############
##################################################	



###################################################################
############ Starting with the 234mPa 1001 keV line ###############
###################################################################
if Pa1001 and Pa:

	c10 = ROOT.TCanvas("c10","Simulated Spectra",870,500)
	c10.SetLogy()
	c10.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_1001keV_"+events+".root"
	
	description1 = fileprefix+" 1001 keV Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h10 = ROOT.TH1F("h10",description1,27000,0,2700)
  
	h10.SetXTitle("Energy [keV]")
	h10.SetYTitle("counts")
	h10.SetStats(0)
	t1.Draw("GeEtot>>h10")
	
	
	################ Calculate the efficiencies ################
	
	#---------- 1001 KeV from 234mPa (no BR)  -----------#
	
	RealEff1001 = h10.GetBinContent(h10.FindBin(1001.026))/NumInEvt
	
	BRxEff1001 = RealEff1001 * BR1001
	
	print "Efficiency x BR (1001): "+str(BRxEff1001)
	print "Real Efficiency (1001): "+str(RealEff1001)
	
	g_effic.SetPoint(point,1001.026,RealEff1001)
	g_effic.SetPointError(point,0.0,0.1*RealEff1001)
	point=point+1
	
	
	c10.SaveAs(workdir+"/"+fileprefix+"_1001keV.C")
	
	t1.Delete()
################################################################
############ Finished with the 234mPa 1001 keV line ############
################################################################
	
	

###################################################################
############ Starting with the 235U 186 keV line ###############
###################################################################
if U186 and U235:

	#This is simulated only if and when it is necessary, otherwise it can be used the 235U isotope decay simulation
	
	c11 = ROOT.TCanvas("c11","Simulated Spectra",870,500)
	c11.SetLogy()
	c11.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_186keV_"+events+".root"
	
	description1 = fileprefix+" 186 keV Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h11 = ROOT.TH1F("h11",description1,27000,0,2700)
  
	h11.SetXTitle("Energy [keV]")
	h11.SetYTitle("counts")
	h11.SetStats(0)
	t1.Draw("GeEtot>>h11")
	

	################ Calculate the efficiencies ################
	
	#----------  185.7 KeV from 235U (no BR)  -----------#
	
	RealEff186 = h11.GetBinContent(h11.FindBin(185.720))/NumInEvt
	
	BRxEff186 = RealEff186 * BR186
	
	print "Efficiency x BR (186): "+str(BRxEff186)
	print "Real Efficiency (186): "+str(RealEff186)
	
	g_effic.SetPoint(point,185.720,RealEff186)
	g_effic.SetPointError(point,0.0,0.1*RealEff186)
	point=point+1
	
	
	c11.SaveAs(workdir+"/"+fileprefix+"_186keV.C")
		
	t1.Delete()
#############################################################
############ Finished with the 235U 186 keV line ############
#############################################################



#########################################################
############ Starting with the 54Mn chain ###############
#########################################################
if Mn54:
	
	c12 = ROOT.TCanvas("c12","Simulated Spectra",870,500)
	c12.SetLogy()
	c12.cd()
	
	FullFN1 = DataPath+"/"+fileprefix+"_54Mn_"+events+".root"
	
	description1 = fileprefix+" ^{54}Mn Simulation"
	
	print "\nOpening <"+FullFN1+">"
	  
	t1 = ROOT.TChain("t1")
	t1.Add(FullFN1)
	
	h12 = ROOT.TH1F("h12",description1,27000,0,2700)
  
	h12.SetXTitle("Energy [keV]")
	h12.SetYTitle("counts")
	h12.SetStats(0)
	t1.Draw("GeEtot>>h12")
	
	
	################ Calculate the efficiencies ################

	#----------  834.838 keV from 54Mn (99.99 % BR)  -----------#
	#NOTE: This is the highest energy line.... no Compton affects it.
	binx = h12.FindBin(834.838)
	BRxEff835 = (h12.GetBinContent(binx) - (h12.Integral(binx-(compbins+1),binx-2)/compbins))/NumInEvt
	
	RealEff835 = BRxEff835 / BR835

	print "Efficiency x BR (835): "+str(BRxEff835)
	print "Real Efficiency (835): "+str(RealEff835)
	
	g_effic.SetPoint(point,834.838,RealEff835)
	g_effic.SetPointError(point,0.0,0.1*RealEff835)
	point=point+1
	

	c12.SaveAs(workdir+"/"+fileprefix+"_186keV.C")
	
	t1.Delete()
###################################################
############ Finished with 54Mn chain ############
###################################################	



effCanv = ROOT.TCanvas("effCanv","",870,500)
effCanv.SetLogy()
effCanv.cd()
g_effic.SetMarkerStyle(20)
g_effic.Draw("AP")
effCanv.SaveAs(workdir+"/"+fileprefix+"_efficplot.root")



###################################################################
########### Print out a file with all the efficiences #############
###################################################################
if writetxt:

	data_out=open(outfilename,"w") #This file will be the input of the sampleanalysis program
	  
	if U238: data_out.write("234\tTh\t92.6\t"+str(BR92)+"\t"+str(RealEff92)+"\t"+str(BRxEff92)+"\n")
	if U235: data_out.write("235\tU\t185.720\t"+str(BR186)+"\t"+str(RealEff186)+"\t"+str(BRxEff186)+"\n")
	if Ra228: data_out.write("212\tPb\t238.632\t"+str(BR239)+"\t"+str(RealEff239)+"\t"+str(BRxEff239)+"\n")
	if U238: data_out.write("214\tPb\t295.224\t"+str(BR295)+"\t"+str(RealEff295)+"\t"+str(BRxEff295)+"\n")
	if Ra228: data_out.write("228\tAc\t338.32\t"+str(BR339)+"\t"+str(RealEff339)+"\t"+str(BRxEff339)+"\n")
	if U238: data_out.write("214\tPb\t351.932\t"+str(BR352)+"\t"+str(RealEff352)+"\t"+str(BRxEff352)+"\n")
	if Ra228: data_out.write("208\tTl\t583.187\t"+str(BR583)+"\t"+str(RealEff583)+"\t"+str(BRxEff583)+"\n")
	if U238: data_out.write("214\tBi\t609.312\t"+str(BR609)+"\t"+str(RealEff609)+"\t"+str(BRxEff609)+"\n")
	if Cs137: data_out.write("137\tCs\t661.657\t"+str(BR662)+"\t"+str(RealEff662)+"\t"+str(BRxEff662)+"\n")
	if Mn54: data_out.write("54\tMn\t834.838\t"+str(BR835)+"\t"+str(RealEff835)+"\t"+str(BRxEff835)+"\n")
	if Sc46: data_out.write("46\tSc\t889.271\t"+str(BR889)+"\t"+str(RealEff889)+"\t"+str(BRxEff889)+"\n")
	if Ra228: data_out.write("228\tAc\t911.196\t"+str(BR911)+"\t"+str(RealEff911)+"\t"+str(BRxEff911)+"\n")
	if Ra228: data_out.write("228\tAc\t968.96\t"+str(BR969)+"\t"+str(RealEff969)+"\t"+str(BRxEff969)+"\n")
	if Pa: data_out.write("234m\tPa\t1001.026\t"+str(BR1001)+"\t"+str(RealEff1001)+"\t"+str(BRxEff1001)+"\n")
	if U238: data_out.write("214\tBi\t1120.287\t"+str(BR1120)+"\t"+str(RealEff1120)+"\t"+str(BRxEff1120)+"\n")
	if Sc46: data_out.write("46\tSc\t1120.537\t"+str(BR1120s)+"\t"+str(RealEff1120s)+"\t"++str(BRxEff1120s)+"\n")
	if Co60: data_out.write("60\tCo\t1173.228\t"+str(BR1173)+"\t"+str(RealEff1173)+"\t"+str(BRxEff1173)+"\n")
	if Co60: data_out.write("60\tCo\t1332.492\t"+str(BR1332)+"\t"+str(RealEff1332)+"\t"+str(BRxEff1332)+"\n")
	if K40: data_out.write("40\tK\t1460.882\t"+str(BR1460)+"\t"+str(RealEff1460)+"\t"+str(BRxEff1460)+"\n")
	if U238: data_out.write("214\tBi\t1764.494\t"+str(BR1765)+"\t"+str(RealEff1765)+"\t"+str(BRxEff1765)+"\n")
	if Ra228: data_out.write("208\tTl\t2614.511\t"+str(BR2615)+"\t"+str(RealEff2615)+"\t"+str(BRxEff2615)+"\n")
	
	data_out.close()
