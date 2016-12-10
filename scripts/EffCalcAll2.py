import ROOT

class isotope_list:
    isotopes = dict()
    
    def AddIsotope(self, isot):
        if isot.name in isotope_list.isotopes:
            print '\nIsotope ' + name + ' already present in the list!\n'
            return
        isotope_list.isotopes[isot.name] = isot
    
class IsotLine:
    def __init__(self,name,mass,element,energy,br):
        self.name = name
        self.mass = mass
        self.element = element
        self.energy = energy
        self.br = br
        self.Eff = self.BrXeff = 0.
        self.combined = False
        self.comblines = [] 
    
    def SetEffic(self, brxeff):
        self.Eff = brxeff/self.br
        self.BrXeff = brxeff
        
        print '\nLine at ' + str(self.energy) + 'keV:'
        print 'BR: ' + str(self.br)
        print 'Eff: ' + str(self.Eff)
        print 'BR x Eff: ' + str(self.BrXeff)
    
    
class isotope:
    def __init__(self, name):
        if name in isotope_list.isotopes:
            print '\nIsotope ' + name + ' already present in the list!\n'
            return None
        self.name = name
        self.lines = {}
        isotope_list().AddIsotope(self)

    
    def AddLine(self, linename, mass, element, energy, br):
        if linename in self.lines:
            print '\nLine "' + linename + '" already present\n'
            return
        self.lines[linename] = IsotLine(linename, mass, element, energy, br)
    
    def CombLines(self, linename, lineslist):
        if len(lineslist) < 2:
            return
        nlines = len(lineslist)
        if linename in lines:
            print '\nLine "' + linename + '" already present\n'
            return
        
        self.lines[linename] = IsotLine(linename, mass, element, energy, br)
        self.lines[linename].combined = True
        self.lines[linename].comblines = lineslist


#Uncomment the line below if you want to run in batch
#ROOT.gROOT.SetBatch()

ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetFuncColor(2)
  
#Here are the setup and initialization stuff to be changed manually
datadir = "/home/atp/fpiastra/g4sims/Screening/xebottle/"
sample="SteelBottle"
events = "1e8"
workdir = datadir + sample + "/effic"
outfilename = "lines.list"
writetxt = False
compbins = 10 #number of bins left and right of lines in which estimate the compton BG


U238 = isotope("238U")
U238.AddLine("92a", "234", "Th", 92.38, 0.0215)
U238.AddLine("92b", "234", "Th", 92.75, 0.0218)
#U238.CombLines("92eff", ["92a", "92b"])
U238.AddLine("1001", "234m", "Pa", 1001.026, 0.00846)

Ra226 = isotope("226Ra")
Ra226.AddLine("295", "214", "Pb", 295.224, 0.184)
Ra226.AddLine("352", "214", "Pb", 351.224, 0.356)
Ra226.AddLine("609", "214", "Bi", 609.312, 0.4549)
Ra226.AddLine("1120", "214", "Bi", 1120.287, 0.1491)
Ra226.AddLine("1764", "214", "Bi", 1764.494, 0.1531)

Ra228 = isotope("228Ra")
Ra228.AddLine("338", "228", "Ac", 338.32, 0.114)
Ra228.AddLine("911", "228", "Ac", 911.196, 0.262)
Ra228.AddLine("969", "228", "Ac", 968.96, 0.159)

Th228 = isotope("228Th")
Th228.AddLine("239", "212", "Pb", 238.632, 0.436)
Th228.AddLine("583", "208", "Tl", 583.187, 0.3054)
Th228.AddLine("2615", "208", "Tl", 2614.511, 0.3584)

K40 = isotope("40K")
K40.AddLine("1461", "40", "K", 1460.88, 0.1055)

Co60 = isotope("60Co")
Co60.AddLine("1173", "60", "Co", 1173.23, 0.9985)
Co60.AddLine("1332", "60", "Co", 1332.49, 0.9998)

Sc46 = isotope("46Sc")
Sc46.AddLine("889", "46", "Sc", 889.271, 0.99984)
Sc46.AddLine("1121", "46", "Sc", 1120.537, 0.99987)

Mn54 = isotope("54Mn")
Mn54.AddLine("835","54","Mn", 834.838 , 0.999746)


nevents = float(events)

point = 0 #Point number of the efficiency TGraphErrors

g_effic = ROOT.TGraphErrors()


for isotname in isotope_list.isotopes:
    isot = isotope_list.isotopes[isotname]
    chainfiles = datadir + sample + "/" + isotname + "/*.root"
    print "\n\nOpening <"+chainfiles+">"
    t1 = ROOT.TChain("t1")
    t1.Add(chainfiles)
    
    print '\n\nIsotope (chain) ' + isotname + ":"
    print '-----------------'
    
    descr = isot.name + "; Energy [keV]; Counts"
    histname = isot.name+"hist"
    histo = ROOT.TH1F(histname, descr, 27000, 0, 2700)
    histo.SetStats(0)
    t1.Draw("GeEtot>>"+histname,"GeEtot>0")
    
    for linename in isot.lines:
        iLine = isot.lines[linename]
        binc = histo.FindBin(iLine.energy)
        binl = binc - 2 #start from here going downward in energy
        binr = binc + 2 #start from here going upward in energy
        
        leftcounts = rightcounts = 0.
        linecounts = histo.GetBinContent(binc)
        
        for iCnt in range(compbins):
            leftcounts += histo.GetBinContent(binl-iCnt)/float(compbins)
            rightcounts += histo.GetBinContent(binr+iCnt)/float(compbins)
        
        if rightcounts > 0.:
            effxbr = (linecounts - (leftcounts+rightcounts)/2.)/nevents
        else:
            effxbr = (linecounts - leftcounts)/nevents
        
        
        iLine.SetEffic(effxbr)
        
        g_effic.SetPoint(point, iLine.energy, iLine.Eff)
        g_effic.SetPointError(point, 0.0, 0.1*iLine.Eff)
        point=point+1


if writetxt:
    data_out=open(outfilename,"w") #This file will be the input of the sampleanalysis program
    
    for isotname in isotope_list.isotopes:
        isot = isotope_list.isotopes[isotname]
        for linename in isot.lines:
            iLine = isot.lines[linename]
            data_out.write(iLine.mass + "\t" + iLine.element + "\t" + str(iLine.br) + "\t" + str(iLine.Eff) + "\t" + str(iLine.BrXeff)  )
