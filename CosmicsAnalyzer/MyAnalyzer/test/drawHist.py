import numpy as np
import ROOT as r

def drawHist(h, xtitle, filename):
    '''
        Feed in an array, nbins, low & up range, x-title, filename
        Plot & return the histgram
    '''
    # r.gStyle.SetOptStat(0)

    can = r.TCanvas()
    can.Draw()
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle("Events")
    h.Draw()
    can.SaveAs(filename+'.pdf')
    return h

def drawHist2d(h, xtitle, ytitle, filename):
    '''
        Feed in an array, nbins, low & up range, x-title, filename
        Plot & return the histgram
    '''
    r.gStyle.SetOptStat(0)

    can = r.TCanvas()
    can.Draw()
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.Draw("colz")
    can.SaveAs(filename+'.pdf')
    return h

def main():
    # load tree
    tree = r.TChain("muonPhiAnalyzer/MuonTree")
    tree.Add("muon_phi_ntuple.root")
    nevents = tree.GetEntries()
    print("nentries = ", nevents)
    
    h_muonPhi = r.TH1D( 'h_muonPhi', '', 50, -3.14, 3.14)
    h_muonPt = r.TH1D( 'h_muonPt', '', 100, 0, 3000)
    h_muonEta = r.TH1D( 'h_muonEta', '', 50, -1.5, 1.5)
    
    h_muonPtvsPhi = r.TH2D( 'h_muonPtvsPhi', '', 50, -3.14, 3.14, 100, 0, 3000)
    h_muonPtvsEta = r.TH2D( 'h_muonPtvsEta', '', 50, -1.5, 1.5, 100, 0, 3000)
    
    event_count = 0
    pfreq = 1000
    
    for event in tree:
        
        if event_count%pfreq == 0:
            print('Processing Event: %s'%(event_count))
        
        for phi in event.muonPhi:
            h_muonPhi.Fill(phi)
        for pt in event.muonPt:
            h_muonPt.Fill(pt)
        for eta in event.muonEta:
            h_muonEta.Fill(eta)
        
        for (pt, phi, eta) in zip(event.muonPt, event.muonPhi, event.muonEta):
            h_muonPtvsPhi.Fill(phi, pt)
            h_muonPtvsEta.Fill(eta, pt)
        
        event_count += 1
    
    drawHist(h_muonPhi, "#phi", "plots/h_muonPhi")
    drawHist(h_muonPt, "p_{T} [GeV]", "plots/h_muonPt")
    drawHist(h_muonEta, "#eta", "plots/h_muonEta")
    drawHist2d(h_muonPtvsPhi, "#phi", "p_{T} [GeV]", "plots/h_muonPtvsPhi")
    drawHist2d(h_muonPtvsEta, "#eta", "p_{T} [GeV]", "plots/h_muonPtvsEta")
    
if __name__=="__main__":
    main()