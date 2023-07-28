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

def drawHistOverlay(h_list, xtitle, filename):
    r.gStyle.SetPalette(r.kRainBow)
    can = r.TCanvas()
    can.Draw()
    
    for h in h_list:
        r.gStyle.SetOptStat(0)
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle("Events")
        h.Draw("SAME PLC PMC")

    r.gPad.BuildLegend()
    can.SaveAs(filename+'.pdf')
    
    
def main():
    muonTypes = ['muon', 'muon1Leg', 'muonEndCapsOnly',\
                'muonNoRPC', 'muonWitht0Correction', 'splitmuon', 'lhcSTAmuon']
    
    h_list_eta = []
    h_list_phi = []
    h_list_pt = []
    h_list_TrackVy_philt0 = []
    h_list_TrackVy_phigt0 = []
    h_list_TrackPy_philt0 = []
    h_list_TrackPy_phigt0 = []
    
    for prcType in muonTypes:
        # load tree
        tree = r.TChain("{}Analyzer/MuonTree".format(prcType))
        tree.Add("muon_ntuple.root")
        nevents = tree.GetEntries()
        print("nentries = ", nevents)
        
        h_muonPhi = r.TH1D( 'h_{}Phi'.format(prcType), prcType, 50, -3.14, 3.14)
        h_muonPt = r.TH1D( 'h_{}Pt'.format(prcType), prcType, 100, 0, 1000)
        h_muonEta = r.TH1D( 'h_{}Eta'.format(prcType), prcType, 50, -1.5, 1.5)
        
        h_muonPtvsPhi = r.TH2D( 'h_{}PtvsPhi'.format(prcType), prcType, 50, -3.14, 3.14, 100, 0, 1000)
        h_muonPtvsEta = r.TH2D( 'h_{}PtvsEta'.format(prcType), prcType, 50, -1.5, 1.5, 100, 0, 1000)
        
        # h_cosmicMuonPhi = r.TH1D( 'h_cosmic{}Phi'.format(prcType), prcType, 50, -3.14, 3.14)
        h_cosmicMuonTrackVy_philt0 = r.TH1D( 'h_cosmic{}TrackVy_philt0'.format(prcType), prcType, 100, -500, 500)
        h_cosmicMuonTrackVy_phigt0 = r.TH1D( 'h_cosmic{}TrackVy_phigt0'.format(prcType), prcType, 100, -500, 500)
        h_cosmicMuonTrackPy_philt0 = r.TH1D( 'h_cosmic{}TrackPy_philt0'.format(prcType), prcType, 110, -500, 50)
        h_cosmicMuonTrackPy_phigt0 = r.TH1D( 'h_cosmic{}TrackPy_phigt0'.format(prcType), prcType, 110, -500, 50)
        
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
        
        drawHist(h_muonPhi, "#phi", "plots/h_{}_Phi".format(prcType))
        drawHist(h_muonPt, "p_{T} [GeV]", "plots/h_{}_Pt".format(prcType))
        drawHist(h_muonEta, "#eta", "plots/h_{}_Eta".format(prcType))
        drawHist2d(h_muonPtvsPhi, "#phi", "p_{T} [GeV]", "plots/h_{}_PtvsPhi".format(prcType))
        drawHist2d(h_muonPtvsEta, "#eta", "p_{T} [GeV]", "plots/h_{}_PtvsEta".format(prcType))
        
        tree = r.TChain("{}Analyzer/CosmicMuonTree".format(prcType))
        tree.Add("muon_ntuple.root")
        
        for event in tree:
            for (track_phi, track_vy, track_py) in zip(event.track_phi, event.track_vy, event.track_py):
                if track_phi < 0:
                    h_cosmicMuonTrackVy_philt0.Fill(track_vy)
                    h_cosmicMuonTrackPy_philt0.Fill(track_py)
                else:
                    h_cosmicMuonTrackVy_phigt0.Fill(track_vy)
                    h_cosmicMuonTrackPy_phigt0.Fill(track_py)
        
        h_list_eta.append(h_muonEta)
        h_list_phi.append(h_muonPhi)
        h_list_pt.append(h_muonPt)
        h_list_TrackVy_philt0.append(h_cosmicMuonTrackVy_philt0)
        h_list_TrackVy_phigt0.append(h_cosmicMuonTrackVy_phigt0)
        h_list_TrackPy_philt0.append(h_cosmicMuonTrackPy_philt0)
        h_list_TrackPy_phigt0.append(h_cosmicMuonTrackPy_phigt0)
    
    drawHistOverlay(h_list_eta, "#eta", "plots/h_Eta")
    drawHistOverlay(h_list_phi, "#phi", "plots/h_Phi")
    drawHistOverlay(h_list_pt, "p_{T} [GeV]", "plots/h_Pt")
    drawHistOverlay(h_list_TrackVy_philt0, "track vy [mm]", "plots/h_TrackVy_philt0")
    drawHistOverlay(h_list_TrackVy_phigt0, "track vy [mm]", "plots/h_TrackVy_phigt0")
    drawHistOverlay(h_list_TrackPy_philt0, "track P_{y} [GeV]", "plots/h_TrackPy_philt0")
    drawHistOverlay(h_list_TrackPy_phigt0, "track P_{y} [GeV]", "plots/h_TrackPy_phigt0")
    
if __name__=="__main__":
    main()