// using namespace std;

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TSystem.h>

// Functions: kinematic quantities
Double_t calcT(P3MVector p, P3MVector pprime);
Double_t calcQ2(P3MVector k, P3MVector kprime);
Double_t calcBjorkenX(P3MVector k, P3MVector kprime, P3MVector p);
Double_t calcM2Miss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);

// Functions: undo afterburner
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k); // Undo procedure AND calculate boost vectors
void undoAfterburn(P3MVector& p); // Undo procedure ONLY

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

// Objects for undoing afterburn boost
Float_t fXAngle{-0.025}; // Crossing angle in radians
Float_t fRotX{};
RotationX rotAboutX;
Float_t fRotY{};
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

void DDIS_TDR(TString fileList){
  std::cout<< " __ __ __ __ __ __ __ __ __ __" <<std::endl;
  std::cout<< "|                             |"<<std::endl;
  std::cout<< "|     ePIC DDIS analysis      |"<<std::endl;
  std::cout<< "|__ __ __ __ __ __ __ __ __ __|"<<std::endl;
  std::cout<< "\nInput filelist: " << fileList <<std::endl;

  ifstream fileListStream;
  fileListStream.open(fileList);
  string fileName;

  //---------------------------------------------------------
  // CREATE TCHAIN FROM INPUT ROOT FILES
  //---------------------------------------------------------
  TChain* events = new TChain("events");
  Int_t nFiles{0};
  while(getline(fileListStream, fileName)){
    TString tmp = fileName;
    auto inputRootFile = TFile::Open(tmp);
    events->Add((TString)fileName);

    nFiles++;
  }
  std::cout<<"\nNo. of files: "<<nFiles<<"; no. of events: "<<events->GetEntries()<<std::endl;
  
  //---------------------------------------------------------
  // DECLARE OUTPUT HISTOGRAMS
  //---------------------------------------------------------
  // Pseudorapidity distributions

  // \eta electron
  TH1D* h_eta_MCe = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_RPe = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  // \eta gamma
  TH1D* h_eta_MCg = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_RPg = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);
  // \eta proton
  TH1D* h_eta_MCp  = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_RPp  = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPPp = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  // t distributions
  TH1D* h_t_MC     = new TH1D("t_MC"   , ";|t|(MC) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_t_RP     = new TH1D("t_RP"    , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_t_RPP    = new TH1D("t_RPP"   , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  // Photon angluar resolutions
  TH1D* h_PhotRes_theta   = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [rad]",600,-1.5,1.5);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [rad]; #delta#theta_{#gamma}",320,0,3.2,600,-1.5,1.5);

  //---------------------------------------------------------
  // DECLARE TTREEREADER AND BRANCHES TO USE
  //---------------------------------------------------------
  TTreeReader tree_reader(events);
  // MC particles
  TTreeReaderArray<float>  mc_px_array          = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float>  mc_py_array          = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float>  mc_pz_array          = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mc_mass_array        = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int>    mc_genStatus_array   = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<int>    mc_pdg_array         = {tree_reader, "MCParticles.PDG"};
  // Reconstructed/MC particle associations - BARREL (using ReconstructedParticles branch)
  TTreeReaderArray<unsigned int> assoc_rec_id   = {tree_reader, "ReconstructedParticleAssociations.recID"};
  TTreeReaderArray<unsigned int> assoc_sim_id   = {tree_reader, "ReconstructedParticleAssociations.simID"};
  // Reconstructed particles - BARREL (using ReconstructedParticles branch)
  TTreeReaderArray<float>  re_px_array          = {tree_reader, "ReconstructedParticles.momentum.x"};
  TTreeReaderArray<float>  re_py_array          = {tree_reader, "ReconstructedParticles.momentum.y"};
  TTreeReaderArray<float>  re_pz_array          = {tree_reader, "ReconstructedParticles.momentum.z"};
  TTreeReaderArray<float>  re_e_array           = {tree_reader, "ReconstructedParticles.energy"};
  TTreeReaderArray<int>    re_pdg_array         = {tree_reader, "ReconstructedParticles.PDG"};
  // Reconstructed/MC particle associations - B0 (using ReconstructedTruthSeededChargedParticles branch)
  TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
  TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
  // Reconstructed particles - B0 (using ReconstructedTruthSeededChargedParticles branch)
  TTreeReaderArray<float>  tsre_px_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
  TTreeReaderArray<float>  tsre_py_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
  TTreeReaderArray<float>  tsre_pz_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
  TTreeReaderArray<float>  tsre_e_array         = {tree_reader, "ReconstructedTruthSeededChargedParticles.energy"};
  TTreeReaderArray<float>  tsre_charge_array    = {tree_reader, "ReconstructedTruthSeededChargedParticles.charge"};
  // RP hits
  TTreeReaderArray<float> global_hit_RP_x       = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.x"};
  TTreeReaderArray<float> global_hit_RP_y       = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.y"};
  TTreeReaderArray<float> global_hit_RP_z       = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.z"};
  TTreeReaderArray<float> rp_px_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
  TTreeReaderArray<float> rp_py_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
  TTreeReaderArray<float> rp_pz_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
  TTreeReaderArray<float> rp_mass_array         = {tree_reader, "ForwardRomanPotRecParticles.mass"};
  TTreeReaderArray<int>   rp_pdg_array          = {tree_reader, "ForwardRomanPotRecParticles.PDG"};
  
  //---------------------------------------------------------
  // RUN OVER TTREEREADER 
  //---------------------------------------------------------
  // Get beams - run over full reader once
  // 4-vectors for beam particles - need these defined outside of file loop
  // Final beams
  P3MVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3MVector beamp4(0,0,0,-1);     // Beam proton (generated)
  // Accumulator variables (need these for calculation)
  P3MVector beame4_acc(0,0,0,-1);     // Beam electron (generated)
  P3MVector beamp4_acc(0,0,0,-1);     // Beam proton (generated)
  
  while (tree_reader.Next()){
    // Beams for each event
    P3MVector beame4_evt(0,0,0,-1);
    P3MVector beamp4_evt(0,0,0,-1);
    TVector3 mctrk;
    // Run over all MCParticles
    // Look for generator status 4 (beam particle), and beams by species
    for(int imc=0;imc<mc_px_array.GetSize();imc++){
      mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
      
      // Beam proton
      if(mc_genStatus_array[imc] == 4 && mc_pdg_array[imc] == 2212) beamp4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), fMass_proton);
      // Beam electron
      if(mc_genStatus_array[imc] == 4 && mc_pdg_array[imc] == 11  ) beame4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), fMass_electron);
    }// END OF MCPARTICLES LOOP
    
    // Add found beams to accumulator
    beame4_acc += beame4_evt;
    beamp4_acc += beamp4_evt;
  } // END OF TREE READER LOOP
    
  // Divide by number of events in file
  beame4.SetCoordinates(beame4_acc.X()/events->GetEntries(), beame4_acc.Y()/events->GetEntries(), beame4_acc.Z()/events->GetEntries(), beame4_acc.M()/events->GetEntries());
  beamp4.SetCoordinates(beamp4_acc.X()/events->GetEntries(), beamp4_acc.Y()/events->GetEntries(), beamp4_acc.Z()/events->GetEntries(), beamp4_acc.M()/events->GetEntries());

  // Run function to undo effect of afterburner and calculate "postburn" variables
  undoAfterburnAndCalc(beamp4,beame4);

  std::cout << "[DEBUG] Found beam energies " << beame4.E() << "x" << beamp4.E() << " GeV" << std::endl;
  /// BEAMS FOUND AND AFTERBURNER REMOVED

  // Main run to find particles and fill histograms
  // Restart TTreeReader
  tree_reader.Restart();
  while (tree_reader.Next()){
    // 4-vectors for MC raw particles
    vector<P3MVector> scate4_gen;   // Scattered electron (generated)
    vector<P3MVector> scatp4_gen;   // Scattered proton   (generated)
    vector<P3MVector> scatg4_gen;   // Scattered photon   (generated)

    // 4-vectors for associated MC particles (ONLY SCATTERED - only need photon for resolution plot)
    vector<P3MVector> scatg4_aso;   // Scattered photon   (associated MC)
    
    // 4-vectors for reconstructed particles (SEPARATE PROTONS FOR B0 AND ROMAN POTS)
    vector<P3MVector> scate4_rec;   // Scattered electron (reconstructed)
    vector<P3MVector> scatp4_rec;   // Scattered proton   (B0 reconstructed)
    vector<P3MVector> scatp4_rom;   // Scattered proton   (Roman Pots reconstructed)
    vector<P3MVector> scatg4_rec;   // Scattered photon   (reconstructed)      

    // Holding 3-vectors
    TVector3 mctrk, assoctrk, recotrk;
      
    //---------------------------------------------------------
    // Fill particle holding arrays
    //---------------------------------------------------------
    // 1. MC generated
    for(int imc=0;imc<mc_px_array.GetSize();imc++){
      mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
      P3MVector q_scat(mctrk.X(),mctrk.Y(),mctrk.Z(),mc_mass_array[imc]);
      
      // Undo afterburner
      undoAfterburn(q_scat);
	
      // Look for scattered particles ==> Generator status 1
      if(mc_genStatus_array[imc] == 1){
      	if(mc_pdg_array[imc] == 2212) scatp4_gen.push_back(q_scat);// proton
	      if(mc_pdg_array[imc] == 11  ) scate4_gen.push_back(q_scat);// electron
	      if(mc_pdg_array[imc] == 22  ) scatg4_gen.push_back(q_scat);// gamma
      } // Found scattered particles
    }// End of generated particles loop
      
    // 2. and 3. Associated MC tracks and their matched reconstucted tracks
    // USING EXPLICIT MATCHING BETWEEN RECO AND GENERATED TRACKS 
    unsigned int mc_assoc_index = -1;
    // LOOK FOR ELECTRONS AND PHOTONS (using ReconstructedParticleAssociations)
    for(unsigned int iAssoc{0};iAssoc<assoc_rec_id.GetSize();iAssoc++){
      mc_assoc_index = assoc_sim_id[iAssoc];
	
      // If reco track isn't associated to an MC track, then skip
      if(mc_assoc_index == -1) continue;
	
      assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]); 
      P3MVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),mc_mass_array[mc_assoc_index]);
      undoAfterburn(q_assoc);

      recotrk.SetXYZ(re_px_array[iAssoc], re_py_array[iAssoc], re_pz_array[iAssoc]);
      P3MVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),mc_mass_array[mc_assoc_index]);
      undoAfterburn(q_reco);
	
      // Fill track vectors based on associated PID
      // Electrons
      if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 11)	scate4_rec.push_back(q_reco); 
      // Photons
      if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 22){ 
	      scatg4_aso.push_back(q_assoc); 
	      scatg4_rec.push_back(q_reco); 
      }
    } // End of associations loop
    
    mc_assoc_index=-1; // Reset association index
    // THEN LOOK FOR PROTONS (using ReconstructedTruthSeededChargedParticleAssociations)
    for(unsigned int iTSAssoc{0};iTSAssoc<tsassoc_rec_id.GetSize();iTSAssoc++){
      mc_assoc_index = tsassoc_sim_id[iTSAssoc];
	
      // Only care about protons here (PID 2212)
      if(mc_assoc_index != -1 && mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 2212){
	      recotrk.SetXYZ(tsre_px_array[iTSAssoc], tsre_py_array[iTSAssoc], tsre_pz_array[iTSAssoc]);
	      P3MVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),mc_mass_array[mc_assoc_index]);
	      undoAfterburn(q_reco);
	      scatp4_rec.push_back(q_reco); 
      }
    } // End of truth-seeded association loop
    
    // Add in RP hits - only looking at protons
    // NO NEED TO UNDO AFTERBURNER FOR FF DETECTORS - NOT APPLIED IN FIRST PLACE
    for(int irpreco{0}; irpreco<rp_px_array.GetSize(); irpreco++){
      recotrk.SetXYZ(rp_px_array[irpreco], rp_py_array[irpreco], rp_pz_array[irpreco]);	
      
      P3MVector q_rpreco(recotrk.X(),recotrk.Y(),recotrk.Z(),rp_mass_array[irpreco]);
      if(rp_pdg_array[irpreco] == 2212){
	      scatp4_rom.push_back(q_rpreco);
      }
    }// End of RP reconstructed particles loop

    //---------------------------------------------------------
    // Fill histograms
    //---------------------------------------------------------
    // Generated particles
    // Need Q2 for electron cuts
    Float_t fQ2{0};
    if(scate4_gen.size() == 0) fQ2 = 0;
    else fQ2 = calcQ2(beame4, scate4_gen[0]);
    // Eta (inclusive tracks)
    if(fQ2 >= 1){
      // e'
      if(scate4_gen.size() == 1) h_eta_MCe->Fill(scate4_gen[0].Eta());
      // gamma
      if(scatg4_gen.size() == 1) h_eta_MCg->Fill(scatg4_gen[0].Eta());
      // p'
      if(scatp4_gen.size() == 1 && scatp4_gen[0].Theta()<0.02){
	      h_eta_MCp->Fill(scatp4_gen[0].Eta());
	      // Add exclusivity cuts for t-distribution
	      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1){
	        // Need to calculate kinematics before cutting on them
	        Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]);
	        Float_t ft = calcT(beamp4, scatp4_gen[0]);
          std:: cout << "[DEBUG]: 279: " << ft << std::endl;
	    
	        // Want MM2 to be close to zero
	        if(TMath::Abs(fM2miss) < 1) h_t_MC->Fill(ft);
	      } // Exclusivity required
      } // Proton tracks done
    } // Q2 limit
    
    // Reconstructed particles
    // Need Q2 for electron cuts
    if(scate4_rec.size() == 0) fQ2 = 0;
    else fQ2 = calcQ2(beame4, scate4_rec[0]);
    // Eta (inclusive tracks)
    if(fQ2 >= 1){
      // e'
      if(scate4_rec.size() == 1) h_eta_RPe->Fill(scate4_rec[0].Eta());
      // gamma
      if(scatg4_rec.size() == 1) h_eta_RPg->Fill(scatg4_rec[0].Eta());
      // p' (B0 - theta between 5.5 and 20 mrad)
      if(scatp4_rec.size() == 1 && scatp4_rec[0].Theta()>0.0055 && scatp4_rec[0].Theta()<0.02){
	      h_eta_RPp->Fill(scatp4_rec[0].Eta());
	      // Add exclusivity cuts for t-distribution
	      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
	        // Need to calculate kinematics before cutting on them
	        Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]);
	        Float_t ft = calcT(beamp4, scatp4_rec[0]);
          std:: cout << "[DEBUG]: 304: " << ft << std::endl;
	    
	        // Want MM2 to be close to zero
	        if(TMath::Abs(fM2miss) < 1) h_t_RP->Fill(ft);
	      } // Exclusivity required
      } // B0 Proton tracks found
      // p' (RP - theta less than 5 mrad)
      if(scatp4_rom.size() == 1 && scatp4_rom[0].Theta()<0.005){
	      h_eta_RPPp->Fill(scatp4_rom[0].Eta());
	      // Add exclusivity cuts for t-distribution
	      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
	        // Need to calculate kinematics before cutting on them
	        Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]);
	        Float_t ft = calcT(beamp4, scatp4_rom[0]);
          std:: cout << "[DEBUG]: 317: " << ft << std::endl;
	    
	        // Want MM2 to be close to zero
	        if(TMath::Abs(fM2miss) < 1 && ft < 0.3) h_t_RPP->Fill(ft);
	      } // Exclusivity required
      } // RP Proton tracks found
    } // Q2 limit
    
    if(scatg4_aso.size()==1 && scatg4_rec.size()==1){
      h_PhotRes_theta->Fill(scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
      h_PhotRes2D_theta->Fill(scatg4_aso[0].Theta(), scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
    }

  } // END OF TREE READER LOOP
  

  //---------------------------------------------------------
  // DRAW HISTOGRAMS
  //---------------------------------------------------------
  gStyle->SetOptStat(00000000);
  // Set header strings
  TString sHead1 = "#it{#bf{ePIC} Simulation}";
  TString sHead2 = "#it{ep 10x100 GeV}";

  // CANVAS 1: ETA DISTRIBUTIONS (separated by species)
  TCanvas* c1 = new TCanvas("c1","",1200,500);
  c1->Divide(3,1,0,0);
  TLatex* tHeadLc1 = new TLatex(0.04, 0.91, sHead1);
  tHeadLc1->SetNDC();
  tHeadLc1->SetTextSize(25);
  tHeadLc1->SetTextFont(43);
  tHeadLc1->SetTextColor(kBlack);
  tHeadLc1->Draw();
  TLatex* tHeadM1c1 = new TLatex(0.29, 0.91, "#bf{EpIC} ep #rightarrow e'p'#gamma");
  tHeadM1c1->SetNDC();
  tHeadM1c1->SetTextSize(25);
  tHeadM1c1->SetTextFont(43);
  tHeadM1c1->SetTextColor(kBlack);
  tHeadM1c1->Draw("same");
  TLatex* tHeadM2c1 = new TLatex(0.59, 0.91, "Q^{2} #geq 1 GeV^{2}");
  tHeadM2c1->SetNDC();
  tHeadM2c1->SetTextSize(25);
  tHeadM2c1->SetTextFont(43);
  tHeadM2c1->SetTextColor(kBlack);
  tHeadM2c1->Draw("same");
  TLatex* tHeadRc1 = new TLatex(0.83, 0.91, sHead2);
  tHeadRc1->SetNDC();
  tHeadRc1->SetTextSize(25);
  tHeadRc1->SetTextFont(43);
  tHeadRc1->SetTextColor(kBlack);
  tHeadRc1->Draw("same");
  c1->cd(1);
  gPad->SetLogy(1);
  h_eta_MCe->SetMaximum(10*h_eta_MCe->GetMaximum());
  h_eta_MCe->GetXaxis()->SetTitle("#eta");
  h_eta_MCe->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCe->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCe->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCe->GetXaxis()->SetLabelOffset(0.01);
  h_eta_MCe->GetXaxis()->CenterTitle();  
  h_eta_MCe->GetYaxis()->SetTitle("Counts");
  h_eta_MCe->GetYaxis()->SetTitleSize(0.05);
  h_eta_MCe->GetYaxis()->SetTitleOffset(0.95);
  h_eta_MCe->GetYaxis()->SetLabelSize(0.04);
  h_eta_MCe->GetYaxis()->SetLabelOffset(0.002);
  h_eta_MCe->SetLineColor(kBlack);
  h_eta_MCe->SetLineWidth(2);
  h_eta_RPe->SetLineColor(kBlack);
  h_eta_RPe->SetMarkerColor(kBlack);
  h_eta_RPe->SetMarkerStyle(20);
  h_eta_MCe->Draw();
  h_eta_RPe->Draw("pesame");
  // Add text
  TLatex* tHeadc1p1 = new TLatex(0.7, 0.92, "e^{-}");
  tHeadc1p1->SetNDC();
  tHeadc1p1->SetTextSize(30);
  tHeadc1p1->SetTextFont(43);
  tHeadc1p1->SetTextColor(kBlack);
  tHeadc1p1->Draw("same");
  // Add legend
  TLegend* lC1p1 = new TLegend(0.7,0.7,0.95,0.9);
  lC1p1->SetLineColorAlpha(kWhite,0);
  lC1p1->SetLineWidth(0);
  lC1p1->SetFillStyle(0);
  lC1p1->SetFillColorAlpha(kWhite,0);
  lC1p1->AddEntry(h_eta_MCe, "MC gen.", "pl");
  lC1p1->AddEntry(h_eta_RPe, "Reco.", "pl");
  lC1p1->Draw();
  c1->cd(2);
  gPad->SetLogy();
  h_eta_MCg->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCg->GetXaxis()->SetTitle("#eta");
  h_eta_MCg->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCg->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCg->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCg->GetXaxis()->SetLabelOffset(0.01); 
  h_eta_MCg->SetLineColor(kRed);
  h_eta_MCg->SetLineWidth(2);
  h_eta_RPg->SetLineColor(kRed);
  h_eta_RPg->SetMarkerColor(kRed);
  h_eta_RPg->SetMarkerStyle(20);
  h_eta_MCg->Draw();
  h_eta_RPg->Draw("pesame");
  // Add text
  TLatex* tHeadc1p2 = new TLatex(0.07, 0.92, "#gamma");
  tHeadc1p2->SetNDC();
  tHeadc1p2->SetTextSize(30);
  tHeadc1p2->SetTextFont(43);
  tHeadc1p2->SetTextColor(kBlack);
  tHeadc1p2->Draw("same");
  c1->cd(3);
  gPad->SetLogy();
  h_eta_MCp->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCp->GetXaxis()->SetTitle("#eta");
  h_eta_MCp->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCp->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCp->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCp->GetXaxis()->SetLabelOffset(0.01);
  h_eta_MCp->SetLineColor(kBlue);
  h_eta_MCp->SetLineWidth(2);
  h_eta_RPp->SetLineColor(kViolet);
  h_eta_RPp->SetMarkerColor(kViolet);
  h_eta_RPp->SetMarkerStyle(20);
  h_eta_RPPp->SetLineColor(kCyan+1);
  h_eta_RPPp->SetMarkerColor(kCyan+1);
  h_eta_RPPp->SetMarkerStyle(20);
  h_eta_MCp->Draw();
  h_eta_RPp->Draw("pesame");
  h_eta_RPPp->Draw("pesame");
  // Add text
  TLatex* tHeadc1p3 = new TLatex(0.07, 0.92, "p");
  tHeadc1p3->SetNDC();
  tHeadc1p3->SetTextSize(30);
  tHeadc1p3->SetTextFont(43);
  tHeadc1p3->SetTextColor(kBlack);
  tHeadc1p3->Draw("same");
  // Add legend
  TLegend* lC1p3 = new TLegend(0.03,0.7,0.7,0.9);
  lC1p3->SetLineColorAlpha(kWhite,0);
  lC1p3->SetLineWidth(0);
  lC1p3->SetFillStyle(0);
  lC1p3->SetFillColorAlpha(kWhite,0);
  lC1p3->AddEntry(h_eta_RPp, "Reco. B0", "pl");
  lC1p3->AddEntry(h_eta_RPPp, "Reco. RP", "pl");
  lC1p3->Draw();

  // CANVAS 2: T-DISTRIBUTION
  // TCanvas* c2 = new TCanvas("c2","",1200,800);
  // h_t_MC->SetMinimum(1);
  // gPad->SetLogy(1);
  // h_t_MC->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  // h_t_MC->GetYaxis()->SetTitle("Counts / 0.02 GeV^{2}");
  // h_t_MC->SetLineColor(kBlack);
  // h_t_MC->SetLineWidth(2);
  // h_t_MC->Draw();
  // h_t_RP->SetLineColor(kBlue);
  // h_t_RP->SetMarkerColor(kBlue);
  // h_t_RP->SetMarkerStyle(20);
  // h_t_RP->Draw("pesame");
  // h_t_RPP->SetLineColor(kCyan+1);
  // h_t_RPP->SetMarkerColor(kCyan+1);
  // h_t_RPP->SetMarkerStyle(20);
  // h_t_RPP->Draw("pesame");
  // // Create header text objects
  // TLatex* tHead1 = new TLatex(0.10, 0.91, sHead1);
  // tHead1->SetNDC();
  // tHead1->SetTextSize(30);
  // tHead1->SetTextFont(43);
  // tHead1->SetTextColor(kBlack);
  // TLatex* tHead2 = new TLatex(0.74, 0.91, sHead2);
  // tHead2->SetNDC();
  // tHead2->SetTextSize(30);
  // tHead2->SetTextFont(43);
  // tHead2->SetTextColor(kBlack);
  // tHead1->Draw("same");
  // tHead2->Draw("same");
  // // Add legend
  // TLatex* tC2 = new TLatex(0.58, 0.83, "#splitline{#bf{EpIC} ep #rightarrow e'p'#gamma, Q^{2} #geq 1 GeV^{2}}{t_{RP} #leq 0.3 GeV^{2}, M_{miss}^{2} < 1 GeV^{2}}");
  // tC2->SetNDC();
  // tC2->SetTextSize(30);
  // tC2->SetTextFont(43);
  // tC2->SetTextColor(kBlack);
  // tC2->Draw("same");
  // TLegend* lC2 = new TLegend(0.57, 0.6, 0.8, 0.77);
  // lC2->SetLineColorAlpha(kWhite,0);
  // lC2->SetFillColorAlpha(kWhite,0);
  // lC2->AddEntry(h_t_MC, "#bf{EpIC} MC gen.", "l");
  // lC2->AddEntry(h_t_RP, "Reco. B0", "lp");
  // lC2->AddEntry(h_t_RPP, "Reco. RP", "lp");
  // lC2->Draw();

  // CANVAS 3A: PHOTON ANGLE RESOLUTION (INSET)
  TCanvas* c3 = new TCanvas("c3","",1200,800);
  // Draw 1D plot as main figure
  h_PhotRes_theta->GetXaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetXaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->GetYaxis()->SetTitle("Counts / 5mrad");
  h_PhotRes_theta->GetYaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetYaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->SetLineColor(kBlack);
  h_PhotRes_theta->SetMarkerColor(kBlack);
  h_PhotRes_theta->SetMarkerStyle(20);
  h_PhotRes_theta->Draw("pe");
  // Create subfigure TPad on top for 2D plot
  TPad* p3 = new TPad("p3","",0.51,0.49,0.89,0.89);
  p3->SetRightMargin(0.01);
  p3->SetLeftMargin(0.15);
  p3->SetTopMargin(0.02);
  p3->SetBottomMargin(0.15);
  p3->Draw("same");
  TLatex* tHeadLc3 = new TLatex(0.1, 0.94, sHead1);
  tHeadLc3->SetNDC();
  tHeadLc3->SetTextSize(30);
  tHeadLc3->SetTextFont(43);
  tHeadLc3->SetTextColor(kBlack);
  tHeadLc3->Draw("same");
  //tHead2->Draw("same");
  p3->cd();
  h_PhotRes2D_theta->RebinX(2);
  h_PhotRes2D_theta->RebinY(8);
  h_PhotRes2D_theta->GetXaxis()->SetTitleSize(0.08);
  h_PhotRes2D_theta->GetXaxis()->SetTitleOffset(0.85);
  h_PhotRes2D_theta->GetXaxis()->SetLabelSize(0.06);
  h_PhotRes2D_theta->GetYaxis()->SetTitle("#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [rad]");
  h_PhotRes2D_theta->GetYaxis()->SetTitleSize(0.08);
  h_PhotRes2D_theta->GetYaxis()->SetTitleOffset(0.85);
  h_PhotRes2D_theta->GetYaxis()->SetLabelSize(0.06);
  h_PhotRes2D_theta->GetZaxis()->SetMaxDigits(2);
  h_PhotRes2D_theta->Draw("same");


  //---------------------------------------------------------
  // SAVE CANVASES AS PNGs AND CLOSE
  //---------------------------------------------------------
  c1->SaveAs("figs/DDIS_eta.png");
  c1->Close();
  // c2->SaveAs("figs/DDIS_t.png");
  // c2->Close();
  c3->SaveAs("figs/DDIS_photres.png");
  c3->Close();

  return;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// Calculate t from beam and scattered proton vectors
// t = (p' - p)^2
Double_t calcT(P3MVector p, P3MVector pprime){
  double t = (pprime - p).M2();
  
  return TMath::Abs(t);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// Calculate Q2 from beam and scattered electron vectors
// Q2 = (k - k')^2
Double_t calcQ2(P3MVector k, P3MVector kprime){
  double q2 = (k - kprime).M2();
  double Q2 = -q2;

  return Q2;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// Calculate Bjorken x from both beam vectors and scattered electron vector
// xB = Q2 / 2(q.p)
Double_t calcBjorkenX(P3MVector k, P3MVector kprime, P3MVector p){
  P3MVector q = k - kprime;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;

  return xB;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// Calculation of missing mass (squared) for fully exclusive e'p'g final state
Double_t calcM2Miss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
// USE BEAM VECTORS
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
  // Holding vectors for beam - undoing crossing angle ONLY
  P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
  P3MVector e_beam(0., 0., -k.E(), k.M());
  
  // Define boost vector to CoM frame
  P3MVector CoM_boost = p_beam+e_beam;
  vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
  
  // Apply boost to beam vectors
  p_beam = boost(p_beam, vBoostToCoM);
  e_beam = boost(e_beam, vBoostToCoM);
  
  // Calculate rotation angles and create rotation objects
  fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

  rotAboutY = RotationY(fRotY);
  rotAboutX = RotationX(fRotX);

  // Apply rotation to beam vectors
  p_beam = rotAboutY(p_beam);
  p_beam = rotAboutX(p_beam);
  e_beam = rotAboutY(e_beam);
  e_beam = rotAboutX(e_beam);

  // Define boost vector back to head-on frame
  P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
  vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

  // Apply boost back to head on frame to beam vectors
  p_beam = boost(p_beam, vBoostToHoF);
  e_beam = boost(e_beam, vBoostToHoF);

  // Make changes to input vectors
  p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
  k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Undo afterburn procedure only
void undoAfterburn(P3MVector& a){
  // Undo AB procedure for single vector, a^{mu}
  a = boost(a, vBoostToCoM); // BOOST TO COM FRAME
  a = rotAboutY(a);          // ROTATE TO Z-AXIS
  a = rotAboutX(a);          // ROTATE TO Z-AXIS
  a = boost(a, vBoostToHoF); // BOOST BACK TO HEAD ON FRAME
}
