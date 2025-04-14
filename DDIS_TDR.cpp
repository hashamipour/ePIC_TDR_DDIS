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
#include <filesystem> // Include the necessary header for file system operations

// Functions: kinematic quantities
Double_t calcT(P3MVector p, P3MVector pprime);
Double_t calcQ2(P3MVector k, P3MVector kprime);
Double_t calcBjorkenX(P3MVector k, P3MVector kprime, P3MVector p);
Double_t calcM2Miss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);

// Functions: undo afterburner
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k); // Undo procedure AND calculate boost vectors
void undoAfterburn(P3MVector& p); // Undo procedure ONLY
double rapidity(const P3MVector& p4); //rapidity by hand
double calculatePseudorapidity(const P3MVector& p4); // pseudo-rapidity by hand

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};
double MASS_PROTON   = fMass_proton;
double MASS_ELECTRON = fMass_electron;

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

    if (!std::filesystem::exists(tmp.Data())) {
      std::cerr << "\nError: File does not exist: " << fileName << std::endl;
      continue; // Skip to the next file
    }

    auto inputRootFile = TFile::Open(tmp);
    events->Add((TString)fileName);
    inputRootFile->Close(); // It's good practice to close the file after adding it to the TChain

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

  double tmin = 1e-3;
  double tmax = 2.0;

  // Define number of bins
  const int nLogBins = 40; // you can adjust this

  // Log step calculation
  double logStep = TMath::Power(tmax / tmin, 1.0 / nLogBins);

  // Array for bin edges
  double t_edges[nLogBins + 1];

  t_edges[0] = tmin;
  for (int i = 1; i <= nLogBins; ++i) {
    t_edges[i] = t_edges[i-1] * logStep;
  }
  TH1D* h_t_MC     = new TH1D("t_MC"    , ";|t|(MC) [(GeV/c^{2})^{2}]"    , nLogBins, t_edges);
  TH1D* h_t_RP     = new TH1D("t_RP"    , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , nLogBins, t_edges);
  TH1D* h_t_RPP    = new TH1D("t_RPP"   , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , nLogBins, t_edges);
  // t resolution
  TH1D* h_PrRes_t        = new TH1D("ProtonRes_t",";t(Reco)-t(MC) [GeV^{2}]",100,-0.5,0.5);
  TH1D* h_PrRelRes_t     = new TH1D("ProtonRelRes_t","BABE;#frac{t(Reco)-t(MC)}{t(MC)}",20,-0.2,0.2);
    // Define bin edges
    Double_t binEdges[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5}; // Example bin edges
    Int_t nBins = sizeof(binEdges) / sizeof(Double_t) - 1; // Number of bins
  TH2D* h_tRelRes_binned = new TH2D("tRelRes_Q2",";Q^{2} [GeV^{2}];#frac{t(Reco)-t(MC)}{t(MC)}",nBins,binEdges, 100,-0.5,0.5);
  // t correlation
  TH2D* h_Corr_t = new TH2D("Corr_t", ";t_{MC};t_{Reco}", 100, 0, 2, 100, 0, 2);
  // Q2 resolution
  TH1D* h_Res_Q2      = new TH1D("Q2_Res",";Q^{2}(Reco)-Q^{2}(MC) [GeV^{2}]",100,-2.5,2.5);
  TH1D* h_RelRes_Q2   = new TH1D("Q2_RelRes","electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",200,-0.2,0.2);
  // Q2 resolution binned
  TH2D* h_RelRes_Q2_binned = new TH2D("tRelRes_binned",";Q^{2} [GeV^{2}];#frac{t(Reco)-t(MC)}{t(MC)}",7,5,40, 100,-0.5,0.5);
  // Q2 correlation
  TH2D* h_Corr_Q2 = new TH2D("Corr_Q2", ";Q^{2}_{MC};Q^{2}_{Reco}", 100, 5, 40, 100, 5, 40);
  // xBj resolution
  TH1D* h_Res_xBj      = new TH1D("xBj_Res",";x_{Bj}(Reco)-x_{Bj}(MC)",200,-0.0002,0.0002);
  TH1D* h_RelRes_xBj   = new TH1D("xBj_RelRes","electron method;#frac{x_{Bj}(Reco)-x_{Bj}(MC)}{x_{Bj}(MC)}",50,-0.03,0.03);
  // xBj resolution binned
  TH2D* h_xBjRelRes_binned = new TH2D("xBjRelRes_binned",";x_{Bj};#frac{x_{Bj}(Reco)-x_{Bj}(MC)}{x_{Bj}(MC)}",5,0.0,0.2,100,-0.5,0.5);
  // xBj correlation
  TH2D* h_Corr_xBj = new TH2D("Corr_xBj", ";#it{x}_{MC};#it{x}_{Reco}", 100, 0, 0.2, 100, 0, 0.2);
  // inelasticity resolution
  TH1D* h_Res_inelasticity     = new TH1D("Res_inelasticity","Inelasticity;y(Reco)-y(MC)",100,-0.05,0.05);
  TH1D* h_RelRes_inelasticity  = new TH1D("RelRes_inelasticity","electron method;#frac{y(Reco)-y(MC)}{y(MC)}",100,-0.2,0.2);
  // inelasticity resolution binned
  TH2D* h_inelasticityRelRes_binned = new TH2D("yRelRes_binned",";y];#frac{y(Reco)-y(MC)}{y(MC)}",20,0.2,0.9,20,0,1.0);
  // inelasticity correlation
  TH2D* h_Corr_inelasticity = new TH2D("Corr_y", ";#it{y}_{MC};#it{y}_{Reco}", 100, 0, 1, 100, 0, 1);

  // Photon angluar resolutions
  TH1D* h_PhotRes_theta   = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [rad]",600,-1.5,1.5);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [rad]; #delta#theta_{#gamma}",320,0,3.2,600,-1.5,1.5);

  // MY HISTOGRAMS
  //events Kinematics MC as I calculated
  TH1D* h_Q2_e     = new TH1D("h_Q2_e",";Q^{2};# of events",20,0,40);
  TH1D* h_y_e      = new TH1D("h_y_e",";y_{e,MC}",100,0,1);
  TH1D* h_x_e      = new TH1D("h_x_e",";x_{e,MC}",30,0,0.20);


  TH1D* h_t_MC_old     = new TH1D("h_t_MC_old",";t_{MC}; counts", nLogBins,t_edges);
  // using the electron method
  TH1D* h_Q2_e_m   = new TH1D("h_Q2_e_m",";Q^{2}",20,0,40);
  TH1D* h_y_e_m    = new TH1D("h_y_e_m",";y_{e,MC}",20,0,1);
  TH1D* h_x_e_m    = new TH1D("h_x_e_m",";x_{e,MC}",30,0,0.20);
  // Kinematic MC Truth from branches 
  TH1D* h_y_truth  = new TH1D("h_y_truth","",20,0,1);
  TH1D* h_x_truth  = new TH1D("h_x_truth","x_{truth}",30,0,0.20);
  TH1D* h_Q2_truth = new TH1D("h_Q2_truth","Q^2;# of events",20,0,40);
  TH1D* h_py_truth = new TH1D("h_py_truth","Q^2",100,-2.5,2.5);
  TH1D* h_px_truth = new TH1D("h_px_truth","Q^2",100,-2.5,2.5);


  // Kinematic using DA
  TH1D* h_y_DA     = new TH1D("h_y_DA","y_{DA}",20,0,1);
  TH1D* h_Q2_DA    = new TH1D("h_Q2_DA",";Q^{2}",20,0,40);
  TH1D* h_x_DA     = new TH1D("h_x_DA",";x_{Bj}",30,0,0.20);



  // Kinematic Reconstructed
  TH1D* h_Q2REC_e        = new TH1D("h_Q2REC_e",";Q^{2}_{e,REC}",100,0,40);
  TH1D* h_yREC_e         = new TH1D("h_yREC_e",";y_{e,REC}",100,0,1);
  TH1D* h_tREC_B0        = new TH1D("h_tREC_B0",";t_{REC}(B0); counts",100,-1,0);

  //---------------------------------------------------------
  // DECLARE TTREEREADER AND BRANCHES TO USE
  //---------------------------------------------------------
  TTreeReader tree_reader(events);
  // MC particles
  TTreeReaderArray<double>  mc_px_array          = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<double>  mc_py_array          = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<double>  mc_pz_array          = {tree_reader, "MCParticles.momentum.z"};
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

  // hadi
  TTreeReaderArray<float> DA_Q2_array         = {tree_reader, "InclusiveKinematicsDA.Q2"};
  
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


  // // Declare variables using TTreeReaderValue
  // TTreeReaderValue<Float_t> electron_Q2(tree_reader, "InclusiveKinematicsElectron.Q2");
  // TTreeReaderValue<Float_t> electron_y(tree_reader, "InclusiveKinematicsElectron.y");
  // TTreeReaderValue<Float_t> electron_x(tree_reader, "InclusiveKinematicsElectron.x");

  // TTreeReaderValue<Float_t> electron_Q2_DA(tree_reader, "InclusiveKinematicsDA.Q2");
  // TTreeReaderValue<Float_t> electron_y_DA(tree_reader, "InclusiveKinematicsDA.y");
  // TTreeReaderValue<Float_t> electron_x_DA(tree_reader, "InclusiveKinematicsDA.x");

  // TTreeReaderValue<Float_t> electron_Q2_truth(tree_reader, "InclusiveKinematicsTruth.Q2");
  // TTreeReaderValue<Float_t> electron_y_truth(tree_reader, "InclusiveKinematicsTruth.y");
  // TTreeReaderValue<Float_t> electron_x_truth(tree_reader, "InclusiveKinematicsTruth.x");

    // Vectors to store hit coordinates
    std::vector<double> mc_px, mc_py;   // MC truth hits
    std::vector<double> reco_px, reco_py; // Reco hits
    // Vectors to store rapidity and azimuthal angle (phi) for protons
    std::vector<double> mc_rapidity, mc_phi;   // MC Truth rapidity-phi
    std::vector<double> rp_rapidity, rp_phi;   // Reconstructed proton rapidity-phi
    

  while (tree_reader.Next()){
    ///////////////////////////////////// stuff from old code ////////////////////////////
    unsigned int mc_elect_index=-1;
    double maxPt=-99.;
    /*
    Beam particles
    */
    TLorentzVector ebeam(0,0,0,0);
    TLorentzVector pbeam(0,0,0,0);
    TLorentzVector scatElectronMC(0,0,0,0);// hh: this should be scattered electron
    TLorentzVector scatProtonMC(0,0,0,0);  // hh: this is scattered proton
  
    for(unsigned int imc=0; imc < mc_px_array.GetSize(); imc++){
      // std::cout << "Filled branch entries: " << std::setw(3)  << rp_rec_px_array.GetSize() << std::endl;
  
      TVector3 mctrk(mc_px_array[imc],mc_py_array[imc],mc_pz_array[imc]);
      if(mc_genStatus_array[imc]==4){//4 is beam particle
        if(mc_pdg_array[imc]==11  ) ebeam.SetVectM(mctrk, MASS_ELECTRON);
        if(mc_pdg_array[imc]==2212) pbeam.SetVectM(mctrk, MASS_PROTON  );
      }
      if(mc_genStatus_array[imc]!=1) continue;
      if(mc_pdg_array[imc]==11&& mctrk.Perp()>maxPt){
        maxPt=mctrk.Perp();
        mc_elect_index=imc;
        scatElectronMC.SetVectM(mctrk,mc_mass_array[imc]);
      }
      if(mc_pdg_array[imc]==2212 && mc_genStatus_array[imc]==1){// hh: this should be scattered proton
        scatProtonMC.SetVectM(mctrk,mc_mass_array[imc]);
      }
    }
    TLorentzVector qbeam       = ebeam - scatElectronMC;
    TLorentzVector deltaProton = scatProtonMC - pbeam;
    double pDotq = pbeam.Dot(qbeam);
    double y     = pDotq/pbeam.Dot(ebeam);
    double Q2    = -(qbeam).Mag2();
    double xBj   = -(qbeam.Dot(qbeam))/(2.0*pDotq);
    double xPom  = -(deltaProton.Dot(qbeam))/(pDotq);
    // std::cout << "xPom = " << xPom << std::endl;
  
  
    h_t_MC_old->Fill(TMath::Abs(deltaProton.Mag2()));
    h_y_e->Fill(y);
    h_x_e->Fill(xBj);
    h_Q2_e->Fill(Q2);	  
    // h_xPom_e->Fill(xPom);
    ///////////////////////////////////// stuff from old code ////////////////////////////
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
    vector<P3MVector> scatp4_gen;   // Scattered proton (generated)
    vector<P3MVector> scatg4_gen;   // Scattered photon (generated)
    // 4-vectors for associated MC particles (ONLY SCATTERED - only need photon for resolution plot)
    vector<P3MVector> scatg4_aso;   // Scattered photon (associated MC)
    // 4-vectors for reconstructed particles (SEPARATE PROTONS FOR B0 AND ROMAN POTS)
    vector<P3MVector> scate4_rec;   // Scattered electron (reconstructed)
    vector<P3MVector> scatp4_rec;   // Scattered proton (B0 reconstructed)
    vector<P3MVector> scatp4_aso;   // Scattered proton (associated MC)
    vector<P3MVector> scate4_aso;   // Scattered electron (associated MC)
    
    vector<P3MVector> scatp4_rom;   // Scattered proton (Roman Pots reconstructed)
    vector<P3MVector> scatg4_rec;   // Scattered photon (reconstructed)      

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
	      if(mc_pdg_array[imc] == 2212){
	        scatp4_gen.push_back(q_scat);
          h_px_truth->Fill(q_scat.Px());
          h_py_truth->Fill(q_scat.Py());
          // mc_px.push_back(q_scat.Px());
          // mc_py.push_back(q_scat.Py());
          // mc_rapidity.push_back(q_scat.Rapidity());
          // mc_phi.push_back(q_scat.Phi());
        } 
	      if(mc_pdg_array[imc] == 11){
	        scate4_gen.push_back(q_scat);
	      }
	      if(mc_pdg_array[imc] == 22){
	        scatg4_gen.push_back(q_scat);
	      }
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
      if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 11){ 
	      scate4_rec.push_back(q_reco); 
        // std::cout << "q_reco_DA = " << DA_Q2_array[mc_assoc_index] << std::endl; 
        P3MVector e_reco(mc_px_array[mc_assoc_index],mc_py_array[mc_assoc_index],mc_pz_array[mc_assoc_index],mc_mass_array[mc_assoc_index]);
        scate4_aso.push_back(e_reco);
      }
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

        P3MVector q_assoc(mc_px_array[mc_assoc_index],mc_py_array[mc_assoc_index],mc_pz_array[mc_assoc_index],mc_mass_array[mc_assoc_index]);
        undoAfterburn(q_assoc);
  
	      scatp4_rec.push_back(q_reco); 
        scatp4_aso.push_back(q_assoc);
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


    // std::cout << "mc_assoc_index = " << mc_assoc_index << std::endl;

    //---------------------------------------------------------
    // Fill histograms
    //---------------------------------------------------------
    // Generated particles
    // Need Q2 for electron cuts
    // Float_t fQ2{0};
    // if(scate4_gen.size() == 0) fQ2 = 0;
    // else fQ2 = calcQ2(beame4, scate4_gen[0]);
    // Eta (inclusive tracks)
    // if(fQ2 >= 1){
      // e'
      // if(scate4_gen.size() == 1) h_eta_MCe->Fill(scate4_gen[0].Eta());
      // gamma
      // if(scatg4_gen.size() == 1) h_eta_MCg->Fill(scatg4_gen[0].Eta());
      // p'
      if(scatp4_gen.size() == 1 && scatp4_gen[0].Theta()<0.02){
	      h_eta_MCp->Fill(scatp4_gen[0].Eta());
        mc_rapidity.push_back(scatp4_gen[0].Eta());
        // std::cout << "pseudo-rapidity check = " << TMath::Nint((scatp4_gen[0].Eta() - calculatePseudorapidity(scatp4_gen[0]))/1e-9)*1e-9 << std::endl;
        // std::cout << "rapidity = " << scatp4_gen[0].Rapidity() << "  " << rapidity(scatp4_gen[0]) << std::endl;
        mc_phi.push_back(scatp4_gen[0].Phi());
        mc_px.push_back(scatp4_gen[0].Px());
        mc_py.push_back(scatp4_gen[0].Py());

	      // Add exclusivity cuts for t-distribution
	      // if(scate4_gen.size() == 1 && scatg4_gen.size() == 1){
	      // Need to calculate kinematics before cutting on them
	      // Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]);
	      Float_t ft = calcT(beamp4, scatp4_gen[0]);
	    
	      // Want MM2 to be close to zero
	      // if(TMath::Abs(fM2miss) < 1) 
        h_t_MC->Fill(ft);
	      // } // Exclusivity required
      } // Proton tracks done
    // } // Q2 limit
    
    // Reconstructed particles
    // Need Q2 for electron cuts
    // if(scate4_rec.size() == 0) fQ2 = 0;
    // else fQ2 = calcQ2(beame4, scate4_rec[0]);
    // Eta (inclusive tracks)
    // if(fQ2 >= 1){
      // e'
      // if(scate4_rec.size() == 1) h_eta_RPe->Fill(scate4_rec[0].Eta());
      // gamma
      // if(scatg4_rec.size() == 1) h_eta_RPg->Fill(scatg4_rec[0].Eta());
      // p' (B0 - theta between 5.5 and 20 mrad)
      if(scatp4_rec.size() == 1 && scatp4_rec[0].Theta()>0.0055 && scatp4_rec[0].Theta()<0.02){
	      h_eta_RPp->Fill(scatp4_rec[0].Eta());
        if(scatp4_aso.size()==1 && scatp4_rec.size()==1){

          double t_rec = calcT(scatp4_rec[0],beamp4);
          double t_aso = calcT(scatp4_aso[0],beamp4);
          h_PrRes_t->Fill(t_rec - t_aso);
          h_PrRelRes_t->Fill((t_rec - t_aso)/t_aso);
          h_tRelRes_binned->Fill(t_aso, (t_rec - t_aso)/t_aso);
          h_Corr_t->Fill(t_aso, t_rec);

          // std::cout << "B0 rec   = " << scatp4_rec[0].Pz() << std::endl;
          // std::cout << "MC truth = " << scatp4_aso[0].Pz() << std::endl;
        }
        // reco_px.push_back(scatp4_rec[0].Px());
        // reco_px.push_back(scatp4_rec[0].Py());
	      // Add exclusivity cuts for t-distribution
	      // if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
	      // Need to calculate kinematics before cutting on them
	      // Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]);
	      Float_t ft = calcT(beamp4, scatp4_rec[0]);
	    
	      // Want MM2 to be close to zero
	      // if(TMath::Abs(fM2miss) < 1) 
        h_t_RP->Fill(ft);
	      // } // Exclusivity required
      } // B0 Proton tracks found
      // p' (RP - theta less than 5 mrad)

      if(scatp4_rom.size() == 1 && scatp4_rom[0].Theta()<0.005){
	      h_eta_RPPp->Fill(scatp4_rom[0].Eta());
        reco_px.push_back(scatp4_rom[0].Px());
        reco_py.push_back(scatp4_rom[0].Py());
        rp_rapidity.push_back(scatp4_rom[0].Eta());
        rp_phi.push_back(scatp4_rom[0].Phi());

        if(scatp4_aso.size()==1 && scatp4_rom.size()==1){

          double t_rec = calcT(scatp4_rom[0],beamp4);
          double t_aso = calcT(scatp4_gen[0],beamp4);
          h_PrRes_t->Fill(t_rec - t_aso);
          // std::cout << scatp4_aso[0] << std::endl;
          // std::cout << "beamp4 =        " << beamp4.px() << " " << beamp4.py() << " " << beamp4.pz() << " " << beamp4.E() << std::endl;
          // std::cout << "scatp4_aso[0] = " << scatp4_gen[0].px() << " " << scatp4_gen[0].py() << " " << scatp4_gen[0].pz() << " " << scatp4_gen[0].E() << std::endl;
          std::cout << "t_rec = " << t_rec << " t_aso = " << t_aso << std::endl;
          h_PrRelRes_t->Fill((t_rec - t_aso)/t_aso);
          h_tRelRes_binned->Fill(t_aso, (t_rec - t_aso)/t_aso);
          h_Corr_t->Fill(t_aso, t_rec);

          // std::cout << "B0 rec   = " << scatp4_rec[0].Pz() << std::endl;
          // std::cout << "MC truth = " << scatp4_aso[0].Pz() << std::endl;
        }

        // if (mc_assoc_index< 4294967295){// this strange no sometimes shows up and the code crashes
        // std::cout << "mc_px_array        = " << mc_pz_array[mc_assoc_index] << std::endl;
        // std::cout << "scatp4_rom[0].Px() = " << scatp4_rom[0].Pz()          << std::endl;
        // }
	      // Add exclusivity cuts for t-distribution
	      // if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
	      // Need to calculate kinematics before cutting on them
	      // Float_t fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]);
	      Float_t ft = calcT(beamp4, scatp4_rom[0]);
	    
	      // Want MM2 to be close to zero
	      // if(TMath::Abs(fM2miss) < 1 && ft < 0.3) 
        h_t_RPP->Fill(ft);
	      // } // Exclusivity required
      } // RP Proton tracks found
    // } // Q2 limit
    
    if(scatg4_aso.size()==1 && scatg4_rec.size()==1){
      h_PhotRes_theta->Fill(scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
      h_PhotRes2D_theta->Fill(scatg4_aso[0].Theta(), scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
    }

    if(scate4_aso.size()==1 && scate4_rec.size()==1){
      // h_PhotRes_theta->Fill(scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
      // h_PhotRes2D_theta->Fill(scatg4_aso[0].Theta(), scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
      double q2_rec = calcQ2(beame4, scate4_rec[0]);
      double q2_aso = calcQ2(beame4, scate4_aso[0]); 
      double pdotq_rec = beamp4.Dot(scate4_rec[0]);
      double pdotq_aso = beamp4.Dot(scate4_aso[0]);
      double s_rec = 4*beame4.E()*beamp4.E(); //(beame4 + beamp4).Dot(beame4 + beamp4);
      double s_aso;
      if (scatp4_aso.size()==1 && scatp4_rec.size()==1){
        s_aso = 4*beame4.E()*scatp4_aso[0].E(); //(beame4 + scate4_aso[0]).Dot(beame4 + scate4_aso[0]);
      } else {s_aso= s_rec;}
      
      /////////////////////////////////////////////////
      h_Res_Q2->Fill(q2_rec - q2_aso);
      h_RelRes_Q2->Fill((q2_rec - q2_aso)/q2_aso);
      // h_Corr_Q2->Fill(q2_aso, q2_rec);
      h_RelRes_Q2_binned->Fill(q2_aso, (q2_rec - q2_aso)/q2_aso);
      ////////////////////////////////////////////////
      double xBj_rec = q2_rec/(2*pdotq_rec);
      double xBj_aso = q2_aso/(2*pdotq_aso);
      h_Res_xBj->Fill(xBj_rec - xBj_aso);
      h_RelRes_xBj->Fill((xBj_rec - xBj_aso)/xBj_aso);
      // h_Corr_xBj->Fill(xBj_rec, xBj_aso);
      h_xBjRelRes_binned->Fill(xBj_aso, (xBj_rec - xBj_aso)/xBj_aso);
      // h_tRelRes_xBj->Fill((xBj_rec - xBj_aso)/xBj_aso, ft);
      ////////////////////////////////////////////////
      double y_rec = q2_rec/(xBj_rec*s_rec);
      double y_aso = q2_aso/(xBj_aso*s_aso);
      h_Res_inelasticity->Fill(y_rec - y_aso);
      h_RelRes_inelasticity->Fill((y_rec - y_aso)/y_aso);
      // h_Corr_inelasticity->Fill(y_aso, y_rec);
      h_inelasticityRelRes_binned->Fill(y_aso, (y_rec - y_aso)/y_aso);
      // std::cout << "y = " << y_rec << " " << y_aso << std::endl;
      // std::cout << "Q2 = " << calcQ2(beame4, scate4_rec[0]) << " " << calcQ2(beame4, scate4_aso[0]) << std::endl;
    }


  }



         // Convert vectors to arrays for TGraph
         int mc_size = mc_rapidity.size();
         int rp_size = rp_rapidity.size();
         double* mc_rapidity_arr = mc_rapidity.data();
         double* mc_phi_arr = mc_phi.data();
         double* rp_rapidity_arr = rp_rapidity.data();
         double* rp_phi_arr = rp_phi.data();
     
         // Create scatter plots
         TGraph* g_mc_hits = new TGraph(mc_size, mc_rapidity_arr, mc_phi_arr);
         TGraph* g_rp_hits = new TGraph(rp_size, rp_rapidity_arr, rp_phi_arr);
     

         g_mc_hits->GetXaxis()->SetLimits(0, 10); // Set rapidity range for MC hits
         g_rp_hits->GetXaxis()->SetLimits(0, 10); // Set rapidity range for Reco hits

         // Style MC truth proton hits (Blue with transparency)
         g_mc_hits->SetMarkerStyle(20);
         g_mc_hits->SetMarkerSize(0.8);
         g_mc_hits->SetMarkerColorAlpha(kBlue, 0.3);  // 30% opacity
         g_mc_hits->GetXaxis()->SetRangeUser(0,10);

     
         // Style Reco proton hits (Red with transparency)
         g_rp_hits->SetMarkerStyle(21);
         g_rp_hits->SetMarkerSize(0.8);
         g_rp_hits->SetMarkerColorAlpha(kRed, 0.3);  // 30% opacity
         g_rp_hits->GetXaxis()->SetRangeUser(0, 10);

     
         // Create a canvas
         TCanvas* c_2 = new TCanvas("c2", "Proton MC vs Reco Hits (#eta-#phi)", 1600, 1200);
         c_2->SetGrid();
     
         // Draw scatter plots
         g_mc_hits->Draw("AP");
         g_mc_hits->SetTitle("Proton MC vs Reco Hits (#eta-#phi); PseudoRapidity ; Azimuthal Angle (#phi)");
         g_rp_hits->Draw("P SAME");
     
         // Add legend
         TLegend* legend_2 = new TLegend(0.7, 0.8, 0.9, 0.9);
         legend_2->AddEntry(g_mc_hits, "MC Truth Protons (Blue, Transparent)", "p");
         legend_2->AddEntry(g_rp_hits, "Reconstructed Protons (Red, Transparent)", "p");
         legend_2->Draw();
     
         // Save the canvas
         c_2->SaveAs("figs/Proton_Eta_Phi.png");



      // Convert vectors to arrays for TGraph
      mc_size = mc_px.size();
      int     reco_size   = reco_px.size();
      double* mc_px_arr   = mc_px.data();
      double* mc_py_arr   = mc_py.data();
      double* reco_px_arr = reco_px.data();
      double* reco_py_arr = reco_py.data();
  
      // Create scatter plots
      TGraph* g_mc_hits_xvy = new TGraph(mc_size, mc_px_arr, mc_py_arr);
      TGraph* g_reco_hits   = new TGraph(reco_size, reco_px_arr, reco_py_arr);
  
      // Style MC truth proton hits (Blue with transparency)
      g_mc_hits_xvy->SetMarkerStyle(20);
      g_mc_hits_xvy->SetMarkerSize(0.8);
      g_mc_hits_xvy->SetMarkerColorAlpha(kBlue, 0.3);  // 30% opacity
  
      // Style Reco proton hits (Red with transparency)
      g_reco_hits->SetMarkerStyle(21);
      g_reco_hits->SetMarkerSize(0.8);
      g_reco_hits->SetMarkerColorAlpha(kRed, 0.3);  // 30% opacity
  
      // Create a canvas
      TCanvas* c_h = new TCanvas("c1", "Proton MC vs Reco Hit Scatter Plot", 1600, 1200);
      c_h->SetGrid();
      g_mc_hits_xvy->GetXaxis()->SetRangeUser(-2, 2);
      g_mc_hits_xvy->GetYaxis()->SetRangeUser(-2, 2);
  
      std::cout << "g_mc_hits points: " << g_mc_hits_xvy->GetN() << std::endl;
      std::cout << "g_reco_hits points: " << g_reco_hits->GetN() << std::endl;
      // Draw scatter plots
      g_mc_hits_xvy->Draw("AP");
      g_mc_hits_xvy->SetTitle("Proton MC vs Reco Hits Scatter Plot; p_{x} [GeV/c]; p_{y} [GeV/c]");
      g_reco_hits->Draw("P SAME");
  
      // Add legend
      TLegend* legend1 = new TLegend(0.7, 0.8, 0.9, 0.9);
      legend1->AddEntry(g_mc_hits_xvy, "MC Truth Protons (Blue, Transparent)", "p");
      legend1->AddEntry(g_reco_hits, "Reco Protons (Red, Transparent)", "p");
      legend1->Draw();
  
      // Save the canvas
      c_h->SaveAs("figs/Proton_MC_vs_Reco_Hits.png");


  

  //---------------------------------------------------------
  // DRAW HISTOGRAMS
  //---------------------------------------------------------
  gStyle->SetOptStat(00000000);
  // Set header strings
  TString sHead1 = "#it{#bf{ePIC} Simulation}";
  TString sHead2 = "#it{ep 10x100 GeV}";

  // CANVAS 1: ETA DISTRIBUTIONS (Proton Only)
TCanvas* c1_proton = new TCanvas("c1_proton", "Proton Eta Distribution", 800, 600); // Adjusted size

TLatex* tHeadLc1_proton = new TLatex(0.04, 0.91, sHead1); // Assuming sHead1 is defined elsewhere
tHeadLc1_proton->SetNDC();
tHeadLc1_proton->SetTextSize(25);
tHeadLc1_proton->SetTextFont(43);
tHeadLc1_proton->SetTextColor(kBlack);
tHeadLc1_proton->Draw();

TLatex* tHeadM1c1_proton = new TLatex(0.29, 0.91, "#bf{EpIC} ep #rightarrow e'p'#gamma");
tHeadM1c1_proton->SetNDC();
tHeadM1c1_proton->SetTextSize(25);
tHeadM1c1_proton->SetTextFont(43);
tHeadM1c1_proton->SetTextColor(kBlack);
tHeadM1c1_proton->Draw("same");

TLatex* tHeadM2c1_proton = new TLatex(0.59, 0.91, "Q^{2} #geq 1 GeV^{2}");
tHeadM2c1_proton->SetNDC();
tHeadM2c1_proton->SetTextSize(25);
tHeadM2c1_proton->SetTextFont(43);
tHeadM2c1_proton->SetTextColor(kBlack);
tHeadM2c1_proton->Draw("same");

TLatex* tHeadRc1_proton = new TLatex(0.83, 0.91, sHead2); // Assuming sHead2 is defined elsewhere
tHeadRc1_proton->SetNDC();
tHeadRc1_proton->SetTextSize(25);
tHeadRc1_proton->SetTextFont(43);
tHeadRc1_proton->SetTextColor(kBlack);
tHeadRc1_proton->Draw("same");

c1_proton->cd(); // No need for Divide, just use the whole canvas
gPad->SetLogy();

h_eta_MCp->SetMaximum(100 * h_eta_MCe->GetMaximum()); // Match maxima across species (if h_eta_MCe is still relevant)
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
h_eta_RPPp->SetLineColor(kCyan + 1);
h_eta_RPPp->SetMarkerColor(kCyan + 1);
h_eta_RPPp->SetMarkerStyle(20);
// Set the y-axis range here:
h_eta_MCp->SetMinimum(0.1); // Set the minimum y-value (adjust as needed)
h_eta_MCp->SetMaximum(1000); // Set the maximum y-value (adjust as needed)
// Set the x-axis range here:
h_eta_MCp->GetXaxis()->SetRangeUser(4, 10); // Set the x-axis range from -5 to 5 (adjust as needed)


h_eta_MCp->Draw();
h_eta_RPp->Draw("pesame");
h_eta_RPPp->Draw("pesame");

// Add text
TLatex* tHeadc1p3_proton = new TLatex(0.07, 0.92, "p");
tHeadc1p3_proton->SetNDC();
tHeadc1p3_proton->SetTextSize(24);
tHeadc1p3_proton->SetTextFont(43);
tHeadc1p3_proton->SetTextColor(kBlack);
tHeadc1p3_proton->Draw("same");

// Add legend
TLegend* lC1p3_proton = new TLegend(0.7, 0.7, 0.9, 0.9);
lC1p3_proton->SetLineColorAlpha(kWhite, 0);
lC1p3_proton->SetLineWidth(0);
lC1p3_proton->SetFillStyle(0);
lC1p3_proton->SetFillColorAlpha(kWhite, 0);
lC1p3_proton->AddEntry(h_eta_RPp, "Reco. B0", "pl");
lC1p3_proton->AddEntry(h_eta_RPPp, "Reco. RP", "pl");
lC1p3_proton->Draw();

// Optionally, save the canvas:
c1_proton->SaveAs("figs/proton_eta_distribution.png"); // or .pdf, .root, etc.


  // CANVAS 2: T-DISTRIBUTION
  TCanvas* c2 = new TCanvas("c2","",1200,800);
  h_t_MC->SetMinimum(0.1);
  gPad->SetLogy(1);
  h_t_MC->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  h_t_MC->GetYaxis()->SetTitle("Counts / 0.02 GeV^{2}");
  h_t_MC->SetLineColor(kBlack);
  h_t_MC->SetLineWidth(2);
  h_t_MC->Draw();
  h_t_RP->SetLineColor(kBlue);
  h_t_RP->SetMarkerColor(kBlue);
  h_t_RP->SetMarkerStyle(20);
  h_t_RP->Draw("pesame");
  h_t_RPP->SetLineColor(kCyan+1);
  h_t_RPP->SetMarkerColor(kCyan+1);
  h_t_RPP->SetMarkerStyle(20);
  h_t_RPP->Draw("pesame");
  // Create header text objects
  TLatex* tHead1 = new TLatex(0.10, 0.91, sHead1);
  tHead1->SetNDC();
  tHead1->SetTextSize(30);
  tHead1->SetTextFont(43);
  tHead1->SetTextColor(kBlack);
  TLatex* tHead2 = new TLatex(0.74, 0.91, sHead2);
  tHead2->SetNDC();
  tHead2->SetTextSize(30);
  tHead2->SetTextFont(43);
  tHead2->SetTextColor(kBlack);
  tHead1->Draw("same");
  tHead2->Draw("same");
  // Add legend
  TLatex* tC2 = new TLatex(0.58, 0.83, "#splitline{#bf{EpIC} ep #rightarrow e'p'#gamma, Q^{2} #geq 1 GeV^{2}}{t_{RP} #leq 0.3 GeV^{2}}");
  tC2->SetNDC();
  tC2->SetTextSize(30);
  tC2->SetTextFont(43);
  tC2->SetTextColor(kBlack);
  // tC2->Draw("same");
  TLegend* lC2 = new TLegend(0.57, 0.6, 0.8, 0.77);
  lC2->SetLineColorAlpha(kWhite,0);
  lC2->SetFillColorAlpha(kWhite,0);
  lC2->AddEntry(h_t_MC, "#bf{RAPGAP} MC gen.", "l");
  lC2->AddEntry(h_t_RP, "Reco. B0", "lp");
  lC2->AddEntry(h_t_RPP, "Reco. RP", "lp");
  lC2->Draw();

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




  // t resulotion plot
    // CANVAS 3A: PHOTON ANGLE RESOLUTION (INSET)
    TCanvas* c4 = new TCanvas("c4","",1200,800);
    // Draw 1D plot as main figure
    h_PrRes_t->GetXaxis()->SetTitleSize(0.05);
    h_PrRes_t->GetXaxis()->SetTitleOffset(0.9);
    h_PrRes_t->GetYaxis()->SetTitle("Counts / 5mrad");
    h_PrRes_t->GetYaxis()->SetTitleSize(0.05);
    h_PrRes_t->GetYaxis()->SetTitleOffset(0.9);
    h_PrRes_t->SetLineColor(kBlue);
    h_PrRes_t->SetMarkerColor(kBlue);
    h_PrRes_t->SetMarkerStyle(20);
    h_PrRes_t->Draw("pe");
        // Fit a Gaussian to the histogram
        TF1* gaussianFit4 = new TF1("gaussianFit", "gaus", h_PrRes_t->GetXaxis()->GetXmin(), h_PrRes_t->GetXaxis()->GetXmax());
        // Set initial parameter values
        gaussianFit4->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_PrRes_t->Fit("gaussianFit", "R");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit4->SetLineColor(kRed);
        gaussianFit4->SetLineWidth(2);
        gaussianFit4->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF
    


    // Q2 resulotion plot
    TCanvas* c5 = new TCanvas("c5","",1200,800);
    // Draw 1D plot as main figure
    h_Res_Q2->GetXaxis()->SetTitleSize(0.05);
    h_Res_Q2->GetXaxis()->SetTitleOffset(0.9);
    h_Res_Q2->GetYaxis()->SetTitle("Counts");
    h_Res_Q2->GetYaxis()->SetTitleSize(0.05);
    h_Res_Q2->GetYaxis()->SetTitleOffset(0.9);
    h_Res_Q2->SetLineColor(kBlue);
    h_Res_Q2->SetMarkerColor(kBlue);
    h_Res_Q2->SetMarkerStyle(20);
    h_Res_Q2->Draw("pe");
        // Fit a Gaussian to the histogram
        double xmin = h_Res_Q2->GetXaxis()->GetXmin();
        double xmax = h_Res_Q2->GetXaxis()->GetXmax();
        TF1* gaussianFit5 = new TF1("gaussianFit", "gaus",xmin, xmax);
        // Set initial parameter values
        gaussianFit5->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_Res_Q2->Fit("gaussianFit", "RQ", "", -0.25, 0.25);
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit5->SetLineColor(kRed);
        gaussianFit5->SetLineWidth(2);
        gaussianFit5->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF

    // xBj resulotion plot
    TCanvas* c6 = new TCanvas("c6","",1200,800);
    // Draw 1D plot as main figure
    // gPad->SetLogy();
    h_Res_xBj->GetXaxis()->SetTitleSize(0.05);
    h_Res_xBj->GetXaxis()->SetTitleOffset(0.9);
    h_Res_xBj->GetYaxis()->SetTitle("Counts");
    h_Res_xBj->GetYaxis()->SetTitleSize(0.05);
    h_Res_xBj->GetYaxis()->SetTitleOffset(0.9);
    h_Res_xBj->SetLineColor(kBlue);
    h_Res_xBj->SetMarkerColor(kBlue);
    h_Res_xBj->SetMarkerStyle(20);
    h_Res_xBj->Draw("pe");
    // Fit a Gaussian to the histogram
    xmin =-2e-5;
    xmax = 1.5e-5;
    TF1* gaussianFit6 = new TF1("gaussianFit", "gaus",xmin, xmax);
    // Set initial parameter values
    gaussianFit6->SetParameters(1, 0, 1e-5); // Amplitude, Mean, Sigma
    // Perform the fit
    h_Res_xBj->Fit("gaussianFit", "RQ", "", xmin, xmax);
    // Get the fit parameters
    // double mean  = gaussianFit->GetParameter(1);
    // double sigma = gaussianFit->GetParameter(2);
    // double meanErr  = gaussianFit->GetParError(1);
    // Draw the fitted function on the same canvas
    gaussianFit6->SetLineColor(kRed);
    gaussianFit6->SetLineWidth(2);
    gaussianFit6->Draw("same");
    gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF
    // gPad->Update();     // force update to render the stats box


    // inelasticity resulotion plot
    TCanvas* c7 = new TCanvas("c7","",1200,800);
    // Draw 1D plot as main figure
    // gPad->SetLogy();
    h_Res_inelasticity->GetXaxis()->SetTitleSize(0.05);
    h_Res_inelasticity->GetXaxis()->SetTitleOffset(0.9);
    h_Res_inelasticity->GetYaxis()->SetTitle("Counts");
    h_Res_inelasticity->GetYaxis()->SetTitleSize(0.05);
    h_Res_inelasticity->GetYaxis()->SetTitleOffset(0.9);
    h_Res_inelasticity->SetLineColor(kBlue);
    h_Res_inelasticity->SetMarkerColor(kBlue);
    h_Res_inelasticity->SetMarkerStyle(20);
    h_Res_inelasticity->Draw("pe");
        // Fit a Gaussian to the histogram
        TF1* gaussianFit7 = new TF1("gaussianFit", "gaus", h_Res_inelasticity->GetXaxis()->GetXmin(), h_Res_inelasticity->GetXaxis()->GetXmax());
        // Set initial parameter values
        gaussianFit7->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_Res_inelasticity->Fit("gaussianFit", "R");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit7->SetLineColor(kRed);
        gaussianFit7->SetLineWidth(2);
        gaussianFit7->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF


    // t relative resulotion plot
    TCanvas* c8 = new TCanvas("c8","",1200,800);
    // Draw 1D plot as main figure
    // gPad->SetLogy();
    c8->SetBottomMargin(0.2); // Increase bottom margin
    h_PrRelRes_t->GetXaxis()->SetTitleSize(0.05);
    h_PrRelRes_t->GetXaxis()->SetTitleOffset(1.3);
    h_PrRelRes_t->GetYaxis()->SetTitle("Counts");
    h_PrRelRes_t->GetYaxis()->SetTitleSize(0.05);
    h_PrRelRes_t->GetYaxis()->SetTitleOffset(0.9);
    h_PrRelRes_t->SetLineColor(kGreen+3);
    h_PrRelRes_t->SetMarkerColor(kGreen+3);
    h_PrRelRes_t->SetMarkerStyle(20);
    h_PrRelRes_t->Draw("pe");
        // Fit a Gaussian to the histogram
        xmin = -0.1;
        xmax = 0.1;
        TF1* gaussianFit8 = new TF1("gaussianFit", "gaus", xmin, xmax);
        // Set initial parameter values
        gaussianFit8->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_PrRelRes_t->Fit("gaussianFit", "RQ");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit8->SetLineColor(kRed);
        gaussianFit8->SetLineWidth(2);
        gaussianFit8->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF

    // Q2 relative resulotion plot
    TCanvas* c9 = new TCanvas("c9","",1200,800);
    // Draw 1D plot as main figure
    c9->SetBottomMargin(0.2); // Increase bottom margin
    h_RelRes_Q2->GetXaxis()->SetTitleSize(0.05);
    h_RelRes_Q2->GetXaxis()->SetTitleOffset(1.3);
    h_RelRes_Q2->GetYaxis()->SetTitle("Counts");
    h_RelRes_Q2->GetYaxis()->SetTitleSize(0.05);
    h_RelRes_Q2->GetYaxis()->SetTitleOffset(0.9);
    h_RelRes_Q2->SetLineColor(kGreen+3);
    h_RelRes_Q2->SetMarkerColor(kGreen+3);
    h_RelRes_Q2->SetMarkerStyle(20);
    h_RelRes_Q2->Draw("pe");
        // Fit a Gaussian to the histogram
        xmin = -0.025;
        xmax = 0.02;
        TF1* gaussianFit9 = new TF1("gaussianFit", "gaus",xmin, xmax);
        // Set initial parameter values
        gaussianFit9->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_RelRes_Q2->Fit("gaussianFit", "RQ");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit9->SetLineColor(kRed);
        gaussianFit9->SetLineWidth(2);
        gaussianFit9->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF


    // xBj relative resulotion plot
    TCanvas* c10 = new TCanvas("c10","",1200,800);
    // Draw 1D plot as main figure
    c10->SetBottomMargin(0.2); // Increase bottom margin
    h_RelRes_xBj->GetXaxis()->SetTitleSize(0.05);
    h_RelRes_xBj->GetXaxis()->SetTitleOffset(1.3);
    h_RelRes_xBj->GetYaxis()->SetTitle("Counts");
    h_RelRes_xBj->GetYaxis()->SetTitleSize(0.05);
    h_RelRes_xBj->GetYaxis()->SetTitleOffset(0.9);
    h_RelRes_xBj->SetLineColor(kGreen+3);
    h_RelRes_xBj->SetMarkerColor(kGreen+3);
    h_RelRes_xBj->SetMarkerStyle(20);
    h_RelRes_xBj->Draw("pe");
        // Fit a Gaussian to the histogram
        xmin = -0.009;// -0.015;
        xmax = 0.006;//0.01;
        TF1* gaussianFit10 = new TF1("gaussianFit", "gaus",xmin, xmax);
        // Set initial parameter values
        gaussianFit10->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_RelRes_xBj->Fit("gaussianFit", "R");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit10->SetLineColor(kRed);
        gaussianFit10->SetLineWidth(2);
        gaussianFit10->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF


    // xBj relative resulotion plot
    TCanvas* c11 = new TCanvas("c11","",1200,800);
    // Draw 1D plot as main figure
    c11->SetBottomMargin(0.2); // Increase bottom margin
    h_RelRes_inelasticity->GetXaxis()->SetTitleSize(0.05);
    h_RelRes_inelasticity->GetXaxis()->SetTitleOffset(1.3);
    h_RelRes_inelasticity->GetYaxis()->SetTitle("Counts");
    h_RelRes_inelasticity->GetYaxis()->SetTitleSize(0.05);
    h_RelRes_inelasticity->GetYaxis()->SetTitleOffset(0.9);
    h_RelRes_inelasticity->SetLineColor(kGreen+3);
    h_RelRes_inelasticity->SetMarkerColor(kGreen+3);
    h_RelRes_inelasticity->SetMarkerStyle(20);
    h_RelRes_inelasticity->Draw("pe");
        // Fit a Gaussian to the histogram
        xmin = -0.03;
        xmax = 0.03;
        TF1* gaussianFit11 = new TF1("gaussianFit", "gaus",xmin, xmax);
        // Set initial parameter values
        gaussianFit11->SetParameters(1, 0, 0.1); // Amplitude, Mean, Sigma
        // Perform the fit
        h_RelRes_inelasticity->Fit("gaussianFit", "RQ");
        // Get the fit parameters
        // double mean  = gaussianFit->GetParameter(1);
        // double sigma = gaussianFit->GetParameter(2);
        // double meanErr  = gaussianFit->GetParError(1);
        // Draw the fitted function on the same canvas
        gaussianFit11->SetLineColor(kRed);
        gaussianFit11->SetLineWidth(2);
        gaussianFit11->Draw("same");
        gStyle->SetOptFit(111); // enables display of fit pars, errors on the pars, 3rd digit (hundreds place) Display chi2/NDF


  

    

    // Create a canvas and draw the histogram
    TCanvas* c13 = new TCanvas("c13","",1200,800);
    //////////////////////////////////////////////
    int nbinsX = h_tRelRes_binned->GetNbinsX();
    TGraphErrors *g = new TGraphErrors(nbinsX);

    for (int j = 1; j <= nbinsX; ++j) {
      // Project X (t_reco - t_mc) for Q² bin j
      TH1D* projY = h_tRelRes_binned->ProjectionY("_py", j, j);
    
      double mean = projY->GetMean();
      double rms  = projY->GetRMS();


      double t_center = h_tRelRes_binned->GetXaxis()->GetBinCenter(j);  // Q² bin center
      double t_width  = h_tRelRes_binned->GetXaxis()->GetBinWidth(j);   // Q² bin width

      std::cout << "t_center: " << t_center << " events: " << projY->GetEntries() << std::endl;
 
      g->SetPoint(j-1, t_center, mean);                   // X=t center, Y=mean of t_diff
      g->SetPointError(j-1,0.0* t_width/2.0, rms);            // X error = half bin width, Y error = RMS

      delete projY;
    }

    // Style and draw
    g->SetTitle("electron method");
    g->GetXaxis()->SetTitle("t_{MC} [GeV^{2}]");
    g->GetYaxis()->SetTitle("#frac{t_{reco} - t_{MC}}{t_{MC}}");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+2);
    g->SetLineColor(kBlue+2);
    g->SetLineWidth(2);
    g->Draw("AP");
    double x1 = g->GetXaxis()->GetXmin();
    double x2 = g->GetXaxis()->GetXmax();
    double y1 = 0;
    double y2 = 0;
    TLine *line = new TLine(x1, y1,x2, y2);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);  // Dashed line
    line->SetLineWidth(2);
    line->Draw();
    ////////////////////////////////////////////////




    // Create a canvas and draw the histogram
    TCanvas *c16 = new TCanvas("c16", "Gaussian 2D", 800, 600);
    // gPad->SetLogy(1);
    // h_Corr_xBj->GetXaxis()->SetLimits(0, 0.2);
    // Create a TLatex object to draw LaTeX
    TLatex *latex = new TLatex();
    h_Corr_t->Draw("COLZ"); // Draw with color palette
    latex->SetNDC();        // Use normalized coordinates (0,1)
    latex->SetTextSize(0.05); // Set text size relative to canvas
    latex->SetTextColor(kBlack);
    latex->SetTextFont(42);  // Set the font (11-61, see TFont)
    latex->DrawLatex(0.15, 0.8, "BABE");
    



    // Create a canvas and draw the histogram
    TCanvas* c17 = new TCanvas("c17","",1200,800);
    c17->SetLeftMargin(0.15); // Increase left margin
    //////////////////////////////////////////////
    nbinsX = h_RelRes_Q2_binned->GetNbinsX();
      TGraphErrors *gQ2 = new TGraphErrors(nbinsX);
  
      for (int j = 1; j <= nbinsX; ++j) {
        // Project X (t_reco - t_mc) for Q² bin j
        TH1D* projY = h_RelRes_Q2_binned->ProjectionY("_pyQ2", j, j);
      
        TF1 *gaus = new TF1("gaus", "gaus");
        projY->Fit(gaus, "Q");  // "Q" = quiet mode, no fit output
        
        double mean  = gaus->GetParameter(1);      // Mean of the fit
        double sigma = gaus->GetParameter(2);      // Sigma of the fit
      
        double t_center = h_RelRes_Q2_binned->GetXaxis()->GetBinCenter(j);  // Q² bin center
        double t_width  = h_RelRes_Q2_binned->GetXaxis()->GetBinWidth(j);   // Q² bin width

        if (t_center < 5) {
          delete projY;
          continue;  // Skip points below x=5
        }
     
        gQ2->SetPoint(j-1, t_center, mean);                   // X=t center, Y=mean of t_diff
        gQ2->SetPointError(j-1, 0.0*t_width/2.0, sigma);            // X error = half bin width, Y error = RMS
    
        delete projY;
        delete gaus;

      }
    
      // Style and draw
      gQ2->SetTitle(";Q^{2}_{MC};#frac{Q^{2}_{reco} - Q^{2}_{MC}}{Q^{2}_{MC}}");
      gQ2->SetMarkerStyle(20);
      gQ2->SetMarkerColor(kBlue+2);
      gQ2->SetLineColor(kBlue+2);
      gQ2->SetLineWidth(2);
      gQ2->GetXaxis()->SetLimits(0, 40);
      gQ2->Draw("AP");
      ////
      x1 = 0;
      x2 = gQ2->GetXaxis()->GetXmax();
      y1 = 0;
      y2 = 0;
      TLine *line_Q2 = new TLine(x1, y1,x2, y2);
      line_Q2->SetLineColor(kRed);
      line_Q2->SetLineStyle(2);  // Dashed line
      line_Q2->SetLineWidth(2);
      line_Q2->Draw();
      ////
      double y_min = gQ2->GetYaxis()->GetXmin();  // Get Y range from the graph
      double y_max = gQ2->GetYaxis()->GetXmax();
      TPave *shade = new TPave(0, y_min, 5, y_max, 0, "br");
      shade->SetFillStyle(3004);  // Hatch pattern (e.g., 3001–3010)
      shade->SetFillColor(kGray+1);
      shade->SetLineColor(kGray+2);  // Optional border
      shade->Draw();
      ////////////////////////////////////////////////


    // Create a canvas and draw the histogram
    TCanvas* c18 = new TCanvas("c18","",1200,800);
    c18->SetLeftMargin(0.15); // Increase left margin
    //////////////////////////////////////////////
    nbinsX = h_xBjRelRes_binned->GetNbinsX();
      TGraphErrors *gxb = new TGraphErrors(nbinsX);
  
      for (int j = 1; j <= nbinsX; ++j) {
        // Project X (t_reco - t_mc) for Q² bin j
        TH1D* projY = h_xBjRelRes_binned->ProjectionY("_pyxBj", j, j);
      
        TF1 *gaus = new TF1("gaus", "gaus");
        projY->Fit(gaus, "Q");  // "Q" = quiet mode, no fit output
        
        double mean  = gaus->GetParameter(1);      // Mean of the fit
        double sigma = gaus->GetParameter(2);      // Sigma of the fit
      
        double t_center = h_xBjRelRes_binned->GetXaxis()->GetBinCenter(j);  // Q² bin center
        double t_width  = h_xBjRelRes_binned->GetXaxis()->GetBinWidth(j);   // Q² bin width

     
        gxb->SetPoint(j-1, t_center, mean);                   // X=t center, Y=mean of t_diff
        gxb->SetPointError(j-1,0.0* t_width/2.0, sigma);            // X error = half bin width, Y error = RMS
    
        delete projY;
      }
    
      // Style and draw
      gxb->SetTitle(";x^{MC}_{Bj};#frac{x^{reco}_{Bj} - x^{MC}_{Bj}}{x^{MC}_{Bj}}");
      gxb->SetMarkerStyle(20);
      gxb->SetMarkerColor(kBlue+2);
      gxb->SetLineColor(kBlue+2);
      gxb->SetLineWidth(2);
      gxb->Draw("AP");
      x1 = gxb->GetXaxis()->GetXmin();
      x2 = gxb->GetXaxis()->GetXmax();
      y1 = 0;
      y2 = 0;
      TLine *line_xBj = new TLine(x1, y1,x2, y2);
      line_xBj->SetLineColor(kRed);
      line_xBj->SetLineStyle(2);  // Dashed line
      line_xBj->SetLineWidth(2);
      line_xBj->Draw();
      ////////////////////////////////////////////////


    // Create a canvas and draw the histogram
    TCanvas* c19 = new TCanvas("c19","",1200,800);
    c19->SetLeftMargin(0.15); // Increase left margin
    //////////////////////////////////////////////
    nbinsX = h_inelasticityRelRes_binned->GetNbinsX();
      TGraphErrors *gy = new TGraphErrors(nbinsX);
  
      for (int j = 1; j <= nbinsX; ++j) {
        // Project X (t_reco - t_mc) for Q² bin j
        TH1D* projY = h_inelasticityRelRes_binned->ProjectionY("_pyY", j, j);
      
        TF1 *gaus = new TF1("gaus", "gaus");
        projY->Fit(gaus, "Q");  // "Q" = quiet mode, no fit output
        
        double mean  = gaus->GetParameter(1);      // Mean of the fit
        double sigma = gaus->GetParameter(2);      // Sigma of the fit
      
        double t_center = h_inelasticityRelRes_binned->GetXaxis()->GetBinCenter(j);  // Q² bin center
        double t_width  = h_inelasticityRelRes_binned->GetXaxis()->GetBinWidth(j);   // Q² bin width
     
        gy->SetPoint(j-1, t_center, mean);                   // X=t center, Y=mean of t_diff
        gy->SetPointError(j-1, 0.0*t_width/2.0, sigma);            // X error = half bin width, Y error = RMS
    
        delete projY;
      }
    
      // Style and draw
      gy->SetTitle(";y_{MC};#frac{y_{reco} - y_{MC}}{y_{MC}}");
      gy->SetMarkerStyle(20);
      gy->SetMarkerColor(kBlue+2);
      gy->SetLineColor(kBlue+2);
      gy->SetLineWidth(2);
      gy->Draw("AP");
      x1 = gy->GetXaxis()->GetXmin();
      x2 = gy->GetXaxis()->GetXmax();
      y1 = 0;
      y2 = 0;
      TLine *line_y = new TLine(x1, y1,x2, y2);
      line_y->SetLineColor(kRed);
      line_y->SetLineStyle(2);  // Dashed line
      line_y->SetLineWidth(2);
      line_y->Draw();
      ////////////////////////////////////////////////

  //---------------------------------------------------------
  // SAVE CANVASES AS PNGs AND CLOSE
  //---------------------------------------------------------
  // c1->SaveAs("figs/DDIS_eta.png");
  // c1->Close();
  // c2->SaveAs("figs/DDIS_t.png");
  c2->Close();
  // c3->SaveAs("figs/DDIS_photres.png");
  c3->Close();
  // c4->SaveAs("figs/DDIS_tPrRes.png");
  c4->Close();
  // c5->SaveAs("figs/DDIS_Q2Res.png");
  c5->Close();
  c6->SaveAs("figs/DDIS_xBjRes.png");
  c6->Close();
  // c7->SaveAs("figs/DDIS_InelasticityRes.png");
  c7->Close();
  c8->SaveAs("figs/DDIS_tPrRelRes.png");
  c8->Close();
  c9->SaveAs("figs/DDIS_Q2RelRes.png");
  c9->Close();
  c10->SaveAs("figs/DDIS_xBjRelRes.png");
  c10->Close();
  c11->SaveAs("figs/DDIS_InelasticityRelRes.png");
  c11->Close();

  c13->SaveAs("figs/DDIS_tRelRes_binned.png");
  c13->Close();

  c17->SaveAs("figs/DDIS_Q2RelRes_binned.png");
  c17->Close();
  c18->SaveAs("figs/DDIS_xBjRelRes_binned.png");
  c18->Close();
  c19->SaveAs("figs/DDIS_InelasticityRelRes_binned.png");
  c19->Close();
  ///////////////////////////////////////////////////////////////////////
  //// the reconstructed variables from the branches









// Declare variables to store branch data
Float_t electron_Q2, electron_y, electron_x;
Float_t electron_Q2_DA, electron_y_DA, electron_x_DA;
Float_t electron_Q2_truth, electron_y_truth, electron_x_truth;

// Set all branch addresses once
events->SetBranchAddress("InclusiveKinematicsElectron.Q2", &electron_Q2);
events->SetBranchAddress("InclusiveKinematicsElectron.y" , &electron_y);
events->SetBranchAddress("InclusiveKinematicsElectron.x" , &electron_x);

events->SetBranchAddress("InclusiveKinematicsDA.Q2", &electron_Q2_DA);
events->SetBranchAddress("InclusiveKinematicsDA.y" , &electron_y_DA);
events->SetBranchAddress("InclusiveKinematicsDA.x" , &electron_x_DA);

events->SetBranchAddress("InclusiveKinematicsTruth.Q2", &electron_Q2_truth);
events->SetBranchAddress("InclusiveKinematicsTruth.y" , &electron_y_truth);
events->SetBranchAddress("InclusiveKinematicsTruth.x" , &electron_x_truth);

// Get total number of entries
Long64_t nentries = events->GetEntries();


// Loop over events **once** and fill all histograms
for (Long64_t i = 0; i < nentries; i++) {
    events->GetEntry(i);

    // Fill histograms for electron method
    h_Q2_e_m->Fill(electron_Q2);
    h_y_e_m->Fill(electron_y);
    h_x_e_m->Fill(electron_x);

    // Fill histograms for DA method
    h_Q2_DA->Fill(electron_Q2_DA);
    h_y_DA->Fill(electron_y_DA);
    h_x_DA->Fill(electron_x_DA);

    // Fill histograms for truth
    h_Q2_truth->Fill(electron_Q2_truth);
    h_y_truth->Fill(electron_y_truth);
    h_x_truth->Fill(electron_x_truth);

    // Fill histograms for resolution
    h_Res_xBj->Fill(electron_x - electron_x_truth);
    h_RelRes_xBj->Fill((electron_x - electron_x_truth)/electron_x_truth);
    // h_Corr_xBj->Fill(xBj_rec, xBj_aso);
    h_xBjRelRes_binned->Fill(electron_x_truth, (electron_x - electron_x_truth)/electron_x_truth);
    /////////////////////////////////////////////////
    h_Res_Q2->Fill(electron_Q2 - electron_Q2_truth);
    h_RelRes_Q2->Fill((electron_Q2 - electron_Q2_truth)/electron_Q2_truth);
    h_RelRes_Q2_binned->Fill(electron_Q2_truth, (electron_Q2 - electron_Q2_truth)/electron_Q2_truth);
    ////////////////////////////////////////////////
    h_Res_xBj->Fill(electron_x - electron_x_truth);
    h_RelRes_xBj->Fill((electron_x - electron_x_truth)/electron_x_truth);
    h_xBjRelRes_binned->Fill(electron_x_truth, (electron_x - electron_x_truth)/electron_x_truth);
    ////////////////////////////////////////////////
    h_Res_inelasticity->Fill(electron_y - electron_y_truth);
    h_RelRes_inelasticity->Fill((electron_y - electron_y_truth)/electron_y_truth);
    h_inelasticityRelRes_binned->Fill(electron_y_truth, (electron_y - electron_y_truth)/electron_y_truth);

    h_Corr_Q2->Fill(electron_Q2_truth, electron_Q2);
    h_Corr_xBj->Fill(electron_x_truth, electron_x);
    h_Corr_inelasticity->Fill(electron_y_truth, electron_y);

}


  /////////////////////////////////////////////////////////////////////
  // Create a canvas fot Q2 hist
  TCanvas* c = new TCanvas("my_canvas", "My Canvas Title", 1200, 900); // Width, height
  // Draw the histogram on the canvas:
  // gPad->SetLogy(1);

  // h_Q2_DA->SetMaximum(100);
  h_Q2_DA->SetLineColor(kBlue);
  h_Q2_DA->SetLineWidth(1);
  h_Q2_DA->SetMarkerStyle(20);
  h_Q2_DA->SetMarkerColor(kBlue);
  h_Q2_DA->Draw("pe");

  // h_Q2_e->SetLineColor(kBlack);
  // h_Q2_e->SetLineWidth(4);
  // h_Q2_e->Draw("SAME");

  h_Q2_e_m->SetLineColor(kRed);
  h_Q2_e_m->SetLineWidth(1);
  h_Q2_e_m->SetMarkerStyle(20);
  h_Q2_e_m->SetMarkerColor(kRed);
  h_Q2_e_m->Draw("peSAME");
  // h_Q2_e_m->Draw("pSAME"); // Overlay on the same canvas

  h_Q2_truth->SetLineColor(kBlack);
  h_Q2_truth->SetLineWidth(2);
  h_Q2_truth->Draw("SAME");

  // Create a legend
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  // legend->AddEntry(h_Q2_e    , "MC: by hand", "l");
  legend->AddEntry(h_Q2_truth, "MC: truth", "l");
  legend->AddEntry(h_Q2_e_m  , "Reco.: electron method", "pe");
  legend->AddEntry(h_Q2_DA   , "Reco. DA", "pe");

  // h_Q2_DA->GetYaxis()->SetRangeUser(0, 900);
  // Draw the legend
  legend->Draw();

  // Set axis titles on the first histogram
  h_Q2_e->GetXaxis()->SetTitle("Q^{2}");
  h_Q2_e->GetYaxis()->SetTitle("# of events");
  // Optional: Adjust sizes and offsets
  h_Q2_e->GetXaxis()->SetTitleSize(0.04);
  h_Q2_e->GetYaxis()->SetTitleSize(0.04);
  h_Q2_e->GetXaxis()->SetTitleOffset(1.0);
  h_Q2_e->GetYaxis()->SetTitleOffset(1.0);

  c->Update();

  // Save the canvas as a JPG:
  c->SaveAs("figs/Q2_hist.png");

    /////////////////////////////////////////////////////////////////////
  // emptying the canvas to plot the x-hist
  c->Clear();
  gPad->SetLogy(1);

  // h_x_e->SetLineColor(kBlack);
  // h_x_e->SetLineWidth(4);
  // h_x_e->Draw();

  h_x_DA->SetLineColor(kBlue);
  h_x_DA->SetLineWidth(1);
  h_x_DA->SetMarkerStyle(20);
  h_x_DA->SetMarkerColor(kBlue);
  h_x_DA->Draw("P E");

  h_x_e_m->SetLineColor(kRed);
  h_x_e_m->SetLineWidth(1);
  h_x_e_m->SetMarkerStyle(20);
  h_x_e_m->SetMarkerColor(kRed);
  h_x_e_m->Draw("P E SAME");


  h_x_truth->SetLineColor(kBlack);
  h_x_truth->SetLineWidth(2);
  h_x_truth->Draw("SAME");

  legend->Clear();  // Removes all entries from the existing legend
  // Set axis titles on the first histogram
  legend->AddEntry(h_x_truth, "MC: truth", "l");
  legend->AddEntry(h_x_DA   , "REC.: DA", "pe");
  // legend->AddEntry(h_x_e    , "MC: by hand", "l");
  legend->AddEntry(h_x_e_m  , "Reco.: electron method", "pe");
  // Draw the legend
  legend->Draw();

  h_x_truth->GetXaxis()->SetTitle("x_{Bj}");
  h_x_truth->GetYaxis()->SetTitle("# of events");
  // Optional: Adjust sizes and offsets
  h_x_truth->GetXaxis()->SetTitleSize(0.04);
  h_x_truth->GetYaxis()->SetTitleSize(0.04);
  h_x_truth->GetXaxis()->SetTitleOffset(1.0);
  h_x_truth->GetYaxis()->SetTitleOffset(1.0);

  c->Update();

  c->SaveAs("figs/x_hist.png");

  /////////////////////////////////////////////////////////////////////
  // emptying the canvas for the Mandelstam t
  c->Clear();
  gPad->SetLogy(1);
  gPad->SetLogx(); // or SetLogy() or SetLogz()
  h_t_MC_old->SetLineColor(kBlack);
  h_t_MC_old->SetLineWidth(2);
  h_t_MC_old->Draw();

  h_t_RPP->SetLineColor(kCyan);
  h_t_RPP->Draw("peSAME");

  h_t_RP->SetLineColor(kBlue);
  h_t_RP->Draw("peSAME");

  legend->Clear();  // Removes all entries from the existing legend
  legend->SetX1NDC(0.15); // left
  legend->SetY1NDC(0.7); // bottom
  legend->SetX2NDC(0.35); // right
  legend->SetY2NDC(0.9); // top
  legend->Draw();
  legend->AddEntry(h_t_MC_old, "MC truth", "l");
  legend->AddEntry(h_t_RPP, "Reco.: RP", "pe");
  legend->AddEntry(h_t_RP , "Reco.: B0", "pe");


  // Set axis titles on the first histogram
  h_t_MC_old->GetXaxis()->SetTitle("|t| (GeV^{2})");
  h_t_MC_old->GetYaxis()->SetTitle("# of events");
  // Optional: Adjust sizes and offsets
  h_t_MC_old->GetXaxis()->SetTitleSize(0.04);
  h_t_MC_old->GetYaxis()->SetTitleSize(0.04);
  h_t_MC_old->GetXaxis()->SetTitleOffset(1.0);
  h_t_MC_old->GetYaxis()->SetTitleOffset(1.0);

  // h_t_RP->Draw("SAME");

  // legend->AddEntry(h_t_RP, "Reco.: B0", "pe");

  c->Update();

  c->SaveAs("figs/t_hist.png");
  
  c->Delete();
  TCanvas* cy = new TCanvas("cy", "My Canvas Title", 1200, 900); // Width, height
  /////////////////////////////////////////////////////////////////////
  // emptying the canvas to plot the y-hist
  cy->Clear();
  // h_y_e->SetLineColor(kBlack);
  // h_y_e->SetLineWidth(4);
  // h_y_e->Draw();
  gPad->SetLogx(0); // or SetLogy() or SetLogz()


  h_y_truth->SetLineColor(kBlack);
  h_y_truth->SetLineWidth(2);
  h_y_truth->Draw();
  // Set axis titles on the first histogram
  h_y_truth->GetXaxis()->SetTitle("y");
  h_y_truth->GetYaxis()->SetTitle("# of events");

  h_y_DA->SetLineColor(kBlue);
  h_y_DA->SetLineWidth(1);
  h_y_DA->SetMarkerStyle(20);
  h_y_DA->SetMarkerColor(kBlue);
  h_y_DA->Draw("e1p SAME");

  h_y_e_m->SetLineColor(kRed);
  h_y_e_m->SetLineWidth(1);
  h_y_e_m->SetMarkerStyle(20);
  h_y_e_m->SetMarkerColor(kRed);
  h_y_e_m->Draw("e1p SAME");

 

  legend->Clear();  // Removes all entries from the existing legend
  // Set axis titles on the first histogram
  // legend->AddEntry(h_y_e    , "MC: by hand", "l");
  TLegend *legend_y = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend_y->SetX1NDC(0.65); // left
  legend_y->SetY1NDC(0.7); // bottom
  legend_y->SetX2NDC(0.85); // right
  legend_y->SetY2NDC(0.9); // top
  legend_y->Draw();
  legend_y->AddEntry(h_y_truth, "MC: truth", "l");
  legend_y->AddEntry(h_y_DA   , "Reco.:DA", "pe");
  legend_y->AddEntry(h_y_e_m  , "Reco.: electron method", "pe");

  // Draw the legend

  h_y_e->GetXaxis()->SetTitle("y");
  h_y_e->GetYaxis()->SetTitle("# of events");
  // Optional: Adjust sizes and offsets
  h_y_e->GetXaxis()->SetTitleSize(0.04);
  h_y_e->GetYaxis()->SetTitleSize(0.04);
  h_y_e->GetXaxis()->SetTitleOffset(1.0);
  h_y_e->GetYaxis()->SetTitleOffset(1.0);

  cy->Update();

  cy->SaveAs("figs/y_hist.png");

    /////////////////////////////////////////////////////////////////////
  // emptying the canvas to plot the px-py-hist
  cy->Clear();
  gPad->SetLogy(1);

  h_py_truth->SetLineColor(kBlue);
  h_py_truth->SetLineWidth(3);
  h_py_truth->Draw();

  // h_x_DA->SetLineColor(kBlue);
  // h_x_DA->SetLineWidth(2);
  // h_x_DA->Draw("SAME");

  h_px_truth->SetLineColor(kGreen);
  h_px_truth->SetLineWidth(3);
  h_px_truth->Draw("SAME");

  // h_x_e_m->SetLineColor(kRed);
  // h_x_e_m->SetLineWidth(2);
  // h_x_e_m->Draw("SAME");

  legend->Clear();  // Removes all entries from the existing legend
  // Set axis titles on the first histogram
  // legend->AddEntry(h_x_DA   , "REC.: DA", "l");
  // legend->AddEntry(h_x_e    , "MC: by hand", "l");
  // legend->AddEntry(h_x_e_m  , "Reco.: electron method", "l");
  legend->AddEntry(h_py_truth, "MC: truth py", "l");
  legend->AddEntry(h_px_truth, "MC: truth px", "l");
  // Draw the legend
  legend->Draw();

  h_py_truth->GetXaxis()->SetTitle("p_{y}");
  h_py_truth->GetYaxis()->SetTitle("# of events");
  // Optional: Adjust sizes and offsets
  h_py_truth->GetXaxis()->SetTitleSize(0.04);
  h_py_truth->GetYaxis()->SetTitleSize(0.04);
  h_py_truth->GetXaxis()->SetTitleOffset(1.0);
  h_py_truth->GetYaxis()->SetTitleOffset(1.0);

  // c->SaveAs("figs/px_py_proton_hist.png");

  /////////////////////////////////////////////////////////////////////

    // Create a canvas and draw the histogram
    TCanvas *c12 = new TCanvas("c12", "Gaussian 2D", 800, 600);

    h_Corr_inelasticity->Draw("COLZ"); // Draw with color palette

    // Create a TLatex object to draw LaTeX
    // TLatex *latex = new TLatex();
    latex->SetNDC();        // Use normalized coordinates (0,1)
    latex->SetTextSize(0.05); // Set text size relative to canvas
    latex->SetTextColor(kBlack);
    latex->SetTextFont(42);  // Set the font (11-61, see TFont)
    latex->DrawLatex(0.15, 0.8, "electron method");
    c12->Update();
    c12->SaveAs("figs/DDIS_Inelasticity_Corr.png");
    c12->Close();


    // Create a canvas and draw the histogram
    TCanvas *c14 = new TCanvas("c14", "Gaussian 2D", 800, 600);
    // gPad->SetLogy(1);
    // h_Corr_xBj->GetXaxis()->SetLimits(0, 0.2);
    h_Corr_xBj->Draw("COLZ"); // Draw with color palette
    // Create a TLatex object to draw LaTeX
    latex->SetNDC();        // Use normalized coordinates (0,1)
    latex->SetTextSize(0.05); // Set text size relative to canvas
    latex->SetTextColor(kBlack);
    latex->SetTextFont(42);  // Set the font (11-61, see TFont)
    latex->DrawLatex(0.15, 0.8, "electron method");
    
        
    // Create a canvas and draw the histogram
    TCanvas *c15 = new TCanvas("c15", "Gaussian 2D", 800, 600);
    // gPad->SetLogy(1);
    // h_Corr_xBj->GetXaxis()->SetLimits(0, 0.2);
    h_Corr_Q2->Draw("COLZ"); // Draw with color palette
    latex->SetNDC();        // Use normalized coordinates (0,1)
    latex->SetTextSize(0.05); // Set text size relative to canvas
    latex->SetTextColor(kBlack);
    latex->SetTextFont(42);  // Set the font (11-61, see TFont)
    latex->DrawLatex(0.15, 0.8, "electron method");





    c14->SaveAs("figs/DDIS_xBj_Corr.png");
    c14->Close();
    c15->SaveAs("figs/DDIS_Q2_Corr.png");
    c15->Close();
    c16->SaveAs("figs/DDIS_t_Corr.png");
    c16->Close();


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



double rapidity(const P3MVector& p4) {
  double E = p4.E();
  double pz = p4.Pz();

  if (E + pz <= 0 || E - pz <= 0) {
    // Handle cases where the logarithm argument is non-positive
    // This could happen due to numerical precision issues or unphysical momenta
    return NAN; // Or handle the error in another way
  }

  return 0.5 * TMath::Log((E + pz) / (E - pz)); // Use TMath::Log
}

/// @brief //////////////////////////////////////////////////
/// @param p4 
/// @return 
double calculatePseudorapidity(const P3MVector& p4) {
  double px = p4.Px();
  double py = p4.Py();
  double pz = p4.Pz();
  double p = TMath::Sqrt(px * px + py * py + pz * pz); // Use TMath::Sqrt

  if (p == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (pz == p) {
      return std::numeric_limits<double>::infinity();
  }

  if (pz == -p) {
      return -std::numeric_limits<double>::infinity();
  }

  return 0.5 * TMath::Log((p + pz) / (p - pz)); // Use TMath::Log
}