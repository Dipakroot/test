     float ditau_m = sortedDiTaus.diTauMass_;


      //cout << "ditau mass=" << ditau_m << endl;
    
      _DitauMass.push_back(ditau_m);
    
    }










    if (dileptonsample){
      //std::cout << "starting tau selection..." << std::endl;

      //added by Dipak
      if (_nTau==2) {
	//std::cout << "starting tau selection......" << std::endl;

	int TauInd1 = selector->Taus.at(0);
	int TauInd2 = selector->Taus.at(1);

	lepVector.SetPtEtaPhiM(tree->tauPt_[TauInd1],
			       tree->tauEta_[TauInd1],
			       tree->tauPhi_[TauInd1],
			       tree->tauMass_[TauInd1]);
	lepVector2.SetPtEtaPhiM(tree->tauPt_[TauInd2],
				tree->tauEta_[TauInd2],
				tree->tauPhi_[TauInd2],
				tree->tauMass_[TauInd2]);

	//_DitauMass = (lepVector+lepVector2).M();
	_DitauDelR = lepVector.DeltaR(lepVector2);
	
	//	nevents= nevents+1;
	//	std::cout << "Number of events with (ntau==2)= "<< nevents << std::endl;


      }

      // std::cout << "Ending tau loop" << std::endl;



	if (_nMu==2) {

	    int muInd1 = selector->Muons.at(0);
	    int muInd2 = selector->Muons.at(1);

	    lepVector.SetPtEtaPhiM(tree->muPt_[muInd1],
				   tree->muEta_[muInd1],
				   tree->muPhi_[muInd1],
				   tree->muMass_[muInd1]);
	    lepVector2.SetPtEtaPhiM(tree->muPt_[muInd2],
				    tree->muEta_[muInd2],
				    tree->muPhi_[muInd2],
				    tree->muMass_[muInd2]);	

	    _DilepMass = (lepVector+lepVector2).M();
	    _DilepDelR = lepVector.DeltaR(lepVector2);
	    
	}
	
	
	if (_nEle==2){
	    int eleInd1 = selector->Electrons.at(0);
	    int eleInd2 = selector->Electrons.at(1);
	    
	    lepVector.SetPtEtaPhiM(tree->elePt_[eleInd1],
				   tree->eleEta_[eleInd1],
				   tree->elePhi_[eleInd1],
				   tree->eleMass_[eleInd1]);
	    lepVector2.SetPtEtaPhiM(tree->elePt_[eleInd2],
				    tree->eleEta_[eleInd2],
				    tree->elePhi_[eleInd2],
				    tree->eleMass_[eleInd2]);
	    
	    _DilepMass = (lepVector+lepVector2).M();
	    _DilepDelR = lepVector.DeltaR(lepVector2);
	    
	}
    }
    
    _passPresel_Ele  = evtPick->passPresel_ele;
    _passPresel_Mu   = evtPick->passPresel_mu;
    _passAll_Ele     = evtPick->passAll_ele;
    _passAll_Mu      = evtPick->passAll_mu;
    //added by Dipak
    //_pass_trigger_tau  = evtPick->Pass_trigger_tau;
    _passPresel_Tau   = evtPick->passPresel_tau;

    int parentPID = -1;
    for (int i_jet = 0; i_jet < _nfwdJet; i_jet++){
	int jetInd = selector->FwdJets.at(i_jet);
	_fwdJetPt.push_back(tree->jetPt_[jetInd]);
	_fwdJetEta.push_back(tree->jetEta_[jetInd]);
	_fwdJetPhi.push_back(tree->jetPhi_[jetInd]);
	_fwdJetMass.push_back(tree->jetMass_[jetInd]);
    }


    jetVectors.clear();
    jetResolutionVectors.clear();
    jetBtagVectors.clear();

    for (int i_jet = 0; i_jet <_nJet; i_jet++){
		
	int jetInd = selector->Jets.at(i_jet);
	_jetPt.push_back(tree->jetPt_[jetInd]);
	_jetEta.push_back(tree->jetEta_[jetInd]);
	_jetPhi.push_back(tree->jetPhi_[jetInd]);
	_jetMass.push_back(tree->jetMass_[jetInd]);

	_jetCMVA.push_back(tree->jetBtagCMVA_[jetInd]);
	_jetCSVV2.push_back(tree->jetBtagCSVV2_[jetInd]);
	_jetDeepB.push_back(tree->jetBtagDeepB_[jetInd]);
	_jetDeepC.push_back(tree->jetBtagDeepC_[jetInd]);
	
	jetVector.SetPtEtaPhiM(tree->jetPt_[jetInd], tree->jetEta_[jetInd], tree->jetPhi_[jetInd], tree->jetMass_[jetInd]);
	
	_jetGenJetIdx.push_back(tree->jetGenJetIdx_[jetInd]);

	double resolution = selector->jet_resolution.at(i_jet);

	jetVectors.push_back(jetVector);
	jetResolutionVectors.push_back(resolution);
	jetBtagVectors.push_back(tree->jetBtagDeepB_[jetInd]);
	
	if (selector->jet_isTagged.at(i_jet)){
	    //	if ( (tree->jetBtagDeepB_[jetInd]) > selector->btag_cut_DeepCSV){
	    bjetVectors.push_back(jetVector);
	    bjetResVectors.push_back(resolution);
	} else {
	    ljetVectors.push_back(jetVector);
	    ljetResVectors.push_back(resolution);
	}
    }	

    //Compute M3
    _M3 = -1.;
    double maxPt = -1;
    if (_nJet>2) {
	TLorentzVector jet1;
	TLorentzVector jet2;
	TLorentzVector jet3;
	for (int i_jet1 = 0; i_jet1 <_nJet-2; i_jet1++){
	    int jetInd1 = selector->Jets.at(i_jet1);
	    jet1.SetPtEtaPhiM(tree->jetPt_[jetInd1],tree->jetEta_[jetInd1],tree->jetPhi_[jetInd1],tree->jetMass_[jetInd1]);
	    
	    for (int i_jet2 = i_jet1+1; i_jet2 <_nJet-1; i_jet2++){
		int jetInd2 = selector->Jets.at(i_jet2);
		jet2.SetPtEtaPhiM(tree->jetPt_[jetInd2],tree->jetEta_[jetInd2],tree->jetPhi_[jetInd2],tree->jetMass_[jetInd2]);

		for (int i_jet3 = i_jet2+1; i_jet3 <_nJet; i_jet3++){
		    int jetInd3 = selector->Jets.at(i_jet3);
		    jet3.SetPtEtaPhiM(tree->jetPt_[jetInd3],tree->jetEta_[jetInd3],tree->jetPhi_[jetInd3],tree->jetMass_[jetInd3]);

		    if ((jet1 + jet2 + jet3).Pt()>maxPt){
			_M3 = (jet1 + jet2 + jet3).M();
			maxPt=(jet1 + jet2 + jet3).Pt();
		    }
		    
		}
	    }
	}
    }

    
    // // // Calculate MET z
    // metZ.SetLepton(lepVector);

    // METVector.SetPtEtaPhiM(tree->MET_pt_,
    // 			   0.,
    // 			   tree->MET_phi_,
    // 			   0.);
	
    // metZ.SetMET(METVector);

    //Calculate transverse mass variables
    //W transverse mass		

    _WtransMass = TMath::Sqrt(2*lepVector.Pt()*tree->MET_pt_*( 1.0 - TMath::Cos( lepVector.DeltaPhi(METVector))));


    // TLorentzVector tempLep;
    // tempLep.SetPtEtaPhiM(lepVector.Pt(),
    // 			 lepVector.Eta(),
    // 			 lepVector.Phi(),
    // 			 0.1056);

    double _met_px = METVector.Px();
    double _met_py = METVector.Py();

    // _nu_pz = metZ.Calculate();
    // _nu_pz_other = metZ.getOther();
    // METVector.SetPz(_nu_pz);
    // cout << _nu_pz << "   ----   " << _nu_pz_other << endl;

    // for (int __j = 0; __j < isBjet.size(); __j++){
    // 	if (isBjet.at(__j)) b_ind.push_back(__j);
    // 	else j_ind.push_back(__j);
    // }
    

    // topEvent.SetBJetVector(bjetVectors);
    // topEvent.SetLJetVector(ljetVectors);
    // topEvent.SetLepton(lepVector);
    // topEvent.SetMET(METVector);
    
    // topEvent.SetBJetResVector(bjetResVectors);
    // topEvent.SetLJetResVector(ljetResVectors);
    // topEvent.SetIgnoreBtag(true);

    topEvent.SetJetVector(jetVectors);
    topEvent.SetJetResVector(jetResolutionVectors);
    topEvent.SetBtagVector(jetBtagVectors);

    topEvent.SetLepton(lepVector);
    topEvent.SetMET(METVector);
    
    topEvent.Calculate();

    if (topEvent.GoodCombination()){
	bhad = jetVectors[topEvent.getBHad()];
	blep = jetVectors[topEvent.getBLep()];
	Wj1 = jetVectors[topEvent.getJ1()];
	Wj2 = jetVectors[topEvent.getJ2()];
	METVector.SetPz(topEvent.getNuPz());

	_chi2 = topEvent.getChi2();
	_M_bjj = ( bhad + Wj1 + Wj2 ).M();
	_M_jj  = ( Wj1 + Wj2 ).M();
	
	_TopHad_pt = ( bhad + Wj1 + Wj2 ).Pt();
	_TopHad_eta = ( bhad + Wj1 + Wj2 ).Eta();
	_TopLep_pt = ( blep + lepVector + METVector ).Pt();
	_TopLep_eta = ( blep + lepVector + METVector ).Eta();
	_TopLep_charge = lepCharge;

	_MassCuts = (_M_bjj > 160 && 
		      _M_bjj < 180 && 
		      _M_jj > 70 &&
		      _M_jj < 90);

	
    }
    
    ljetVectors.clear();
    bjetVectors.clear();
    
    ljetResVectors.clear();
    bjetResVectors.clear();
    
    
    if (isMC){

	// Float_t LHE scale variation weights (w_var / w_nominal); 
	// [0] is mur=0.5 muf=0.5 ; 
	// [1] is mur=0.5 muf=1 ; 
	// [2] is mur=0.5 muf=2 ; 
	// [3] is mur=1 muf=0.5 ; 
	// [4] is mur=1 muf=1 ; 
	// [5] is mur=1 muf=2 ; 
	// [6] is mur=2 muf=0.5 ; 
	// [7] is mur=2 muf=1 ; 
	// [8] is mur=2 muf=2 

	_q2weight_Up = 1.;
	_q2weight_Do = 1.;

	if (tree->nLHEScaleWeight_==9){
	    for (int i = 0; i < 9; i++){
		if(i==2||i==6){continue;}
		_genScaleSystWeights.push_back(tree->LHEScaleWeight_[i]);
	    }
            double nomWeight=tree->LHEScaleWeight_[4];
            if (nomWeight!=0){
                _q2weight_Up = *max_element(_genScaleSystWeights.begin(), _genScaleSystWeights.end())/nomWeight;
                _q2weight_Do = *min_element(_genScaleSystWeights.begin(), _genScaleSystWeights.end())/nomWeight;
            }
	}

	if (tree->nLHEScaleWeight_==44){
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[0]);
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[5]);
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[15]);
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[24]);
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[34]);
	    _genScaleSystWeights.push_back(tree->LHEScaleWeight_[39]);

	    _q2weight_Up = *max_element(_genScaleSystWeights.begin(), _genScaleSystWeights.end());
	    _q2weight_Do = *min_element(_genScaleSystWeights.begin(), _genScaleSystWeights.end());
	}

	double pdfMean = 0.;
	for (int j=0; j < tree->nLHEPdfWeight_; j++ ){
	    _pdfSystWeight.push_back(tree->LHEPdfWeight_[j]);
	    pdfMean += tree->LHEPdfWeight_[j];
	}
	pdfMean = pdfMean/_pdfSystWeight.size();
	    
	double pdfVariance = 0.;
	for (int j=0; j < _pdfSystWeight.size(); j++){
	    pdfVariance += pow((_pdfSystWeight[j]-pdfMean),2.);
	}
        if (pdfMean=0) pdfMean=1;
	_pdfuncer = sqrt(pdfVariance/_pdfSystWeight.size())/pdfMean;
	_pdfweight_Up = (1. + _pdfuncer);
	_pdfweight_Do = (1. - _pdfuncer);

	_ISRweight_Up = 1.;
	_ISRweight_Do = 1.;

	_FSRweight_Up = 1.;
	_FSRweight_Do = 1.;
	
	if (tree->nPSWeight_==4){
            if (tree->genWeight_ != 0){
                _ISRweight_Up = tree->PSWeight_[2];
                _ISRweight_Do = tree->PSWeight_[0];

                _FSRweight_Up = tree->PSWeight_[3];
                _FSRweight_Do = tree->PSWeight_[1];
            }
        }

    }

    
    
    for (int i_mc = 0; i_mc <_nGenPart; i_mc++){
	_genPt.push_back(tree->GenPart_pt_[i_mc]);
	_genPhi.push_back(tree->GenPart_phi_[i_mc]);
	_genEta.push_back(tree->GenPart_eta_[i_mc]);
	_genMass.push_back(tree->GenPart_mass_[i_mc]);
	_genStatus.push_back(tree->GenPart_status_[i_mc]);
	_genStatusFlag.push_back(tree->GenPart_statusFlags_[i_mc]);
	_genPDGID.push_back(tree->GenPart_pdgId_[i_mc]);
	_genMomIdx.push_back(tree->GenPart_genPartIdxMother_[i_mc]);
    }

    for (int i_genJet = 0; i_genJet < _nGenJet; i_genJet++){
	_genJetPt.push_back(tree->GenJet_pt_[i_genJet]);
	_genJetEta.push_back(tree->GenJet_eta_[i_genJet]);
	_genJetPhi.push_back(tree->GenJet_phi_[i_genJet]);
	_genJetMass.push_back(tree->GenJet_mass_[i_genJet]);
    }
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
double makeRecoNtuple::SFtop(double pt){
    return exp(0.0615 - 0.0005*pt);
}

double makeRecoNtuple::topPtWeight(){
    double toppt=0.0;
    double antitoppt=0.0;
    double weight = 1.0;

    // TODO needs to be reimplemented with NANOAOD
    for(int mcInd=0; mcInd<tree->nGenPart_; ++mcInd){
    	if(tree->GenPart_pdgId_[mcInd]==6  && tree->GenPart_statusFlags_[mcInd]>>13&1) toppt = tree->GenPart_pt_[mcInd];
    	if(tree->GenPart_pdgId_[mcInd]==-6 && tree->GenPart_statusFlags_[mcInd]>>13&1) antitoppt = tree->GenPart_pt_[mcInd];
    }
    if(toppt > 0.001 && antitoppt > 0.001)
	weight = sqrt( SFtop(toppt) * SFtop(antitoppt) );
    
    //This has been changed, the new prescription is to not use the top pt reweighting, and the syst is using it
    return weight;
    
}

void makeRecoNtuple::loadBtagEff(string sampleName, string year){
    std::string fName = "weight/BtagSF/btag_efficiencies_"+year+".root";
    std::string effType = "Other";
    if (sampleType.find("TTGamma") != std::string::npos){
	effType = "Top";
    }
    if (sampleType.find("TTbar") != std::string::npos){
	effType = "Top";
    }
    
    std::string leffName = effType+"_l_efficiency";
    std::string ceffName = effType+"_c_efficiency";
    std::string beffName = effType+"_b_efficiency";

    TFile* inputFile = TFile::Open(fName.c_str(),"read");
    l_eff = (TH2D*) inputFile->Get(leffName.c_str());
    c_eff = (TH2D*) inputFile->Get(ceffName.c_str());
    b_eff = (TH2D*) inputFile->Get(beffName.c_str());
}				   

float makeRecoNtuple::getBtagSF_1a(string sysType, BTagCalibrationReader reader, bool verbose){

    double weight = 1.0;

    double jetPt;
    double jetEta;
    double jetBtag;
    int jetFlavor;
    double SFb;
    double Eff;

    double pMC=1;
    double pData=1;
	
    string b_sysType = "central";
    string l_sysType = "central";
    if (sysType=="b_up"){
	b_sysType = "up";
    } else if (sysType=="b_down"){
	b_sysType = "down";
    } else if (sysType=="l_up"){
	l_sysType = "up";
    } else if (sysType=="l_down"){
	l_sysType = "down";
    }	
    if (verbose){
	cout << "Btagging Scale Factors"<<endl;
    }

    for(std::vector<int>::const_iterator jetInd = selector->Jets.begin(); jetInd != selector->Jets.end(); jetInd++){

	jetPt = tree->jetPt_[*jetInd];
	jetEta = fabs(tree->jetEta_[*jetInd]);
	jetFlavor = abs(tree->jetHadFlvr_[*jetInd]);
	jetBtag = tree->jetBtagDeepB_[*jetInd];

	if (jetFlavor == 5){
	    SFb = reader.eval_auto_bounds(b_sysType, BTagEntry::FLAV_B, jetEta, jetPt); 
	    int xbin = b_eff->GetXaxis()->FindBin(min(jetPt,799.));
	    int ybin = b_eff->GetYaxis()->FindBin(abs(jetEta));
	    Eff = b_eff->GetBinContent(xbin,ybin);
	}
	else if(jetFlavor == 4){
	    SFb = reader.eval_auto_bounds(b_sysType, BTagEntry::FLAV_C, jetEta, jetPt); 
	    int xbin = c_eff->GetXaxis()->FindBin(min(jetPt,799.));
	    int ybin = c_eff->GetYaxis()->FindBin(abs(jetEta));
	    Eff = c_eff->GetBinContent(xbin,ybin);
	}
	else {
	    SFb = reader.eval_auto_bounds(l_sysType, BTagEntry::FLAV_UDSG, jetEta, jetPt); 
	    int xbin = l_eff->GetXaxis()->FindBin(min(jetPt,799.));
	    int ybin = l_eff->GetYaxis()->FindBin(abs(jetEta));
	    Eff = l_eff->GetBinContent(xbin,ybin);
	}

	if (jetBtag>selector->btag_cut_DeepCSV){
	    pMC *= Eff;
	    pData *= Eff*SFb;
	} else {
	    pMC *= 1. - Eff;
	    pData *= 1. - (Eff*SFb);
	}
	if (verbose){
	    cout << "    jetPt="<<jetPt<<"  jetEta="<<jetEta<<"  jetFlavor="<<jetFlavor<<"  jetBtag="<<jetBtag<<"  Tagged="<<(jetBtag>selector->btag_cut_DeepCSV)<<"  Eff="<<Eff<<"  SF="<<SFb<<endl;
	    cout << "          --p(MC)="<<pMC<<"  --p(Data)="<<pData << endl;
	}
    }

    //    weight = pData/pMC;
    if (pMC==0){
	//      cout << "Inf weight" << endl;
	//	cout << pData << " / " << pMC << endl;
	weight = 0.;
    } else {
	weight = pData/pMC;
    }
    if (verbose){
	cout << "  FinalWeight="<<weight<<endl;
    }

    return weight;

}


vector<float> makeRecoNtuple::getBtagSF_1c(string sysType, BTagCalibrationReader reader, vector<float> &btagSF){

    // Saving weights w(0|n), w(1|n), w(2|n)
    vector<float> btagWeights;

    double weight0tag = 1.0; 		//w(0|n)
    double weight1tag = 0.0;		//w(1|n)

    double jetPt;
    double jetEta;
    int jetFlavor;
    double SFb;
	
    string b_sysType = "central";
    string l_sysType = "central";
    if (sysType=="b_up"){
	b_sysType = "up";
    } else if (sysType=="b_down"){
	b_sysType = "down";
    } else if (sysType=="l_up"){
	l_sysType = "up";
    } else if (sysType=="l_down"){
	l_sysType = "down";
    }	


    for(std::vector<int>::const_iterator bjetInd = selector->bJets.begin(); bjetInd != selector->bJets.end(); bjetInd++){
	jetPt = tree->jetPt_[*bjetInd];
	jetEta = fabs(tree->jetEta_[*bjetInd]);
	jetFlavor = abs(tree->jetHadFlvr_[*bjetInd]);
		
	if (jetFlavor == 5) SFb = reader.eval_auto_bounds(b_sysType, BTagEntry::FLAV_B, jetEta, jetPt); 
	else if(jetFlavor == 4) SFb = reader.eval_auto_bounds(b_sysType, BTagEntry::FLAV_C, jetEta, jetPt); 
	else {
	    SFb = reader.eval_auto_bounds(l_sysType, BTagEntry::FLAV_UDSG, jetEta, jetPt); 
	}

	btagSF.push_back(SFb);
    }

    if(selector->bJets.size() == 0) {
	btagWeights.push_back(1.0);
	btagWeights.push_back(0.0);
	btagWeights.push_back(0.0);

	return btagWeights;

    } else if (selector->bJets.size() == 1) {
	btagWeights.push_back(1-btagSF.at(0));
	btagWeights.push_back(btagSF.at(0));
	btagWeights.push_back(0.0);
		
	return btagWeights;

    } else {

	// We are following the method 1SFc from the twiki
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1c_Event_reweighting_using_scale
	for (int i = 0; i < selector->bJets.size(); i++){
	    SFb = btagSF.at(i);
	    weight0tag *= 1.0 - SFb;
	    double prod = SFb;
	    for (int j = 0; j < selector->bJets.size(); j++){
		if (j==i) {continue;}
		prod *= (1.-btagSF.at(j));
	    }
	    weight1tag += prod;
	}
	btagWeights.push_back(weight0tag);
	btagWeights.push_back(weight1tag);
	btagWeights.push_back(1.0 - weight0tag - weight1tag);
	return btagWeights;
    }
}


#endif

int main(int ac, char** av){
  makeRecoNtuple(ac, av);
  return 0;
}
