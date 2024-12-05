#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8Plugins/FastJet3.h"
#include "eicBeamShape.hh"

#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "eicJetMaker.hh"

#include <iostream>
#include <memory>
#include <filesystem>

using namespace Pythia8;

void initTree(TTree* tree){
    //Tree Branches take the ownership of the variables defined in jetMaker.hh
    tree->Branch("iEvent", &iEvent);
    tree->Branch("x", &x);
    tree->Branch("y", &y);
    tree->Branch("Q2", &Q2);
    tree->Branch("W2", &W2);
    tree->Branch("partonPDGID", &partonPDGID);
    tree->Branch("partonPt", &partonPt);
    tree->Branch("partonEta", &partonEta);
    tree->Branch("partonPhi", &partonPhi);
    tree->Branch("partonE", &partonE);
    tree->Branch("iJet", &iJet);
    tree->Branch("jetPt", &jetPt);
    tree->Branch("jetEta", &jetEta);
    tree->Branch("jetPhi", &jetPhi);
    tree->Branch("jetE", &jetE);
    tree->Branch("nConstituents", &nConstituents);
    tree->Branch("constituentPDGID", &constituentPDGID);
    tree->Branch("constituentPt", &constituentPt);
    tree->Branch("constituentEta", &constituentEta);
    tree->Branch("constituentPhi", &constituentPhi);
    tree->Branch("constituentE", &constituentE);
}

void clear(){
    //Clear the vectors before filling them again
    constituentPDGID.clear();
    constituentPt.clear();
    constituentEta.clear();
    constituentPhi.clear();
    constituentE.clear();
}

int main(){
    // Create Pythia instance. Shorthands for event and info.
    Pythia pythia("/home/tanmaypani/heplibs/pythia8310/share/Pythia8/xmldoc");
    auto& event = pythia.event;
    auto& info = pythia.info; 

    const unsigned char configFlag = 2; //configFlag = 1 is for hiDiv, configFlag = 2 is for hiAcc, and configFlag = 3 is for when running eA energies (18x110, 10x110, 5x110, and 5x41).
    const double hadE = 275.0, lepE = 18.0, xing = 0.025; //hadE = hadron beam energy, lepE = lepton beam energy, xing = crossing angle
    // A class to generate beam parameters according to custom parametrization.
    BeamShapePtr myBeamShape = std::make_shared<eicBeamShape>(configFlag,hadE,lepE,xing);
    pythia.setBeamShapePtr(myBeamShape);// Hand pointer to Pythia.

    std::string config = "";
    switch(configFlag){
        case 1: config = "hiDiv"; break; //High divergence configuration
        case 2: config = "hiAcc"; break; //High acceptance configuration
        case 3: config = "eA"; break; //eA configuration
    }

    std::string cardFile = Form("cardFiles/dis_eicBeam_%s_%ix%i.cmnd", config.c_str(), unsigned(lepE), unsigned(hadE));
    // Read in Steering File & Make Other Settings
    pythia.readFile(cardFile);
    pythia.readString("Main:timesAllowErrors = 10000"); // allow more errors, eP is brittle
    pythia.readString("Random:setSeed = on"); //Indicates whether a user-set seed should be used every time the Pythia::init routine is called
    pythia.readString("Random:seed = 0"); //The seed to be used, A negative value gives the default seed, a value 0 gives a random seed based on the time, and a value between 1 and 900,000,000 a unique different random number sequence.

    // Initialize.
    pythia.init();

    double Q2Min = pythia.settings.parm("PhaseSpace:Q2Min");

    double jetR = 1.0;
    double jetPtMin = 1.0;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    //fastjet::Selector constituentSelector = fastjet::SelectorEtaMax(3.5)&& fastjet::SelectorPtMin(0.2) && fastjet::SelectorPtMax(100);
    fastjet::Selector jetSelector = fastjet::SelectorEtaMax(2.5) && fastjet::SelectorPtMin(5);

    unsigned long nEventsNeeded = 1000000;
    // Create ROOT TFile and TTree
    std::string outDir = "output";
    if(!std::filesystem::exists(outDir)){
        std::filesystem::create_directory(outDir);
    }
    std::string fileName = Form("%s/pythia8NCDIS_%ix%i_xAngle%0.3f_%s_Q2Min%0.1f_jetR%0.1f_nEv%li.root", outDir.c_str(), unsigned(lepE), unsigned(hadE), xing, config.c_str(), Q2Min, jetR, nEventsNeeded);
    std::unique_ptr<TFile> outFile(TFile::Open(fileName.c_str(), "RECREATE"));
    auto outTree = std::make_unique<TTree>("jet", "NC-DIS jets with matched partons");
    initTree(outTree.get());
    
    while(iEvent < nEventsNeeded){ //begin event loop
        if(!pythia.next()) continue; // Generate event; skip if generation failed.

        auto& incomingNucleus = event[1]; //2nd entry in event record is the incoming nucleus (proton/A, status code = -12)
        auto& incomingLepton = event[4]; //5th entry in event record is the incoming lepton (status code = -12)
        auto& scatteredParton = event[5]; //6th entry in event record is the scattered parton (status code = -23)
        auto& scatteredLepton = event[6]; //7th entry in event record is the scattered lepton (status code = -23)

        // Q2, W2, Bjorken x, y.
        auto photonKinematics = incomingLepton.p() - scatteredLepton.p();
        W2 = (incomingNucleus.p() + photonKinematics).mCalc();
        Q2 = -photonKinematics.m2Calc();
        x = Q2/(2.*incomingNucleus.p()*photonKinematics);
        y = (incomingNucleus.p()*photonKinematics)/(incomingNucleus.p()*incomingLepton.p());

        //event level cuts for DIS
        if(y < 0.01 || y > 0.95) continue;

        auto particle = event.begin(); //pointer to the begining of the particle list in the event record
        std::vector<fastjet::PseudoJet> allConstituents; //all particles in the event for jet clustering
        while(particle != event.end()){ //loop over all particles in the event record; stop when pointer reaches the end of the list
            unsigned int index = particle - event.begin();

            //particle level cuts
            if(index <= 7){particle++; continue;}
            if(!(particle->isFinal())){particle++; continue;}
            if(abs(particle->eta()) > 4.0){particle++; continue;}

            allConstituents.emplace_back(*particle);
            particle++; //move to the next particle in the event record, by incrementing the pointer
        }//end of loop over all particles in the event record

        //Jet clustering, ClusterSequence object is created containing history of the clustering process
        fastjet::ClusterSequence cs(allConstituents, jetDef);
        //Get the jets from the ClusterSequence object, sorted by pt and jetSelector applied to get the jets that pass the selection criteria
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(jetSelector(cs.inclusive_jets(jetPtMin)));
        auto nJets = jets.size();

        //No jets in the event, so it doesnt count
        if(nJets == 0){continue;}

        iEvent++;

        fastjet::PseudoJet disParton(scatteredParton); //create a PseudoJet object from the scattered parton for consistency with the jet objects
        iJet = 0;
        for(auto& jet : jets){ //loop over all jets clustered in the event
            clear(); //clear the vectors for the next jet

            jetPt = jet.pt();
            jetEta = jet.eta();
            jetPhi = jet.phi();
            jetE = jet.e();

            //Matching the jet to the DIS parton
            double dEta = std::abs(jet.eta() - disParton.eta());
            double dPhi = std::abs(jet.phi() - disParton.phi());
            if(dPhi > fastjet::pi) dPhi = fastjet::twopi - dPhi;
            double dR = std::sqrt(dEta*dEta + dPhi*dPhi);
            if(dR > jetR){continue;}

            iJet++;
            partonPDGID = disParton.user_info<Particle>().id();
            partonPt = disParton.pt();
            partonEta = disParton.eta();
            partonPhi = disParton.phi();
            partonE = disParton.e();

            nConstituents = jet.constituents().size();
            for(auto& constituent : jet.constituents()){//loop over all constituents of the jet
                constituentPDGID.emplace_back(constituent.user_info<Particle>().id());
                constituentPt.emplace_back(constituent.pt());
                constituentEta.emplace_back(constituent.eta());
                constituentPhi.emplace_back(constituent.phi());
                constituentE.emplace_back(constituent.e());
            }//end of loop over all constituents of the jet

            outTree->Fill(); //fill the tree with the current jet and matched parton information
        }//end of loop over all jets clustered in the event
    }//end of event loop

    // Write the tree to the file
    outTree->Write();
    // Done, print statistics.
    pythia.stat();
    return 0;
}