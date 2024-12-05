#pragma once

#include <vector>

unsigned long long iEvent = 0;
double x = 0;
double y = 0;
double Q2 = 0;
double W2 = 0;

int partonPDGID = 0;
double partonPt = 0;
double partonEta = 0;
double partonPhi = 0;
double partonE = 0;

int iJet = 0;
double jetPt = 0;
double jetEta = 0;
double jetPhi = 0;
double jetE = 0;

unsigned int nConstituents = 0;
std::vector<int> constituentPDGID = {};
std::vector<double> constituentPt = {};
std::vector<double> constituentEta = {};
std::vector<double> constituentPhi = {};
std::vector<double> constituentE = {};

