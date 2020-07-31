#pragma once
#include <stdio.h>
#include <math.h>

//Same to lepton
double M2ScalarToLepton(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2VectorTauToElectronTVCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);

//Vector electron coupling
double M2VectorElectronRLCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);

//Vecotr muon coupling
double M2VectorMuonRLCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2VectorMuonLCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2VectorElectronLCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2VectorMuonRCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2VectorElectronRCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);

//Scalar to muon
double M2ScalarMuonCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2ScalarElectronCoupling2(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);


double M2ScalarYukawa(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);

//Scalar to electron
double M2ScalarElectronCoupling(double E1, double E2, double E3, double ct1, double ct2, double ct3, double ph2, double ph3, double Gf, double m, double me, double m2, double g);
