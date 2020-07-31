#pragma once
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

class Spectrum {
public:
  vector<vector<double>> d_spectrum;
  int orderX;
  int orderY;
  vector<double> X_values;
  vector<double> Y_values;
  double totalWidth;
  Spectrum();
  Spectrum(int _orderX, int _orderY, vector<vector<double>> _d_spectrum,
           vector<double> _X_values, vector<double> _Y_values,
           double _totalWidth);
};
