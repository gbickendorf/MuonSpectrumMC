#include "Spectrum.h"

Spectrum::Spectrum() {}

Spectrum::Spectrum(int _orderX, int _orderY, vector<vector<double>> _d_spectrum,
         vector<double> _X_values, vector<double> _Y_values,
         double _totalWidth) {
  orderX = _orderX;
  orderY = _orderY;
  X_values = _X_values;
  Y_values = _Y_values;
  d_spectrum = _d_spectrum;
  totalWidth = _totalWidth;
}
