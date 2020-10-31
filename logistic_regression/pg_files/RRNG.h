// Copyright 2012 Jesse Windle - jesse.windle@gmail.com

// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////////////////////

// YOU MUST ALWAYS CALL GetRNGSeed() and PutRNGSeed() WHEN USING THESE FUNCTIONS!!!

//////////////////////////////////////////////////////////////////////

#ifndef __BASICRNG__
#define __BASICRNG__

#include "R.h"
#include "Rmath.h"
// #include "Matrix.h"

class BasicRNG {

 public:

  // Random variates.
  double unif  ();                             // Uniform
  double expon_mean(double mean);                  // Exponential
  double expon_rate(double rate);                  // Exponential
  double chisq (double df);                    // Chisq
  double norm  (double sd);                    // Normal
  double norm  (double mean , double sd);      // Normal
  double gamma_scale (double shape, double scale); // Gamma_Scale
  double gamma_rate  (double shape, double rate);  // Gamma_Rate
  double igamma(double shape, double scale);   // Inv-Gamma
  double flat  (double a=0  , double b=1  );   // Flat
  double beta  (double a=1.0, double b=1.0);   // Beta

  int bern  (double p);                     // Bernoulli

  // CDF
  static double p_norm (double x, int use_log=0);
  static double p_gamma_rate(double x, double shape, double rate, int use_log=0);

  // Density
  static double d_beta(double x, double a, double b);

  // Utility
  static double Gamma (double x, int use_log=0);

}; // BasicRNG

#endif
