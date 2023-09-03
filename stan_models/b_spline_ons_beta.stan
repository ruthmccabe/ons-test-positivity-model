// Copyright 2023 Ruth McCabe
// Copyright 2021 Oliver Eales
// Copyright 2017 Milad Kharratzadeh
//
  // Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
  //
  // 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
  // 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
  // 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}


data {
  int num_data;              // number of data points
  int num_knots;             // num of knots
  vector[num_knots] knots;   // the sequence of knots
  int spline_degree;         // the degree of spline (is equal to order - 1)
  real X[num_data];
  real alpha_ons[num_data];  // parameter of beta distribution
  real beta_ons[num_data];   // parameter of beta distribution
  //real Y_mean[num_data];
  //real Y_var[num_data];      // variance of the ONS data 
}

transformed data {
  real Y[num_data];
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  for (i in 1:num_data){
    Y[i] = logit(beta_rng(alpha_ons[i],beta_ons[i]));
  };
}


parameters {
  row_vector[num_basis] a_raw;
  real<lower=0> tau;
  real<lower=0> gamma;
  // real a0;  // intercept
  //real Y_est[num_data];
  // alpha_est[num_data];
  //real beta_est[num_data];
}


transformed parameters {
  row_vector[num_basis] a;
  vector[num_data] Y_hat;
  //vector[num_data] Y_hat_logit;
  //vector[num_data] alpha_est;
  //vector[num_data] beta_est;
  a[1] = a_raw[1];
  //for (i in 2:num_basis)
    //  a[i] = a[i-1] + a_raw[i]*tau;
    a[2] = a_raw[2];
    for (i in 3:num_basis)
      //a[i] = 2*a[i-1] - a[i-2] + a_raw[i]*tau;
      a[i] = a[i-1] + a_raw[i]*tau;
    Y_hat =  to_vector(a*B);
    //for (i in 1:num_data){
      //  Y_hat_logit[i] = exp(Y_hat[i])/(1+exp(Y_hat[i]));
      //  alpha_est[i] = Y_hat_logit[i]*(((Y_hat_logit[i] * (1- Y_hat_logit[i]))/Y_var[i]) - 1);
      //  beta_est[i] = alpha_est[i]*((1 - Y_hat_logit[i])/Y_hat_logit[i]);
      //}
    
}

model {
  // Priors
  a_raw[1:num_basis] ~ normal(0,1);
  tau ~ inv_gamma(0.0001, 0.0001);
  //gamma ~ normal(0,1);
  gamma ~ inv_gamma(0.0001,0.0001);
  //Y_est ~ beta(alpha_ons,beta_ons);
  //a_raw ~ normal(0, 1);
  //a0 ~ normal(0, 1);
  //tau ~ normal(0, 1);
  
  //Likelihood
  Y ~ normal(Y_hat,gamma);
  //Y ~ normal(Y_hat,Y_var);
  //Y ~ binomial_logit(N,Y_hat);
} 

