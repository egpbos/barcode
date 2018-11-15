/*
 * Barcode
 * Copyright E.G.P. Bos
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */


////////////
// unused //
////////////

// Model: Gaussian likelihood //
void gaussian_likelihood_partial_f_delta_x_log_like_interpolate(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  // EGP: partial -log L/partial delta_x = (Lambda-nobs)/sigma**2

  fftw_array<real_prec> deltaXQ(n->N), nobsXQ(n->N), noiseXQ(n->N), windowXQ(n->N);

  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, deltaX, n->N, deltaXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->nobs, n->N, nobsXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->noise, n->N, noiseXQ);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, hd->window, n->N, windowXQ);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif // MULTITHREAD
  for(ULONG i=0;i<n->N;i++)
  {
    real_prec Lambda = windowXQ[i] * hd->rho_c * pow(num_1 + hd->biasP * deltaXQ[i], hd->biasE);
    if ((windowXQ[i]>0.) && (Lambda > 0.0))
    {
      real_prec resid = nobsXQ[i]-Lambda;
      out[i] = resid/(noiseXQ[i]*noiseXQ[i]);
    }
    else
      out[i] = 0.0;
  }
}
// N.B.: grad_f_delta_x_comp should no longer be used. Only necessary for calc_h != 2, which is bogus!
void gaussian_likelihood_grad_f_delta_x_comp_interpolate(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out,
                                                         unsigned int component)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> deltaXQ(n->N);
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, deltaX, n->N, deltaXQ);

  gradfft(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaXQ, out, component);
  //gradfindif(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, deltaX, out, component);
}

// End model: Gaussian likelihood  //


void likelihood_calc_h_old(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N), dummy(n->N), dummy2(n->N);
  fftw_array<complex_prec> outC(n->Nhalf), dummyC(n->Nhalf);

  fillZero(outC, n->Nhalf);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  for (unsigned i = 1; i <= 3; i++){
    // compute i-th component of x-gradient of f(deltaX) -> dummy
    hd->grad_f_delta_x_comp(hd, deltaX, dummy, i);
    //if (n->iGibbs == 1)
    //fill_one(dummy, n->N); // testing: when delta_1 == 0 (which it is with zero initial guess), then deltaX and its derivative will be zero, breaking the whole algorithm for this first step

    // multiply with partLike -> dummy
    multiplyArrays(partLike, dummy, dummy, n->N);
    // interpolate to x(q_i) => g_i(q)
    //{
    //switch (n->mk)
    //{
    //case 0:
    //throw runtime_error("likelihood_calc_h: NGP interpolation not yet implemented!");
    //break;
    //case 1:
    interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz,dummy,n->N,dummy2);
    //break;
    //case 2:
    //throw runtime_error("likelihood_calc_h: TSC interpolation not yet implemented!");
    //break;
    //}
    //}
    // transform to FS => g_i(k) -> AUX
    fftR2C(n->N1, n->N2, n->N3, dummy2, dummyC);
    // multiply by -ik_i/k^2
    bool rfft = true;
    grad_inv_lap_FS(n->N1,n->N2,n->N3, n->L1,n->L2,n->L3, dummyC, dummyC, i, rfft);
    // sum the resulting h(k) term to outC
    add_to_array(dummyC, outC, n->Nhalf);
  }
  // transform h(k) to real space -> h(q)
  fftC2R(n->N1, n->N2, n->N3, outC, out);
}

void likelihood_calc_h_simple(struct HAMIL_DATA *hd, real_prec *deltaX, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> partLike(n->N);

  hd->partial_f_delta_x_log_like(hd, deltaX, partLike);

  // interpolate to x(q_i) => g_i(q)
  //switch (n->mk)
  //{
  //case 0:
  //throw runtime_error("likelihood_calc_h: NGP interpolation not yet implemented!");
  //break;
  //case 1:
  interpolate_CIC(n->N1,n->N2,n->N3,n->L1,n->L2,n->L3,n->d1,n->d2,n->d3, hd->posx, hd->posy, hd->posz, partLike, n->N, out);
  //break;
  //case 2:
  //throw runtime_error("likelihood_calc_h: TSC interpolation not yet implemented!");
  //break;
  //}
}

real_prec SPH_kernel_radius(struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  return(SPH_kernel_radius(n->particle_kernel, n->particle_kernel_h));
}

void grad_SPH_kernel_3D(real_prec x, real_prec y, real_prec z, real_prec h_inv, real_prec h_sq_inv, real_prec norm,
                        real_prec &out_x, real_prec &out_y, real_prec &out_z) {
  // N.B.: when using this for density estimation, you need to normalize the
  // result afterwards with V/N! See e.g. getDensity_SPH. Same goes for
  // SPH_kernel_3D.

  // norm must be 1/(PI * h^4)

  // derivative of Monaghan kernel W_4
  real_prec r_sq = x*x + y*y + z*z;
  real_prec partial;
  // if (r_sq > 4*h_sq)
  //   partial = 0.;
  // else if (r_sq > h_sq) {
  //   real_prec r = sqrt(r_sq);
  //   real_prec q = r * h_inv;
  //   real_prec qmin2 = q - 2;
  //   partial = -0.75*qmin2*qmin2 * norm / r;
  // } else {
  //   real_prec r = sqrt(r_sq);
  //   partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
  // }
  auto q_sq_i = static_cast<int>(r_sq * h_sq_inv);
  switch (q_sq_i) {
    case 0:
    {
      real_prec r = sqrt(r_sq);
      partial = (2.25*r*h_inv*h_inv - 3*h_inv) * norm;
    }
      break;
    case 1:
    case 2:
    case 3:
    {
      real_prec r = sqrt(r_sq);
      real_prec q = r * h_inv;
      real_prec qmin2 = q - 2;
      partial = -0.75*qmin2*qmin2 * norm / r;
    }
      break;
    default:
      partial = 0.;
      break;
  }

  out_x = partial * x;
  out_y = partial * y;
  out_z = partial * z;
}

int SPH_kernel_3D_cells_count(struct HAMIL_DATA *hd) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  return(SPH_kernel_3D_cells_count(n->particle_kernel, n->particle_kernel_h,
                                   n->d1, n->d2, n->d3));
}

// std::vector<int> SPH_kernel_3D_cells(struct HAMIL_DATA *hd)
void SPH_kernel_3D_cells(struct HAMIL_DATA *hd, vector<int> &out_i,
                         vector<int> &out_j, vector<int> &out_k) {
  struct HAMIL_NUMERICAL *n = hd->numerical;
  SPH_kernel_3D_cells(n->particle_kernel, n->particle_kernel_h,
                      n->d1, n->d2, n->d3, out_i, out_j, out_k);
}

void likelihood_grad_log_like_old(struct HAMIL_DATA *hd, real_prec *delta, real_prec *out)
{
  struct HAMIL_NUMERICAL *n = hd->numerical;
  fftw_array<real_prec> dummy_h(n->N);

  // calculate deltaX and pos
  {
    unsigned facL=1;
    bool reggrid=true;
    gsl_rng *seed = nullptr; // empty: reggrid is true anyway
    real_prec kernel_scale = SPH_kernel_scale(hd);
    Lag2Eul(delta, hd->deltaX, hd->posx, hd->posy, hd->posz, n->N1, n->N2, n->N3, n->L1, n->L2, n->L3, n->d1, n->d2,
            n->d3, n->min1, n->min2, n->min3, hd->D1, hd->D2, hd->ascale, hd->OM, hd->OL, hd->sfmodel, n->mk, n->kth, facL,
            reggrid, seed, "", kernel_scale, n->R2Cplan, n->C2Rplan);
  }

  likelihood_calc_h_old(hd, hd->deltaX, dummy_h);

  // EGP: for Zel'dovich, -h is the entire term. For other models more has to be done.
  multiply_factor_array(-1., dummy_h, out, n->N); // EGP: Zel'dovich: -gradLogLike = -h
}
