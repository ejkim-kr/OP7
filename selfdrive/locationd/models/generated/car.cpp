#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3733011136508428965) {
   out_3733011136508428965[0] = delta_x[0] + nom_x[0];
   out_3733011136508428965[1] = delta_x[1] + nom_x[1];
   out_3733011136508428965[2] = delta_x[2] + nom_x[2];
   out_3733011136508428965[3] = delta_x[3] + nom_x[3];
   out_3733011136508428965[4] = delta_x[4] + nom_x[4];
   out_3733011136508428965[5] = delta_x[5] + nom_x[5];
   out_3733011136508428965[6] = delta_x[6] + nom_x[6];
   out_3733011136508428965[7] = delta_x[7] + nom_x[7];
   out_3733011136508428965[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7214253961443358426) {
   out_7214253961443358426[0] = -nom_x[0] + true_x[0];
   out_7214253961443358426[1] = -nom_x[1] + true_x[1];
   out_7214253961443358426[2] = -nom_x[2] + true_x[2];
   out_7214253961443358426[3] = -nom_x[3] + true_x[3];
   out_7214253961443358426[4] = -nom_x[4] + true_x[4];
   out_7214253961443358426[5] = -nom_x[5] + true_x[5];
   out_7214253961443358426[6] = -nom_x[6] + true_x[6];
   out_7214253961443358426[7] = -nom_x[7] + true_x[7];
   out_7214253961443358426[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2366253611531111470) {
   out_2366253611531111470[0] = 1.0;
   out_2366253611531111470[1] = 0;
   out_2366253611531111470[2] = 0;
   out_2366253611531111470[3] = 0;
   out_2366253611531111470[4] = 0;
   out_2366253611531111470[5] = 0;
   out_2366253611531111470[6] = 0;
   out_2366253611531111470[7] = 0;
   out_2366253611531111470[8] = 0;
   out_2366253611531111470[9] = 0;
   out_2366253611531111470[10] = 1.0;
   out_2366253611531111470[11] = 0;
   out_2366253611531111470[12] = 0;
   out_2366253611531111470[13] = 0;
   out_2366253611531111470[14] = 0;
   out_2366253611531111470[15] = 0;
   out_2366253611531111470[16] = 0;
   out_2366253611531111470[17] = 0;
   out_2366253611531111470[18] = 0;
   out_2366253611531111470[19] = 0;
   out_2366253611531111470[20] = 1.0;
   out_2366253611531111470[21] = 0;
   out_2366253611531111470[22] = 0;
   out_2366253611531111470[23] = 0;
   out_2366253611531111470[24] = 0;
   out_2366253611531111470[25] = 0;
   out_2366253611531111470[26] = 0;
   out_2366253611531111470[27] = 0;
   out_2366253611531111470[28] = 0;
   out_2366253611531111470[29] = 0;
   out_2366253611531111470[30] = 1.0;
   out_2366253611531111470[31] = 0;
   out_2366253611531111470[32] = 0;
   out_2366253611531111470[33] = 0;
   out_2366253611531111470[34] = 0;
   out_2366253611531111470[35] = 0;
   out_2366253611531111470[36] = 0;
   out_2366253611531111470[37] = 0;
   out_2366253611531111470[38] = 0;
   out_2366253611531111470[39] = 0;
   out_2366253611531111470[40] = 1.0;
   out_2366253611531111470[41] = 0;
   out_2366253611531111470[42] = 0;
   out_2366253611531111470[43] = 0;
   out_2366253611531111470[44] = 0;
   out_2366253611531111470[45] = 0;
   out_2366253611531111470[46] = 0;
   out_2366253611531111470[47] = 0;
   out_2366253611531111470[48] = 0;
   out_2366253611531111470[49] = 0;
   out_2366253611531111470[50] = 1.0;
   out_2366253611531111470[51] = 0;
   out_2366253611531111470[52] = 0;
   out_2366253611531111470[53] = 0;
   out_2366253611531111470[54] = 0;
   out_2366253611531111470[55] = 0;
   out_2366253611531111470[56] = 0;
   out_2366253611531111470[57] = 0;
   out_2366253611531111470[58] = 0;
   out_2366253611531111470[59] = 0;
   out_2366253611531111470[60] = 1.0;
   out_2366253611531111470[61] = 0;
   out_2366253611531111470[62] = 0;
   out_2366253611531111470[63] = 0;
   out_2366253611531111470[64] = 0;
   out_2366253611531111470[65] = 0;
   out_2366253611531111470[66] = 0;
   out_2366253611531111470[67] = 0;
   out_2366253611531111470[68] = 0;
   out_2366253611531111470[69] = 0;
   out_2366253611531111470[70] = 1.0;
   out_2366253611531111470[71] = 0;
   out_2366253611531111470[72] = 0;
   out_2366253611531111470[73] = 0;
   out_2366253611531111470[74] = 0;
   out_2366253611531111470[75] = 0;
   out_2366253611531111470[76] = 0;
   out_2366253611531111470[77] = 0;
   out_2366253611531111470[78] = 0;
   out_2366253611531111470[79] = 0;
   out_2366253611531111470[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3749368701677080985) {
   out_3749368701677080985[0] = state[0];
   out_3749368701677080985[1] = state[1];
   out_3749368701677080985[2] = state[2];
   out_3749368701677080985[3] = state[3];
   out_3749368701677080985[4] = state[4];
   out_3749368701677080985[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3749368701677080985[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3749368701677080985[7] = state[7];
   out_3749368701677080985[8] = state[8];
}
void F_fun(double *state, double dt, double *out_6714729038286529007) {
   out_6714729038286529007[0] = 1;
   out_6714729038286529007[1] = 0;
   out_6714729038286529007[2] = 0;
   out_6714729038286529007[3] = 0;
   out_6714729038286529007[4] = 0;
   out_6714729038286529007[5] = 0;
   out_6714729038286529007[6] = 0;
   out_6714729038286529007[7] = 0;
   out_6714729038286529007[8] = 0;
   out_6714729038286529007[9] = 0;
   out_6714729038286529007[10] = 1;
   out_6714729038286529007[11] = 0;
   out_6714729038286529007[12] = 0;
   out_6714729038286529007[13] = 0;
   out_6714729038286529007[14] = 0;
   out_6714729038286529007[15] = 0;
   out_6714729038286529007[16] = 0;
   out_6714729038286529007[17] = 0;
   out_6714729038286529007[18] = 0;
   out_6714729038286529007[19] = 0;
   out_6714729038286529007[20] = 1;
   out_6714729038286529007[21] = 0;
   out_6714729038286529007[22] = 0;
   out_6714729038286529007[23] = 0;
   out_6714729038286529007[24] = 0;
   out_6714729038286529007[25] = 0;
   out_6714729038286529007[26] = 0;
   out_6714729038286529007[27] = 0;
   out_6714729038286529007[28] = 0;
   out_6714729038286529007[29] = 0;
   out_6714729038286529007[30] = 1;
   out_6714729038286529007[31] = 0;
   out_6714729038286529007[32] = 0;
   out_6714729038286529007[33] = 0;
   out_6714729038286529007[34] = 0;
   out_6714729038286529007[35] = 0;
   out_6714729038286529007[36] = 0;
   out_6714729038286529007[37] = 0;
   out_6714729038286529007[38] = 0;
   out_6714729038286529007[39] = 0;
   out_6714729038286529007[40] = 1;
   out_6714729038286529007[41] = 0;
   out_6714729038286529007[42] = 0;
   out_6714729038286529007[43] = 0;
   out_6714729038286529007[44] = 0;
   out_6714729038286529007[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6714729038286529007[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6714729038286529007[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6714729038286529007[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6714729038286529007[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6714729038286529007[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6714729038286529007[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6714729038286529007[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6714729038286529007[53] = -9.8000000000000007*dt;
   out_6714729038286529007[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6714729038286529007[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6714729038286529007[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6714729038286529007[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6714729038286529007[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6714729038286529007[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6714729038286529007[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6714729038286529007[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6714729038286529007[62] = 0;
   out_6714729038286529007[63] = 0;
   out_6714729038286529007[64] = 0;
   out_6714729038286529007[65] = 0;
   out_6714729038286529007[66] = 0;
   out_6714729038286529007[67] = 0;
   out_6714729038286529007[68] = 0;
   out_6714729038286529007[69] = 0;
   out_6714729038286529007[70] = 1;
   out_6714729038286529007[71] = 0;
   out_6714729038286529007[72] = 0;
   out_6714729038286529007[73] = 0;
   out_6714729038286529007[74] = 0;
   out_6714729038286529007[75] = 0;
   out_6714729038286529007[76] = 0;
   out_6714729038286529007[77] = 0;
   out_6714729038286529007[78] = 0;
   out_6714729038286529007[79] = 0;
   out_6714729038286529007[80] = 1;
}
void h_25(double *state, double *unused, double *out_7582225353975460504) {
   out_7582225353975460504[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4929349576647638877) {
   out_4929349576647638877[0] = 0;
   out_4929349576647638877[1] = 0;
   out_4929349576647638877[2] = 0;
   out_4929349576647638877[3] = 0;
   out_4929349576647638877[4] = 0;
   out_4929349576647638877[5] = 0;
   out_4929349576647638877[6] = 1;
   out_4929349576647638877[7] = 0;
   out_4929349576647638877[8] = 0;
}
void h_24(double *state, double *unused, double *out_450546908299936874) {
   out_450546908299936874[0] = state[4];
   out_450546908299936874[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8243900214148900145) {
   out_8243900214148900145[0] = 0;
   out_8243900214148900145[1] = 0;
   out_8243900214148900145[2] = 0;
   out_8243900214148900145[3] = 0;
   out_8243900214148900145[4] = 1;
   out_8243900214148900145[5] = 0;
   out_8243900214148900145[6] = 0;
   out_8243900214148900145[7] = 0;
   out_8243900214148900145[8] = 0;
   out_8243900214148900145[9] = 0;
   out_8243900214148900145[10] = 0;
   out_8243900214148900145[11] = 0;
   out_8243900214148900145[12] = 0;
   out_8243900214148900145[13] = 0;
   out_8243900214148900145[14] = 1;
   out_8243900214148900145[15] = 0;
   out_8243900214148900145[16] = 0;
   out_8243900214148900145[17] = 0;
}
void h_30(double *state, double *unused, double *out_7425038685624331359) {
   out_7425038685624331359[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7447682535154887504) {
   out_7447682535154887504[0] = 0;
   out_7447682535154887504[1] = 0;
   out_7447682535154887504[2] = 0;
   out_7447682535154887504[3] = 0;
   out_7447682535154887504[4] = 1;
   out_7447682535154887504[5] = 0;
   out_7447682535154887504[6] = 0;
   out_7447682535154887504[7] = 0;
   out_7447682535154887504[8] = 0;
}
void h_26(double *state, double *unused, double *out_2039736216078379492) {
   out_2039736216078379492[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1187846257773582653) {
   out_1187846257773582653[0] = 0;
   out_1187846257773582653[1] = 0;
   out_1187846257773582653[2] = 0;
   out_1187846257773582653[3] = 0;
   out_1187846257773582653[4] = 0;
   out_1187846257773582653[5] = 0;
   out_1187846257773582653[6] = 0;
   out_1187846257773582653[7] = 1;
   out_1187846257773582653[8] = 0;
}
void h_27(double *state, double *unused, double *out_5862845452073245467) {
   out_5862845452073245467[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1773110065280394232) {
   out_1773110065280394232[0] = 0;
   out_1773110065280394232[1] = 0;
   out_1773110065280394232[2] = 0;
   out_1773110065280394232[3] = 1;
   out_1773110065280394232[4] = 0;
   out_1773110065280394232[5] = 0;
   out_1773110065280394232[6] = 0;
   out_1773110065280394232[7] = 0;
   out_1773110065280394232[8] = 0;
}
void h_29(double *state, double *unused, double *out_4820633364815633708) {
   out_4820633364815633708[0] = state[1];
}
void H_29(double *state, double *unused, double *out_911884590834422863) {
   out_911884590834422863[0] = 0;
   out_911884590834422863[1] = 1;
   out_911884590834422863[2] = 0;
   out_911884590834422863[3] = 0;
   out_911884590834422863[4] = 0;
   out_911884590834422863[5] = 0;
   out_911884590834422863[6] = 0;
   out_911884590834422863[7] = 0;
   out_911884590834422863[8] = 0;
}
void h_28(double *state, double *unused, double *out_939644728809704636) {
   out_939644728809704636[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2875514862399749114) {
   out_2875514862399749114[0] = 1;
   out_2875514862399749114[1] = 0;
   out_2875514862399749114[2] = 0;
   out_2875514862399749114[3] = 0;
   out_2875514862399749114[4] = 0;
   out_2875514862399749114[5] = 0;
   out_2875514862399749114[6] = 0;
   out_2875514862399749114[7] = 0;
   out_2875514862399749114[8] = 0;
}
void h_31(double *state, double *unused, double *out_2367858038470351120) {
   out_2367858038470351120[0] = state[8];
}
void H_31(double *state, double *unused, double *out_561638155540231177) {
   out_561638155540231177[0] = 0;
   out_561638155540231177[1] = 0;
   out_561638155540231177[2] = 0;
   out_561638155540231177[3] = 0;
   out_561638155540231177[4] = 0;
   out_561638155540231177[5] = 0;
   out_561638155540231177[6] = 0;
   out_561638155540231177[7] = 0;
   out_561638155540231177[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3733011136508428965) {
  err_fun(nom_x, delta_x, out_3733011136508428965);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7214253961443358426) {
  inv_err_fun(nom_x, true_x, out_7214253961443358426);
}
void car_H_mod_fun(double *state, double *out_2366253611531111470) {
  H_mod_fun(state, out_2366253611531111470);
}
void car_f_fun(double *state, double dt, double *out_3749368701677080985) {
  f_fun(state,  dt, out_3749368701677080985);
}
void car_F_fun(double *state, double dt, double *out_6714729038286529007) {
  F_fun(state,  dt, out_6714729038286529007);
}
void car_h_25(double *state, double *unused, double *out_7582225353975460504) {
  h_25(state, unused, out_7582225353975460504);
}
void car_H_25(double *state, double *unused, double *out_4929349576647638877) {
  H_25(state, unused, out_4929349576647638877);
}
void car_h_24(double *state, double *unused, double *out_450546908299936874) {
  h_24(state, unused, out_450546908299936874);
}
void car_H_24(double *state, double *unused, double *out_8243900214148900145) {
  H_24(state, unused, out_8243900214148900145);
}
void car_h_30(double *state, double *unused, double *out_7425038685624331359) {
  h_30(state, unused, out_7425038685624331359);
}
void car_H_30(double *state, double *unused, double *out_7447682535154887504) {
  H_30(state, unused, out_7447682535154887504);
}
void car_h_26(double *state, double *unused, double *out_2039736216078379492) {
  h_26(state, unused, out_2039736216078379492);
}
void car_H_26(double *state, double *unused, double *out_1187846257773582653) {
  H_26(state, unused, out_1187846257773582653);
}
void car_h_27(double *state, double *unused, double *out_5862845452073245467) {
  h_27(state, unused, out_5862845452073245467);
}
void car_H_27(double *state, double *unused, double *out_1773110065280394232) {
  H_27(state, unused, out_1773110065280394232);
}
void car_h_29(double *state, double *unused, double *out_4820633364815633708) {
  h_29(state, unused, out_4820633364815633708);
}
void car_H_29(double *state, double *unused, double *out_911884590834422863) {
  H_29(state, unused, out_911884590834422863);
}
void car_h_28(double *state, double *unused, double *out_939644728809704636) {
  h_28(state, unused, out_939644728809704636);
}
void car_H_28(double *state, double *unused, double *out_2875514862399749114) {
  H_28(state, unused, out_2875514862399749114);
}
void car_h_31(double *state, double *unused, double *out_2367858038470351120) {
  h_31(state, unused, out_2367858038470351120);
}
void car_H_31(double *state, double *unused, double *out_561638155540231177) {
  H_31(state, unused, out_561638155540231177);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
