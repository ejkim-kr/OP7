#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3733011136508428965);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7214253961443358426);
void car_H_mod_fun(double *state, double *out_2366253611531111470);
void car_f_fun(double *state, double dt, double *out_3749368701677080985);
void car_F_fun(double *state, double dt, double *out_6714729038286529007);
void car_h_25(double *state, double *unused, double *out_7582225353975460504);
void car_H_25(double *state, double *unused, double *out_4929349576647638877);
void car_h_24(double *state, double *unused, double *out_450546908299936874);
void car_H_24(double *state, double *unused, double *out_8243900214148900145);
void car_h_30(double *state, double *unused, double *out_7425038685624331359);
void car_H_30(double *state, double *unused, double *out_7447682535154887504);
void car_h_26(double *state, double *unused, double *out_2039736216078379492);
void car_H_26(double *state, double *unused, double *out_1187846257773582653);
void car_h_27(double *state, double *unused, double *out_5862845452073245467);
void car_H_27(double *state, double *unused, double *out_1773110065280394232);
void car_h_29(double *state, double *unused, double *out_4820633364815633708);
void car_H_29(double *state, double *unused, double *out_911884590834422863);
void car_h_28(double *state, double *unused, double *out_939644728809704636);
void car_H_28(double *state, double *unused, double *out_2875514862399749114);
void car_h_31(double *state, double *unused, double *out_2367858038470351120);
void car_H_31(double *state, double *unused, double *out_561638155540231177);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}