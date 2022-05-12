#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_7883823041686261173);
void live_err_fun(double *nom_x, double *delta_x, double *out_4326014548724243586);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8804837556079097701);
void live_H_mod_fun(double *state, double *out_5951945452246248407);
void live_f_fun(double *state, double dt, double *out_7045241452713039968);
void live_F_fun(double *state, double dt, double *out_4744992820336334836);
void live_h_4(double *state, double *unused, double *out_4063494398576480882);
void live_H_4(double *state, double *unused, double *out_3574752001031558963);
void live_h_9(double *state, double *unused, double *out_1440340947536770573);
void live_H_9(double *state, double *unused, double *out_3712466934232888507);
void live_h_10(double *state, double *unused, double *out_4493613389207261310);
void live_H_10(double *state, double *unused, double *out_323923986299902565);
void live_h_12(double *state, double *unused, double *out_947592367452026341);
void live_H_12(double *state, double *unused, double *out_1444704407000402832);
void live_h_31(double *state, double *unused, double *out_7272778378623806227);
void live_H_31(double *state, double *unused, double *out_208089943658951587);
void live_h_32(double *state, double *unused, double *out_517681623910710321);
void live_H_32(double *state, double *unused, double *out_7738111779312575270);
void live_h_13(double *state, double *unused, double *out_1893351072680475);
void live_H_13(double *state, double *unused, double *out_2286584468871788719);
void live_h_14(double *state, double *unused, double *out_1440340947536770573);
void live_H_14(double *state, double *unused, double *out_3712466934232888507);
void live_h_33(double *state, double *unused, double *out_7428219120733597680);
void live_H_33(double *state, double *unused, double *out_2942467060979906017);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}