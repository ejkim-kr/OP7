#!/usr/bin/env python3

from cereal import car
from common.numpy_fast import interp
from selfdrive.config import Conversions as CV
from selfdrive.car.hyundai.values import CAR, Buttons, CarControllerParams
from selfdrive.car import STD_CARGO_KG, scale_rot_inertia, scale_tire_stiffness, gen_empty_fingerprint, get_safety_config
from selfdrive.car.interfaces import CarInterfaceBase
from common.params import Params
from selfdrive.controls.lib.desire_helper import LANE_CHANGE_SPEED_MIN

GearShifter = car.CarState.GearShifter
EventName = car.CarEvent.EventName
ButtonType = car.CarState.ButtonEvent.Type

class CarInterface(CarInterfaceBase):
  def __init__(self, CP, CarController, CarState):
    super().__init__(CP, CarController, CarState)
    self.cp2 = self.CS.get_can2_parser(CP)
    self.mad_mode_enabled = Params().get_bool('MadModeEnabled')

  @staticmethod
  def get_pid_accel_limits(CP, current_speed, cruise_speed):

    v_current_kph = current_speed * CV.MS_TO_KPH

    gas_max_bp = [0., 10., 20., 50., 70., 130., 150.]
    gas_max_v = [CarControllerParams.ACCEL_MAX, 1.4, 0.82, 0.65, 0.47, 0.18, 0.1]

    return CarControllerParams.ACCEL_MIN, interp(v_current_kph, gas_max_bp, gas_max_v)

  @staticmethod
  def get_params(candidate, fingerprint=gen_empty_fingerprint(), car_fw=[]):  # pylint: disable=dangerous-default-value
    ret = CarInterfaceBase.get_std_params(candidate, fingerprint)

    ret.openpilotLongitudinalControl = Params().get_bool('LongControlEnabled')

    ret.carName = "hyundai"
    ret.safetyConfigs = [get_safety_config(car.CarParams.SafetyModel.hyundaiLegacy, 0)]

    tire_stiffness_factor = 1.
    ret.maxSteeringAngleDeg = 1000.

    # lateral
    ret.lateralTuning.init('pid')
#  NIRO_EV Only:
    ret.mass = 1737. + STD_CARGO_KG
    ret.wheelbase = 2.7

    tire_stiffness_factor = 1.0 #0.95 #0.95 #반비례
    ret.steerActuatorDelay = 0.15 #0.1
    ret.steerRateCost = 1.0
    ret.steerLimitTimer = 0.3 #0.1 #0.3
    ret.steerRatio = 15.4
      
    ret.lateralTuning.pid.kf = 0.00004
    ret.lateralTuning.pid.kpBP = [0., 0.]
    ret.lateralTuning.pid.kpV = [0.05, 0.05]
    ret.lateralTuning.pid.kiBP = [0., 0.]
    ret.lateralTuning.pid.kiV = [0.015,  0.015]

    ret.radarTimeStep = 0.05
    ret.centerToFront = ret.wheelbase * 0.4

    ret.steerMaxBP = [0.]
    ret.steerMaxV = [1.5]


    # longitudinal
    ret.longitudinalTuning.kpBP = [0., 5.*CV.KPH_TO_MS, 10.*CV.KPH_TO_MS, 20.*CV.KPH_TO_MS, 130.*CV.KPH_TO_MS]
    ret.longitudinalTuning.kpV = [1.6, 1.18, 0.9, 0.78, 0.48]
    ret.longitudinalTuning.kiBP = [0., 130. * CV.KPH_TO_MS]
    ret.longitudinalTuning.kiV = [0.1, 0.06]
    ret.longitudinalActuatorDelayLowerBound = 0.3
    ret.longitudinalActuatorDelayUpperBound = 0.3

    ret.stopAccel = -2.0
    ret.stoppingDecelRate = 0.5  # brake_travel/s while trying to stop
    ret.vEgoStopping = 0.5
    ret.vEgoStarting = 0.5  # needs to be >= vEgoStopping to avoid state transition oscillation


    # TODO: get actual value, for now starting with reasonable value for
    # civic and scaling by mass and wheelbase
    ret.rotationalInertia = scale_rot_inertia(ret.mass, ret.wheelbase)

    # TODO: start from empirically derived lateral slip stiffness for the civic and scale by
    # mass and CG position, so all cars will have approximately similar dyn behaviors
    ret.tireStiffnessFront, ret.tireStiffnessRear = scale_tire_stiffness(ret.mass, ret.wheelbase, ret.centerToFront,
                                                                         tire_stiffness_factor=tire_stiffness_factor)

    # no rear steering, at least on the listed cars above
    ret.steerRatioRear = 0.
    ret.steerControlType = car.CarParams.SteerControlType.torque

    ret.stoppingControl = True

    ret.enableBsm = 0x58b in fingerprint[0]
    ret.enableAutoHold = 1151 in fingerprint[0]

    # ignore CAN2 address if L-CAN on the same BUS
    ret.mdpsBus = 1 if 593 in fingerprint[1] and 1296 not in fingerprint[1] else 0
    ret.sasBus = 1 if 688 in fingerprint[1] and 1296 not in fingerprint[1] else 0
    ret.sccBus = 0 if 1056 in fingerprint[0] else 1 if 1056 in fingerprint[1] and 1296 not in fingerprint[1] \
                                                                     else 2 if 1056 in fingerprint[2] else -1

    if ret.sccBus >= 0:
      ret.hasScc13 = 1290 in fingerprint[ret.sccBus]
      ret.hasScc14 = 905 in fingerprint[ret.sccBus]

    ret.hasEms = 608 in fingerprint[0] and 809 in fingerprint[0]
    ret.hasLfaHda = 1157 in fingerprint[0]

    ret.radarOffCan = ret.sccBus == -1
    ret.pcmCruise = not ret.radarOffCan

    # set safety_hyundai_community only for non-SCC, MDPS harrness or SCC harrness cars or cars that have unknown issue
    if ret.radarOffCan or ret.mdpsBus == 1 or ret.openpilotLongitudinalControl or ret.sccBus == 1 or Params().get_bool('MadModeEnabled'):
      ret.safetyConfigs = [get_safety_config(car.CarParams.SafetyModel.hyundaiCommunity, 0)]
    return ret

  def update(self, c, can_strings):
    self.cp.update_strings(can_strings)
    self.cp2.update_strings(can_strings)
    self.cp_cam.update_strings(can_strings)

    ret = self.CS.update(self.cp, self.cp2, self.cp_cam)
    ret.canValid = self.cp.can_valid and self.cp2.can_valid and self.cp_cam.can_valid

    if self.CP.pcmCruise and not self.CC.scc_live:
      self.CP.pcmCruise = False
    elif self.CC.scc_live and not self.CP.pcmCruise:
      self.CP.pcmCruise = True

    # most HKG cars has no long control, it is safer and easier to engage by main on

    if self.mad_mode_enabled:
      ret.cruiseState.enabled = ret.cruiseState.available

    # turning indicator alert logic
    if not self.CC.keep_steering_turn_signals and (ret.leftBlinker or ret.rightBlinker or self.CC.turning_signal_timer) and ret.vEgo < LANE_CHANGE_SPEED_MIN - 1.2:
      self.CC.turning_indicator_alert = True
    else:
      self.CC.turning_indicator_alert = False

    # low speed steer alert hysteresis logic (only for cars with steer cut off above 10 m/s)
    if ret.vEgo < (self.CP.minSteerSpeed + 0.2) and self.CP.minSteerSpeed > 10.:
      self.low_speed_alert = True
    if ret.vEgo > (self.CP.minSteerSpeed + 0.7):
      self.low_speed_alert = False

    buttonEvents = []
    if self.CS.cruise_buttons != self.CS.prev_cruise_buttons:
      be = car.CarState.ButtonEvent.new_message()
      be.pressed = self.CS.cruise_buttons != 0
      but = self.CS.cruise_buttons if be.pressed else self.CS.prev_cruise_buttons
      if but == Buttons.RES_ACCEL:
        be.type = ButtonType.accelCruise
      elif but == Buttons.SET_DECEL:
        be.type = ButtonType.decelCruise
      elif but == Buttons.GAP_DIST:
        be.type = ButtonType.gapAdjustCruise
      #elif but == Buttons.CANCEL:
      #  be.type = ButtonType.cancel
      else:
        be.type = ButtonType.unknown
      buttonEvents.append(be)
    if self.CS.cruise_main_button != self.CS.prev_cruise_main_button:
      be = car.CarState.ButtonEvent.new_message()
      be.type = ButtonType.altButton3
      be.pressed = bool(self.CS.cruise_main_button)
      buttonEvents.append(be)
    ret.buttonEvents = buttonEvents

    events = self.create_common_events(ret)

    if self.CC.longcontrol and self.CS.cruise_unavail:
      events.add(EventName.brakeUnavailable)
    #if abs(ret.steeringAngleDeg) > 90. and EventName.steerTempUnavailable not in events.events:
    #  events.add(EventName.steerTempUnavailable)
    if self.low_speed_alert and not self.CS.mdps_bus:
      events.add(EventName.belowSteerSpeed)
    if self.CC.turning_indicator_alert:
      events.add(EventName.turningIndicatorOn)
    if self.mad_mode_enabled and EventName.pedalPressed in events.events:
      events.events.remove(EventName.pedalPressed)

  # handle button presses
    for b in ret.buttonEvents:
      # do disable on button down
      if b.type == ButtonType.cancel and b.pressed:
        events.add(EventName.buttonCancel)
      if self.CC.longcontrol and not self.CC.scc_live:
        # do enable on both accel and decel buttons
        if b.type in [ButtonType.accelCruise, ButtonType.decelCruise] and not b.pressed:
          events.add(EventName.buttonEnable)
        if EventName.wrongCarMode in events.events:
          events.events.remove(EventName.wrongCarMode)
        if EventName.pcmDisable in events.events:
          events.events.remove(EventName.pcmDisable)
      elif not self.CC.longcontrol and ret.cruiseState.enabled:
        # do enable on decel button only
        if b.type == ButtonType.decelCruise and not b.pressed:
          events.add(EventName.buttonEnable)

    # scc smoother
    if self.CC.scc_smoother is not None:
      self.CC.scc_smoother.inject_events(events)

    ret.events = events.to_msg()

    self.CS.out = ret.as_reader()
    return self.CS.out

  # scc smoother - hyundai only
  def apply(self, c, controls):
    ret = self.CC.update(c, c.enabled, self.CS, self.frame, c, c.actuators,
                               c.cruiseControl.cancel, c.hudControl.visualAlert, c.hudControl.leftLaneVisible,
                               c.hudControl.rightLaneVisible, c.hudControl.leftLaneDepart, c.hudControl.rightLaneDepart,
                               c.hudControl.setSpeed, c.hudControl.leadVisible, controls)
    self.frame += 1
    return ret