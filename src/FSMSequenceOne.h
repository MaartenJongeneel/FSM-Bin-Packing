#pragma once

#include <mc_control/fsm/Controller.h>

#include "api.h"

struct FSMSequenceOne_DLLAPI FSMSequenceOne : public mc_control::fsm::Controller
{
  FSMSequenceOne(mc_rbdyn::RobotModulePtr rm, double dt, const mc_rtc::Configuration & config);

  bool run() override;

  void reset(const mc_control::ControllerResetData & reset_data) override;

  void nextPose();

  mc_rtc::Configuration config_;
  
  std::shared_ptr<mc_tasks::EndEffectorTask> efTask_;
  std::shared_ptr<mc_solver::DynamicsConstraint> boxDynamics;
  std::vector<Eigen::Vector6d> Poses;
  std::vector<Eigen::Vector6d> IO;
  std::vector<double> Precision;
  std::vector<double> v_max_;
  std::vector<double> a_max_;
  std::vector<double> w_max_;
  std::vector<double> dw_max_;

  int i = -1;
  int StateUpdateCounter = 0;
  bool Sleep{false};
  bool isVacuum{false};
  bool isBlowOff{false};

  bool initdone{false};
  Eigen::Vector3d p_err;
  Eigen::AngleAxisd ang_err;
  double Task_err = 0;

  double tolerance;
  double a_max = 1.5;     //  maximum linear acceleration                               [m/s^2]
  double v_max = 0.3;     //  maximum linear velocity                                   [m/s]
  double a_ = 0.0;        //  current linear acceleration                               [m/s^2]
  double v_ = 0.0;        //  current linear velocity                                   [m/s]  
  double last_linVel = 0; //  linear velocity reached at the end of acceleration phase  [m/s] 
  double t_b;             //  blending time                                             [s]
  double L_p;             //  norm of the linear path to be travelled      

  double dw_max = 1.5;    //  maximum angular acceleration                              [rad/s^2]
  double w_max = 1.5;     //  maximum angular velocity                                  [rad/s]
  double dw_ = 0.0;       //  current angular acceleration                              [rad/s^2]
  double w_ = 0.0;        //  current angular velocity                                  [rad/s]
  double t_b_w;           //  blending time                                             [s]
  double last_angVel = 0; //  angular velocity reached at the end of acceleration phase [m/s] 
  double L_w;             //  norm of the angular path to be travelled 

  double t_ff;            //  total duration of the motion to be generated that satisfies dynamic constraints

  Eigen::Matrix3d WR_B_desired;    
  Eigen::Vector3d Wo_B_desired;

  // Using time point and system_clock
  std::chrono::time_point<std::chrono::system_clock> t_start;
  std::chrono::time_point<std::chrono::system_clock> t_final;
  std::chrono::time_point<std::chrono::system_clock> current_time;
  std::chrono::duration<double> t_elapsed ;

  double t_sim{0.0};

  Eigen::Vector3d Wo_B_init;
  Eigen::Matrix3d WR_B_init;

  std::vector<sva::PTransformd> interim_poses;

  //For the FT sensor
  Eigen::Vector6d ft;

  // logger
  mc_rtc::Logger log{mc_rtc::Logger::Policy::NON_THREADED, "", ""}; // these dummy parameters are overwritten in .cpp

private:
  std::vector<double> jointTorquesAGX_;
  Eigen::VectorXd jointTorquesAGX;
  Eigen::VectorXd jointTorquesMCRTC;

  
};
