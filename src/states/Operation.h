#pragma once

#include <mc_control/fsm/State.h>
#include <mc_tasks/SurfaceTransformTask.h>
#include <mc_tasks/EndEffectorTask.h>
#include <mc_iam/devices/IO.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<unistd.h>
#include<cstdlib>

#include <math.h>
#include <chrono>

struct Operation : mc_control::fsm::State
{
    mc_rtc::Configuration state_conf_;

    void configure(const mc_rtc::Configuration & config) override;
    void start(mc_control::fsm::Controller & ctl) override;
    bool run(mc_control::fsm::Controller & ctl) override;
    void teardown(mc_control::fsm::Controller & ctl) override;   
    
private:
    std::shared_ptr<mc_tasks::SurfaceTransformTask> surfTransTask_;
    std::string robot_;
    double t_start_sim{0.0};
    double current_time_sim{0.0};
    int real{0};
};
