#include "Operation.h"
#include "../FSMSequenceOne.h"

#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<unistd.h>

   

void Operation::configure(const mc_rtc::Configuration & config)
{
}

void Operation::start(mc_control::fsm::Controller & ctl_)
{
    auto & ctl = static_cast<FSMSequenceOne &>(ctl_);   
    real = int(ctl.config_("real"));
    auto io = mc_iam::IO::get(ctl.robot("ur10_FTSensor_sr_gripper"));
    // io->stop_blowoff();
    io->stop_vacuum();
    // io->start_optitrack();

    // initialize the time variables
    ctl.t_final = std::chrono::system_clock::now();
    
    // ctl.efTask_ = std::make_shared<mc_tasks::EndEffectorTask>("control_link", ctl.robots(), 0, 15, 1e6); // Seq 2
    ctl.efTask_ = std::make_shared<mc_tasks::EndEffectorTask>("control_link", ctl.robots(), 0, 10, 1e6); // Seq 1 3 4 5 
    ctl.efTask_->positionTask->weight(5e6);
    ctl.efTask_->orientationTask->weight(5e6);
    ctl.solver().addTask(ctl.efTask_);

    Eigen::Vector3d Wo_B_real(ctl.realRobot().bodyPosW("control_link").translation());
    Eigen::Matrix3d WR_B_real(ctl.realRobot().bodyPosW("control_link").rotation());
    sva::PTransformd WH_B(WR_B_real.transpose(), Wo_B_real);
    ctl.efTask_->set_ef_pose(WH_B);
    ctl.Task_err = 0;

    // ctl.a_ = ctl.a_max;
    // ctl.v_ = ctl.v_max;
    // ctl.w_ = ctl.w_max; 
    // ctl.dw_ = ctl.dw_max;
}
bool Operation::run(mc_control::fsm::Controller & ctl_) // Called every iteration until it returns true
{
    auto & ctl = static_cast<FSMSequenceOne &>(ctl_);    
    auto io = mc_iam::IO::get(ctl.realRobot("ur10_FTSensor_sr_gripper"));
    // ctl.StateUpdateCounter++;
    
    if(ctl.initdone && ctl.Task_err > ctl.Precision[ctl.i]  && !ctl.Sleep){
        // here we create the parameterization of the time trajectory

        // compute the current time for the trajectory
        
        double t_{0.0};
        if(real == 1){
            ctl.t_elapsed =  std::chrono::system_clock::now() - ctl.t_start;
            t_ = ctl.t_elapsed.count();
        }else{
            t_ = (ctl.t_sim - t_start_sim);
        }
        
        double sig, dsig, ddsig, sig_w, dsig_w, ddsig_w;
        
        if(t_ >= 0 && t_ <= ctl.t_b ){
            sig = 0.5*ctl.a_*pow(t_,2);
            dsig = ctl.a_*t_;
            ddsig = ctl.a_;
            ctl.last_linVel = dsig;
        }
        else{
            if( t_ > ctl.t_b && t_ <= ctl.t_ff - ctl.t_b  ){
                sig = ctl.last_linVel*t_ + 0.5 * (ctl.L_p - ctl.last_linVel*ctl.t_ff);
                dsig = ctl.last_linVel; 
                ddsig = 0;
            }
            else
            {
                if(t_ > ctl.t_ff - ctl.t_b && t_ <= ctl.t_ff  ){
                    sig = ctl.L_p - 0.5*ctl.a_*pow(ctl.t_ff - t_,2);
                    dsig = ctl.last_linVel -ctl.a_*(t_ - (ctl.t_ff - ctl.t_b) );
                    ddsig =  -ctl.a_;
                }
                else{
                    sig = ctl.L_p;
                    dsig = 0;
                    ddsig = 0;
                }
            }
        }

        if(t_ >= 0 && t_ <= ctl.t_b_w ){ 
            sig_w = 0.5*ctl.dw_*pow(t_,2);
            dsig_w = ctl.dw_*t_;
            ddsig_w = ctl.dw_;
            ctl.last_angVel = dsig_w;
        }
        else{
            if( t_ > ctl.t_b_w && t_ <= ctl.t_ff - ctl.t_b_w  ){
                sig_w = ctl.last_angVel*t_ + 0.5 * (ctl.L_w - ctl.last_angVel*ctl.t_ff);
                dsig_w = ctl.last_angVel; 
                ddsig_w = 0;
            }
            else{
                if(t_ > ctl.t_ff - ctl.t_b_w && t_ <= ctl.t_ff  ){
                    sig_w = ctl.L_w - 0.5*ctl.dw_*pow(ctl.t_ff - t_,2);
                    dsig_w = ctl.last_angVel -ctl.dw_*(t_ - (ctl.t_ff - ctl.t_b_w) );
                    ddsig_w =  -ctl.dw_;
                }
                else{
                    sig_w = ctl.L_w;
                    dsig_w = 0;
                    ddsig_w = 0;
                } 
            }
        }

        Eigen::Vector3d Wo_B_cmd; Eigen::Matrix3d WR_B_cmd;

        Wo_B_cmd = ctl.Wo_B_init + (sig/ctl.L_p)*(ctl.Poses[ctl.i].head(3) - ctl.Wo_B_init);
        WR_B_cmd = ctl.WR_B_init*Eigen::Matrix3d( Eigen::AngleAxisd((sig_w/ctl.L_w)*ctl.ang_err.angle(),ctl.ang_err.axis()) );

        // //--------------------------------------------------------------------------------------
        // std::cout << "Current Joint config: [ ";
        // for(auto e_val: ctl.realRobot("ur10_FTSensor_sr_gripper").encoderValues()) {
        //     std::cout << e_val << " " ;
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "Desired pos: " << ctl.Poses[ctl.i].head(3).transpose() << std::endl;
        // std::cout << "Current pos: " << ctl.realRobot().bodyPosW("control_link").translation().transpose() << std::endl;
        // std::cout << "Desired orien: \n" << ctl.WR_B_init  << std::endl;
        // std::cout << "Current orien: " << ctl.realRobot().bodyPosW("control_link").rotation() << std::endl;
        //--------------------------------------------------------------------------------------

        ctl.efTask_->set_ef_pose(sva::PTransformd(WR_B_cmd.transpose(), Wo_B_cmd));
        // ctl.efTask_->positionTask->refVel( (dsig/ctl.L_p)*(ctl.Poses[ctl.i].head(3) - ctl.Wo_B_init) ); 
        // ctl.efTask_->positionTask->refAccel( (ddsig/ctl.L_p)*(ctl.Poses[ctl.i].head(3) - ctl.Wo_B_init) );
        // ctl.efTask_->orientationTask->refVel(ctl.WR_B_init *(dsig_w*ctl.ang_err.axis())); 
        // ctl.efTask_->orientationTask->refAccel((ddsig_w/ctl.L_w)*ctl.ang_err.axis()*ctl.ang_err.angle());


        // compute the task error
        Eigen::AngleAxisd ang; 
        ang.fromRotationMatrix( ctl.realRobot().bodyPosW("control_link").rotation() * ctl.WR_B_desired );
 
        ctl.Task_err = (ctl.Poses[ctl.i].head(3) - ctl.realRobot().bodyPosW("control_link").translation()).norm() + 180*ang.angle()/ (500*M_PI) ;
        // std::cout << "Precision : " << ctl.Precision[ctl.i] << " , task error: " << ctl.Task_err << std::endl;
        // std::cout << "In command phase\n";
            
    }else{ // once converged to the desired pose, generate the trajectory parameters to go to the next one

        if(!ctl.initdone){
            ctl.initdone = true;
            mc_rtc::log::success("Initialization done!");
            ctl.nextPose();

            double t_f, t_f_w;
                
            // check which trajectory parameterization will last longer: translation or rotation
            if( ctl.L_p > std::pow(ctl.v_max_[ctl.i],2)/ctl.a_max_[ctl.i] ){
                ctl.t_b = ctl.v_max_[ctl.i]/ctl.a_max_[ctl.i];
                t_f = (ctl.L_p - ctl.a_max_[ctl.i]*std::pow(ctl.t_b,2))/ctl.v_max_[ctl.i] + 2*ctl.t_b;
            }else{
                ctl.t_b = std::sqrt(ctl.L_p/ctl.a_max_[ctl.i]);
                t_f = 2*ctl.t_b;
            }

            if( ctl.L_w > std::pow(ctl.w_max_[ctl.i],2)/ctl.dw_max_[ctl.i] ){
                ctl.t_b_w = ctl.w_max_[ctl.i]/ctl.dw_max_[ctl.i];
                t_f_w = (ctl.L_w - ctl.dw_max_[ctl.i]*std::pow(ctl.t_b_w,2))/ctl.w_max_[ctl.i] + 2*ctl.t_b_w;
            }else{
                ctl.t_b_w = std::sqrt(ctl.L_w/ctl.dw_max_[ctl.i]);
                t_f_w = 2*ctl.t_b_w;
            }

            if( t_f >= t_f_w){ // rotational part is faster, recompute it at slower pace
                ctl.t_ff = t_f;

                // can I do it in 2 phases?
                if( ctl.w_max_[ctl.i] >= 2*ctl.L_w/t_f && ctl.dw_max_[ctl.i] >= 4*ctl.L_w/std::pow(t_f,2) ){
                    ctl.dw_ = 4*ctl.L_w/std::pow(t_f,2);
                    ctl.w_  = 2*ctl.L_w/t_f;
                    ctl.t_b_w = ctl.t_ff/2;
                }else{ // otherwise, made it in 3 phases at max vel
                    // trajectory exceed w_max constraints, so must be done in 3 phases
                    ctl.w_ = ctl.w_max_[ctl.i];
                    ctl.t_b_w = ctl.t_ff - ctl.L_w/ctl.w_max_[ctl.i];
                    ctl.dw_ = ctl.L_w/(ctl.t_b_w*ctl.t_ff - std::pow(ctl.t_b_w,2));
                    
                }

            }else{             // transtional part is faster, recompute it at slower pace
                ctl.t_ff = t_f_w;

                // can I do it in 2 phases?
                if( ctl.v_max_[ctl.i] >= 2*ctl.L_p/ctl.t_ff && ctl.a_max_[ctl.i] >= 4*ctl.L_p/std::pow(ctl.t_ff,2) ){
                    ctl.a_ = 4*ctl.L_p/std::pow(ctl.t_ff,2);
                    ctl.v_  = 2*ctl.L_p/ctl.t_ff;
                    ctl.t_b = ctl.t_ff/2;
                }else{ // otherwise, made it in 3 phases at max vel
                    // trajectory exceed w_max constraints, so must be done in 3 phases
                    ctl.v_ = ctl.v_max_[ctl.i];
                    ctl.t_b = ctl.t_ff - ctl.L_p/ctl.v_max_[ctl.i];
                    ctl.a_ = ctl.L_p/(ctl.t_b*ctl.t_ff - std::pow(ctl.t_b,2));
                }
            }
        }
        else{
            if(!ctl.Sleep){
                mc_rtc::log::success("Completed task " + std::to_string(ctl.i + 1));
                // go to the next pose
                ctl.nextPose();
                if (ctl.IO[ctl.i][0] == 1){ 
                    ctl.current_time = std::chrono::system_clock::now();
                    current_time_sim = ctl.t_sim;
                    if(ctl.i >= 1 && ctl.IO[ctl.i-1][0] != ctl.IO[ctl.i][0]){
                        ctl.Sleep = true;
                    }
                    
                    io->start_vacuum();
                    std::cout << "Vacuum ON \n";
                   
                }else{
                    io->stop_vacuum();
                    ctl.current_time = std::chrono::system_clock::now();
                    current_time_sim = ctl.t_sim;
                    if(ctl.i >= 1 && ctl.IO[ctl.i-1][0] != ctl.IO[ctl.i][0]){
                        ctl.Sleep = true;
                    }
                    std::cout << "Vacuum OFF \n";
                }

                if (ctl.IO[ctl.i][1] == 1){
                    ctl.current_time = std::chrono::system_clock::now();
                    current_time_sim = ctl.t_sim;
                    if(ctl.i >= 1 && ctl.IO[ctl.i-1][1] != ctl.IO[ctl.i][1]){
                        ctl.Sleep = true;
                    }
                    // io->start_blowoff();
                    std::cout << "BlowOff ON \n";
                   
                }else{
                    // io->stop_blowoff();
                    std::cout << "BlowOff OFF \n";
                }

                double t_f, t_f_w;
                
                // check which trajectory parameterization will last longer: translation or rotation
                if( ctl.L_p > std::pow(ctl.v_max_[ctl.i],2)/ctl.a_max_[ctl.i] ){
                    ctl.t_b = ctl.v_max_[ctl.i]/ctl.a_max_[ctl.i];
                    t_f = (ctl.L_p - ctl.a_max_[ctl.i]*std::pow(ctl.t_b,2))/ctl.v_max_[ctl.i] + 2*ctl.t_b;
                }else{
                    ctl.t_b = std::sqrt(ctl.L_p/ctl.a_max_[ctl.i]);
                    t_f = 2*ctl.t_b;
                }

                if( ctl.L_w > std::pow(ctl.w_max_[ctl.i],2)/ctl.dw_max_[ctl.i] ){
                    ctl.t_b_w = ctl.w_max_[ctl.i]/ctl.dw_max_[ctl.i];
                    t_f_w = (ctl.L_w - ctl.dw_max_[ctl.i]*std::pow(ctl.t_b_w,2))/ctl.w_max_[ctl.i] + 2*ctl.t_b_w;
                }else{
                    ctl.t_b_w = std::sqrt(ctl.L_w/ctl.dw_max_[ctl.i]);
                    t_f_w = 2*ctl.t_b_w;
                }

                if( t_f >= t_f_w){ // rotational part is faster, recompute it at slower pace
                    ctl.t_ff = t_f;

                    // can I do it in 2 phases?
                    if( ctl.w_max_[ctl.i] >= 2*ctl.L_w/t_f && ctl.dw_max_[ctl.i] >= 4*ctl.L_w/std::pow(t_f,2) ){
                        ctl.dw_ = 4*ctl.L_w/std::pow(t_f,2);
                        ctl.w_  = 2*ctl.L_w/t_f;
                        ctl.t_b_w = ctl.t_ff/2;
                    }else{ // otherwise, made it in 3 phases at max vel
                        // trajectory exceed w_max constraints, so must be done in 3 phases
                        ctl.w_ = ctl.w_max_[ctl.i];
                        ctl.t_b_w = ctl.t_ff - ctl.L_w/ctl.w_max_[ctl.i];
                        ctl.dw_ = ctl.L_w/(ctl.t_b_w*ctl.t_ff - std::pow(ctl.t_b_w,2));
                        
                    }

                }else{             // transtional part is faster, recompute it at slower pace
                    ctl.t_ff = t_f_w;

                    // can I do it in 2 phases?
                    if( ctl.v_max_[ctl.i] >= 2*ctl.L_p/ctl.t_ff && ctl.a_max_[ctl.i] >= 4*ctl.L_p/std::pow(ctl.t_ff,2) ){
                        ctl.a_ = 4*ctl.L_p/std::pow(ctl.t_ff,2);
                        ctl.v_  = 2*ctl.L_p/ctl.t_ff;
                        ctl.t_b = ctl.t_ff/2;
                    }else{ // otherwise, made it in 3 phases at max vel
                        // trajectory exceed w_max constraints, so must be done in 3 phases
                        ctl.v_ = ctl.v_max_[ctl.i];
                        ctl.t_b = ctl.t_ff - ctl.L_p/ctl.v_max_[ctl.i];
                        ctl.a_ = ctl.L_p/(ctl.t_b*ctl.t_ff - std::pow(ctl.t_b,2));
                    }
                }
                                
                // Check if we're done running through states
                if(ctl.i > ctl.Poses.size() -1 ){
                    ctl.efTask_->set_ef_pose(sva::PTransformd(ctl.WR_B_init.transpose(), ctl.Wo_B_init));
                    // ctl.efTask_->positionTask->refVel( Eigen::Vector3d::Zero() );
                    // ctl.efTask_->orientationTask->refVel( Eigen::Vector3d::Zero() );
                    output("OK");
                    return true;
                }
                
            }
            else{// you enter here for sleeping
                // std::cout << "idle position: " << Eigen::Vector3d( ctl.Poses[ctl.i-1].head(3) ).transpose() << std::endl;
                // std::cout << "idle position: " << ctl.Wo_B_init.transpose() << std::endl;
                // std::cout << "Current pos  :" << ctl.realRobot().bodyPosW("control_link").translation().transpose() << std::endl;
                ctl.efTask_->set_ef_pose(sva::PTransformd(ctl.WR_B_init.transpose(), ctl.Wo_B_init));
                // ctl.efTask_->positionTask->refVel( Eigen::Vector3d::Zero() );
                // ctl.efTask_->orientationTask->refVel( Eigen::Vector3d::Zero() );
                // std::cout << "I am Sleeping in Vacuum ON phase\n";
                if(real == 1 && ((double)(std::chrono::system_clock::now() - ctl.current_time).count())/1000000000 >= 0.5){
                        ctl.Sleep = false;
                        std::cout << "Stop Sleeping in Vacuum ON phase\n";
                }
                if(real == 0 && ((ctl.t_sim - current_time_sim) >= 0.5)){
                        ctl.Sleep = false;
                        std::cout << "Stop Sleeping in Vacuum ON phase\n";
                }
            }
        }// end initDone     
        
        // initialize the clock for the next trajectory
        if(real == 1){
            ctl.t_start = std::chrono::system_clock::now();
        }
        if(real == 0){
            t_start_sim = ctl.t_sim;
        }     
    }

         
    return false; 
}

void Operation::teardown(mc_control::fsm::Controller & ctl_)
{
    auto & ctl = static_cast<FSMSequenceOne &>(ctl_);
    ctl_.solver().removeTask(ctl.efTask_);
    auto io = mc_iam::IO::get(ctl.realRobot("ur10_FTSensor_sr_gripper"));

    // io->stop_optitrack();
    // io->stop_blowoff();
    io->stop_vacuum();
    
    output("OK");
    exit(EXIT_SUCCESS);
}

EXPORT_SINGLE_STATE("Operation", Operation)
