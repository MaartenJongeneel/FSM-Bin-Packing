#include "FSMSequenceOne.h"

#include <mc_iam/devices/IO.h>
#include <string>
#include <stdexcept>


FSMSequenceOne::FSMSequenceOne(mc_rbdyn::RobotModulePtr rm, double dt, const mc_rtc::Configuration & config)
: mc_control::fsm::Controller(rm, dt, config)
{
    config_.load(config);
    // READ GRAB BOX POSES FROM CSV FILE
    std::string csvFile("/home/administrator/devel/controllers/fsm-sequence-1/src/states/BoxPoses/poses.csv");
    std::fstream Posesfile (csvFile, std::ios::in);
    Eigen::Vector6d row;
    Eigen::Vector6d IOs;
    std::string line, word;
    std::vector<std::string> states{"PrePick", "Box", "PostPick", "PrePlace",   "Place" , "PreTwist", "Twist", "PreRepick", "PreRePick", "RePick", "PostPlace" ,  "PrePlace018", "home"  , "PostRePick", "midpoint", "PrePlace017_seq3" };
    std::vector<double> precisions{   0.01   , 0.005,    0.05   ,   0.005   ,    0.005  ,    0.05  ,   0.005 ,      0.005 ,      0.005 ,   0.005 ,    0.01     ,    0.005      ,  0.005  ,    0.01     ,      0.2  ,     0.005          };
    std::vector<double> A_max{        3.0    , 4.0  ,    4.0    ,   3.0     ,    3.0    ,    3.0    ,   3.0  ,      4.0   ,      4.0   ,   3.0   ,    4.0      ,    0.8        ,    3.0  ,    3.0      ,      2.0  ,     0.8            };
    std::vector<double> V_max{        0.35   , 0.6  ,    0.6    ,   0.5     ,    0.4    ,    0.4    ,   0.4  ,      0.6   ,      0.6   ,   0.5   ,    0.6      ,    0.2        ,    0.5  ,    0.5      ,      0.5  ,     0.2            };
    std::vector<double> W_max{        2.0    , 1.0  ,    1.5    ,   2.0     ,    1.0    ,    2.0    ,   2.0  ,      2.0   ,      2.0   ,   2.0   ,    1.5      ,    1.0        ,    1.5  ,    1.5      ,      1.5  ,     1.0            };
    std::vector<double> dW_max{       2.0    , 1.5  ,    2.0    ,   2.0     ,    1.0    ,    1.5    ,   1.5  ,      1.5   ,      1.5   ,   1.5   ,    2.0      ,    1.5        ,    0.5  ,    1.5      ,      2.0  ,     1.5            };
    if(Posesfile.is_open()) 
    {
      while(getline(Posesfile, line)) 
      {
          std::stringstream str(line);
          int ii = 0;
          while(getline(str, word, ','))
          {
              //First 6 elements of the csv file define the pose (x,y,z pos and x,y,z rot vector)
              if(ii<6){
                row[ii] = std::stod(word);
              }
              //Seventh element of the csv file is the precision we use for the error tolerance
              // this is provided as a string which mentions the current phase we are (e.g. pre-pick or pick pose) 
              // from which we can infer the dynamic limits and tollerances to be used. 
              else if(ii==6)
              {
                // std::cout << "current state: " << word << std::endl; 
                int idx = -1;
                for(int idy=0;idy < states.size(); idy++)
                {
                  if (word.rfind(states[idy], 0) == 0)
                  {
                    idx = idy;
                  }
                }
                if(idx == -1){
                  // std::cout << "Value not detected... for word " << word  << std::endl;
                  throw std::invalid_argument("Value not detected... for word " + word);
                }
                if(idx == 1 || idx == 9){
                  row[2] = row[2] - 0.02; //picking and re-picking lower due to non-modeled non-holding phase
                }
                if(idx == 7 || idx == 8 || idx ==9){
                  row[0] = row[0] + 0.03; row[1] = row[1] + 0.01; //Sequence 1 re-picking at different location
                  // row[0] = row[0] + 0.05; row[1] = row[1] + 0.01; //Sequence 2 to 5 re-picking at different location
                }
                if(idx==4){
                  row[2] = row[2] - 0.005; // placing 5mm lower due to non-modeled release phase
                }
                // std::cout << "found   state: " << states[idx] << std::endl; 
                Precision.push_back(precisions[idx]);
                // set max vel max acc based on the phase
                v_max_.push_back(V_max[idx]);
                a_max_.push_back(A_max[idx]);
                w_max_.push_back(W_max[idx]);
                dw_max_.push_back(dW_max[idx]);

                // Precision.push_back(std::stod(word));
              }
              //Next 6 elements are reserved for the I/O of the UR10 (vacuum,blowoff,convcw,convccw,optitrack,reserved)
              else{
                IOs[ii-7] = std::stod(word);
              }
              ii++;
          };
          IO.push_back(IOs);
          Poses.push_back(row);            
      }
      std::cout << "length of precision: " << Precision.size() << std::endl;
      // std::cout << "length of v_max: " << v_max_.size() << std::endl;
      // std::cout << "length of a_max: " << a_max_.size() << std::endl;
      // std::cout << "length of w_max: " << w_max_.size() << std::endl;
      // std::cout << "length of dw_max: " << dw_max_.size() << std::endl;

    } else { std::cout<<"Could not open the file\n"; }

  log.addLogEntry("UR10_ft",  [this]() { return ft; });

  mc_rtc::log::success("FSMSequenceOne init done ");
}








bool FSMSequenceOne::run()
{
  // TO TEST IF VELOCITIES ARE TRACKED
  // auto vel = realRobot("box013").bodyVelW("box1_body");
  // mc_rtc::log::success(vel);

  // OPEN-LOOP CONTROL (can run on real system)
  // return mc_control::fsm::Controller::run(); // OPEN-LOOP CONTROL

  // CLOSED-LOOP(ish) CONTROL (can run on real system, watch out)
  // return mc_control::fsm::Controller::run(mc_solver::FeedbackType::ObservedRobots); // CLOSED-LOOP(ish) CONTROL

  // REAL CLOSED-LOOP CONTROL (only in sim)
  // return mc_control::fsm::Controller::run(mc_solver::FeedbackType::ClosedLoopIntegrateReal); // REAL CLOSED-LOOP CONTROL (don't do this on real setup)

  // REAL CLOSED-LOOP CONTROL (only in sim)
  // Reset qi = qr and dqi = dqr (internal state to real state) in each time step (do this only in simulation, due to noisy vel. signals of real robot)
  // robot().mbc().q = realRobot().mbc().q;
  // robot().mbc().alpha = realRobot().mbc().alpha;
  // return mc_control::fsm::Controller::run();

  // REAL CLOSED-LOOP CONTROL (can run on real system, watch out)
  // // for every 10 (or so) time steps{
  //   robot().mbc().q = realRobot().mbc().q;
  //   robot().mbc().alpha = realRobot().mbc().alpha;
  // // }
  // return mc_control::fsm::Controller::run(); //or use the closed-loop-ish control here (Jari says might be better)
  // if(StateUpdateCounter == 10){
  //   robot().mbc().q = realRobot().mbc().q;
  //   robot().mbc().alpha = realRobot().mbc().alpha;
    
  //   // std::cout << "resetting Internal counter.. \n" ;
  //   StateUpdateCounter = 0;
  // }
  // // std::cout << "counter value: " << StateUpdateCounter << std::endl ;
  // // return mc_control::fsm::Controller::run(mc_solver::FeedbackType::ObservedRobots);
  // return mc_control::fsm::Controller::run();

  // FOR SIMULATION
  t_sim = t_sim + solver().dt();
  // ft = realRobot().forceSensor("UR10FTSensor").wrenchWithoutGravity(realRobot()).vector();
  
  return mc_control::fsm::Controller::run(mc_solver::FeedbackType::ObservedRobots); // CLOSED-LOOP(ish) CONTROL

}







void FSMSequenceOne::reset(const mc_control::ControllerResetData & reset_data)
{
  mc_control::fsm::Controller::reset(reset_data);
  
  auto io = mc_iam::IO::get(robot());
  if(!io)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[URBasicTossSceneSRGripper] {} has no I.AM IO attached", robot().name());
  }
  auto gui = this->gui();
  if(!gui)
  {
    return;
  }
  // gui->addElement({"AGX"}, mc_rtc::gui::Robot("box013", [this]() -> const mc_rbdyn::Robot & { return realRobots().robot("box013"); }));
}

void FSMSequenceOne::nextPose(){

  Wo_B_init = this->realRobot().bodyPosW("control_link").translation();
  WR_B_init = this->realRobot().bodyPosW("control_link").rotation().transpose();

  if(std::abs(i) < Poses.size())
  {
    i++;
  }
  
  if(i < Poses.size()){
    Eigen::Vector3d u(Poses[i].segment(3,3)); 
    double theta = u.norm();    
    u[0] = u[0]/theta;
    u[1] = u[1]/theta; 
    u[2] = u[2]/theta; 
    WR_B_desired = Eigen::Matrix3d( Eigen::AngleAxisd(theta,u) );    
    Wo_B_desired = Eigen::Vector3d( Poses[i].head(3) );  

    p_err = Poses[i].head(3) - Wo_B_init;
    ang_err.fromRotationMatrix( WR_B_init.transpose() * WR_B_desired );

    L_p = p_err.norm();
    L_w = ang_err.angle();

    Task_err = L_p + 180*L_w/(1000*M_PI);

    // reset values
    v_  = v_max = v_max_[i]; 
    a_  = a_max = a_max_[i]; 
    w_  = w_max = w_max_[i]; 
    dw_ = dw_max = dw_max_[i]; 
  }
  
}