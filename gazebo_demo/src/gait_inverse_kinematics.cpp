#include <ros/ros.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDynUrdf/Reader.h>
#include <geometry_msgs/PoseArray.h>
#include <std_msgs/Float64.h>
#include <sensor_msgs/JointState.h>
#include "task_define.h"

rbdyn_urdf::Urdf robot_data;
MultiTaskPtr tasks;
InverseMethodPtr method;
std::string ee_names[4] = {"link1_4", "link2_4", "link3_4", "link4_4"};
//std::string ee_names[4] = {"l_wrist", "r_wrist", "r_ankle", "l_ankle"};
ros::Subscriber sub;
ros::Publisher pub;
ros::Publisher command_pub[16];

void JointStateFromMBC(rbd::MultiBody mb, rbd::MultiBodyConfig mbc, sensor_msgs::JointState& msg)
{
    int count = 0;
    for (auto itr = mbc.q.begin(); itr != mbc.q.end(); ++itr) {
        if (mb.joint(count).type() == rbd::Joint::Type::Rev ||
            mb.joint(count).type() == rbd::Joint::Type::Prism) {// 1dof joint
            msg.name.push_back(mb.joint(count).name());
            msg.position.push_back(*(itr->begin()));
        } 
        count++;
    }
    msg.header.stamp = ros::Time::now();
}


void gaitCallback(const geometry_msgs::PoseArray::ConstPtr& msg, rbd::MultiBody mb, rbd::MultiBodyConfig &mbc)
{
  // geometry_msgs/PoseArray has data in the order of legs. 
  for (int i = 0; i < 4; i++)
  {
    auto pos = msg->poses[i].position;
    auto X_O_Ti = sva::PTransformd(Vector3d(pos.x, pos.y, pos.z));
    boost::static_pointer_cast<BodyTask>(tasks[i].second)->_X_O_T = X_O_Ti;
  }

  TaskMin(mb, mbc, tasks, method, 1.0, 200, 1e-6);
  for (auto&& var : mbc.q)
  {
    for (auto&& var2 : var)
    {
      while(var2 < -M_PI) var2+=2*M_PI;
      while(+M_PI < var2) var2-=2*M_PI;
    }
  }
  sensor_msgs::JointState js_msg;
  JointStateFromMBC(mb, mbc, js_msg);
  pub.publish(js_msg);

  for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
          std::string name = "joint" + std::to_string(i+1) + "_" + std::to_string(j) + std::to_string(j+1);
          std_msgs::Float64 v;
          v.data = mbc.q[mb.jointIndexByName(name)][0];
          command_pub[i*4+ j].publish(v);
      }
  }
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "gait_ik");
  ros::NodeHandle nh;
  std::string urdf_string;
  if (nh.getParam("/robot_description", urdf_string))
  {
    robot_data = rbdyn_urdf::readUrdf(urdf_string);
  }
  else ROS_ERROR("could't get param");
  rbd::MultiBody mb = robot_data.mbg.makeMultiBody("base_link", rbd::Joint::Fixed);
  rbd::MultiBodyConfig mbc(mb);
  mbc.zero(mb);
  rbd::forwardKinematics(mb, mbc);
  rbd::forwardVelocity(mb, mbc);
  TaskPtr posturetask(new PostureTask(mb, mbc));
  sva::PTransformd X_O_T = sva::PTransformd(Eigen::Vector3d(0,0,0)); //Temporal value
  sva::PTransformd x_b_p = sva::PTransformd(Eigen::Vector3d(0.140,0,0)); //To end-effector point
  TaskPtr bodytask1(new BodyTask(mb, ee_names[0], X_O_T, x_b_p, "position", "bodytask1"));
  TaskPtr bodytask2(new BodyTask(mb, ee_names[1], X_O_T, x_b_p, "position", "bodytask2"));
  TaskPtr bodytask3(new BodyTask(mb, ee_names[2], X_O_T, x_b_p, "position", "bodytask3"));
  TaskPtr bodytask4(new BodyTask(mb, ee_names[3], X_O_T, x_b_p, "position", "bodytask4"));
  tasks.push_back(std::pair<double, TaskPtr>(100, bodytask1));
  tasks.push_back(std::pair<double, TaskPtr>(100, bodytask2));
  tasks.push_back(std::pair<double, TaskPtr>(100, bodytask3));
  tasks.push_back(std::pair<double, TaskPtr>(100, bodytask4));
  //tasks.push_back(std::pair<double, TaskPtr>(1.0, posturetask));

  int row = 0;
  int col = mb.nrDof();
  for (auto&& var : tasks) {
      auto j = var.second->J(mb, mbc);
      row+=j.rows();
  }
  MatrixXd Wl = MatrixXd::Identity(col, col);
  MatrixXd We = MatrixXd::Identity(row, row)*100;
  //Set limits for LMInvConsideredSolvalityWithLimit
  VectorXd HiLimit = VectorXd::Ones(col) * (+150.0) * M_PI/180.0;
  VectorXd LoLimit = VectorXd::Ones(col) * (-150.0) * M_PI/180.0;
  rbd::MultiBodyConfig mbchi(mb), mbclo(mb);
  mbchi.q = rbd::vectorToDof(mb, HiLimit);
  mbclo.q = rbd::vectorToDof(mb, LoLimit);
  //Each value setting
  mbchi.q[mb.jointIndexByName("joint1_01")][0] = (+30.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint2_01")][0] = (-30.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint3_01")][0] = (-30.0) * M_PI/180.0;
  mbchi.q[mb.jointIndexByName("joint4_01")][0] = (+30.0) * M_PI/180.0;

  mbchi.q[mb.jointIndexByName("joint1_23")][0] = (+10.0) * M_PI/180.0;
  mbchi.q[mb.jointIndexByName("joint2_23")][0] = (+10.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint3_23")][0] = (-10.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint4_23")][0] = (-10.0) * M_PI/180.0;

  mbchi.q[mb.jointIndexByName("joint1_34")][0] = (+10.0) * M_PI/180.0;
  mbchi.q[mb.jointIndexByName("joint2_34")][0] = (+10.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint3_34")][0] = (-10.0) * M_PI/180.0;
  mbclo.q[mb.jointIndexByName("joint4_34")][0] = (-10.0) * M_PI/180.0;
  HiLimit = rbd::dofToVector(mb, mbchi.q);
  LoLimit = rbd::dofToVector(mb, mbclo.q);

  //initial joint should not be exceeded from limits
  mbc.q[mb.jointIndexByName("joint1_01")][0] = (-90.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint2_01")][0] = (+90.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint3_01")][0] = (+90.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint4_01")][0] = (-90.0) * M_PI/180.0;
  
  mbc.q[mb.jointIndexByName("joint1_12")][0] = (-80.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint2_12")][0] = (+80.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint3_12")][0] = (-80.0) * M_PI/180.0;
  mbc.q[mb.jointIndexByName("joint4_12")][0] = (+80.0) * M_PI/180.0;

  //method = InverseMethodPtr(new LMInvConsideredSolvality(Wl, We));
  method = InverseMethodPtr(new LMInvConsideredSolvalityWithLimit(Wl, We, HiLimit, LoLimit));

  pub = nh.advertise<sensor_msgs::JointState>("joint_states", 1);
  sub = nh.subscribe<geometry_msgs::PoseArray>("gait_pose", 1, boost::bind(gaitCallback, _1, mb, mbc));

  for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
          std::string name = "joint" + std::to_string(i+1) + "_" + std::to_string(j) + std::to_string(j+1) + "_position_controller/command";
          command_pub[i*4 + j] = nh.advertise<std_msgs::Float64>(name, 1);
      }
  }

  ROS_INFO("start ik");
  ros::spin();

  return 0;
}


