#include <ros/ros.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDynUrdf/Reader.h>
#include <geometry_msgs/PoseArray.h>
#include <sensor_msgs/JointState.h>
#include "task_define.h"

rbdyn_urdf::Urdf robot_data;
MultiTaskPtr tasks;
//std::string ee_names[4] = {"link1_4", "link2_4", "link4_4", "link3_4"};
std::string ee_names[4] = {"l_wrist", "r_wrist", "r_ankle", "l_ankle"};
ros::Subscriber sub;
ros::Publisher pub;

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
  // geometry_msgs/PoseArray has data in the order of quadraits. 
  boost::shared_ptr<BodyTask> btask;
  auto m1 = msg->poses[0].position;
  sva::PTransformd X_O_T1 = sva::PTransformd(Vector3d(m1.x, m1.y, m1.z));
  btask = boost::dynamic_pointer_cast<BodyTask>(tasks[0].second);
  btask->_X_O_T = X_O_T1;
  auto m2 = msg->poses[1].position;
  sva::PTransformd X_O_T2 = sva::PTransformd(Vector3d(m2.x, m2.y, m2.z));
  btask = boost::dynamic_pointer_cast<BodyTask>(tasks[1].second);
  btask->_X_O_T = X_O_T2;
  auto m3 = msg->poses[2].position;
  sva::PTransformd X_O_T3 = sva::PTransformd(Vector3d(m3.x, m3.y, m2.z));
  btask = boost::dynamic_pointer_cast<BodyTask>(tasks[2].second);
  btask->_X_O_T = X_O_T3;
  auto m4 = msg->poses[3].position;
  sva::PTransformd X_O_T4 = sva::PTransformd(Vector3d(m4.x, m4.y, m4.z));
  //btask = boost::dynamic_pointer_cast<BodyTask>(tasks[3].second);
  //btask->_X_O_T = X_O_T4;
  boost::dynamic_pointer_cast<BodyTask>(tasks[3].second)->_X_O_T = X_O_T4;

  manyTaskMin(mb, mbc, tasks, 1.0, 200);
  sensor_msgs::JointState js_msg;
  JointStateFromMBC(mb, mbc, js_msg);
  pub.publish(js_msg);
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
  sva::PTransformd X_O_T = sva::PTransformd(Eigen::Vector3d(0,0,0));
  sva::PTransformd x_b_p = sva::PTransformd(Eigen::Vector3d(0.100,0,0));
  TaskPtr bodytask1(new BodyTask(mb, ee_names[0], X_O_T, x_b_p, "bodytask1"));
  TaskPtr bodytask2(new BodyTask(mb, ee_names[1], X_O_T, x_b_p, "bodytask2"));
  TaskPtr bodytask3(new BodyTask(mb, ee_names[2], X_O_T, x_b_p, "bodytask3"));
  TaskPtr bodytask4(new BodyTask(mb, ee_names[3], X_O_T, x_b_p, "bodytask4"));
  tasks.push_back(std::pair<double, TaskPtr>(10000000, bodytask1));
  tasks.push_back(std::pair<double, TaskPtr>(10000000, bodytask2));
  tasks.push_back(std::pair<double, TaskPtr>(10000000, bodytask3));
  tasks.push_back(std::pair<double, TaskPtr>(10000000, bodytask4));
  //tasks.push_back(std::pair<double, TaskPtr>(1.0, posturetask));

  pub = nh.advertise<sensor_msgs::JointState>("joint_states", 1);
  sub = nh.subscribe<geometry_msgs::PoseArray>("gait_pose", 1, boost::bind(gaitCallback, _1, mb, mbc));

  ROS_INFO("start ik");
  ros::spin();

  return 0;
}


