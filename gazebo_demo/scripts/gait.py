#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rospy
import numpy as np
from geometry_msgs.msg import PoseArray, Pose

class GaitData():
    def __init__(self, hz=0.1, height=0.050):
        self.cycle_hz = hz
        self.cycle_time = 1.0/self.cycle_hz
        self.current_time = 0
        self.swing_height = height
        self.frame_id = 'base_link'
        self.stride = 0.100
        # base_link, yz-plane, the order of leg
        self.position_refs = [[+0.280, +0.155, +0.150 ], #limb1
                              [+0.280, -0.155, +0.150 ], 
                              [+0.280, +0.155, -0.150 ], 
                              [+0.280, -0.155, -0.150 ] ]

    # Call this at each sampling time
    def trot_phases(self, dt):
        if self.current_time >= self.cycle_time:
            self.current_time = 0
        self.current_time += dt
        phases = [ self.current_time ,\
                   self.current_time -1*self.cycle_time/4., \
                   self.current_time -2*self.cycle_time/4., \
                   self.current_time -3*self.cycle_time/4. ]
        for i in range(len(phases)):
            if phases[i] < 0:
                phases[i] += self.cycle_time

        return phases


def trot_gait(gait_data, pub_hz):
    tps = gait_data.trot_phases(1.0/pub_hz)
    poses = PoseArray()
    poses.poses = [Pose() for i in range(4)]
    # decide leg order
    leg_order = [0, 2, 3, 1] # TODO: Decide it based on target
    stamp = rospy.get_rostime()
    for i in range(len(tps)):
        leg_idx = leg_order[i]
        if 0<=tps[i]<=gait_data.cycle_time/4.: #swing
            poses.poses[leg_idx] = swing_pos(gait_data, leg_idx, tps[i], gait_data.cycle_time/4.0)
        else: #returning
            poses.poses[leg_idx] = stance_pos(gait_data, leg_idx, tps[i]-gait_data.cycle_time/4.0, gait_data.cycle_time*3.0/4.0)
        poses.header.stamp = stamp
        poses.header.frame_id = gait_data.frame_id

    return poses

# 0 <= phase <= end_phase
def swing_pos(gait_data, leg_idx, phase, end_phase):
    pose = Pose()

    ref_p = gait_data.position_refs[leg_idx]
    # TODO: Decide trajectory, temp go straight
    pose.position.x = ref_p[0] - 0.5*gait_data.swing_height*(1 - np.cos(2*np.pi/(end_phase)*phase))
    pose.position.y = ref_p[1]
    pose.position.z = ref_p[2] - 0.5*gait_data.stride*np.cos(np.pi/(end_phase)*phase)
    pose.orientation.w = 1

    return pose


# phase <= end_phase
def stance_pos(gait_data, leg_idx, phase, end_phase):
    pose = Pose()

    ref_p = gait_data.position_refs[leg_idx]
    # TODO: Decide trajectory, temp go straight
    pose.position.x = ref_p[0]
    pose.position.y = ref_p[1]
    pose.position.z = ref_p[2] + 0.5*gait_data.stride*np.cos(np.pi/(end_phase)*phase)
    pose.orientation.w = 1

    return pose

if __name__ == '__main__':
    try:
        rospy.init_node('gait', anonymous=True)
        rospy.loginfo("gait node start")
        pub = rospy.Publisher('gait_pose', PoseArray, queue_size=1)
        gd = GaitData()
        pub_hz = 100.0
        r = rospy.Rate(pub_hz)
        while not rospy.is_shutdown():
            pub.publish(trot_gait(gd, pub_hz))
            r.sleep()

    except rospy.ROSInterruptException: pass
