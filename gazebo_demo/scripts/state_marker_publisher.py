#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rospy
from geometry_msgs.msg import TwistStamped
from geometry_msgs.msg import Twist
from geometry_msgs.msg import PointStamped
from std_msgs.msg import Float32

## subscribe geometry_msgs/Twist
def twistCommandCallback(data, args):
    frame_id = args[0]
    pub = args[1]

    twistSt = TwistStamped()
    twistSt.header.frame_id = frame_id
    twistSt.header.stamp = rospy.get_rostime()
    twistSt.twist = data
    pub.publish(twistSt)

## subscrive geometry_msgs/PointStamped
def pressureHeightCallback(data, args):
    pub = args

    pub.publish(data.point.z)

if __name__ == '__main__':
    try:
        rospy.init_node('state_marker_publisher', anonymous=True)
        rospy.loginfo("marker publishing start")
        kameTwPub = rospy.Publisher('kame1/twist_stamped', TwistStamped, queue_size=1)
        toriTwPub = rospy.Publisher('tori1/twist_stamped', TwistStamped, queue_size=1)
        toriPressureHeightPub = rospy.Publisher('tori1/pressure_height_z', Float32, queue_size=1)

        kameTwSub = rospy.Subscriber("kame1/teleop_velocity_smoother/raw_cmd_vel", Twist, twistCommandCallback, ("kame1/base_link", kameTwPub))
        toriTwSub = rospy.Subscriber("tori1/cmd_vel", Twist, twistCommandCallback, ("tori1/base_link", toriTwPub))
        toriHeightSub = rospy.Subscriber("tori1/pressure_height", PointStamped, pressureHeightCallback, (toriPressureHeightPub))
        rospy.spin()

    except rospy.ROSInterruptException: pass
