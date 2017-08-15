#!/usr/bin/env python

import rospy
import matplotlib.pyplot as plt
from pointcloud_feature.msg import *

def callback(feature):
    plt.ion()
    plt.clf()
    plt.plot(feature.data)
    print max(feature.data), min(feature.data)
    plt.axis([0,128,-1,6])
    plt.show()
    plt.pause(0.5)

def listener():
    rospy.init_node('feature_plot', anonymous=True)
    rospy.Subscriber("/soft_object_manipulator_gpr_FeatureNet_node/feature_output", Feature, callback, queue_size = 1)
    rospy.spin()

if __name__ == "__main__":
    listener()
