#!/usr/bin/env python

import rospy
from pointcloud_feature.msg import *
import numpy as np

def talker():
    pub = rospy.Publisher('feature_output', Feature, queue_size=1)
    rospy.init_node('talker', anonymous=True)
    rate = rospy.Rate(1)
    while not rospy.is_shutdown():
        f = Feature()
        f.data = list(np.random.rand(128))
        pub.publish(f)
        rate.sleep();

if __name__ == '__main__':
  try:
      talker()
  except rospy.ROSInterruptException:
      pass
