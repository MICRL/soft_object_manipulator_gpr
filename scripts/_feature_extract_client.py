#!/usr/bin/env python

import rospy
import numpy as np

from pointcloud_feature.srv import *

def feature_extract_client():
    rospy.wait_for_service('feature_extract')
    try:
        feature_extract = rospy.ServiceProxy('feature_extract', FeatureExtract)
        resp1 = feature_extract(tuple(np.random.rand(1,1,32,32,32).reshape(32*32*32)))
        for f in resp1.feature:
            print f
        return resp1.feature
    except rospy.ServiceException, e:
        print "Service call failed: %s" %e

if __name__ == "__main__":
    print "Requesting service"
    feature_extract_client()
