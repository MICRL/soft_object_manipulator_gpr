//system and self
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "onlineGPR.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "onlineGPR_types.h"
#include "onlineGPR_terminate.h"
#include "onlineGPR_initialize.h"

//vision process library
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/passthrough.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/region_growing_rgb.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/point_types_conversion.h>
#include <pcl/point_types.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

//ros
#include <ros/ros.h>
#include <ros/package.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Int8.h>
#include <sensor_msgs/JointState.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/PointCloud2.h>
#include <cv_bridge/cv_bridge.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/chain.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <pcl_ros/point_cloud.h>

//stateRecorder
#include "StateRecorder.h"



#define N 4
#define M 2
#define F 2
#define FPS 10
#define MAX_VEL 0.04  //max allowed velocity
#define TOLERANCE 0.005
typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image,
                                                        sensor_msgs::PointCloud2> syncPolicy;
cv_bridge::CvImagePtr cv_ptr_RGB;
pcl::PointCloud<pcl::PointXYZRGB >::Ptr cldC(new pcl::PointCloud<pcl::PointXYZRGB>);
pcl::PointCloud<pcl::PointXYZRGB >::Ptr cldROIC(new pcl::PointCloud<pcl::PointXYZRGB>);
pcl::visualization::CloudViewer viewer ("Simple Cloud Viewer");

//real cartesian velocity of ee reciving from yumi topic
double cart_v_real[3*M] = {0};
//global variable feedback points
double p[3*F] = {0};

bool received = false,getFeedbackPoint = false,moved = false;


double norm(const double x[3])
{
  double y;
  double scale;
  int k;
  double absxk;
  double t;
  y = 0.0;
  scale = 2.2250738585072014E-308;
  for (k = 0; k < 3; k++) {
    absxk = std::abs(x[k]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0 + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * std::sqrt(y);
}
static void getFeatures(const double p[6],double c_p[4])
{
    int i,i0;
    double b_p[3];
    for (i = 0; i < 3; i++) {
        b_p[i] = p[i] - p[i + 3];
    }

    for (i0 = 0; i0 < 3; i0++) {
        c_p[i0] = (p[i0] + p[3 + i0]) / 2.0;
    }

    c_p[3] = norm(b_p);

}

//return the distance between two points
double dist(pcl::PointXYZRGB p1,pcl::PointXYZRGB p2)
{
    return (std::sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z)));
}

//return value higher than threshold
void highPass(double *pointer,int size, double threshold)
{
    for(int i=0;i<size;i++){
        if(std::abs(pointer[i])<threshold)
            pointer[i] = 0;
    }
}

//return soomth value
void lowPass(double *pointer,pcl::PointXYZRGB centroid11, pcl::PointXYZRGB centroid22)
{
    double lambda = 0.2;
    pointer[0] = lambda*pointer[0] + (1-lambda)*centroid11.x;
    pointer[1] = lambda*pointer[1] + (1-lambda)*centroid11.y;
    pointer[2] = lambda*pointer[2] + (1-lambda)*centroid11.z;
    pointer[3] = lambda*pointer[3] + (1-lambda)*centroid22.x;
    pointer[4] = lambda*pointer[4] + (1-lambda)*centroid22.y;
    pointer[5] = lambda*pointer[5] + (1-lambda)*centroid22.z;
}

//return the centroid of given points
pcl::PointXYZRGB getCentroid(std::vector<pcl::PointXYZRGB> p)
{
    pcl::PointXYZRGB centroid;
    for(int i=0;i<p.size();i++){
        centroid.x += p[i].x;
        centroid.y += p[i].y;
        centroid.z += p[i].z;
    }
    centroid.x /= p.size();
    centroid.y /= p.size();
    centroid.z /= p.size();
    return centroid;

}

//callback to get joint state
void cartVelCallback(const std_msgs::Float64MultiArray& msg)
{
    for (int i = 0; i < msg.data.size(); i++) {
        cart_v_real[i] = msg.data[i];
    }
    //high pass to filter noise
    highPass(cart_v_real,msg.data.size(),1e-4);

    for(int i=0;i<6;i++)
        if(cart_v_real[i] != 0)
            moved = true;
    received = true;

    // std::cout << "received real cartesian velocity from yumi topic\n";

}

bool start = false;
void startCallback(const std_msgs::Bool& msg)
{
    start = msg.data;
}


//time to start tracking stage
bool confirmed = false;
pcl::ExtractIndices<pcl::PointXYZRGB> eifilter (true); // Initializing with true will allow us to extract the removed indices
pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());
pcl::PointXYZRGB centroid11;
pcl::PointXYZRGB centroid22;
void pointCallback(const sensor_msgs::PointCloud2ConstPtr& pclPts){
    pcl::fromROSMsg<pcl::PointXYZRGB >(*pclPts,*cldC);

    pcl::PassThrough<pcl::PointXYZRGB> pass;
	pass.setInputCloud (cldC);

    pass.setFilterFieldName ("x");
    if(confirmed)
        pass.setFilterLimits (-0.3, 0.3);
    else
        pass.setFilterLimits (-0.1, 0.1);


	pass.filter (*cldROIC);


    pcl::PointXYZHSV temp;
    double xCenter = 0,yCenter = 0,zCenter = 0;
	for(int i=0;i<cldROIC->size();i++){
		pcl::PointXYZRGBtoXYZHSV(cldROIC->points[i],temp);
		if(temp.h>160 && temp.h<210 && temp.s>0.4){
			inliers->indices.push_back(i);
			xCenter +=cldROIC->points[i].x;

// 			yCenter +=cldROIC->points[i].y;
// 			zCenter +=cldROIC->points[i].z;
		}
	}
	xCenter /= inliers->indices.size();

	eifilter.setInputCloud (cldROIC);
	eifilter.setIndices (inliers);
	eifilter.filter (*cldROIC);
    inliers->indices.clear();


    std::vector<pcl::PointXYZRGB> selectedPoint1;
    std::vector<pcl::PointXYZRGB> selectedPoint2;
    std::vector<pcl::PointXYZRGB> selectedPoint11;
    std::vector<pcl::PointXYZRGB> selectedPoint22;
    for(int i=0;i<cldROIC->size();i++){
        if(confirmed){
            double dist1 = dist(cldROIC->points[i],centroid11);
            double dist2 = dist(cldROIC->points[i],centroid22);
            if(dist1 < std::min(dist2,0.02)){
                selectedPoint1.push_back(cldROIC->points[i]);
                cldROIC->points[i].g = 255;
            }
            else if(dist2 < std::min(dist1,0.02)){
                selectedPoint2.push_back(cldROIC->points[i]);
                cldROIC->points[i].g = 255;
            }

        }
        else{
            if(cldROIC->points[i].x < xCenter){
                selectedPoint1.push_back(cldROIC->points[i]);
            }
            else{
                selectedPoint2.push_back(cldROIC->points[i]);
            }
        }

    }


//     std::cout<<"size1: "<<sizes[0]<<"  size2:"<<sizes[1]<<std::endl;
    if(selectedPoint1.size() > 5 && selectedPoint2.size() > 5){
        pcl::PointXYZRGB centroid1,centroid2;
        //if confirmed, use last centroid to estimate current centroid.Otherwise compute the centroid use all points

        centroid1 = getCentroid(selectedPoint1);
        centroid2 = getCentroid(selectedPoint2);
        if(confirmed){
            selectedPoint11 = selectedPoint1;
            selectedPoint22 = selectedPoint2;
        }
        else{
            for(int i=0;i<selectedPoint1.size();i++){
                if(dist(selectedPoint1[i],centroid1)<0.1)
                    selectedPoint11.push_back(selectedPoint1[i]);
                }
            for(int i=0;i<selectedPoint2.size();i++){
                if(dist(selectedPoint2[i],centroid2)<0.1)
                    selectedPoint22.push_back(selectedPoint2[i]);
                }
        }
        if(selectedPoint11.size() > 5 && selectedPoint22.size() > 5){
            centroid11 = getCentroid(selectedPoint11);
            centroid22 = getCentroid(selectedPoint22);
            centroid11.r = 255;
            centroid22.r = 255;

            //high pass to filter noise
            highPass(p,6,1e-3);
            //low pass to filter noise
            lowPass(p,centroid11,centroid22);
            cldROIC->points.push_back(centroid11);
            cldROIC->points.push_back(centroid22);

            getFeedbackPoint = true;
            confirmed = true;
        }
    }




	viewer.showCloud(cldROIC);

}



void main_soft_object_manipulator(int argc,char ** argv)
{

    //onlineGPR input variable
    double x[4] = {0};    // value of the input of GPR;
    double x_v[4] = {0};  // real input of GPR; the difference between of the last frame and this frame
    double x_old[4] = {0};// last frame of the input
    double y[6] = {0};    // label of GPR
    double x_star[4] = {0}; // target position - current position
    static double gram_matrix_data[10000] = {0};  // inverse of K
    int gram_matrix_size[2] = {0};
    double X_data[400] = {0}; // aggregation of x_v
    int X_size[2] = {0};
    double Y_data[600] = {0}; // aggregation of y
    int Y_size[2] = {0};
    static double K_data[10000] = {0};
    int K_size[2] = {0};
    static double gram_matrix_update_data[10201] = {0};
    int gram_matrix_update_size[2] = {0};
    double y_star[6] = {0};
    double cov_star;
    double X_update_data[404] = {0};
    int X_update_size[2] = {0};
    double Y_update_data[606] = {0};
    int Y_update_size[2] = {0};
    static double K_update_data[10201] = {0};
    int K_update_size[2] = {0};

    //others
    double xd[4] = {0}; // target location
    //record state for debug
    StateRecorder stateRecorder;

    ros::init(argc,argv, "soft_object_manipulator", ros::init_options::NoSigintHandler);
    ros::NodeHandle node;
    //subscrib topic /yumi/joint_state to get joint state
    ros::Subscriber subCartVel = node.subscribe("/yumi/ikSloverVel_controller/state", 1000, cartVelCallback);
    //subscrib the start signal
    ros::Subscriber subStart = node.subscribe("/yumi/ikSloverVel_controller/go", 1000, startCallback);
    //publish velocity control command
    ros::Publisher cartVelPublisher = node.advertise<std_msgs::Float64MultiArray> ("/yumi/ikSloverVel_controller/command", 1);
    //publish pcd record command: 1 for record and 2 for replay
    ros::Publisher pcdRecorderPublisher = node.advertise<std_msgs::Int8> ("/pcd_recorder/command", 2);

//     message_filters::Subscriber<sensor_msgs::Image> imageRGB_sub(node, "/camera/rgb/image_color", 2);
// 	message_filters::Subscriber<sensor_msgs::PointCloud2> imagePoint_sub(node, "/camera/depth_registered/points", 2);
// 	message_filters::Synchronizer<syncPolicy> sync(syncPolicy(10),imageRGB_sub, imagePoint_sub);
// 	sync.registerCallback(boost::bind(&chatterCallbackCVC, _1, _2));
    ros::Subscriber imagePoint_sub = node.subscribe("/camera/depth_registered/points_SR300_611205001943", 2, pointCallback);

    //command
    std_msgs::Float64MultiArray armCommand;
    armCommand.data.resize(3*M);


    ros::Rate rate(FPS);
    bool first_time = true;
    bool stopFlag = false;
    int count = 0;
    while(ros::ok())
    {
        if(viewer.wasStopped())
            break;
        //reached the target
        if(stopFlag){
            if(start){
                for(int i=0;i<6;i++)
                    armCommand.data[i] = 0;
                cartVelPublisher.publish(armCommand);
            }
            std::cout << "Reached the target SUCCESSFULLY\n";
            break;
        }

        if(received && getFeedbackPoint){
            for(int i=0;i<4;i++)
                x_old[i] = x[i];
            getFeatures(p,x);
            for(int i=0;i<4;i++)
                x_v[i] = (x[i] - x_old[i])*FPS;
            for(int i=0;i<4;i++)
                x_star[i] = xd[i] - x[i];
            for(int i=0;i<6;i++)
                y[i] = cart_v_real[i];
            highPass(y,6,1e-4);
            if(first_time){
                double C[4] = {0,-0.05,0,0};
                for(int i=0;i<4;i++)
                    xd[i] = x[i] + C[i];
            }

            if(moved){
                if(first_time){
                    /*Initialize*/
                    //X
                    for(int i=0;i<4;i++)
                        X_data[i] = x_v[i];
                    X_size[0] = 4;
                    X_size[1] = 1;
                    //Y
                    for(int i=0;i<6;i++)
                        Y_data[i] = y[i];
                    Y_size[0] = 1;
                    Y_size[1] = 6;
                    //K
                    K_data[0] = 1;
                    K_size[0] = 1;
                    K_size[1] = 1;
                    //gram_matrix
                    double keps = 1e-4;
                    gram_matrix_data[0] = 1/(1+keps);
                    gram_matrix_size[0] = 1;
                    gram_matrix_size[1] = 1;

                    stateRecorder.record(x_v,y,x_star,y_star,X_data,Y_data,K_data,gram_matrix_data);
                    //record the target point yd
                    ofstream target_file;
                    target_file.open(ros::package::getPath("soft_object_manipulator")+"/record_data/yd.txt");
                    for(int i=0;i<4;i++)
                        target_file<<xd[i]<<' ';

                    first_time = false;
                }

                else {
                    //for debug
                    std::cout<<"K: ";
                    for(int i=0;i<9;i++)
                        std::cout<<K_data[i]<<" ";
                    std::cout<<"\n";
                    std::cout<<"gram_matrix: ";
                    for(int i=0;i<9;i++)
                        std::cout<<gram_matrix_data[i]<<" ";
                    std::cout<<"\n";



                    // Call the entry-point 'onlineGPR'.
                    onlineGPR(x_v, y, x_star, gram_matrix_data, gram_matrix_size, X_data, X_size,
                        Y_data, Y_size, K_data, K_size, gram_matrix_update_data,
                        gram_matrix_update_size, y_star, &cov_star, X_update_data,
                        X_update_size, Y_update_data, Y_update_size, K_update_data,
                        K_update_size);
                    stateRecorder.record(x_v,y,x_star,y_star,X_data,Y_data,K_data,gram_matrix_data);
                    //update K and gram_matrix
                    K_size[0] = K_update_size[0];
                    K_size[1] = K_update_size[1];
                    gram_matrix_size[0] = gram_matrix_update_size[0];
                    gram_matrix_size[1] = gram_matrix_update_size[1];
                    for(int i=0;i<K_size[1];i++){
                        for(int j=0;j<K_size[0];j++){
                            K_data[i*K_size[0]+j] = K_update_data[i*K_size[0]+j];
                            gram_matrix_data[i*K_size[0]+j] = gram_matrix_update_data[i*K_size[0]+j];
                        }
                    }

                    //update X
                    X_size[0] = X_update_size[0];
                    X_size[1] = X_update_size[1];
                    for(int i=0;i<X_size[1];i++)
                        for(int j=0;j<X_size[0];j++)
                            X_data[i*X_size[0]+j] = X_update_data[i*X_size[0]+j];

                    //update Y
                    Y_size[0] = Y_update_size[0];
                    Y_size[1] = Y_update_size[1];
                    for(int i=0;i<Y_size[1];i++)
                        for(int j=0;j<Y_size[0];j++)
                            Y_data[i*Y_size[0]+j] = Y_update_data[i*Y_size[0]+j];

                }
            }
//             std::cout<<"command v: ";
//             for(int i=0;i<6;i++)
//                 std::cout<< y_star[i]<<" ";

            received = false;
            getFeedbackPoint = false;
            moved = false;

            //check if reach the desired vaule
            stopFlag = true;
            std::cout<<"deltX: ";
            for(int i=0;i<4;i++){
                if(std::abs(xd[i] - x[i]) > TOLERANCE)
                    stopFlag = false;
                std::cout<<xd[i] - x[i]<<' ';
            }
            std::cout << "\n";

            std::cout << "Velocity: ";
            for(int i=0;i<6;i++){
                y_star[i] = y_star[i]<MAX_VEL?y_star[i]:MAX_VEL;
                y_star[i] = y_star[i]>-MAX_VEL?y_star[i]:-MAX_VEL;
                if(stopFlag)
                    armCommand.data[i] = 0;
                else
                    armCommand.data[i] = y_star[i];
                std::cout<<armCommand.data[i]<<' ';
            }
            std::cout << "\n";


            if(start){
                //record pcd and rgb image
                std_msgs::Int8 pcdRecorderCommand;
                pcdRecorderCommand.data = 1;
                pcdRecorderPublisher.publish(pcdRecorderCommand);

                //send arm velocity
                cartVelPublisher.publish(armCommand);
                std::cout << "sended velocity command\n";
            }
            std::cout << "\n-----------------------\n";

            }
            else{
                for(int i=0;i<3*M;i++)
//                     cart_v[i] = 0;
                std::cout << "Didn't receive the cart_v_real or get feedback points\n";
            }

            ros::spinOnce();
            rate.sleep();
        }
}


int main(int argc, char ** argv)
{
    // Initialize the application.
    // You do not need to do this more than one time.
    onlineGPR_initialize();

    // Invoke the entry-point functions.
    // You can call entry-point functions multiple times.
    main_soft_object_manipulator(argc,argv);

    // Terminate the application.
    // You do not need to do this more than one time.
    onlineGPR_terminate();
    return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
