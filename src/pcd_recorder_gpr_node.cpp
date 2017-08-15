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


#define FPS 5

typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, 
                                                        sensor_msgs::PointCloud2> syncPolicy;
cv_bridge::CvImagePtr cv_ptr_RGB;
pcl::PointCloud<pcl::PointXYZRGB >::Ptr cldC(new pcl::PointCloud<pcl::PointXYZRGB>);
pcl::PointCloud<pcl::PointXYZRGB >::Ptr cldROIC(new pcl::PointCloud<pcl::PointXYZRGB>);

boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));


bool received = false,getFeedbackPoint = false;

bool start = false;
void startCallback(const std_msgs::Bool& msg)
{
    start = msg.data;
}


pcl::ExtractIndices<pcl::PointXYZRGB> eifilter (true); // Initializing with true will allow us to extract the removed indices
pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());
pcl::PointXYZRGB centroid11;
pcl::PointXYZRGB centroid22;

std::string data_path;
int count = 0;
int command_number = 0;
float yd[4] = {0};
bool target_received = false;
bool writer_have_open = false;
bool reader_have_open = false;
//callback to get pointcloud
pcl::PointCloud<pcl::PointXYZRGB >::Ptr pcd(new pcl::PointCloud<pcl::PointXYZRGB>);
void chatterCallbackCVC(const sensor_msgs::ImageConstPtr& imgRGB,const sensor_msgs::PointCloud2ConstPtr& pclPts){
    pcl::fromROSMsg<pcl::PointXYZRGB >(*pclPts,*cldC);
    cv_ptr_RGB = cv_bridge::toCvCopy(imgRGB,sensor_msgs::image_encodings::BGR8);
    
    
    //filtering to speed up
    pcl::VoxelGrid<pcl::PointXYZRGB> sor;
    sor.setInputCloud (cldC);
    sor.setLeafSize (0.001f, 0.001f, 0.001f);
    sor.filter (*cldC);

    //start mode
    if(command_number == 0) return;
    
    //record mode
    if(command_number == 1){
        viewer->removePointCloud("cloud");
        pcl::io::savePCDFile(data_path+std::to_string(count++)+".pcd",*cldC);
        viewer->addPointCloud(cldC);
        received = true;
    }
    
}


//callback to get command to decide to record or replay
void commandCallback(const std_msgs::Int8& msg)
{
    if(command_number != msg.data)
        count = 0;
    command_number = msg.data;
}






int main(int argc, char ** argv)
{
    ros::init(argc,argv, "pcd_recorder_node", ros::init_options::NoSigintHandler);
    ros::NodeHandle node;
    ros::Subscriber command_sub = node.subscribe("/pcd_recorder/command", 2, commandCallback);
    
    message_filters::Subscriber<sensor_msgs::Image> imageRGB_sub(node, "/camera/image/rgb_611205001943", 2);
	message_filters::Subscriber<sensor_msgs::PointCloud2> imagePoint_sub(node, "/camera/depth_registered/points_SR300_611205001943", 2);
	message_filters::Synchronizer<syncPolicy> sync(syncPolicy(10),imageRGB_sub, imagePoint_sub);
	sync.registerCallback(boost::bind(&chatterCallbackCVC, _1, _2));
    
    data_path = ros::package::getPath("soft_object_manipulator")+"/record_data/";
    
    //opencv
    cv::namedWindow("Live");
    cv::namedWindow("Record");
    cv::VideoWriter writer;
    cv::VideoCapture reader;
    
    
    //pcl
    Eigen::Affine3f cs;
    cs =  Eigen::AngleAxisf(M_PI, Eigen::Vector3f::UnitZ());
    viewer->setBackgroundColor (0, 0, 0);
//     viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    viewer->setCameraPosition(0,0,0,0,0,0.5,0,-1,0);


    ros::Rate rate(FPS);
    while(ros::ok() && !viewer->wasStopped())
    {
        if(received){
            if(!writer_have_open){
                writer.open(data_path+"video.avi", CV_FOURCC('M', 'J', 'P', 'G'), FPS, cv::Size(640,480));
                writer_have_open = true;
            }
            writer.write(cv_ptr_RGB->image);
            cv::imshow("Live", cv_ptr_RGB->image);
            cv::waitKey(5);
            received = false;
        }
        
        //replay mode
        if(command_number == 2){
            //show image
            if(!reader_have_open){
                reader.open(data_path+"video.avi");
                reader_have_open = true;
            }
            cv::Mat record_image;
            reader>>record_image;
            cv::imshow("Record",record_image);
            cv::waitKey(5);
            
            //show pcd
            viewer->removePointCloud("cloud");
            int ret = pcl::io::loadPCDFile(data_path+std::to_string(count++)+".pcd",*pcd);
            if(!target_received){
                viewer->removeAllShapes();
                ifstream target_file;
                target_file.open(data_path+"yd.txt");
                for(int i=0;i<4;i++)
                    target_file>>yd[i];
                pcl::PointXYZRGB center;
                center.x = yd[0];
                center.y = yd[1];
                center.z = yd[2];
                viewer->addSphere(center,yd[3]/2,1,0,0,"sphere");
                target_file.close();
                target_received = true;
            }
            if(ret==0){
                viewer->addPointCloud(pcd);
                
            }
            else{
                command_number = 0;
                target_received = false;
            }
                
        }
        viewer->spinOnce();
        ros::spinOnce();
        rate.sleep();
    }
    writer.release();
        
    return 0;
}

