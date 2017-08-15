//system and self
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <boost/format.hpp>
#include "rt_nonfinite.h"
#include "onlineGPR.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "onlineGPR_types.h"
#include "onlineGPR_terminate.h"
#include "onlineGPR_initialize.h"
#include <pointcloud_feature/FeatureExtract.h>
#include <pointcloud_feature/Feature.h>

#include <omp.h>

//vision process library
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>
// #include <pcl/segmentation/lccp_segmentation.h>
#include <pcl/segmentation/supervoxel_clustering.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/point_types_conversion.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

//ros
#include <ros/ros.h>
#include <ros/package.h>
#include <std_msgs/Float32.h>
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
#include <pcl_conversions/pcl_conversions.h>

//stateRecorder
#include "StateRecorder.h"

#define USETIME 1

#if (USETIME)
  #include <time.h>
  struct timeval starttime, endtime;
  long mtime, seconds, useconds;
#endif

double RADIUS; //0.015f
int K; //10
int TOP; // 10
double DOWN_SAMPLE_SIZE_1; // 0.005f
double DOWN_SAMPLE_SIZE_2; // 0.01f

int GRID_SIZE ; // 32
int GRID_SIZE_2;
int GRID_SIZE_3;
int INTEREST_SIZE; // 24
int FEATURE_SIZE; // 128
int STEP; // 2
const int N = 4; //4
const int M = 2; //2
const int F = 1; //1
int FPS; //5
double MAX_VEL; // 0.04, max allowed velocity
double TOLERANCE; //0.005

typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, sensor_msgs::PointCloud2> syncPolicy;

int counter = 0;

ros::Publisher filtered_pub;
ros::Publisher final_pub;
ros::ServiceClient client;
ros::Publisher feature_pub;

cv_bridge::CvImagePtr cv_ptr_RGB;
pcl::PointCloud <pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud <pcl::PointXYZRGB>);
pcl::PointCloud <pcl::PointXYZRGB>::Ptr cloud_filtered (new pcl::PointCloud <pcl::PointXYZRGB>);
pcl::PointCloud <pcl::PointXYZRGB>::Ptr cloud_final (new pcl::PointCloud <pcl::PointXYZRGB>);
pcl::PointCloud <pcl::PointXYZRGB>::Ptr cloud_temp (new pcl::PointCloud <pcl::PointXYZRGB>);
pcl::PointCloud <pcl::PointXYZRGB>::Ptr cloud_downsampled (new pcl::PointCloud <pcl::PointXYZRGB>);
pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());
pcl::ExtractIndices<pcl::PointXYZRGB> eifilter (true); // Initializing with true will allow us to extract the removed indices
pcl::visualization::CloudViewer viewer ("Simple Cloud Viewer");
std::vector<pcl::PointXYZRGB> center_point;
std::vector<double> grid;
std::vector<double> feature;

//real cartesian velocity of ee reciving from yumi topic
double cart_v_real[3*M] = {0};
//global variable feedback points
double p[3*F] = {0};
double surface_curvature = 0;

bool received = false,getFeedbackPoint = false,moved = false;

int restrict_N(const int &a)
{
  if (a < 0) return 0;
  if (a > INTEREST_SIZE - 1) return INTEREST_SIZE - 1;
  return a;
}

static void getFeatures(std::vector<double> p,double *c)
{
  c = &p[0];
}

//return the distance between two points

double dist(pcl::PointXYZRGB p1,pcl::PointXYZRGB p2)
{
  return (std::sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z)));
}

//return value higher than threshold
void highPass(double *pointer,int size, double threshold)
{
  for(int i=0;i<size;i++)
  {
    if(std::abs(pointer[i])<threshold)
      pointer[i] = 0;
  }
}

//return soomth value
void lowPass(double *pointer, std::vector<pcl::PointXYZRGB> centroids)
{
  double lambda = 0.2;
  for (int i=0; i<F; i++)
  {
    pointer[0+i*3] = lambda*pointer[0+i*3] + (1-lambda)*centroids[i].x;
    pointer[1+i*3] = lambda*pointer[1+i*3] + (1-lambda)*centroids[i].y;
    pointer[2+i*3] = lambda*pointer[2+i*3] + (1-lambda)*centroids[i].z;
  }
}

//return the centroid of given points
pcl::PointXYZRGB getCentroid(pcl::PointCloud <pcl::PointXYZL> p)
{
  pcl::PointXYZRGB centroid;
  for(int i=0;i<p.size();i++)
  {
    centroid.x += p[i].x;
    centroid.y += p[i].y;
    centroid.z += p[i].z;
  }
  centroid.x /= p.size();
  centroid.y /= p.size();
  centroid.z /= p.size();
  return centroid;
}

//return the centroid of given points
pcl::PointXYZRGB getCentroid(pcl::PointCloud <pcl::PointXYZRGB> p)
{
  pcl::PointXYZRGB centroid;
  for(int i=0;i<p.size();i++)
  {
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
  for (int i = 0; i < msg.data.size(); i++)
  {
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
void pointCallback(const sensor_msgs::PointCloud2ConstPtr& cloud_msg)
{
  // Read the cloud from the ROS msg
  pcl::fromROSMsg(*cloud_msg, *cloud);

  // Perform the down-sampling
  // std::cout << cloud->size() << std::endl;
  pcl::VoxelGrid<pcl::PointXYZRGB> ds;
  ds.setInputCloud (cloud);
  ds.setLeafSize (DOWN_SAMPLE_SIZE_1, DOWN_SAMPLE_SIZE_1, DOWN_SAMPLE_SIZE_1);
  ds.filter (*cloud_downsampled);
  // std::cout << cloud_downsampled->size() << std::endl;

  // Cut the depth of the cloud
  pcl::PassThrough<pcl::PointXYZRGB> pass;
  pass.setInputCloud (cloud_downsampled);
  pass.setFilterFieldName ("z");
  pass.setFilterLimits (0.0, 2.0);
  pass.filter (*cloud_filtered);

  // Extract the points with the aprropriate color
  pcl::PointXYZHSV temp;
  for(int i=0;i<cloud_filtered->size();i++)
  {
		pcl::PointXYZRGBtoXYZHSV(cloud_filtered->points[i],temp);
		if(temp.h>160 && temp.h<210 && temp.s>0.6 && temp.v>0.25)
    {
			inliers->indices.push_back(i);
      cloud_filtered->points[i].g = 255;
		}
	}
  eifilter.setInputCloud (cloud_filtered);
	eifilter.setIndices (inliers);
	eifilter.filter (*cloud_filtered);

  inliers->indices.clear();

  if (cloud_filtered->size() > 5)
  {
    ds.setInputCloud (cloud_filtered);
    ds.setLeafSize (DOWN_SAMPLE_SIZE_2, DOWN_SAMPLE_SIZE_2, DOWN_SAMPLE_SIZE_2);
    ds.filter (*cloud_filtered);

    //Calculate the curvature and output the max and avg
    //for normal computing
    pcl::NormalEstimationOMP<pcl::PointXYZRGB, pcl::Normal> ne;
    ne.setInputCloud (cloud_filtered);
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB> ());
    ne.setSearchMethod (tree);
    ne.setKSearch (K);
    pcl::PointCloud<pcl::Normal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::Normal>);
    ne.compute (*cloud_with_normals);

    pcl::RegionGrowing<pcl::PointXYZRGB, pcl::Normal> reg;
    reg.setMinClusterSize (1);
    reg.setMaxClusterSize (1000000);
    reg.setSearchMethod (tree);
    reg.setNumberOfNeighbours (15);
    reg.setInputCloud (cloud_filtered);
    reg.setInputNormals (cloud_with_normals);
    reg.setSmoothnessThreshold (60.0 / 180.0 * M_PI);
    reg.setCurvatureThreshold (5.0);

    std::vector <pcl::PointIndices> clusters;
    reg.extract (clusters);

    if (center_point.size() == 0)
    {
      center_point.push_back(getCentroid(*cloud_filtered));
    }

    double min_distance = -1.0;
    for (int i = 0; i< clusters.size(); i++)
    {
      pcl::PointIndices::Ptr temp (new pcl::PointIndices ());
      *temp = clusters[i];
      eifilter.setInputCloud (cloud_filtered);
    	eifilter.setIndices (temp);
    	eifilter.filter (*cloud_temp);

      if (min_distance < 0)
      {
        min_distance = dist(getCentroid(*cloud_temp), center_point[0]);
        *cloud_final = *cloud_temp;
      }
      else
      {
        if (dist(getCentroid(*cloud_temp), center_point[0]) < min_distance)
        {
          *cloud_final = *cloud_temp;
        }
      }
    }

    // recompute the Normal
    ne.setInputCloud (cloud_final);
    ne.setSearchMethod (tree);
    ne.setKSearch (K);
    ne.compute (*cloud_with_normals);

    int point_cloud_size = cloud_final->size();

    if (point_cloud_size > 0)
    {
      // Compute the AABB
      float x_min = (*cloud_final)[0].x,
            x_max = (*cloud_final)[0].x,
            y_min = (*cloud_final)[0].y,
            y_max = (*cloud_final)[0].y,
            z_min = (*cloud_final)[0].z,
            z_max = (*cloud_final)[0].z;

      #pragma omp parallel for reduction(min:x_min, y_min, z_min), reduction(max:x_max, y_max, z_max)
      for (int i = 1; i < point_cloud_size; i++)
      {
        if ((*cloud_final)[i].x < x_min) x_min = (*cloud_final)[i].x;
        if ((*cloud_final)[i].y < y_min) y_min = (*cloud_final)[i].y;
        if ((*cloud_final)[i].z < z_min) z_min = (*cloud_final)[i].z;

        if ((*cloud_final)[i].x > x_max) x_max = (*cloud_final)[i].x;
        if ((*cloud_final)[i].y > y_max) y_max = (*cloud_final)[i].y;
        if ((*cloud_final)[i].z > z_max) z_max = (*cloud_final)[i].z;
      }

      float x_length = x_max - x_min,
            y_length = y_max - y_min,
            z_length = z_max - z_min;
      float fixed_length = max(max(x_length,y_length),z_length);
      x_min = (x_max + x_min - fixed_length) / 2;
      y_min = (y_max + y_min - fixed_length) / 2;
      z_min = (z_max + z_min - fixed_length) / 2;

      // Do the Voxel Gridlization
      grid.clear();
      grid.resize(GRID_SIZE_3);

      #pragma omp parallel for
      for (int i = 0; i < point_cloud_size; i++)
      {
        int x = restrict_N(int(((*cloud_final)[i].x - x_min) * INTEREST_SIZE / fixed_length)) + STEP;
        int y = restrict_N(int(((*cloud_final)[i].y - y_min) * INTEREST_SIZE / fixed_length)) + STEP;
        int z = restrict_N(int(((*cloud_final)[i].z - z_min) * INTEREST_SIZE / fixed_length)) + STEP;
        grid[z * GRID_SIZE_2 + y * GRID_SIZE + x] = 1;
      }

      pointcloud_feature::FeatureExtract srv;
      srv.request.grid = grid;

      // if (client)
      // {
      //   if (client.call(srv))
      //   {
      //     feature = srv.response.feature;
      //   }
      //   else
      //   {
      //     ROS_ERROR("Failed to call serive");
      //   }
      // }
      // else
      // {
      //   ROS_ERROR("Failed to connect service");
      // }

      if (ros::service::call("feature_extract", srv))
      {
        feature = srv.response.feature;
        std::cout << feature.size() << std::endl;
      }
      else
      {
        ROS_ERROR("Failed to call serive");
      }

      center_point.clear();
      center_point.push_back(getCentroid(*cloud_final));

      center_point[0].r = 255;
      cloud_final->push_back(center_point[0]);
      cloud_final->push_back(center_point[0]);

      //high pass to filter noise
      highPass(p,6,1e-3);
      //low pass to filter noise
      // lowPass(p,center_point);

      getFeedbackPoint = true;
      confirmed = true;
    }
  }

  // Convert to ROS data type
  sensor_msgs::PointCloud2 filtered_output;
  sensor_msgs::PointCloud2 final_output;
  pcl::toROSMsg(*cloud_filtered, filtered_output);
  pcl::toROSMsg(*cloud_final, final_output);
  filtered_output.header = cloud_msg->header;
  final_output.header = cloud_msg->header;

  // Publish the data
  filtered_pub.publish (filtered_output);
  final_pub.publish (final_output);

  pointcloud_feature::Feature feature_msg;
  feature_msg.data = feature;
  std::cout << "Here: " << feature_msg.data.size() << std::endl;
  feature_pub.publish (feature_msg);

  viewer.showCloud(cloud_final);
}



void main_soft_object_manipulator(int argc,char ** argv)
{
  //onlineGPR input variable
  double x[128] = {0};    // value of the input of GPR;
  double x_v[128] = {0};  // real input of GPR; the difference between of the last frame and this frame
  double x_old[128] = {0};// last frame of the input
  double y[6] = {0};    // label of GPR
  double x_star[128] = {0}; // target position - current position
  static double gram_matrix_data[10000] = {0};  // inverse of K
  int gram_matrix_size[2] = {0};
  double X_data[12800] = {0}; // aggregation of x_v
  int X_size[2] = {0};
  double Y_data[600] = {0}; // aggregation of y
  int Y_size[2] = {0};
  static double K_data[10000] = {0};
  int K_size[2] = {0};
  static double gram_matrix_update_data[10201] = {0};
  int gram_matrix_update_size[2] = {0};
  double y_star[6] = {0};
  double cov_star;
  double X_update_data[12928] = {0};
  int X_update_size[2] = {0};
  double Y_update_data[606] = {0};
  int Y_update_size[2] = {0};
  static double K_update_data[10201] = {0};
  int K_update_size[2] = {0};

  //others
  double xd[128] = {0}; // target location
  //record state for debug
  StateRecorder stateRecorder;

  ros::init(argc,argv, "soft_object_manipulator", ros::init_options::NoSigintHandler);
  ros::NodeHandle node("~");
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

  filtered_pub = node.advertise<sensor_msgs::PointCloud2> ("filtered_output", 1);
  final_pub = node.advertise<sensor_msgs::PointCloud2> ("final_output", 1);
  feature_pub = node.advertise<pointcloud_feature::Feature> ("feature_output", 1);

  client = node.serviceClient <pointcloud_feature::FeatureExtract> ("feature_extract", true);

  //command
  std_msgs::Float64MultiArray armCommand;
  armCommand.data.resize(3*M);

  #if (USETIME)
    gettimeofday(&starttime, NULL);
  #endif

  node.param("RADIUS", RADIUS, 0.015);
  node.param("K", K, 10);
  node.param("TOP", TOP, 10);
  node.param("DOWN_SAMPLE_SIZE_1", DOWN_SAMPLE_SIZE_1, 0.005);
  node.param("DOWN_SAMPLE_SIZE_2", DOWN_SAMPLE_SIZE_2, 0.01);

  node.param("FPS", FPS, 5);
  node.param("MAX_VEL", MAX_VEL, 0.04);
  node.param("TOLERANCE", TOLERANCE, 0.005);

  node.param("GRID_SIZE", GRID_SIZE, 32);
  GRID_SIZE_2 = GRID_SIZE*GRID_SIZE;
  GRID_SIZE_3 = GRID_SIZE_2*GRID_SIZE;
  node.param("INTEREST_SIZE", INTEREST_SIZE, 24);
  STEP = (GRID_SIZE - INTEREST_SIZE) / 2;
  node.param("FEATURE_SIZE", FEATURE_SIZE, 128);

  grid.resize(GRID_SIZE_3);

  ros::Rate rate(FPS);
  bool first_time = true;
  bool stopFlag = false;
  int count = 0;
  while(ros::ok())
  {

    #if (USETIME)
      gettimeofday(&endtime, NULL);
      seconds  = endtime.tv_sec  - starttime.tv_sec;
      useconds = endtime.tv_usec - starttime.tv_usec;
      mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
      std::cout << "total:\t" << mtime << std::endl << std::endl;
      gettimeofday(&starttime, NULL);
    #endif


    if(viewer.wasStopped()) break;
    //reached the target
    if(stopFlag)
    {
      if(start)
      {
        for(int i=0;i<6;i++)
            armCommand.data[i] = 0;
        cartVelPublisher.publish(armCommand);
      }
      std::cout << "Reached the target SUCCESSFULLY\n";
      break;
    }

    if(received && getFeedbackPoint)
    {
      for(int i=0;i<FEATURE_SIZE;i++)
        x_old[i] = x[i];
      getFeatures(feature,x);
      for(int i=0;i<FEATURE_SIZE;i++)
        x_v[i] = (x[i] - x_old[i])*FPS;
      for(int i=0;i<FEATURE_SIZE;i++)
        x_star[i] = xd[i] - x[i];
      for(int i=0;i<6;i++)
        y[i] = cart_v_real[i];
      highPass(y,6,1e-4);
      if(first_time)
      {
        double C[FEATURE_SIZE] = {0};
        for(int i=0;i<FEATURE_SIZE;i++)
          xd[i] = x[i] + C[i];
      }

      if(moved)
      {
        if(first_time)
        {
          /*Initialize*/
          //X
          for(int i=0;i<FEATURE_SIZE;i++)
            X_data[i] = x_v[i];
          X_size[0] = FEATURE_SIZE;
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
          for(int i=0;i<FEATURE_SIZE;i++)
            target_file<<xd[i]<<' ';

          first_time = false;
        }

        else
        {
          //for debug
          // std::cout<<"K: ";
          // for(int i=0;i<9;i++)
          //   std::cout<<K_data[i]<<" ";
          // std::cout<<"\n";
          // std::cout<<"gram_matrix: ";
          // for(int i=0;i<9;i++)
          //   std::cout<<gram_matrix_data[i]<<" ";
          // std::cout<<"\n";

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
          for(int i=0;i<K_size[1];i++)
          {
            for(int j=0;j<K_size[0];j++)
            {
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
      // std::cout<<"command v: ";
      // for(int i=0;i<6;i++)
      //     std::cout<< y_star[i]<<" ";

      received = false;
      getFeedbackPoint = false;
      moved = false;

      //check if reach the desired vaule
      stopFlag = true;
      std::cout<<"deltX: ";
      for(int i=0;i<FEATURE_SIZE;i++)
      {
        if(std::abs(xd[i] - x[i]) > TOLERANCE)
            stopFlag = false;
        std::cout<<xd[i] - x[i]<<' ';
      }
      std::cout << "\n";

      std::cout << "Velocity: ";
      for(int i=0;i<6;i++)
      {
        y_star[i] = y_star[i]<MAX_VEL?y_star[i]:MAX_VEL;
        y_star[i] = y_star[i]>-MAX_VEL?y_star[i]:-MAX_VEL;
        if(stopFlag)
            armCommand.data[i] = 0;
        else
            armCommand.data[i] = y_star[i];
        std::cout<<armCommand.data[i]<<' ';
      }
      std::cout << "\n";

      if(start)
      {
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
    else
    {
      if (!received && !getFeedbackPoint)
      {
        std::cout << "Didn't receive the cart_v_real and get feedback points\n";
      }
      else if (!getFeedbackPoint)
      {
        std::cout << "Didn't get feedback points\n";
      }
      else
      {
        std::cout << "Didn't receive the cart_v_real\n";
      }
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
