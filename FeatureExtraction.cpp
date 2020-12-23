//基于法向量的点云特征提取
// FeaturePointExtraction.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/radius_outlier_removal.h>

using namespace std;
using namespace pcl;
using namespace Eigen;

//计算两点距离的平方
int PointDistance(PointXYZ P1, PointXYZ P2)
{
	int ans = ((P1.x - P2.x)*(P1.x - P2.x) + (P1.y - P2.y)*(P1.y - P2.y) + (P1.z - P2.z)*(P1.z - P2.z));
	return ans;
}
//计算某点法向量
Vector3d GetNormal(PointCloud<PointXYZ>::Ptr Cloud, int idx, int neibor_distance, vector<int> &neibor_point_idx)
{
	//vector<int> nNeiborPointIdx;//邻域点索引
	Vector3d vCurPoint = { Cloud->at(idx).x,Cloud->at(idx).y, Cloud->at(idx).z };
	for (int i = 0; i < Cloud->width; i++)
	{
		if (PointDistance(Cloud->at(idx), Cloud->at(i)) < neibor_distance)//邻域范围
			neibor_point_idx.push_back(i);
	}
	PointXYZ pNeiborHeart;
	//计算邻域点质心
	double dSumX(0.0), dSumY(0.0), dSumZ(0.0);
	for (int i = 0; i < neibor_point_idx.size(); i++)
	{
		dSumX += Cloud->at(neibor_point_idx.at(i)).x;
		dSumY += Cloud->at(neibor_point_idx.at(i)).y;
		dSumZ += Cloud->at(neibor_point_idx.at(i)).z;
	}
	pNeiborHeart.x = dSumX / neibor_point_idx.size();
	pNeiborHeart.y = dSumY / neibor_point_idx.size();
	pNeiborHeart.z = dSumZ / neibor_point_idx.size();
	//计算协方差矩阵
	Matrix3d Conv = Matrix3d::Zero();
	Vector3d vNeiborHeart = { pNeiborHeart.x, pNeiborHeart.y, pNeiborHeart.z };
	Vector3d vNeiborPoint;
	for (int i = 0; i < neibor_point_idx.size(); i++)
	{
		vNeiborPoint = { Cloud->at(neibor_point_idx.at(i)).x,Cloud->at(neibor_point_idx.at(i)).y, Cloud->at(neibor_point_idx.at(i)).z };
		Conv += (vNeiborPoint - vNeiborHeart)*(vNeiborPoint - vNeiborHeart).transpose();
	}
	Conv /= neibor_point_idx.size();
	//找到协方差矩阵的最小特征值
	SelfAdjointEigenSolver<Matrix3d> EigenSolver(Conv);
	Vector3d vEigenValues = EigenSolver.eigenvalues();//特征值
	Matrix3d mEigenVectors = EigenSolver.eigenvectors();//特征向量
	MatrixXf::Index MinEigenValueIdx;//最小特征值位置
	vEigenValues.real().rowwise().sum().minCoeff(&MinEigenValueIdx);
	Vector3d vMinEigenVector;
	vMinEigenVector << mEigenVectors.real()(0, MinEigenValueIdx), mEigenVectors.real()(1, MinEigenValueIdx), mEigenVectors.real()(2, MinEigenValueIdx);
	if (vMinEigenVector.dot(-vCurPoint) < 0)//法向量和坐标向量的点积
		vMinEigenVector = -vMinEigenVector;
	return vMinEigenVector;
}
//计算邻域点法向量
Vector3d GetNeiborNormal(PointCloud<PointXYZ>::Ptr Cloud, int idx, int neibor_distance)
{
	vector<int> neibor_point_idx;//邻域点索引
	Vector3d vCurPoint = { Cloud->at(idx).x,Cloud->at(idx).y, Cloud->at(idx).z };
	for (int i = 0; i < Cloud->width; i++)
	{
		if (PointDistance(Cloud->at(idx), Cloud->at(i)) < neibor_distance)//邻域范围
			neibor_point_idx.push_back(i);
	}
	PointXYZ pNeiborHeart;
	//计算邻域点质心
	double dSumX(0.0), dSumY(0.0), dSumZ(0.0);
	for (int i = 0; i < neibor_point_idx.size(); i++)
	{
		dSumX += Cloud->at(neibor_point_idx.at(i)).x;
		dSumY += Cloud->at(neibor_point_idx.at(i)).y;
		dSumZ += Cloud->at(neibor_point_idx.at(i)).z;
	}
	pNeiborHeart.x = dSumX / neibor_point_idx.size();
	pNeiborHeart.y = dSumY / neibor_point_idx.size();
	pNeiborHeart.z = dSumZ / neibor_point_idx.size();
	//计算协方差矩阵
	Matrix3d Conv = Matrix3d::Zero();
	Vector3d vNeiborHeart = { pNeiborHeart.x, pNeiborHeart.y, pNeiborHeart.z };
	Vector3d vNeiborPoint;
	for (int i = 0; i < neibor_point_idx.size(); i++)
	{
		vNeiborPoint = { Cloud->at(neibor_point_idx.at(i)).x,Cloud->at(neibor_point_idx.at(i)).y, Cloud->at(neibor_point_idx.at(i)).z };
		Conv += (vNeiborPoint - vNeiborHeart)*(vNeiborPoint - vNeiborHeart).transpose();
	}
	Conv /= neibor_point_idx.size();
	//找到协方差矩阵的最小特征值
	SelfAdjointEigenSolver<Matrix3d> EigenSolver(Conv);
	Vector3d vEigenValues = EigenSolver.eigenvalues();//特征值
	Matrix3d mEigenVectors = EigenSolver.eigenvectors();//特征向量
	MatrixXf::Index MinEigenValueIdx;//最小特征值位置
	vEigenValues.real().rowwise().sum().minCoeff(&MinEigenValueIdx);
	Vector3d vMinEigenVector;
	vMinEigenVector << mEigenVectors.real()(0, MinEigenValueIdx), mEigenVectors.real()(1, MinEigenValueIdx), mEigenVectors.real()(2, MinEigenValueIdx);
	if (vMinEigenVector.dot(-vCurPoint) < 0)//法向量和坐标向量的点积
		vMinEigenVector = -vMinEigenVector;
	return vMinEigenVector;
}
//点云可视化
void visualize_pcd(PointCloud<PointXYZ>::Ptr pcd_src, PointCloud<PointXYZ>::Ptr pcd_tgt, PointCloud<PointXYZ>::Ptr pcd_final)
{
	//int vp_1, vp_2;
	// Create a PCLVisualizer object
	pcl::visualization::PCLVisualizer viewer("registration Viewer");
	//viewer.createViewPort (0.0, 0, 0.5, 1.0, vp_1);
	// viewer.createViewPort (0.5, 0, 1.0, 1.0, vp_2);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_h(pcd_src, 0, 255, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> tgt_h(pcd_tgt, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> final_h(pcd_final, 0, 0, 255);
	viewer.addPointCloud(pcd_src, src_h, "source cloud");
	viewer.addPointCloud(pcd_tgt, tgt_h, "tgt cloud");
	viewer.addPointCloud(pcd_final, final_h, "final cloud");
	//viewer.addCoordinateSystem(1.0);
	while (!viewer.wasStopped())
	{
		viewer.spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}
int main()
{
	//读取点云
	char* strFile = new char[80];
	strcpy(strFile, "e:\\data\\office1.pcd");
	PointCloud<PointXYZ>::Ptr Cloud(new PointCloud<PointXYZ>);
	PointCloud<PointXYZ>::Ptr CloudVoxeled(new PointCloud<PointXYZ>);
	PointCloud<PointXYZ>::Ptr CloudFiltered(new PointCloud<PointXYZ>);
	io::loadPCDFile(strFile, *Cloud);
	//体素滤波
	VoxelGrid<PointXYZ> sor;
	sor.setInputCloud(Cloud);
	sor.setLeafSize(10, 10, 10);
	sor.filter(*CloudVoxeled);
	//离群点滤波
	pcl::RadiusOutlierRemoval<pcl::PointXYZ> pcFilter;
	pcFilter.setInputCloud(CloudVoxeled);
	pcFilter.setRadiusSearch(40);	// 搜索半径
	pcFilter.setMinNeighborsInRadius(3);	// 最少的邻居数目
	pcFilter.filter(*CloudFiltered);
	int nNeiborDistance = 400;//邻域范围
	PointCloud<PointXYZ>::Ptr FeatureCloud(new PointCloud<PointXYZ>);
	for (int iter = 0; iter < CloudFiltered->size(); iter++)
	{
		vector<int> nNeiborPointIdx;
		Vector3d vNormal = GetNormal(CloudFiltered, iter, nNeiborDistance, nNeiborPointIdx);//获得法向量
		cout << iter << " ";

		vector<Vector3d> vNeiborEigenNormals;
		for (int i = 0; i < nNeiborPointIdx.size(); i++)//计算邻域点法向量
			vNeiborEigenNormals.push_back(GetNeiborNormal(CloudFiltered, nNeiborPointIdx.at(i), nNeiborDistance));
		double nThreshold = 0.86;//阈值 判断是否为特征点
		int bIsFeaturePoint = 0;
		for (int i = 0; i < nNeiborPointIdx.size(); i++)
		{
			if (abs(vNormal.dot(vNeiborEigenNormals.at(i))) < nThreshold)
			{
				//FeatureCloud->push_back(CloudFiltered->at(iter));
				//break;
				bIsFeaturePoint++;
			}
		}
		if (bIsFeaturePoint > 2)
			FeatureCloud->push_back(CloudFiltered->at(iter));
	}
	visualize_pcd(FeatureCloud, FeatureCloud, FeatureCloud);//绿 红 蓝
	getchar();
	return 0;
}
