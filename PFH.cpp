//计算点云PFH特征
//一点的PFH数组
struct PFH
{
	float fHistogram[125];
	//static int descriptorSize() { return 125; }
	//friend std::ostream& operator << (std::ostream& os, const PFHSignature125& p);
}

//pfh特征直方图计算
	vector<PFH> FinalPFH;
	int iSearchRadius = 400;//搜寻邻域点阈值 距离平方
	float fPFHTuples[4] = { 0.0,0.0,0.0,0.0 };//点对PFH四元组
	int iFeatureBinIndex[3];//点对四要素标准化的索引
	int iBinNum = 5;//特征值范围的子区间bin个数
	int iHistogramIndex = 0;//特征值归入直方图的索引
	int iMulti = 1;//特征值归入直方图的计算乘数
	for (int iter = 0; iter < Cloud->size(); ++iter)
	{
		PFH PointPFH;
		for (int i = 0; i < 125; ++i) 
			PointPFH.fHistogram[i] = 0.0;
		//查找邻域点
		vector<int> vNeiborPointIdx;
		for (int i = 0; i < FeatureCloud->size(); ++i)
		{
			if (PointDistance(FeatureCloud->at(iter), FeatureCloud->at(i)) < iSearchRadius)
				vNeiborPointIdx.push_back(i);
		}
		float fHistogramIncrease = 100.0f / static_cast<float> (vNeiborPointIdx.size() * (vNeiborPointIdx.size() - 1) / 2);
		//计算该点的PFH
		for (int i = 0; i < vNeiborPointIdx.size(); ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				Vector4f p1(FeatureCloud->at(vNeiborPointIdx.at(i)).x, FeatureCloud->at(vNeiborPointIdx.at(i)).y, FeatureCloud->at(vNeiborPointIdx.at(i)).z, 1.0f);
				Vector4f p2(FeatureCloud->at(vNeiborPointIdx.at(j)).x, FeatureCloud->at(vNeiborPointIdx.at(j)).y, FeatureCloud->at(vNeiborPointIdx.at(j)).z, 1.0f);
				Vector4f n1(FeatureNormals->at(vNeiborPointIdx.at(i)).normal_x, FeatureNormals->at(vNeiborPointIdx.at(i)).normal_y, FeatureNormals->at(vNeiborPointIdx.at(i)).normal_z, 0.0f);
				Vector4f n2(FeatureNormals->at(vNeiborPointIdx.at(j)).normal_x, FeatureNormals->at(vNeiborPointIdx.at(j)).normal_y, FeatureNormals->at(vNeiborPointIdx.at(j)).normal_z, 0.0f);
				computePairFeatures(p1, n1, p2, n2, fPFHTuples[0], fPFHTuples[1], fPFHTuples[2], fPFHTuples[3]);
				//标准化点对四要素到0 1 2 3 4
				iFeatureBinIndex[0] = static_cast<int> (floor(iBinNum * ((fPFHTuples[0] + M_PI) * (1.0/(2.0*M_PI)))));
				if (iFeatureBinIndex[0] < 0)        iFeatureBinIndex[0] = 0;
				if (iFeatureBinIndex[0] >= iBinNum) iFeatureBinIndex[0] = iBinNum - 1;
				iFeatureBinIndex[1] = static_cast<int> (floor(iBinNum * ((fPFHTuples[1] + 1.0) * 0.5)));
				if (iFeatureBinIndex[1] < 0)        iFeatureBinIndex[1] = 0;
				if (iFeatureBinIndex[1] >= iBinNum) iFeatureBinIndex[1] = iBinNum - 1;
				iFeatureBinIndex[2] = static_cast<int> (floor(iBinNum * ((fPFHTuples[2] + 1.0) * 0.5)));
				if (iFeatureBinIndex[2] < 0)        iFeatureBinIndex[2] = 0;
				if (iFeatureBinIndex[2] >= iBinNum) iFeatureBinIndex[2] = iBinNum - 1;
				// 归入直方图
				iHistogramIndex = 0;
				iMulti = 1;
				for (int d = 0; d < 3; ++d)
				{
					iHistogramIndex += iMulti * iFeatureBinIndex[d];
					iMulti *= iBinNum;
				}
				PointPFH.fHistogram[iHistogramIndex] += fHistogramIncrease;
			}
		}
	FinalPFH.push_back(PointPFH);
	}
