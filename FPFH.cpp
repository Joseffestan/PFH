//计算点云FPFH特征
//计算SPFH
	vector<FPFH> FinalFPFH;
	MatrixXf fSPFHHistF1, fSPFHHistF2, fSPFHHistF3;
	float fSearchRadius = 400;//邻域点搜寻范围（半径的平方）
	int iBinsNumF1 = 11, iBinsNumF2 = 11, iBinsNumF3 = 11;//直方图分区间数
	fSPFHHistF1.setZero(FeatureCloud->size(), iBinsNumF1);
	fSPFHHistF2.setZero(FeatureCloud->size(), iBinsNumF2);
	fSPFHHistF3.setZero(FeatureCloud->size(), iBinsNumF3);
	vector<int> vNearestPointsIdx;
	vector<float> vNearestPointsDist;
	float fPointsDist = 0.0;
	Vector4f fFPFHTuples = { 0.0,0.0,0.0,0.0 };
	for (int iter = 0; iter < FeatureCloud->size(); ++iter)
	{
		vNearestPointsIdx.clear();
		for (int i = 0; i < FeatureCloud->size(); ++i)
		{//寻找近邻点
			fPointsDist = PointDistance(FeatureCloud->at(iter), FeatureCloud->at(i));
			if (fPointsDist < fSearchRadius)
				vNearestPointsIdx.push_back(i);
		}
		float fHistogramIncrease = 100.0f /(vNearestPointsIdx.size() - 1);
		for (size_t idx = 0; idx < vNearestPointsIdx.size(); ++idx)
		{		
			if (0 == PointDistance(FeatureCloud->at(iter), FeatureCloud->at(vNearestPointsIdx.at(idx))))
				continue;//避免和自己进行计算
			// 计算点对四元组
			Vector4f p1(FeatureCloud->at(iter).x, FeatureCloud->at(iter).y, FeatureCloud->at(iter).z, 1.0f);
			Vector4f p2(FeatureCloud->at(vNearestPointsIdx.at(idx)).x, FeatureCloud->at(vNearestPointsIdx.at(idx)).y, FeatureCloud->at(vNearestPointsIdx.at(idx)).z, 1.0f);
			Vector4f n1(FeatureNormals->at(iter).normal_x, FeatureNormals->at(iter).normal_y, FeatureNormals->at(iter).normal_z, 0.0f);
			Vector4f n2(FeatureNormals->at(vNearestPointsIdx.at(idx)).normal_x, FeatureNormals->at(vNearestPointsIdx.at(idx)).normal_y, FeatureNormals->at(vNearestPointsIdx.at(idx)).normal_z, 0.0f);
			computePointPairFeatures(p1, n1, p2, n2, fFPFHTuples[0], fFPFHTuples[1], fFPFHTuples[2], fFPFHTuples[3]);
			//标准化 f1, f2, f3 并写入直方图
			int iHistBinIndex = floor(iBinsNumF1 * ((fFPFHTuples[0] + M_PI) * (1.0 / (2.0*M_PI))));
			if (iHistBinIndex < 0)           iHistBinIndex = 0;
			if (iHistBinIndex >= iBinsNumF1) iHistBinIndex = iBinsNumF1 - 1;
			fSPFHHistF1(iter, iHistBinIndex) += fHistogramIncrease;
			iHistBinIndex = floor(iBinsNumF2 * ((fFPFHTuples[1] + 1.0) * 0.5));
			if (iHistBinIndex < 0)           iHistBinIndex = 0;
			if (iHistBinIndex >= iBinsNumF2) iHistBinIndex = iBinsNumF2 - 1;
			fSPFHHistF2(iter, iHistBinIndex) += fHistogramIncrease;
			iHistBinIndex = floor(iBinsNumF3 * ((fFPFHTuples[2] + 1.0) * 0.5));
			if (iHistBinIndex < 0)           iHistBinIndex = 0;
			if (iHistBinIndex >= iBinsNumF3) iHistBinIndex = iBinsNumF3 - 1;
			fSPFHHistF3(iter, iHistBinIndex) += fHistogramIncrease;
		}
	}
	//对SPFH加权得到FPFH
	FPFH PointFPFH;//某点的SPFH
	for (int iter = 0; iter < FeatureCloud->size(); ++iter)
	{
		vNearestPointsIdx.clear();
		vNearestPointsDist.clear();
		for (int i = 0; i < 33; ++i)
			PointFPFH.fHistogram[i] = 0.0;
		for (int i = 0; i < FeatureCloud->size(); ++i)
		{//寻找近邻点
			fPointsDist = PointDistance(FeatureCloud->at(iter), FeatureCloud->at(i));
			if (fPointsDist < fSearchRadius)
			{
				vNearestPointsIdx.push_back(i);
				vNearestPointsDist.push_back(fPointsDist);
			}
		}
		float fSumF1 = 0.0, fSumF2 = 0.0, fSumF3 = 0.0;
		float fWeight = 0.0, fValueF1, fValueF2, fValueF3;
		// 加权该点的所有邻域点的SPFH
		for (size_t idx = 0; idx < vNearestPointsIdx.size(); ++idx)
		{
			if (vNearestPointsDist[idx] == 0)
				continue;
			fWeight = 1.0f / vNearestPointsDist[idx];//该邻域点的权重
			for (int f1_i = 0; f1_i < iBinsNumF1; ++f1_i)//加权计算FPFH
			{
				fValueF1 = fSPFHHistF1(vNearestPointsIdx.at(idx), f1_i) * fWeight;
				fSumF1 += fValueF1;
				PointFPFH.fHistogram[f1_i] += fValueF1;
			}
			for (int f2_i = 0; f2_i < iBinsNumF2; ++f2_i)
			{
				fValueF2 = fSPFHHistF2(vNearestPointsIdx.at(idx), f2_i) * fWeight;
				fSumF2 += fValueF2;
				PointFPFH.fHistogram[f2_i + iBinsNumF1] += fValueF2;
			}
			for (int f3_i = 0; f3_i < iBinsNumF3; ++f3_i)
			{
				fValueF3 = fSPFHHistF3(vNearestPointsIdx.at(idx), f3_i) * fWeight;
				fSumF3 += fValueF3;
				PointFPFH.fHistogram[f3_i + iBinsNumF1 + iBinsNumF2] += fValueF3;
			}
		}
		if (fSumF1 != 0) fSumF1 = 100.0 / fSumF1;//点的直方图求和标准化为100
		if (fSumF2 != 0) fSumF2 = 100.0 / fSumF2;
		if (fSumF3 != 0) fSumF3 = 100.0 / fSumF3;
		for (int f1_i = 0; f1_i < iBinsNumF1; ++f1_i)
			PointFPFH.fHistogram[f1_i] *= fSumF1;
		for (int f2_i = 0; f2_i < iBinsNumF2; ++f2_i)
			PointFPFH.fHistogram[f2_i + iBinsNumF1] *= fSumF2;
		for (int f3_i = 0; f3_i < iBinsNumF3; ++f3_i)
			PointFPFH.fHistogram[f3_i + iBinsNumF1 + iBinsNumF2] *= fSumF3;
		FinalFPFH.push_back(PointFPFH);		
	}
	FILE* pf2;
	fopen_s(&pf2, "e:\\my_fpfh.txt", "w");
	for (int i = 0; i < FinalFPFH.size(); ++i)
	{
		for (int j = 0; j < 33; ++j)
		{
			fprintf(pf2, "%f ", FinalFPFH.at(i).fHistogram[j]);
		}
		fprintf(pf2, "\n");
	}
	fclose(pf2);
