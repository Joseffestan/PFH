//求解点对四元组
bool computePairFeatures11(const Eigen::Vector4f &p1, const Eigen::Vector4f &n1,
	const Eigen::Vector4f &p2, const Eigen::Vector4f &n2,
	float &f1, float &f2, float &f3, float &f4)
{ 
	Eigen::Vector4f dp2p1 = p2 - p1;
	dp2p1[3] = 0.0f;
	// f4是两相邻点的距离
	f4 = dp2p1.norm();

	// 在进行数学运算前后，需要检查输出错误信息
	if (f4 == 0.0f)
	{
		PCL_DEBUG("[pcl::computePairFeatures] Euclidean distance between points is 0!\n");
		f1 = f2 = f3 = f4 = 0.0f;
		return (false);
	}

	// 计算f3，首先先确定谁是source，谁是target
	Eigen::Vector4f n1_copy = n1,
		n2_copy = n2;
	n1_copy[3] = n2_copy[3] = 0.0f;
	float angle1 = n1_copy.dot(dp2p1) / f4;

	// Make sure the same point is selected as 1 and 2 for each pair
	float angle2 = n2_copy.dot(dp2p1) / f4;
	if (std::acos(std::fabs(angle1)) > std::acos(std::fabs(angle2)))
	{
		// switch p1 and p2
		n1_copy = n2;
		n2_copy = n1;
		n1_copy[3] = n2_copy[3] = 0.0f;
		dp2p1 *= (-1);
		f3 = -angle2;
	}
	else
		f3 = angle1;

	// Create a Darboux frame coordinate system u-v-w
	// u = n1; v = (p_idx - q_idx) x u / || (p_idx - q_idx) x u ||; w = u x v
	// 齐次坐标相乘，第四维不需要管，只需要计算完叉乘后赋值0即可，因为叉乘是计算与现有向量都垂直的向量，因此前三项肯定满足
	Eigen::Vector4f v = dp2p1.cross3(n1_copy);
	v[3] = 0.0f;
	float v_norm = v.norm();
	// 计算完叉乘后，检查结果模长
	if (v_norm == 0.0f)
	{
		PCL_DEBUG("[pcl::computePairFeatures] Norm of Delta x U is 0!\n");
		f1 = f2 = f3 = f4 = 0.0f;
		return (false);
	}
	// Normalize v
	v /= v_norm;

	Eigen::Vector4f w = n1_copy.cross3(v);
	// Do not have to normalize w - it is a unit vector by construction

	v[3] = 0.0f;
	f2 = v.dot(n2_copy);
	w[3] = 0.0f;
	// Compute f1 = arctan (w * n2, u * n2) i.e. angle of n2 in the x=u, y=w coordinate system
	f1 = std::atan2(w.dot(n2_copy), n1_copy.dot(n2_copy)); // @todo optimize this

	return (true);
}
