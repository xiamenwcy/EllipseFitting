#include "DirectEllipseFit.h"

int main()
{
	double X[]={185,182.854,181.225,180.086,178.973,177.868,176.755,175.622,174.296,173.13,171.969,170.798,169.279,168.079,166.893,164.909,163.696,162.138,160.956,159.762,158.209,157.017,155.826,154.192,152.988,150.933,148.857,147.078,144.993,142.906,140.81,139.01,136.911,135.068,132.962,130.85,128.978,127.072,124.962,123.033,120.929,118.985,117.026,114.933,112.97,110.997,109.017,107.03,104.968,102.983,100.992,98.997,96.9998,95.0016,93.0035,91.007,89.0142,86.9663,84.9676,82.9798,81.0034,79.0373,76.9419,74.9862,73.043,70.9357,69.0133,66.9003,64.9952,62.8779,60.9968,58.8842,57.0397,54.9433,52.8446,51.0558,48.9744,46.8964,44.8226,43.1216,41.9258,39.8931,38.2492,37.0425,35.8484,34.2504,33.055,31.8722,30.6703,29.1848,28.0179,26.8576,25.6873,24.2945,23.136,21.9922,20.8515,19.7037,18.3912,17.2477};
	double Y[]={61.198,61.5732,62.3475,63.1181,63.9648,64.8266,65.6698,66.469,67.4408,68.2048,68.9483,69.6458,70.5159,71.1517,71.7898,72.8201,73.4095,74.264,74.9162,75.5511,76.3976,77.0326,77.6426,78.4249,78.9721,79.8223,80.5938,81.2261,81.9788,82.7216,83.4149,84.0314,84.6993,85.2384,85.8604,86.4056,86.9039,87.3347,87.8113,88.1805,88.5699,88.8985,89.1848,89.4768,89.7476,89.9753,90.1773,90.3597,90.5298,90.663,90.7535,90.804,90.82,90.8109,90.7843,90.7399,90.6681,90.5499,90.3878,90.1938,89.9716,89.7166,89.3934,89.0862,88.7556,88.3275,87.9361,87.4582,87.0212,86.503,86.0123,85.4124,84.8682,84.1757,83.4614,82.8385,82.0713,81.272,80.4356,79.7146,79.1696,78.2374,77.4562,76.9088,76.3161,75.5034,74.8949,74.2351,73.5835,72.6877,71.9707,71.2272,70.4893,69.5496,68.7953,68.0116,67.2182,66.4283,65.4503,64.6586};
	//������ת��Ϊ����
	Eigen::Map<Eigen::VectorXd> xVec(X,100);
	Eigen::Map<Eigen::VectorXd> yVec(Y,100);
	//yVec��xVec���Ե�������ʹ�ã����߿�����ʽת��Ϊ����
	DirectEllipseFit ellipFit(xVec, yVec);
 	Ellipse ellip = ellipFit.doEllipseFit();

	return 0;
}