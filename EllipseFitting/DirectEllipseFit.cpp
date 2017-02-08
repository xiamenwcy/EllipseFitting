

#define _USE_MATH_DEFINES
#include <math.h>  //M_PI_2

#include "DirectEllipseFit.h"

Ellipse::Ellipse()
{
    algeFlag = false;
    a = b = c = d = e = f = 0;

    geomFlag = false;
    cx = cy = 0;
    rl = rs = 0;
    phi = 0;
}

void Ellipse::alge2geom()
{
	 if(!algeFlag)
	return;

	double tmp1 = b*b - 4*a*c;
	double tmp2 = sqrt((a-c)*(a-c)+b*b);
	double tmp3 = a*e*e + c*d*d - b*d*e + tmp1*f;

	double r1 = -sqrt(2*tmp3*(a+c+tmp2)) / tmp1;
	double r2 = -sqrt(2*tmp3*(a+c-tmp2)) / tmp1;
	rl = r1>=r2 ? r1 : r2;
	rs = r1<=r2 ? r1 : r2;

	cx = (2*c*d - b*e) / tmp1;
	cy = (2*a*e - b*d) / tmp1;

	phi = 0.5 * atan2(b, a-c);
	if(r1>r2)
	phi += M_PI_2;

	geomFlag = true;
}

void Ellipse::geom2alge()
{
	if(!geomFlag)
		return;

	a = rl*rl*sin(phi)*sin(phi) + rs*rs*cos(phi)*cos(phi);
	b = 2 * (rs*rs - rl*rl) * sin(phi) * cos(phi);
	c = rl*rl*cos(phi)*cos(phi) + rs*rs*sin(phi)*sin(phi);
	d = -2*a*cx - b*cy;
	e = -b*cx - 2*c*cy;
	f = a*cx*cx + b*cx*cy + c*cy*cy - rl*rl*rs*rs;

	algeFlag = true;
}

/*******************************************************************************
 * Template Class Defination
 ******************************************************************************/
DirectEllipseFit::DirectEllipseFit( const Eigen::VectorXd &xData,  const Eigen::VectorXd &yData)
{
	m_xData = xData;
	m_yData = yData;
}


Ellipse DirectEllipseFit::doEllipseFit()
{
    //Data preparation: normalize data
    Eigen::VectorXd xData = symmetricNormalize(m_xData);
    Eigen::VectorXd yData = symmetricNormalize(m_yData);

    //Bulid n*6 design matrix, n is size of xData or yData
   Eigen::MatrixXd dMtrx = getDesignMatrix(xData, yData);

    //Bulid 6*6 scatter matrix
    Eigen::MatrixXd sMtrx = getScatterMatrix(dMtrx);

    //Build 6*6 constraint matrix
    Eigen::MatrixXd cMtrx = getConstraintMatrix();

    //Solve eigensystem
	Eigen::VectorXd eigVal;
	Eigen::MatrixXd eigVec;
    bool flag = solveGeneralEigens(sMtrx, cMtrx,eigVal, eigVec);
    if(!flag)
      {
		  std::cout<<"Eigensystem solving failure!"<<std::endl;
		  return  Ellipse();
	  }

    Ellipse ellip = calcEllipsePara(eigVal, eigVec);

    return ellip;
}

double DirectEllipseFit::getMeanValue(const Eigen::VectorXd &data)
{
	/* T mean=0;
	for(int i=0; i<data.size(); ++i)
	mean += data(i); 

	return mean/data.size();*/
	return data.mean();
}

double DirectEllipseFit::getMaxValue(const Eigen::VectorXd &data)
{
	/*  T max = data[0];
	for(int i=1; i<data.size(); ++i)
	if(data.at(i)>max)
	max = data(i);

	return max;*/
	return data.maxCoeff();
	
}


double DirectEllipseFit::getMinValue(const Eigen::VectorXd &data)
{
	/* T min = data[0];
	for(int i=1; i<data.size(); ++i)
	if(data.at(i)<min)
	min = data(i);

	return min;*/
	return data.minCoeff();
}

double DirectEllipseFit::getScaleValue(const Eigen::VectorXd &data)
{

    return (0.5 * (getMaxValue(data) - getMinValue(data)));
}

Eigen::VectorXd DirectEllipseFit::symmetricNormalize(const Eigen::VectorXd &data)
{
	double mean = getMeanValue(data);
	double normScale = getScaleValue(data);

	int N=data.size();
	Eigen::VectorXd symData(N);
	for(int i=0; i<N; ++i)
		symData(i)=(data(i) - mean) / normScale;

	return symData;
}

Eigen::VectorXd DirectEllipseFit::dotMultiply(const Eigen::VectorXd &xData,
	const Eigen::VectorXd &yData)
{
	int  N=xData.size();
	Eigen::VectorXd product(N);

	for(int i=0; i<N; ++i)
		product(i)=xData(i)*yData(i);

	return product;
}


Eigen::MatrixXd  DirectEllipseFit::getDesignMatrix(const Eigen::VectorXd &xData, const Eigen::VectorXd &yData)
{
	int N=xData.size();
	Eigen::MatrixXd designMtrx(N,6);
	designMtrx.col(0)=dotMultiply(xData, xData);
	designMtrx.col(1)=dotMultiply(xData, yData);
	designMtrx.col(2)=dotMultiply(yData, yData);
	designMtrx.col(3)=xData;
	designMtrx.col(4)=yData;
	designMtrx.col(5)= Eigen::VectorXd::Ones(N);

	return designMtrx;
}

Eigen::MatrixXd DirectEllipseFit::getConstraintMatrix()
{
	Eigen::MatrixXd consMtrx= Eigen::MatrixXd::Zero(6,6);

	consMtrx(1,1) = 1;
	consMtrx(0,2) = -2;
	consMtrx(2,0) = -2;

	return consMtrx;
}


Eigen::MatrixXd DirectEllipseFit::getScatterMatrix(const Eigen::MatrixXd &dMtrx)
{
	return dMtrx.transpose()*dMtrx ;
}


bool DirectEllipseFit::solveGeneralEigens(const Eigen::MatrixXd &sMtrx,const Eigen::MatrixXd &cMtrx,Eigen::VectorXd &eigVal,Eigen::MatrixXd &eigVec)                                    
{
	integer N = sMtrx.cols(); // Number of columns of A and B. Number of rows of v.
	if (sMtrx.cols() != N  || sMtrx.rows()!=N || cMtrx.rows()!=N)
		return false;

	eigVec.resize(N,N);
	eigVal.resize(N);
	// lambda are the eigenvalues and are stored in lambda.
	// The eigenvalues are stored as: (lambda(:, 1) + lambda(:, 2)*i)./lambda(:, 3)
	Eigen::MatrixXd lambda;
	lambda.resize(N, 3);
	
	//每列的个数，即行数
	integer LDA = sMtrx.outerStride();
	integer LDB = cMtrx.outerStride();
	integer LDV = eigVec.outerStride();

	double WORKDUMMY;
	integer LWORK = -1; // Request optimum work size.
	integer INFO = 0;

	double * alphar = const_cast<double*>(lambda.col(0).data());
	double * alphai = const_cast<double*>(lambda.col(1).data());
	double * beta   = const_cast<double*>(lambda.col(2).data());

	// Get the optimum work size.
	dggev_("N", "V", &N, const_cast<double*>(sMtrx.data()), &LDA, const_cast<double*>(cMtrx.data()), &LDB, alphar, alphai, beta, 0, &LDV, eigVec.data(), &LDV, &WORKDUMMY, &LWORK, &INFO);

	LWORK = int(WORKDUMMY) + 32;
	Eigen::VectorXd WORK(LWORK);

	dggev_("N", "V", &N, const_cast<double*>(sMtrx.data()), &LDA, const_cast<double*>(cMtrx.data()), &LDB, alphar, alphai, beta, 0, &LDV, eigVec.data(), &LDV, WORK.data(), &LWORK, &INFO);

	

    for(int i=0; i<N; ++i)
	{
		eigVal(i)=alphar[i]/beta[i];
    }

   	return INFO==0;


}



Ellipse DirectEllipseFit::calcEllipsePara(const Eigen::VectorXd &eigVal,const Eigen::MatrixXd &eigVec)
{
    //Extract eigenvector corresponding to negative eigenvalue
    int eigIdx = -1;
    for(int i=0; i<eigVal.size(); ++i)
	{
        if(eigVal(i)<1e-6 && _finite(eigVal(i)))
		{
            eigIdx = i;
            break;
        }
    }
    if(eigIdx<0)
        return Ellipse();

    //Unnormalize and get coefficients of conic section
    double tA = eigVec(0,eigIdx);
    double tB = eigVec(1,eigIdx);
    double tC = eigVec(2,eigIdx);
    double tD = eigVec(3,eigIdx);
    double tE = eigVec(4,eigIdx);
    double tF = eigVec(5,eigIdx);

    double mx = getMeanValue(m_xData);
    double my = getMeanValue(m_yData);
    double sx = getScaleValue(m_xData);
    double sy = getScaleValue(m_yData);

    Ellipse ellip;
    ellip.a = tA * sy * sy;
    ellip.b = tB * sx * sy;
    ellip.c = tC * sx * sx;
    ellip.d = -2*tA*sy*sy*mx - tB*sx*sy*my + tD*sx*sy*sy;
    ellip.e = -tB*sx*sy*mx - 2*tC*sx*sx*my + tE*sx*sx*sy;
    ellip.f = tA*sy*sy*mx*mx + tB*sx*sy*mx*my + tC*sx*sx*my*my
            - tD*sx*sy*sy*mx - tE*sx*sx*sy*my + tF*sx*sx*sy*sy;
    ellip.algeFlag = true;

    ellip.alge2geom();

    return ellip;
}