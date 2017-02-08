/*******************************************************************************
 * ProName: DirectEllipseFit.h
 * Author:  Caiyong wang
 * Data:    2017/1/10
 * -----------------------------------------------------------------------------
 * INSTRUCTION: Perform ellipse fitting by direct method
 * DEPENDANCE:  clapack and Eigen
 * REFERENCE:
 *      (1) Fitzgibbon, A., et al. (1999). "Direct least square fitting of ellipses."
 *          IEEE Transactions on pattern analysis and machine intelligence 21(5):
 *          476-480.
 *      (2) http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/
 * CONVENTION:
 *      (1) Matrix expressed as QVector<VectorXd >, fast direction is rows,
 *          slow direction is columns, that is, internal QVector indicates column
 *          vector of T, external QVector indicates row vector of column vectors
 *      (2) Matrix expressed as 1-order array, arranged in rows, that is, row by
 *          row.
 ******************************************************************************/
#ifndef DIRECTELLIPSEFIT_H
#define DIRECTELLIPSEFIT_H


//Eigen
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>

//clapack,必须放在Eigen后面

#include <f2c.h>
#include <clapack.h>


class  Ellipse
{
public:
    Ellipse();
    /**
     * @brief alge2geom:    algebraic parameters to geometric parameters
     * @ref:    https://en.wikipedia.org/wiki/Ellipse#In_analytic_geometry
     *          http:homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/
     * @note:   The calculation of phi refer to wikipedia is not correct,
     *          refer to Bob Fisher's matlab program.
     *          What's more, the calculated geometric parameters can't back to
     *          initial algebraic parameters from geom2alge();
     */
    void alge2geom();
    /**
     * @brief geom2alge:    geometric parameters to algebraic parameters
     * @ref:    https://en.wikipedia.org/wiki/Ellipse#In_analytic_geometry
     */
    void geom2alge();

public:
    //algebraic parameters as coefficients of conic section
    float a, b, c, d, e, f;
    bool algeFlag;

    //geometric parameters
    float cx;   //centor in x coordinate
    float cy;   //centor in y coordinate
    float rl;   //semimajor: large radius
    float rs;   //semiminor: small radius
    float phi;  //azimuth angel in radian unit
    bool geomFlag;
};


class DirectEllipseFit
{
public:
	//Eigen Matrix

     DirectEllipseFit(const Eigen::VectorXd &xData,  const Eigen::VectorXd &yData);

    Ellipse doEllipseFit();

private:
    double getMeanValue(const Eigen::VectorXd &data);
    double getMaxValue(const Eigen::VectorXd &data);
    double getMinValue(const Eigen::VectorXd &data);
    double getScaleValue(const Eigen::VectorXd &data);
    Eigen::VectorXd symmetricNormalize(const Eigen::VectorXd &data);
    //Make sure xData and yData are of same size
    Eigen::VectorXd dotMultiply(const Eigen::VectorXd &xData, const Eigen::VectorXd &yData);
    //Get n*6 design matrix D, make sure xData and yData are of same size
	  Eigen::MatrixXd getDesignMatrix(const Eigen::VectorXd &xData,	const Eigen::VectorXd &yData);
    //Get 6*6 constraint matrix C
     Eigen::MatrixXd getConstraintMatrix();
    //Get 6*6 scatter matrix S from design matrix
    Eigen::MatrixXd getScatterMatrix(const Eigen::MatrixXd &dMtrx);

    /**
     * @brief solveGeneralEigens:   Solve generalized eigensystem
	 *  Solve:
	 *  A * v(j) = lambda(j) * B * v(j).
     * @note        For real eiginsystem solving.
     * @param sMtrx:    6*6 square matrix in this application
     * @param cMtrx:    6*6 square matrix in this application
     * @param eigVal:    eigenvalues, 6*1 matrix
	 * @param eigVec:     eigenvectors, 6*6 matrix
     * @return  success or failure status
     */

    bool solveGeneralEigens(const Eigen::MatrixXd &sMtrx,const Eigen::MatrixXd &cMtrx,Eigen::VectorXd &eigVal, Eigen::MatrixXd &eigVec);
                            
    /**
     * @brief calcEllipsePara:  calculate ellipse parameter form eigen information
	 * @param eigVal:    eigenvalues, 6*1 matrix
	 * @param eigVec:     eigenvectors, 6*6 matrix
     * @return ellipse parameter
     */
   Ellipse calcEllipsePara(const Eigen::VectorXd &eigVal,const Eigen::MatrixXd &eigVec);

private:
	Eigen::VectorXd m_xData;
	Eigen::VectorXd m_yData;

};


#endif // DIRECTELLIPSEFIT_H
