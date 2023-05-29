/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2021 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef COMMON_MATH_UTILS_H
#define COMMON_MATH_UTILS_H

#include "common/types.h"
#include <limits>
#include "Eigen/Dense"

namespace celib
{

    const double kExpNormTolerance = 1e-14;
    const double kLogTraceTolerance = 3.0 - kExpNormTolerance;


    const double kSqrt2 = std::sqrt(2.0);
    const double kSqrtPi = std::sqrt(M_PI);


    inline celib::Row9 mat3ToRow(celib::Mat3 R)
    {
        celib::Row9 output;
        output = Eigen::Map<celib::Row9>(R.data());
        return output;
    }


    inline celib::Mat3_9 jacobianYXv(const celib::Mat3& Y, const celib::Vec3& v)
    {
        celib::Mat3_9 output;
        output <<
            mat3ToRow((v * Y.row(0)).transpose()) ,
            mat3ToRow((v * Y.row(1)).transpose()) ,
            mat3ToRow((v * Y.row(2)).transpose());
        return output;
    }

    inline celib::Row3 jacobianNorm(const celib::Vec3& v)
    {
        celib::Row3 output;
        output = v.transpose()/(v.norm());
        return output;
    }

    inline celib::Mat3 jacobianNormalisation(const celib::Vec3& v)
    {
        celib::Mat3 output;
        output = ( (Mat3::Identity()*(v.norm())) - (v*jacobianNorm(v) ) )/(v.norm()*v.norm());
        return output;
    }

    inline celib::Mat9 jacobianYXW(const celib::Mat3& Y, const celib::Mat3& W)
    {
        celib::Mat9 output;
        output <<
            mat3ToRow((W.col(0) * Y.row(0)).transpose()) ,
            mat3ToRow((W.col(0) * Y.row(1)).transpose()) ,
            mat3ToRow((W.col(0) * Y.row(2)).transpose()) ,
            mat3ToRow((W.col(1) * Y.row(0)).transpose()) ,
            mat3ToRow((W.col(1) * Y.row(1)).transpose()) ,
            mat3ToRow((W.col(1) * Y.row(2)).transpose()) ,
            mat3ToRow((W.col(2) * Y.row(0)).transpose()) ,
            mat3ToRow((W.col(2) * Y.row(1)).transpose()) ,
            mat3ToRow((W.col(2) * Y.row(2)).transpose());
        return output;
    }

    inline celib::Mat9 jacobianXtYX(celib::Mat3 X, celib::Mat3 Y)
    {
        celib::Mat9 output;
        celib::Vec3 q1(1.0, 0.0, 0.0);
        celib::Vec3 q2(0.0, 1.0, 0.0);
        celib::Vec3 q3(0.0, 0.0, 1.0);
        output <<
            mat3ToRow(Y.transpose()*X*q1*(q1.transpose()) + Y*X*q1*(q1.transpose())), 
            mat3ToRow(Y.transpose()*X*q2*(q1.transpose()) + Y*X*q1*(q2.transpose())), 
            mat3ToRow(Y.transpose()*X*q3*(q1.transpose()) + Y*X*q1*(q3.transpose())), 
            mat3ToRow(Y.transpose()*X*q1*(q2.transpose()) + Y*X*q2*(q1.transpose())), 
            mat3ToRow(Y.transpose()*X*q2*(q2.transpose()) + Y*X*q2*(q2.transpose())), 
            mat3ToRow(Y.transpose()*X*q3*(q2.transpose()) + Y*X*q2*(q3.transpose())), 
            mat3ToRow(Y.transpose()*X*q1*(q3.transpose()) + Y*X*q3*(q1.transpose())), 
            mat3ToRow(Y.transpose()*X*q2*(q3.transpose()) + Y*X*q3*(q2.transpose())), 
            mat3ToRow(Y.transpose()*X*q3*(q3.transpose()) + Y*X*q3*(q3.transpose()));
        return output;
    }



    inline celib::Mat9 jacobianYX(const celib::Mat3& Y)
    {
        celib::Mat9 output;
        output = celib::Mat9::Zero();
        output.block<3,3>(0,0) = Y;
        output.block<3,3>(3,3) = Y;
        output.block<3,3>(6,6) = Y;
        return output;
    }
    inline celib::Mat3_9 jacobianXv(const celib::Vec3& v)
    {
        celib::Mat3_9 output;
        output <<
            v[0]   , 0.0    , 0.0    , v[1]   , 0.0    , 0.0    , v[2]   , 0.0    , 0.0,
            0.0    , v[0]   , 0.0    , 0.0    , v[1]   , 0.0    , 0.0    , v[2]   , 0.0,
            0.0    , 0.0    , v[0]   , 0.0    , 0.0    , v[1]   , 0.0    , 0.0    , v[2];
        return output;
    }

    inline celib::Mat9 jacobianTranspose()
    {
        celib::Mat9 output;
        output <<
            1.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 1.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 1.0    , 0.0    , 0.0,
            0.0    , 1.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 0.0    , 1.0    , 0.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 1.0    , 0.0,
            0.0    , 0.0    , 1.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 1.0    , 0.0    , 0.0    , 0.0,
            0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 0.0    , 1.0;
        return output;
    }

    inline celib::Mat9 jacobianXY(const celib::Mat3& Y)
    {
        celib::Mat9 output;
        output <<
            Y(0,0) , 0.0    , 0.0    , Y(1,0) , 0.0    , 0.0    , Y(2,0) , 0.0    , 0.0,
            0.0    , Y(0,0) , 0.0    , 0.0    , Y(1,0) , 0.0    , 0.0    , Y(2,0) , 0.0,
            0.0    , 0.0    , Y(0,0) , 0.0    , 0.0    , Y(1,0) , 0.0    , 0.0    , Y(2,0),
            Y(0,1) , 0.0    , 0.0    , Y(1,1) , 0.0    , 0.0    , Y(2,1) , 0.0    , 0.0,
            0.0    , Y(0,1) , 0.0    , 0.0    , Y(1,1) , 0.0    , 0.0    , Y(2,1) , 0.0,
            0.0    , 0.0    , Y(0,1) , 0.0    , 0.0    , Y(1,1) , 0.0    , 0.0    , Y(2,1),
            Y(0,2) , 0.0    , 0.0    , Y(1,2) , 0.0    , 0.0    , Y(2,2) , 0.0    , 0.0,
            0.0    , Y(0,2) , 0.0    , 0.0    , Y(1,2) , 0.0    , 0.0    , Y(2,2) , 0.0,
            0.0    , 0.0    , Y(0,2) , 0.0    , 0.0    , Y(1,2) , 0.0    , 0.0    , Y(2,2);
        return output;
    }


    inline Eigen::MatrixXd seKernel(const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const double l2, const double sf2)
    {
        Eigen::MatrixXd output(x1.size(), x2.size());
        Eigen::MatrixXd X1(x1.size(), x2.size());
        X1 = x1.replicate(1,x2.size());
        Eigen::MatrixXd X2(x1.size(), x2.size());
        X2 = x2.transpose().replicate(x1.size(),1);
        output = (((X1-X2).array().pow(2) * (-0.5/l2)).exp() * sf2).matrix();
        return output;
    }




    inline MatX cholesky(const MatX& K)
    {
        MatX L(K.rows(), K.cols());
        Eigen::LLT<Eigen::MatrixXd> lltOfA(K);
        L = lltOfA.matrixL();

        return L;
    }

    inline VecX solveKinvY(const MatX& K, const VecX& Y)
    {
        MatX L = cholesky(K);
        VecX alpha;
        alpha = L.triangularView<Eigen::Lower>().solve(Y);
        L.triangularView<Eigen::Lower>().transpose().solveInPlace(alpha);

        return alpha;
    }

    inline VecX sparseSolveKinvY(const MatX& K, const VecX& Y, const double epsilon)
    {
        sMatX sK = K.sparseView(0, epsilon);
        Eigen::SimplicialCholesky<sMatX> chol(sK);
        VecX alpha = chol.solve(Y);
        return alpha;
    }

    inline VecX iterSolveKinvY(const MatX& K, const VecX& Y, const double epsilon)
    {
        sMatX sK = K.sparseView(0, epsilon);
        Eigen::ConjugateGradient<sMatX, Eigen::Lower|Eigen::Upper> cg;
        cg.compute(sK);
        VecX alpha = cg.solve(Y);
        return alpha;
    }

    inline Mat2 angleToRotMat(const double theta)
    {
        Mat2 output;
        output << std::cos(theta), (-std::sin(theta)),
            std::sin(theta), std::cos(theta);
        return output;
    }

    inline Mat3 rotPosToHomogeneous(Mat2 rot, Vec2 pos)
    {
        Mat3 output = Mat3::Zero();
        output(2,2) = 1.0;
        output.block<2,2>(0,0) = rot;
        output.block<2,1>(0,2) = pos;
        return output;
    }

    inline Mat3 invHomogeneous(Mat3 H)
    {
        Mat3 output = Mat3::Zero();
        output(2,2) = 1.0;
        output.block<2,2>(0,0) = H.block<2,2>(0,0).transpose();
        output.block<2,1>(0,2) = -(H.block<2,2>(0,0).transpose())*(H.block<2,1>(0,2));
        return output;
    }

    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> seKernelTCov(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X, const T l2, const T sf2)
    {
        int nb_data = X.cols();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> K(nb_data, nb_data);
        for(int c = 0; c < nb_data; ++c)
        {
            for(int r = c; r < nb_data; ++r)
            {
                Eigen::Matrix<T, Eigen::Dynamic, 1> temp = X.col(c) - X.col(r);
                K(r,c) = temp.squaredNorm();
                if(r!=c) K(c,r) = K(r,c);
            }
        }
        return ((K*(-0.5/l2)).array().exp() *sf2).matrix();
    }

    template <typename T>
    inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> seKernelXD(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X_r, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& X_c, const T l2, const T sf2)
    {
        int nb_row = X_r.cols();
        int nb_col = X_c.cols();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> K(nb_row, nb_col);
        for(int r = 0; r < nb_row; ++r)
        {
            for(int c = 0; c < nb_col; ++c)
            {
                Eigen::Matrix<T, Eigen::Dynamic, 1> temp = X_r.col(r) - X_c.col(c);
                K(r,c) = temp.squaredNorm();
            }
        }
        return ((K*(-0.5/l2)).array().exp() *sf2).matrix();
    }

    // Implementation of the conjugate gradient method for solving linear systems
    inline VecX conjugateGradient(const MatX& A, const VecX& b, const VecX& x0 = VecX(0), const double epsilon=1e-6)
    {
        // Set x to x0 if it is not empty, otherwise set it to zero
        VecX x = x0.size() == 0 ? VecX::Zero(b.size()) : x0;
        VecX r = b - A*x;
        VecX p = r;
        double rTr = r.dot(r);
        double rTr_old = rTr;
        int iter = 0;
        while(rTr > epsilon)
        {
            VecX Ap = A*p;
            double alpha = rTr/(p.dot(Ap));
            x = x + alpha*p;
            r = r - alpha*Ap;
            rTr_old = rTr;
            rTr = r.dot(r);
            p = r + (rTr/rTr_old)*p;
            iter++;
        }
        std::cout << "Conjugate gradient iterations: " << iter << std::endl;
        return x;
    }

    
}

#endif