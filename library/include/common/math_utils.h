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


    //inline celib::Mat3 eulToRotMat(double eul_z, double eul_y, double eul_x)
    //{

    //    celib::Mat3 transform;
    //    double c1 = std::cos(eul_x);
    //    double c2 = std::cos(eul_y);
    //    double c3 = std::cos(eul_z);
    //    double s1 = std::sin(eul_x);
    //    double s2 = std::sin(eul_y);
    //    double s3 = std::sin(eul_z);
    //    transform << c1*c2, c1*s2*s3 - c3*s1, s1*s3 + c1*c3*s2,
    //              c2*s1, c1*c3 + s1*s2*s3, c3*s1*s2 - c1*s3,
    //            -s2, c2*s3, c2*c3;
    
    //    return transform;
    //}

    //inline celib::Mat3 eulToRotMat(std::vector<double> eul)
    //{
    //    if(eul.size() != 3) throw std::range_error("Wrong vector size for Euler to Rotation matrix conversion");
    //    return eulToRotMat(eul[2], eul[1], eul[0]);
    //}


    //// SO3 Log mapping
    //inline Vec3 logMap(const Mat3& rot_mat){
    //    double trace_mat = rot_mat.trace();
    //    if(trace_mat < kLogTraceTolerance){
    //        double phi = std::acos( (trace_mat-1) * 0.5);
    //        Mat3 skew_mat = (phi/(2.0*std::sin(phi))) * (rot_mat - rot_mat.transpose());
    //        Vec3 output;
    //        output << skew_mat(2,1), skew_mat(0,2), skew_mat(1,0);
    //        return output;
    //    }else{
    //        return Vec3::Zero();
    //    }
    //}   
    
    //// SO3 Exp mapping
    //inline Mat3 expMap(const Vec3& vec){
    //    double vec_norm = vec.norm();
    //    if(vec_norm > kExpNormTolerance){
    //        Mat3 skew_mat;
    //        skew_mat << 0.0, -vec(2), vec(1),
    //                    vec(2), 0.0, -vec(0),
    //                    -vec(1), vec(0), 0.0;
    //        return  Mat3::Identity()
    //            + ( (std::sin(vec_norm)/vec_norm) * skew_mat)
    //            + ( ( (1 - std::cos(vec_norm))/(vec_norm * vec_norm)) * skew_mat * skew_mat);


    //    }else{
    //        return Mat3::Identity();
    //    }
    //}

    //inline Mat3 toSkewSymMat(const Vec3& rot_vec)
    //{
    //    Mat3 skew_mat;
    //    skew_mat << 0.0, -rot_vec(2), rot_vec(1),
    //                rot_vec(2), 0.0, -rot_vec(0),
    //                -rot_vec(1), rot_vec(0), 0.0;
    //    return skew_mat;

    //}



    //// Righthand Jacobian of SO3 Exp mapping
    //template<typename T>
    //inline Eigen::Matrix<T, 3, 3> jacobianRighthandSO3( Eigen::Matrix<T, 3, 1> rot_vec)
    //{
    //    Eigen::Matrix<T, 3, 3> output = Eigen::Matrix<T, 3, 3>::Identity();
    //    T vec_norm = rot_vec.norm();

    //    
    //    Eigen::Matrix<T, 3, 1> vec = rot_vec;

    //    if(vec_norm > kExpNormTolerance)
    //    {

    //        Eigen::Matrix<T, 3, 3> skew_mat;
    //        skew_mat << T(0.0), T(-vec(2)), T(vec(1)),
    //                    T(vec(2)), T(0.0), T(-vec(0)),
    //                    T(-vec(1)), T(vec(0)), T(0.0);
    //        
    //        output += ( (vec_norm - sin(vec_norm)) / (vec_norm*vec_norm*vec_norm) )*skew_mat*skew_mat  - ( (1.0 - cos(vec_norm))/(vec_norm*vec_norm) )*skew_mat;
    //    }
    //    return output;
    //}

    //// Inverse Righthand Jacobian of SO3 Exp mapping
    //template<typename T>
    //inline Eigen::Matrix<T, 3, 3> inverseJacobianRighthandSO3( Eigen::Matrix<T, 3, 1> rot_vec)
    //{
    //    Eigen::Matrix<T, 3, 3> output = Eigen::Matrix<T, 3, 3>::Identity();
    //    T vec_norm = rot_vec.norm();

    //    
    //    Eigen::Matrix<T, 3, 1> vec = rot_vec;

    //    if(vec_norm > kExpNormTolerance)
    //    {

    //        Eigen::Matrix<T, 3, 3> skew_mat;
    //        skew_mat << T(0.0), T(-vec(2)), T(vec(1)),
    //                    T(vec(2)), T(0.0), T(-vec(0)),
    //                    T(-vec(1)), T(vec(0)), T(0.0);
    //        
    //        output += 0.5*skew_mat + ( ( (1.0/(vec_norm*vec_norm)) - ((1+std::cos(vec_norm))/(2.0*vec_norm*std::sin(vec_norm))) )*skew_mat*skew_mat);
    //    }
    //    return output;
    //}

    //// Inverse Lefthand Jacobian of SO3 Exp mapping
    //template<typename T>
    //inline Eigen::Matrix<T, 3, 3> inverseJacobianLefthandSO3( Eigen::Matrix<T, 3, 1> rot_vec)
    //{
    //    return inverseJacobianRighthandSO3<T>(-rot_vec);
    //}

    //// Lefthand Jacobian of SO3 Exp mapping
    //inline celib::Mat3 jacobianLefthandSO3( celib::Vec3 rot_vec)
    //{
    //    return jacobianRighthandSO3<double>(-rot_vec);
    //}

   //// Lefthand Jacobian of SE3 Exp mapping
    //inline celib::Mat6 jacobianLefthandSE3( celib::Vec6 vec)
    //{
    //    celib::Mat6 output = Mat6::Zero();
    //    output.block<3,3>(0,0) = jacobianLefthandSO3(vec.segment<3>(0));
    //    output.block<3,3>(3,3) = output.block<3,3>(0,0);

    //    Mat3 rho_skew = toSkewSymMat(vec.segment<3>(3));
    //    Mat3 phi_skew = toSkewSymMat(vec.segment<3>(0));
    //    double phi = vec.segment<3>(0).norm();

    //    Mat3 q = 0.5 * rho_skew
    //            + ((phi - std::sin(phi))/(std::pow(phi,3)))*(phi_skew*rho_skew + rho_skew*phi_skew + phi_skew*rho_skew*phi_skew)
    //            + ((std::pow(phi,2)+2.0*std::cos(phi)-2.0)/(2.0*std::pow(phi,4)))*(phi_skew*phi_skew*rho_skew + rho_skew*phi_skew*phi_skew - 3.0*phi_skew*rho_skew*phi_skew)
    //            + ((2.0*phi-3.9*std::sin(phi)+phi*std::cos(phi))/(2.0*std::pow(phi,5)))*(phi_skew*rho_skew*phi_skew*phi_skew + phi_skew*phi_skew*rho_skew*phi_skew);

    //    output.block<3,3>(0,3) = q;
    //    
    //    return output;
    //}

    //// Right hand Jacobian of SE3 Exp mapping
    //inline celib::Mat6 jacobianRighthandSE3( celib::Vec6 vec)
    //{   
    //    return jacobianLefthandSE3( - vec);
    //}

   //// Lefthand Jacobian of SE3 Exp mapping
    //inline celib::Mat6 inverseJacobianLefthandSE3( celib::Vec6 vec)
    //{
    //    celib::Mat6 output = Mat6::Zero();
    //    celib::Mat3 j_l_inv = jacobianLefthandSO3(vec.segment<3>(0));
    //    output.block<3,3>(0,0) = j_l_inv;
    //    output.block<3,3>(3,3) = j_l_inv;

    //    celib::Mat3 rho_skew = toSkewSymMat(vec.segment<3>(3));
    //    celib::Mat3 phi_skew = toSkewSymMat(vec.segment<3>(0));
    //    double phi = vec.segment<3>(0).norm();

    //    celib::Mat3 q = 0.5 * rho_skew
    //            + ((phi - std::sin(phi))/(std::pow(phi,3)))*(phi_skew*rho_skew + rho_skew*phi_skew + phi_skew*rho_skew*phi_skew)
    //            + ((std::pow(phi,2)+2.0*std::cos(phi)-2.0)/(2.0*std::pow(phi,4)))*(phi_skew*phi_skew*rho_skew + rho_skew*phi_skew*phi_skew - 3.0*phi_skew*rho_skew*phi_skew)
    //            + ((2.0*phi-3.9*std::sin(phi)+phi*std::cos(phi))/(2.0*std::pow(phi,5)))*(phi_skew*rho_skew*phi_skew*phi_skew + phi_skew*phi_skew*rho_skew*phi_skew);

    //    output.block<3,3>(0,3) = -j_l_inv*q*j_l_inv;
    //    
    //    return output;
    //}

    //// Right hand Jacobian of SE3 Exp mapping
    //inline celib::Mat6 inverseJacobianRighthandSE3( celib::Vec6 vec)
    //{   
    //    return inverseJacobianLefthandSE3( - vec);
    //}

    //// Adjoint element of SE3
    //inline celib::Mat6 adjointSE3(const celib::Mat3 rot, const celib::Vec3 pos)
    //{
    //    celib::Mat6 output = celib::Mat6::Zero();
    //    output.block<3,3>(0,0) = rot;
    //    output.block<3,3>(3,3) = rot;
    //    output.block<3,3>(0,3) = toSkewSymMat(pos)*rot;
    //    return output;
    //}

    //// Exponential map on SE3
    //inline std::pair<celib::Mat3, celib::Vec3> expMapSE3(celib::Vec6 vec)
    //{
    //    return {expMap(vec.segment<3>(0)), jacobianLefthandSO3(vec.segment<3>(0))*vec.segment<3>(3)};
    //}


    //// Exponential map on SE3
    //inline celib::Vec6 expMapSE3(celib::Mat3 rot, celib::Vec3 pos)
    //{
    //    celib::Vec6 output;
    //    output.segment<3>(0) = logMap(rot);
    //    output.segment<3>(3) = inverseJacobianLefthandSO3<double>(output.segment<3>(0))*pos;
    //    return output;
    //}

    //// Jacobian of the SO3 Exp mapping
    //inline celib::Mat9_3 jacobianExpMap(const celib::Vec3& rot_vec){

    //    celib::Mat9_3 output;
    //    double vec_norm = rot_vec.norm();

    //    if(vec_norm > kExpNormTolerance)
    //    {                                                                           
    //                                                                                                                
    //        double r1_2 = rot_vec(0) * rot_vec(0);
    //        double r2_2 = rot_vec(1) * rot_vec(1);
    //        double r3_2 = rot_vec(2) * rot_vec(2);
    //        double r12_22_23_15 = std::pow(r1_2+r2_2+r3_2, 1.5);
    //        double r12_22_23_2 = std::pow(r1_2+r2_2+r3_2, 2);
    //        double r_norm = std::sqrt(r1_2+r2_2+r3_2);
    //        double r1 = rot_vec(0);
    //        double r2 = rot_vec(1);
    //        double r3 = rot_vec(2);
    //                                                                                                                
    //        // Equation from MATLAB symbolic toolbox (might have a better formualtion, to inspect later)            
    //        output(0,0) = - (r1*std::sin(r_norm)*(r2_2 + r3_2))/r12_22_23_15 - (2*r1*(r2_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(0,1) = (2*r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r2*std::sin(r_norm)*(r2_2 + r3_2))/r12_22_23_15 - (2*r2*(r2_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(0,2) = (2*r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r3*std::sin(r_norm)*(r2_2 + r3_2))/r12_22_23_15 - (2*r3*(r2_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(1,0) = (r1*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*r3*std::sin(r_norm))/r12_22_23_15 + (r1_2*r2*std::sin(r_norm))/r12_22_23_15 + (2*r1_2*r2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(1,1) = (r2*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r2*r3*std::sin(r_norm))/r12_22_23_15 + (r1*r2_2*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(1,2) = std::sin(r_norm)/r_norm - (r3_2*std::sin(r_norm))/r12_22_23_15 + (r3_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(2,0) = (r1*r2*std::sin(r_norm))/r12_22_23_15 - (r1*r2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r1_2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1_2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(2,1) = (r2_2*std::sin(r_norm))/r12_22_23_15 - std::sin(r_norm)/r_norm - (r2_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(2,2) = (r2*r3*std::sin(r_norm))/r12_22_23_15 - (r2*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r1*r3_2*std::sin(r_norm))/r12_22_23_15 + (2*r1*r3_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(3,0) = (r1*r3*std::sin(r_norm))/r12_22_23_15 - (r1*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r1_2*r2*std::sin(r_norm))/r12_22_23_15 + (2*r1_2*r2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(3,1) = (r2*r3*std::sin(r_norm))/r12_22_23_15 - (r2*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r1*r2_2*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(3,2) = (r3_2*std::sin(r_norm))/r12_22_23_15 - std::sin(r_norm)/r_norm - (r3_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(4,0) = (2*r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*std::sin(r_norm)*(r1_2 + r3_2))/r12_22_23_15 - (2*r1*(r1_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(4,1) = - (r2*std::sin(r_norm)*(r1_2 + r3_2))/r12_22_23_15 - (2*r2*(r1_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(4,2) = (2*r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r3*std::sin(r_norm)*(r1_2 + r3_2))/r12_22_23_15 - (2*r3*(r1_2 + r3_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(5,0) = std::sin(r_norm)/r_norm - (r1_2*std::sin(r_norm))/r12_22_23_15 + (r1_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(5,1) = (r1*r2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*r2*std::sin(r_norm))/r12_22_23_15 + (r2_2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r2_2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(5,2) = (r1*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*r3*std::sin(r_norm))/r12_22_23_15 + (r2*r3_2*std::sin(r_norm))/r12_22_23_15 + (2*r2*r3_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(6,0) = (r1*r2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*r2*std::sin(r_norm))/r12_22_23_15 + (r1_2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1_2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(6,1) = std::sin(r_norm)/r_norm - (r2_2*std::sin(r_norm))/r12_22_23_15 + (r2_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(6,2) = (r2*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r2*r3*std::sin(r_norm))/r12_22_23_15 + (r1*r3_2*std::sin(r_norm))/r12_22_23_15 + (2*r1*r3_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(7,0) = (r1_2*std::sin(r_norm))/r12_22_23_15 - std::sin(r_norm)/r_norm - (r1_2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) + (r1*r2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r1*r2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(7,1) = (r1*r2*std::sin(r_norm))/r12_22_23_15 - (r1*r2*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r3*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r2_2*r3*std::sin(r_norm))/r12_22_23_15 + (2*r2_2*r3*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(7,2) = (r1*r3*std::sin(r_norm))/r12_22_23_15 - (r1*r3*std::cos(r_norm))/(r1_2 + r2_2 + r3_2) - (r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) + (r2*r3_2*std::sin(r_norm))/r12_22_23_15 + (2*r2*r3_2*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(8,0) = (2*r1*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r1*std::sin(r_norm)*(r1_2 + r2_2))/r12_22_23_15 - (2*r1*(r1_2 + r2_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(8,1) = (2*r2*(std::cos(r_norm) - 1))/(r1_2 + r2_2 + r3_2) - (r2*std::sin(r_norm)*(r1_2 + r2_2))/r12_22_23_15 - (2*r2*(r1_2 + r2_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //        output(8,2) = - (r3*std::sin(r_norm)*(r1_2 + r2_2))/r12_22_23_15 - (2*r3*(r1_2 + r2_2)*(std::cos(r_norm) - 1))/r12_22_23_2;
    //    }else{
    //        output <<   0,   0,   0,
    //                    0,   0,   1,
    //                    0,  -1,   0,
    //                    0,   0,  -1,
    //                    0,   0,   0,
    //                    1,   0,   0,
    //                    0,   1,   0,
    //                   -1,   0,   0,
    //                    0,   0,   0;
    //    }
    //    return output;
    //}


    //inline celib::Mat3_9 jacobianLogMap(const celib::Mat3& rot_mat)
    //{
    //    celib::Mat3_9 output;
    //    double trace_mat = rot_mat.trace();
    //    
    //    if(trace_mat < kLogTraceTolerance){

    //        // Equation from MATLAB symbolic toolbox (might have a better formualtion, to inspect later)
    //        output(0,0) = - (rot_mat(1,2) - rot_mat(2,1))/(4*(std::pow(rot_mat(0,0)/2.0
    //                        + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) - 
    //                        (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*
    //                         (rot_mat(1,2) - rot_mat(2,1))*(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0
    //                             + rot_mat(2,2)/2.0 - 0.5))/(4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 
    //                                 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(0,1) = 0.0;
    //        output(0,2) = 0.0;
    //        output(0,3) = 0.0;
    //        output(0,4) = - (rot_mat(1,2) - rot_mat(2,1))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) -
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*
    //             (rot_mat(1,2) - rot_mat(2,1))*(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(0,5) = std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)
    //            /(2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(0,6) = 0.0;
    //        output(0,7) = -std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)/
    //            (2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(0,8) = - (rot_mat(1,2) - rot_mat(2,1))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) -
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*
    //             (rot_mat(1,2) - rot_mat(2,1))*(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(1,0) = (rot_mat(0,2) - rot_mat(2,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) + 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*
    //             (rot_mat(0,2) - rot_mat(2,0))*(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(1,1) = 0.0;
    //        output(1,2) = -std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)/
    //            (2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(1,3) = 0.0;
    //        output(1,4) = (rot_mat(0,2) - rot_mat(2,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) + 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*
    //             (rot_mat(0,2) - rot_mat(2,0))*(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/(
    //             4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(1,5) = 0.0;
    //        output(1,6) = std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)/
    //            (2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(1,7) = 0.0;
    //        output(1,8) = (rot_mat(0,2) - rot_mat(2,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) + 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*(rot_mat(0,2) - rot_mat(2,0))*
    //             (rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(2,0) = - (rot_mat(0,1) - rot_mat(1,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) - 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*(rot_mat(0,1) - rot_mat(1,0))*
    //             (rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(2,1) = std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)/
    //            (2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(2,2) = 0.0;
    //        output(2,3) = -std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)/
    //            (2*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),0.5));
    //        output(2,4) = - (rot_mat(0,1) - rot_mat(1,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) - 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*(rot_mat(0,1) - rot_mat(1,0))*
    //             (rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //        output(2,5) = 0.0;
    //        output(2,6) = 0.0;
    //        output(2,7) = 0.0;
    //        output(2,8) = - (rot_mat(0,1) - rot_mat(1,0))/
    //            (4*(std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2) - 1)) - 
    //            (std::acos(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5)*(rot_mat(0,1) - rot_mat(1,0))*
    //             (rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5))/
    //            (4*std::pow(1 - std::pow(rot_mat(0,0)/2.0 + rot_mat(1,1)/2.0 + rot_mat(2,2)/2.0 - 0.5,2),1.5));
    //    }else{
    //        output <<   0,   0,    0,    0, 0, 0.5,   0, -0.5, 0,
    //                    0,   0, -0.5,    0, 0,   0, 0.5,    0, 0,
    //                    0, 0.5,    0, -0.5, 0,   0,   0,    0, 0;
    //    }
    //    return output;
    //}

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

    //inline void testJacobianXtYX()
    //{
    //    celib::Mat3 X = celib::Mat3::Random();
    //    celib::Mat3 Y = celib::Mat3::Random();

    //    celib::Mat3 res = X.transpose() * Y * X;
    //    double quantum = 0.0001;

    //    auto jacob = jacobianXtYX(X,Y);
    //    celib::Mat9 num_jacob;

    //    for(int r = 0; r <3; ++r)
    //    {
    //        for(int c = 0; c < 3; ++c)
    //        {
    //            celib::Mat3 temp = X;
    //            temp(r,c) += quantum;
    //            auto temp_res = temp.transpose() * Y * temp;

    //            num_jacob.col(c*3+r) = mat3ToRow((temp_res-res)/quantum).transpose();

    //        }
    //    }
    //    std::cout << "Test jacobian XtYX" << std::endl;
    //    std::cout << "Analitical\n" << jacob << std::endl;
    //    std::cout << "Numerical\n" << num_jacob << std::endl;

    //}


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



    //Eigen::MatrixXd seKernelIntegral(const double a, const Eigen::VectorXd& b, const Eigen::VectorXd& x2, const double l2, const double sf2);



    //Eigen::MatrixXd seKernelIntegralDt(const double a, const Eigen::VectorXd& b, const Eigen::VectorXd& x2, const double l2, const double sf2);



    //Eigen::MatrixXd seKernelIntegral2(const double a, const Eigen::VectorXd& b, const Eigen::VectorXd& x2, const double l2, const double sf2);


    //Eigen::MatrixXd seKernelIntegral2Dt(const double a, const Eigen::VectorXd& b, const Eigen::VectorXd& x2, const double l2, const double sf2);


    //inline double seKernelIntegralBothSides(const double a, const double b, const double l2, const double sf2)
    //{
    //    double sqrt2 = std::sqrt(2);
    //    double sqrt_inv_l2 = std::sqrt(1.0/l2);
    //    return 2.0*l2*sf2*std::exp(-std::pow(a - b,2)/(2.0*l2)) - 2.0*l2*sf2 + (sqrt2*sf2*std::sqrt(M_PI)*erf((sqrt2*(a - b)*sqrt_inv_l2)/2.0)*(a - b))/sqrt_inv_l2;

    //}
    //inline double seKernelIntegral2BothSides(const double a, const double b, const double l2, const double sf2)
    //{
    //    return 1.0;
    //}


    
    //inline void preintPredict(const PreintMeas& preint, const Mat3& state_rot, const Vec3& state_pos, const Vec3& state_vel, const Vec3& state_bf, const Vec3& state_bw, const double& state_dt, const Vec3& gravity, Mat3& out_rot, Vec3& out_pos, Vec3& out_vel)
    //{
    //    Mat3 d_r = preint.delta_R * expMap(preint.d_delta_R_d_bw*state_bw + preint.d_delta_R_d_t * state_dt);
    //    Vec3 d_p = preint.delta_p + preint.d_delta_p_d_bf*state_bf + preint.d_delta_p_d_bw*state_bw + preint.d_delta_p_d_t * state_dt;
    //    Vec3 d_v = preint.delta_v + preint.d_delta_v_d_bf*state_bf + preint.d_delta_v_d_bw*state_bw + preint.d_delta_v_d_t * state_dt;
    //    out_rot = state_rot * d_r;
    //    out_vel = state_vel + (gravity*preint.dt) + (state_rot*d_v);
    //    out_pos = state_pos + (state_vel*preint.dt) + (gravity*preint.dt_sq_half) + state_rot*d_p;
    //}

    //inline void preintPredict(const PreintMeas& preint, const FrameState& state, const Vec3& gravity, Mat3& out_rot, Vec3& out_pos, Vec3& out_vel)
    //{
    //    const Mat3 state_rot = state.rot();
    //    const Vec3 state_pos = state.pos();
    //    const Vec3 state_vel = state.vel();
    //    const Vec3 state_bf = state.accBias();
    //    const Vec3 state_bw = state.gyrBias();
    //    const double state_dt = state.d_t[0];
    //    preintPredict(preint, state_rot, state_pos, state_vel, state_bf, state_bw, state_dt, gravity, out_rot, out_pos, out_vel);
    //}

    //// Output the Jacobian as
    //// d_rot_d_rot, d_rot_d_bw, d_rot_d_dt, d_vel_d_rot, d_vel_d_vel, d_vel_d_bf, d_vel_d_bw, d_vel_d_dt, d_vel_d_g, d_p_d_rot, d_p_d_vel, d_p_d_p, d_p_d_bf, d_p_d_bw, d_p_d_dt, d_p_d_g
    //inline std::tuple<Mat9, Mat9_3, Vec9, Mat3_9, Mat3, Mat3, Mat3, Vec3, Mat3, Mat3_9, Mat3, Mat3, Mat3, Mat3, Vec3, Mat3> preintPredictAndJacobian(const PreintMeas& preint,  const Mat3& state_rot, const Vec3& state_pos, const Vec3& state_vel, const Vec3& state_bf, const Vec3& state_bw, const double& state_dt, const Vec3& gravity, Mat3& out_rot, Vec3& out_pos, Vec3& out_vel)
    //{

    //    Vec3 delta_R_correction = preint.d_delta_R_d_bw*state_bw + preint.d_delta_R_d_t * state_dt;
    //    Mat3 d_r = preint.delta_R * expMap(delta_R_correction);
    //    Vec3 d_p = preint.delta_p + preint.d_delta_p_d_bf*state_bf + preint.d_delta_p_d_bw*state_bw + preint.d_delta_p_d_t * state_dt;
    //    Vec3 d_v = preint.delta_v + preint.d_delta_v_d_bf*state_bf + preint.d_delta_v_d_bw*state_bw + preint.d_delta_v_d_t * state_dt;
    //    out_rot = state_rot * d_r;
    //    out_vel = state_vel + (gravity*preint.dt) + (state_rot*d_v);
    //    out_pos = state_pos + (state_vel*preint.dt) + (gravity*preint.dt_sq_half) + state_rot*d_p;


    //    Mat9 temp_R_1 = jacobianYX(state_rot*preint.delta_R);
    //    Mat9_3 j_exp = jacobianExpMap(delta_R_correction);
    //    Mat9_3 temp_R_2 = temp_R_1*j_exp;
    //    Mat9_3 d_rot_d_bw = temp_R_2*preint.d_delta_R_d_bw;
    //    Vec9 d_rot_d_dt = temp_R_2*preint.d_delta_R_d_t;
    //    Mat9 d_rot_d_rot = jacobianXY(d_r);

    //    
    //    Mat3_9 d_vel_d_rot = jacobianXv(d_v);
    //    Mat3 d_vel_d_vel = Mat3::Identity();
    //    Mat3 d_vel_d_bf = state_rot*preint.d_delta_v_d_bf;
    //    Mat3 d_vel_d_bw = state_rot*preint.d_delta_v_d_bw;
    //    Vec3 d_vel_d_dt = state_rot*preint.d_delta_v_d_t;
    //    Mat3 d_vel_d_g = Mat3::Identity()*preint.dt;
    //    Mat3_9 d_p_d_rot = jacobianXv(d_p);
    //    Mat3 d_p_d_vel = Mat3::Identity()*preint.dt;
    //    Mat3 d_p_d_p = Mat3::Identity();
    //    Mat3 d_p_d_bf = state_rot*preint.d_delta_p_d_bf;
    //    Mat3 d_p_d_bw = state_rot*preint.d_delta_p_d_bw;
    //    Vec3 d_p_d_dt = state_rot*preint.d_delta_p_d_t;
    //    Mat3 d_p_d_g = Mat3::Identity()*preint.dt_sq_half;



    //    return {d_rot_d_rot, d_rot_d_bw, d_rot_d_dt, d_vel_d_rot, d_vel_d_vel, d_vel_d_bf, d_vel_d_bw, d_vel_d_dt, d_vel_d_g, d_p_d_rot, d_p_d_vel, d_p_d_p, d_p_d_bf, d_p_d_bw, d_p_d_dt, d_p_d_g};
    //}

    //// Output the Jacobian as
    //// d_rot_d_rot, d_rot_d_bw, d_rot_d_dt, d_vel_d_rot, d_vel_d_vel, d_vel_d_bf, d_vel_d_bw, d_vel_d_dt, d_vel_d_g, d_p_d_rot, d_p_d_vel, d_p_d_p, d_p_d_bf, d_p_d_bw, d_p_d_dt, d_p_d_g
    //inline std::tuple<Mat9, Mat9_3, Vec9, Mat3_9, Mat3, Mat3, Mat3, Vec3, Mat3, Mat3_9, Mat3, Mat3, Mat3, Mat3, Vec3, Mat3> preintPredictAndJacobian(const PreintMeas& preint, const FrameState& state, const Vec3& gravity, Mat3& out_rot, Vec3& out_pos, Vec3& out_vel)
    //{
    //    const Mat3 state_rot = state.rot();
    //    const Vec3 state_pos = state.pos();
    //    const Vec3 state_vel = state.vel();
    //    const Vec3 state_bf = state.accBias();
    //    const Vec3 state_bw = state.gyrBias();
    //    const double state_dt = state.d_t[0];
    //    return preintPredictAndJacobian(preint, state_rot, state_pos, state_vel, state_bf, state_bw, state_dt, gravity, out_rot, out_pos, out_vel);
    //}

    //inline double distToPlane(const Vec3& point, const Vec3& p_0, const Vec3& p_1, const Vec3& p_2)
    //{

    //    Vec3 cross_prod = (p_0 - p_1).cross(p_0 - p_2);
    //    double cross_norm = cross_prod.norm();
    //    return (point - p_0).dot(cross_prod) / cross_norm;
    //}


    //inline double distToLine(const Vec3& point, const Vec3& line_pt_1, const Vec3& line_pt_2)
    //{
    //    Vec3 line = line_pt_2 - line_pt_1;

    //    return ((point - line_pt_1).cross(point - line_pt_2)).norm() / (line.norm());
    //    
    //}

    //inline std::tuple<double, Row3, Row3, Row3> distToLineAndJacobian(const Vec3& point, const Vec3& line_pt_1, const Vec3& line_pt_2)
    //{
    //    Vec3 line = line_pt_2 - line_pt_1;

    //    double distance = ((point - line_pt_1).cross(point - line_pt_2)).norm() / (line.norm());

    //    double line_norm = line.norm();
    //    double line_norm_sq = line_norm*line_norm;
    //    Vec3 cross_prod = (point - line_pt_1).cross(point - line_pt_2);
    //    double cross_norm = cross_prod.norm();
    //    cross_prod = cross_prod/cross_norm;

    //    Row3 j_point = (cross_prod.transpose()*toSkewSymMat(line)*line_norm)/line_norm_sq;
    //    Row3 j_l_1 = (cross_prod.transpose()*toSkewSymMat(point-line_pt_2)*line_norm + (line.transpose()/line_norm)*cross_norm)/line_norm_sq;
    //    Row3 j_l_2 = (cross_prod.transpose()*toSkewSymMat(line_pt_1 - point)*line_norm - (line.transpose()/line_norm)*cross_norm)/line_norm_sq;


    //    return {distance, j_point, j_l_1, j_l_2};
    //}

    //// Transform a xyz point into polar coordinate: azimuth, elevation, range
    //inline std::tuple<double, double, double> xyzToPolar(Vec3 p)
    //{
    //    double r = p.norm();
    //    double az = std::atan2(p(1), p(0));
    //    double el = std::asin(p(2)/r);
    //    return {az, el, r};
    //}

    //// Propagate noise on polar coordinates to xyz
    //inline Mat3 polarCovToXyzCov(Mat3 polar_cov, double az, double el, double r)
    //{
    //    Mat3 j_polar_to_xyz;
    //    j_polar_to_xyz <<  
    //        -r*cos(el)*std::sin(az), -r*cos(az)*std::sin(el),  cos(az)*cos(el),
    //        r*std::cos(az)*std::cos(el), -r*std::sin(az)*std::sin(el), std::cos(el)*std::sin(az),
    //        0,               r*std::cos(el),                   sin(el);

    //    return j_polar_to_xyz*polar_cov*j_polar_to_xyz.transpose();
    //}    

    //inline std::tuple<Vec3, Mat3, Vec3> eigenValueAndVectorsAndMean(const Mat3X& pts, const int id_start, const int id_finish)
    //{
    //    Vec3 mean = Vec3::Zero();
    //    for(int j = id_start; j <= id_finish; ++j)
    //    {
    //        mean += pts.col(j);
    //    }
    //    mean = (1.0/(double)(id_finish-id_start+1)) * mean;

    //    Mat3 cov = Mat3::Zero();
    //    for(int j = id_start; j <= id_finish; ++j)
    //    {
    //        Vec3 temp = pts.col(j) - mean;
    //        cov += temp*(temp.transpose());
    //    }
    //    cov = (1.0/(double)(id_finish-id_start))*cov;
    //    Eigen::SelfAdjointEigenSolver<Mat3> eig(cov);
    //    
    //    return {eig.eigenvalues(), eig.eigenvectors(), mean};
    //}
    //inline std::pair<Vec3, Mat3> eigenValueAndVectors(const Mat3X& pts, const int id_start, const int id_finish)
    //{
    //    auto [eig_val, eig_vec, mean] = eigenValueAndVectorsAndMean(pts, id_start, id_finish);
    //    return {eig_val, eig_vec};
    //}



    //inline celib::Mat3 orientationPrior(const PreintMeas& preint)
    //{
    //    Vec3 mean_acc = preint.delta_v / preint.delta_v.norm();
    //    Vec3 ideal_g(0,0,1);
    //    Mat3 new_rot;
    //    if(ideal_g == mean_acc){
    //        new_rot = Mat3::Identity();
    //    }else if(ideal_g == -mean_acc){
    //        new_rot = -Mat3::Identity();
    //    }else{
    //        Vec3 cross_vec = mean_acc.cross(ideal_g);
    //        double dot_prod = mean_acc.dot(ideal_g);
    //        Mat3 temp_mat = toSkewSymMat(cross_vec);
    //        new_rot = Mat3::Identity() + temp_mat + (temp_mat * temp_mat / (1 + dot_prod));
    //    }
    //    return new_rot;
    //}


    //inline celib::Vec3 addN2Pi(const celib::Vec3& r, const int n)
    //{
    //    double norm_r = r.norm();
    //    if(norm_r != 0)
    //    {
    //        celib::Vec3 unit_r = r/norm_r;
    //        return unit_r*(2.0*M_PI*n + norm_r);
    //    }
    //    else
    //    {
    //        return r;
    //    }
    //}

    //inline std::pair<celib::Vec3, int> getClosest(const celib::Vec3& t, const std::vector<celib::Vec3> s)
    //{
    //    int id_min = 0;
    //    double dist_min = std::numeric_limits<double>::max();
    //    for(int i = 0; i < s.size(); ++i)
    //    {
    //        if((t - s[i]).norm() < dist_min)
    //        {
    //            dist_min = (t - s[i]).norm();
    //            id_min = i;
    //        }
    //    }
    //    return {s[id_min], id_min};
    //}



// F//OR DEBUG: To be removed eventyally !
    //inline void testpreintPredicJacobian()
    //{
    //    PreintMeas preint;
    //    preint.delta_R = expMap(Vec3::Random());
    //    preint.delta_v = Vec3::Random();
    //    preint.delta_p = Vec3::Random();
    //    preint.d_delta_R_d_bw = Mat3::Random();
    //    preint.d_delta_R_d_t = Vec3::Random();
    //    preint.d_delta_v_d_bw = Mat3::Random();
    //    preint.d_delta_v_d_bf = Mat3::Random();
    //    preint.d_delta_v_d_t = Vec3::Random();
    //    preint.d_delta_p_d_bw = Mat3::Random();
    //    preint.d_delta_p_d_bf = Mat3::Random();
    //    preint.d_delta_p_d_t = Vec3::Random();
    //    preint.dt = Vec3::Random()[0];
    //    preint.dt_sq_half = preint.dt*preint.dt*0.5;

    //    Vec3 gravity = Vec3::Random();

    //    FrameState state;
    //    Eigen::Map<Mat3> state_rot(state.pose.data());
    //    Eigen::Map<Vec3> state_pos(state.pose.data()+9);
    //    Eigen::Map<Vec3> state_vel(state.v.data());
    //    Eigen::Map<Vec3> state_bf(state.b_f.data());
    //    Eigen::Map<Vec3> state_bw(state.b_w.data());
    //    state.d_t[0] = (Vec3::Random())[0];
    //    state_rot = Mat3::Random();
    //    state_pos = Vec3::Random();
    //    state_vel = Vec3::Random();
    //    state_bf = Vec3::Random();
    //    state_bw = Vec3::Random();


    //    Mat3 rot;
    //    Vec3 vel, pos;


    //    auto [d_rot_d_rot, d_rot_d_bw, d_rot_d_dt, d_vel_d_rot, d_vel_d_vel, d_vel_d_bf, d_vel_d_bw, d_vel_d_dt, d_vel_d_g, d_p_d_rot, d_p_d_vel, d_p_d_p, d_p_d_bf, d_p_d_bw, d_p_d_dt, d_p_d_g] = preintPredictAndJacobian(preint, state, gravity, rot, pos, vel);

    //    Mat9 num_d_rot_d_rot;
    //    Mat9_3 num_d_rot_d_bw;
    //    Vec9 num_d_rot_d_dt;
    //    Mat3_9 num_d_vel_d_rot;
    //    Mat3 num_d_vel_d_vel;
    //    Mat3 num_d_vel_d_bf;
    //    Mat3 num_d_vel_d_bw;
    //    Vec3 num_d_vel_d_dt;
    //    Mat3 num_d_vel_d_g;
    //    Mat3_9 num_d_p_d_rot;
    //    Mat3 num_d_p_d_vel;
    //    Mat3 num_d_p_d_p;
    //    Mat3 num_d_p_d_bf;
    //    Mat3 num_d_p_d_bw;
    //    Vec3 num_d_p_d_dt;
    //    Mat3 num_d_p_d_g;

    //    double quantum = 0.001;

    //    for(int i = 0; i < 9; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.pose[i] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_rot_d_rot.col(i) = (mat3ToRow(rot_dx) - mat3ToRow(rot)).transpose()/quantum;
    //        num_d_vel_d_rot.col(i) = (vel_dx - vel)/quantum;
    //        num_d_p_d_rot.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.pose[i+9] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_p_d_p.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.v[i] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_vel_d_vel.col(i) = (vel_dx - vel)/quantum;
    //        num_d_p_d_vel.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.b_f[i] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_vel_d_bf.col(i) = (vel_dx - vel)/quantum;
    //        num_d_p_d_bf.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.b_w[i] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_rot_d_bw.col(i) = (mat3ToRow(rot_dx) - mat3ToRow(rot)).transpose()/quantum;
    //        num_d_vel_d_bw.col(i) = (vel_dx - vel)/quantum;
    //        num_d_p_d_bw.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        Vec3 g_dx = gravity;
    //        g_dx[i] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state, g_dx, rot_dx, pos_dx, vel_dx);
    //        num_d_vel_d_g.col(i) = (vel_dx - vel)/quantum;
    //        num_d_p_d_g.col(i) = (pos_dx - pos)/quantum;
    //    }
    //    {
    //        FrameState state_dx = state;
    //        state_dx.d_t[0] += quantum;
    //        Mat3 rot_dx;
    //        Vec3 vel_dx, pos_dx;

    //        preintPredict(preint, state_dx, gravity, rot_dx, pos_dx, vel_dx);
    //        num_d_rot_d_dt = (mat3ToRow(rot_dx) - mat3ToRow(rot)).transpose()/quantum;
    //        num_d_vel_d_dt = (vel_dx - vel)/quantum;
    //        num_d_p_d_dt = (pos_dx - pos)/quantum;
    //    }

    //    std::cout << "------------------------ Testing Jacobians for preint prediction" << std::endl;
    //    std::cout << "d_rot_d_rot"
    //        << std::endl << d_rot_d_rot << std::endl
    //        << std::endl << num_d_rot_d_rot << std::endl << std::endl;
    //    std::cout << "d_rot_d_bw"
    //        << std::endl << d_rot_d_bw << std::endl
    //        << std::endl << num_d_rot_d_bw << std::endl << std::endl;
    //    std::cout << "d_rot_d_dt"
    //        << std::endl << d_rot_d_dt << std::endl
    //        << std::endl << num_d_rot_d_dt << std::endl << std::endl;
    //    std::cout << "d_vel_d_rot"
    //        << std::endl << d_vel_d_rot << std::endl
    //        << std::endl << num_d_vel_d_rot << std::endl << std::endl;
    //    std::cout << "d_vel_d_vel"
    //        << std::endl << d_vel_d_vel << std::endl
    //        << std::endl << num_d_vel_d_vel << std::endl << std::endl;
    //    std::cout << "d_vel_d_bf"
    //        << std::endl << d_vel_d_bf << std::endl
    //        << std::endl << num_d_vel_d_bf << std::endl << std::endl;
    //    std::cout << "d_vel_d_bw"
    //        << std::endl << d_vel_d_bw << std::endl
    //        << std::endl << num_d_vel_d_bw << std::endl << std::endl;
    //    std::cout << "d_vel_d_dt"
    //        << std::endl << d_vel_d_dt << std::endl
    //        << std::endl << num_d_vel_d_dt << std::endl << std::endl;
    //    std::cout << "d_vel_d_g"
    //        << std::endl << d_vel_d_g << std::endl
    //        << std::endl << num_d_vel_d_g << std::endl << std::endl;
    //    std::cout << "d_p_d_rot"
    //        << std::endl << d_p_d_rot << std::endl
    //        << std::endl << num_d_p_d_rot << std::endl << std::endl;
    //    std::cout << "d_p_d_vel"
    //        << std::endl << d_p_d_vel << std::endl
    //        << std::endl << num_d_p_d_vel << std::endl << std::endl;
    //    std::cout << "d_p_d_p"
    //        << std::endl << d_p_d_p << std::endl
    //        << std::endl << num_d_p_d_p << std::endl << std::endl;
    //    std::cout << "d_p_d_bf"
    //        << std::endl << d_p_d_bf << std::endl
    //        << std::endl << num_d_p_d_bf << std::endl << std::endl;
    //    std::cout << "d_p_d_bw"
    //        << std::endl << d_p_d_bw << std::endl
    //        << std::endl << num_d_p_d_bw << std::endl << std::endl;
    //    std::cout << "d_p_d_dt"
    //        << std::endl << d_p_d_dt << std::endl
    //        << std::endl << num_d_p_d_dt << std::endl << std::endl;
    //    std::cout << "d_p_d_g"
    //        << std::endl << d_p_d_g << std::endl
    //        << std::endl << num_d_p_d_g << std::endl << std::endl;

    //}


// F//OR DEBUG: To be removed eventyally !
    //inline void testPointPreintTransformJacobian()
    //{
    //    PreintMeas preint;
    //    preint.delta_R = expMap(Vec3::Random());
    //    preint.delta_v = Vec3::Random();
    //    preint.delta_p = Vec3::Random();
    //    preint.d_delta_R_d_bw = Mat3::Random();
    //    preint.d_delta_R_d_t = Vec3::Random();
    //    preint.d_delta_v_d_bw = Mat3::Random();
    //    preint.d_delta_v_d_bf = Mat3::Random();
    //    preint.d_delta_v_d_t = Vec3::Random();
    //    preint.d_delta_p_d_bw = Mat3::Random();
    //    preint.d_delta_p_d_bf = Mat3::Random();
    //    preint.d_delta_p_d_t = Vec3::Random();
    //    preint.dt = Vec3::Random()[0];
    //    preint.dt_sq_half = preint.dt*preint.dt*0.5;

    //    std::array<double,3> gravity;
    //    Eigen::Map<Vec3> temp_g(gravity.data());
    //    temp_g = Vec3::Random();

    //    FrameState state;
    //    Eigen::Map<Mat3> state_rot(state.pose.data());
    //    Eigen::Map<Vec3> state_pos(state.pose.data()+9);
    //    Eigen::Map<Vec3> state_vel(state.v.data());
    //    Eigen::Map<Vec3> state_bf(state.b_f.data());
    //    Eigen::Map<Vec3> state_bw(state.b_w.data());
    //    state.d_t[0] = (Vec3::Random())[0];
    //    state_rot = Mat3::Random();
    //    state_pos = Vec3::Random();
    //    state_vel = Vec3::Random();
    //    state_bf = Vec3::Random();
    //    state_bw = Vec3::Random();
    //    std::array<double, 12> calib;
    //    Eigen::Map<Mat3> state_calib_rot(calib.data());
    //    Eigen::Map<Vec3> state_calib_pos(calib.data()+9);
    //    state_calib_rot = expMap(Vec3::Random());
    //    state_calib_pos = Vec3::Random();


    //    PointPreint point_preint;
    //    point_preint.p = Vec3::Identity();
    //    point_preint.preint = std::shared_ptr<PreintMeas>(new PreintMeas);
    //    *point_preint.preint = preint;


    //    Mat3 rot;
    //    Vec3 vel, pos;

    //    auto [point, d_x_d_R, d_x_d_v, d_x_d_p, d_x_d_bf, d_x_d_bw, d_x_d_dt, d_x_d_gravity, d_x_d_R_calib, d_x_d_p_calib] = point_preint.transformAndJacobian(state, gravity, calib);

    //    Mat3_9 num_d_x_d_R;
    //    Mat3 num_d_x_d_v;
    //    Mat3 num_d_x_d_p;
    //    Mat3 num_d_x_d_bf;
    //    Mat3 num_d_x_d_bw;
    //    Vec3 num_d_x_d_dt;
    //    Mat3 num_d_x_d_gravity;
    //    Mat3_9 num_d_x_d_R_calib;
    //    Mat3 num_d_x_d_p_calib;


    //    double quantum = 0.001;

    //    for(int i = 0; i < 9; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.pose[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_R.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.pose[9+i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_p.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.v[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_v.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.b_f[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_bf.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        FrameState state_dx = state;
    //        state_dx.b_w[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_bw.col(i) = (point_dx - point)/quantum;
    //    }
    //    {
    //        FrameState state_dx = state;
    //        state_dx.d_t[0] += quantum;

    //        Vec3 point_dx = point_preint.transform(state_dx, gravity, calib); 
    //        num_d_x_d_dt = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        std::array<double, 3> g_dx = gravity;
    //        g_dx[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state, g_dx, calib); 
    //        num_d_x_d_gravity.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 9; ++i)
    //    {
    //        std::array<double, 12> calib_dx = calib;
    //        calib_dx[i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state, gravity, calib_dx); 
    //        num_d_x_d_R_calib.col(i) = (point_dx - point)/quantum;
    //    }
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        std::array<double, 12> calib_dx = calib;
    //        calib_dx[9+i] += quantum;

    //        Vec3 point_dx = point_preint.transform(state, gravity, calib_dx); 
    //        num_d_x_d_p_calib.col(i) = (point_dx - point)/quantum;
    //    }

    //    std::cout << "------------------------ Testing Jacobians for point transform" << std::endl;
    //    std::cout << "num_d_x_d_R"
    //            << std::endl << d_x_d_R << std::endl
    //            << std::endl << num_d_x_d_R << std::endl << std::endl;
    //    std::cout << "num_d_x_d_v"
    //            << std::endl << d_x_d_v << std::endl
    //            << std::endl << num_d_x_d_v << std::endl << std::endl;
    //    std::cout << "num_d_x_d_p"
    //            << std::endl << d_x_d_p << std::endl
    //            << std::endl << num_d_x_d_p << std::endl << std::endl;
    //    std::cout << "num_d_x_d_bf"
    //            << std::endl << d_x_d_bf << std::endl
    //            << std::endl << num_d_x_d_bf << std::endl << std::endl;
    //    std::cout << "num_d_x_d_bw"
    //            << std::endl << d_x_d_bw << std::endl
    //            << std::endl << num_d_x_d_bw << std::endl << std::endl;
    //    std::cout << "num_d_x_d_dt"
    //            << std::endl << d_x_d_dt << std::endl
    //            << std::endl << num_d_x_d_dt << std::endl << std::endl;
    //    std::cout << "num_d_x_d_gravity"
    //            << std::endl << d_x_d_gravity << std::endl
    //            << std::endl << num_d_x_d_gravity << std::endl << std::endl;
    //    std::cout << "9 num_d_x_d_R_calib"
    //            << std::endl << d_x_d_R_calib << std::endl
    //            << std::endl << num_d_x_d_R_calib << std::endl << std::endl;
    //    std::cout << "num_d_x_d_p_calib"
    //            << std::endl << d_x_d_p_calib << std::endl
    //            << std::endl << num_d_x_d_p_calib << std::endl << std::endl;

    //}

// F//OR DEBUG: To be removed eventyally !
    //inline void testResidualJacobians()
    //{
    //    // Case of plane
    //    {
    //        std::vector<Vec3, Eigen::aligned_allocator<Vec3> > current_point(1);
    //        current_point[0] = Vec3::Random();
    //        std::vector<std::vector<Vec3, Eigen::aligned_allocator<Vec3> > > points(3, std::vector<Vec3, Eigen::aligned_allocator<Vec3> >(1));
    //        points[0][0] = Vec3::Random();
    //        points[1][0] = Vec3::Random();
    //        points[2][0] = Vec3::Random();
    //        
    //        LidarMatch match;
    //        match.type = LidarFeatureType::kPlane;

    //        auto [ residual, j_current, j_points] = match.residualAndJacobian(current_point, points);

    //        double quantum = 0.001;

    //        std::cout << "----------- Plane residual jacobians" << std::endl;
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto current_dx = current_point;
    //                current_dx[0][i] += quantum;
    //                double residual_dx = match.residual(current_dx, points);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Current " << std::endl << j_current[0] << std::endl << num_j << std::endl << std::endl;
    //        }

    //        for(int k = 0; k < 3; ++k)
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto points_dx = points;
    //                points_dx[k][0][i] += quantum;
    //                double residual_dx = match.residual(current_point, points_dx);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Points " << k << std::endl << j_points[k][0] << std::endl << num_j << std::endl << std::endl;
    //        }
    //    }
    //    // Case of outward edge
    //    {
    //        std::vector<Vec3, Eigen::aligned_allocator<Vec3> > current_point(1);
    //        current_point[0] = Vec3::Random();
    //        std::vector<std::vector<Vec3, Eigen::aligned_allocator<Vec3> > > points(2, std::vector<Vec3, Eigen::aligned_allocator<Vec3> >(1));
    //        points[0][0] = Vec3::Random();
    //        points[1][0] = Vec3::Random();
    //        
    //        LidarMatch match;
    //        match.type = LidarFeatureType::kEdgeOutward;

    //        auto [ residual, j_current, j_points] = match.residualAndJacobian(current_point, points);

    //        double quantum = 0.001;

    //        std::cout << "----------- Outward edge residual jacobians" << std::endl;
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto current_dx = current_point;
    //                current_dx[0][i] += quantum;
    //                double residual_dx = match.residual(current_dx, points);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Current " << std::endl << j_current[0] << std::endl << num_j << std::endl << std::endl;
    //        }

    //        for(int k = 0; k < 2; ++k)
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto points_dx = points;
    //                points_dx[k][0][i] += quantum;
    //                double residual_dx = match.residual(current_point, points_dx);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Points " << k << std::endl << j_points[k][0] << std::endl << num_j << std::endl << std::endl;
    //        }
    //    }
    //    // Case of inward edge
    //    {
    //        std::vector<Vec3, Eigen::aligned_allocator<Vec3> > current_point(1);
    //        current_point[0] = Vec3::Random();
    //        std::vector<std::vector<Vec3, Eigen::aligned_allocator<Vec3> > > points(2, std::vector<Vec3, Eigen::aligned_allocator<Vec3> >(1));
    //        points[0][0] = Vec3::Random();
    //        points[1][0] = Vec3::Random();
    //        
    //        LidarMatch match;
    //        match.type = LidarFeatureType::kEdgeInward;

    //        auto [ residual, j_current, j_points] = match.residualAndJacobian(current_point, points);

    //        double quantum = 0.00001;

    //        std::cout << "----------- Inward edge residual jacobians" << std::endl;
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto current_dx = current_point;
    //                current_dx[0][i] += quantum;
    //                double residual_dx = match.residual(current_dx, points);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Current " << std::endl << j_current[0] << std::endl << num_j << std::endl << std::endl;
    //        }

    //        for(int k = 0; k < 2; ++k)
    //        {
    //            Row3 num_j;
    //            for(int i = 0; i < 3; ++i)
    //            {
    //                auto points_dx = points;
    //                points_dx[k][0][i] += quantum;
    //                double residual_dx = match.residual(current_point, points_dx);
    //                num_j[i] = (residual_dx - residual)/quantum;
    //            }
    //            std::cout << "Points " << k << std::endl << j_points[k][0] << std::endl << num_j << std::endl << std::endl;
    //        }
    //    }
    //    // Case of pillar 
    //    {
    //        std::vector<Vec3, Eigen::aligned_allocator<Vec3> > current_point(2);
    //        current_point[0] = Vec3::Random();
    //        current_point[1] = Vec3::Random();
    //        std::vector<std::vector<Vec3, Eigen::aligned_allocator<Vec3> > > points(3, std::vector<Vec3, Eigen::aligned_allocator<Vec3> >(2));
    //        points[0][0] = Vec3::Random();
    //        points[0][1] = Vec3::Random();
    //        points[1][0] = Vec3::Random();
    //        points[1][1] = Vec3::Random();
    //        points[2][0] = Vec3::Random();
    //        points[2][1] = Vec3::Random();
    //        
    //        LidarMatch match;
    //        match.type = LidarFeatureType::kPillar;

    //        auto [ residual, j_current, j_points] = match.residualAndJacobian(current_point, points);

    //        double quantum = 0.001;

    //        std::cout << "----------- Pillar residual jacobians" << std::endl;
    //        for(int f = 0; f < 2; ++f)
    //        {
    //            {
    //                Row3 num_j;
    //                for(int i = 0; i < 3; ++i)
    //                {
    //                    auto current_dx = current_point;
    //                    current_dx[f][i] += quantum;
    //                    double residual_dx = match.residual(current_dx, points);
    //                    num_j[i] = (residual_dx - residual)/quantum;
    //                }
    //                std::cout << "Current " << f << std::endl << j_current[f] << std::endl << num_j << std::endl << std::endl;
    //            }

    //            for(int k = 0; k < 2; ++k)
    //            {
    //                Row3 num_j;
    //                for(int i = 0; i < 3; ++i)
    //                {
    //                    auto points_dx = points;
    //                    points_dx[k][f][i] += quantum;
    //                    double residual_dx = match.residual(current_point, points_dx);
    //                    num_j[i] = (residual_dx - residual)/quantum;
    //                }
    //                std::cout << "Points " << k << "  " << f << std::endl << j_points[k][f] << std::endl << num_j << std::endl << std::endl;
    //            }
    //        }
    //    }
    //}


    //inline void testNormalisationJacobian()
    //{
    //    Vec3 vec = Vec3::Random();
    //    Vec3 unit_vec = vec / (vec.norm());

    //    Mat3 jacobian = jacobianNormalisation(vec);

    //    Mat3 num_jacobian = Mat3::Zero();
    //    double quantum = 0.00001;
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        Vec3 temp_vec = vec;
    //        temp_vec(i) += quantum;
    //        temp_vec = temp_vec / (temp_vec.norm());
    //        num_jacobian.col(i) = (temp_vec - unit_vec)/quantum;
    //    }

    //    std::cout << "Numerical jacobian of Normalisation " << std::endl << num_jacobian << std::endl;
    //    std::cout << "Analitical jacobian of Normalisation " << std::endl << jacobian << std::endl;

    //}
 

    //inline void testNormJacobian()
    //{
    //    Vec3 vec = Vec3::Random();
    //    

    //    Row3 jacobian = jacobianNorm(vec);

    //    Row3 num_jacobian;
    //    double quantum = 0.00001;
    //    for(int i = 0; i < 3; ++i)
    //    {
    //        Vec3 temp_vec = vec;
    //        temp_vec(i) += quantum;
    //        
    //        num_jacobian(i) = ((temp_vec.norm()) - (vec.norm()) )/quantum;
    //    }

    //    std::cout << "Numerical jacobian of Norm " << std::endl << num_jacobian << std::endl;
    //    std::cout << "Analitical jacobian of Norm " << std::endl << jacobian << std::endl;

    //}

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