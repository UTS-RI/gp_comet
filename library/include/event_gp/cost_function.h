/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef EVENT_GP_COST_FUNCTIONS_H
#define EVENT_GP_COST_FUNCTIONS_H

#include <ceres/ceres.h>


#include "common/types.h"
#include "common/gp_state_manager.h"
#include "common/math_utils.h"

namespace celib 
{

    typedef std::shared_ptr<Mat2X> Mat2XPtr;
    typedef std::shared_ptr<VecX> VecXPtr;



    // TODO To make more efficient:
    // . Investigate the use of KISS GP
    class CostFunctionSe2FrontEnd: public ceres::CostFunction
    {
        private:
            Mat2X events_xy_;
            Mat2X fixed_events_;
            std::vector<IndexVec> rot_ids_;
            std::vector<IndexVec> pos_ids_;
            std::vector<VecX> ks_K_inv_;

            std::vector<IndexVec> state_to_event_;
            std::vector<int> state_size_;
        
            HyperParam xy_hyper_;


        public:
            CostFunctionSe2FrontEnd(
                    const std::vector<EventPtr>& events,
                    const std::vector<int>& event_ids,
                    const std::vector<EventPtr>& fixed_events,
                    GpStateManager& state_manager,
                    const HyperParam& hyper,
                    std::vector<double*>& state_ptrs,
                    const bool filter = false,
                    const int downsample = 0
                    );

            virtual bool Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const;

    };



    class CostFunctionSe2ApproxFrontEnd: public ceres::CostFunction
    {
        private:
            Mat2X events_xy_;
            std::vector<IndexVec> rot_ids_;
            std::vector<IndexVec> pos_ids_;
            std::vector<VecX> ks_K_inv_;

            std::vector<IndexVec> state_to_event_;
            std::vector<int> state_size_;
        
            HyperParam xy_hyper_;


        public:
            CostFunctionSe2ApproxFrontEnd(
                    const std::vector<EventPtr>& events,
                    const std::vector<int>& event_ids,
                    GpStateManager& state_manager,
                    const HyperParam& hyper,
                    std::vector<double*>& state_ptrs,
                    const bool filter = false,
                    const int downsample = 0
                    );

            virtual bool Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const;

            bool testJacobian();

    };





    class CostFunctionDiscreteHomographyFrontEnd: public ceres::CostFunction
    {
        private:
            const Vec2 event_xy_;
            const Mat2XPtr template_;
            const VecXPtr alpha_;
            const double weight_;
        
            HyperParam xy_hyper_;


        public:
            CostFunctionDiscreteHomographyFrontEnd(
                    const EventPtr& event,
                    const Mat2XPtr& template_events,
                    const VecXPtr& alpha,
                    const HyperParam& hyper,
                    const double weight = 1.0
                    );

            virtual bool Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const;

    };


    class CostFunctionDiscreteInvHomographyFrontEnd: public ceres::CostFunction
    {
        private:
            const Vec2 event_xy_;
            const Mat2XPtr template_;
            const VecXPtr alpha_;
            const double weight_;
        
            HyperParam xy_hyper_;


        public:
            CostFunctionDiscreteInvHomographyFrontEnd(
                    const EventPtr& event,
                    const Mat2XPtr& template_events,
                    const VecXPtr& alpha,
                    const HyperParam& hyper,
                    const double weight = 1.0
                    );

            virtual bool Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const;

    };


} //namespace
#endif