/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/
#ifndef GP_STATE_MANAGER_H
#define GP_STATE_MANAGER_H


#include "common/types.h"
#include "ceres/ceres.h"

namespace celib
{
    const int kLengthScaleFactor = 10;
    const int kBiasLengthScaleFactor = 20;
    const int kCutoffFactor = 3;

    enum GpStateType {kSe2, kHomography};

    typedef std::vector<int> IndexVec;

    typedef std::vector<double> SubState;
    typedef std::shared_ptr<SubState> SubStatePtr;


    class CostFunctionGPNormaliser: public ceres::CostFunction
    {
        private:
            int state_id_;
            int nb_states_;
            int size_;
            VecX jacobian_;
            VecX ks_Kinv_;

        public:
            CostFunctionGPNormaliser(const double* state, const std::vector<double*>& normalise_state, const std::vector<int> sizes, const VecX& ksKinv);
                


            // Cost function and jacobians
            virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const;
    };

    class GpStateManager
    {
        public:
            GpStateManager( const GpStateType type, const double state_frequency);

            void setStartingTime(const double t);


            // For poses, return the pointers to the state variables, the id mappings, and the inference vectors
            std::tuple<
                std::vector<double*>,   // Pointers to the states
                std::vector<int>,       // State sizes
                std::vector<IndexVec>,  // Pos ID mapping of the timestamps to the ptrs
                std::vector<IndexVec>,  // Pos ID mapping of the timestamps to the ptrs
                std::vector<VecX>       // ksKinv
                > getStatePosePtrs(const std::vector<double>& time, bool query_only = false);


            // For homography, return the pointers to the state variables, the id mappings, and the inference vectors
            std::tuple<
                std::vector<double*>,   // Pointers to the states
                std::vector<IndexVec>,  // Id of the states for each of the timestamps  
                std::vector<VecX>       // ksKinv
                > getStateHomographyPtrs(const std::vector<double>& time, bool query_only = false);

            // Returns pointers to the closest set of states
            std::vector<double*> getClosestStatePtrs(const double time);
            std::vector<double*> getStatePtrsUpTo(const double time);

            std::pair<VecX, VecX> getPoseRotVec(const double time);
            Vec9 getHomographyVec(const double time);
            Mat3 getHomography(const double time);

            void plotPose(const std::string window_prefix);
            void plotHomography(const std::string window_prefix);


            void addNormalisersToProblem(ceres::Problem& problem);
            void resetUsedFlags();

            void setConstantTill(ceres::Problem& problem, const double time);
            void setConstantCloseTo(ceres::Problem& problem, const double time) const;

        private:
            GpStateType type_;

            HyperParam hyper_;
            HyperParam bias_hyper_;
            double state_period_;
            double cutoff_;
            double cutoff_state_;

            std::map<int, MatX> K_inv_;

            std::vector<SubStatePtr> homography_;
            std::vector<SubStatePtr> rot_;
            std::vector<SubStatePtr> pos_;
            std::vector<SubStatePtr> acc_bias_;
            std::vector<SubStatePtr> gyr_bias_;
            std::vector<double> time_;
            Buffer<int> time_buffer_;

            std::vector<bool> used_;


            // Will add states and timestamps so that the last state is in the cutoff from the last timestamp provided.
            // times need to be sorted
            void updateTime( const std::vector<double>& times);
            
    };



    inline VecX ceresParamToEigen(double const* const* parameters, const std::vector<int>& ids, const VecX& ks_K_inv, const int size)
    {
        if(ids.size() != ks_K_inv.size()) throw std::invalid_argument("ceresParamToEigen(gp_state_manager.h): The number of ids do not match the size of the ks_K_inv vector");

        MatX state(ids.size(), size);
        for(int s = 0; s < ids.size(); ++s)
        {
            const Eigen::Map<const RowX> line(parameters[ids.at(s)], size);
            state.row(s) = line;
        }

        return (ks_K_inv.transpose() * state).transpose();
    }

    inline VecX ceresParamToEigen(double const* const* parameters, const VecX& ks_K_inv, const int size)
    {
        MatX state(ks_K_inv.size(), size);
        for(int s = 0; s < ks_K_inv.size(); ++s)
        {
            const Eigen::Map<const RowX> line(parameters[s], size);
            state.row(s) = line;
        }

        return (ks_K_inv.transpose() * state).transpose();
    }




} // namespace celib





#endif