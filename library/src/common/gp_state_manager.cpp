/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#include "common/gp_state_manager.h"
#include "common/math_utils.h"

#include "plot/plotter.hpp"


namespace celib
{


    CostFunctionGPNormaliser::CostFunctionGPNormaliser(const double* state, const std::vector<double*>& normalise_state, const std::vector<int> sizes, const VecX& ksKinv)
        :
        size_(sizes.at(0))
        , nb_states_(normalise_state.size())
    {

        // Find the id of the state to be "normalised"
        bool found = false;
        for(int i = 0; (!found) && (i < nb_states_); ++i)
        {
            if(normalise_state.at(i) == state)
            {
                found = true;
                state_id_ = i;
            }
        }
        if(!found) throw std::invalid_argument("CostFunctionGPNormaliser (Constructor): The state does not seem to be in the list of states");


        // Precompute the jacobian and ksKinv
        ks_Kinv_ = ksKinv;
        jacobian_ = ksKinv;
        jacobian_(state_id_) -= 1.0;

        // Fill up the information for ceres
        set_num_residuals(size_);
        std::vector<int>* block_sizes = mutable_parameter_block_sizes();
        for(const auto& s: sizes) block_sizes->push_back(s);
    }


    // Cost function and jacobians
    bool CostFunctionGPNormaliser::Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
    {
        // Map a VecX to the ceres residuals
        Eigen::Map<VecX> res(residuals, size_, 1);

        // Create a MatX form the ceres state variables (parameters)
        MatX states(nb_states_,size_);
        for(int r = 0; r < nb_states_; ++r)
            for(int c = 0; c < size_; ++c)
                states(r,c) = parameters[r][c];

        // Get the state to normalise
        VecX state(size_);
        for(int r = 0; r < size_; ++r)
            state(r) = parameters[state_id_][r];

        
        // Compute the inference
        VecX inf = (ks_Kinv_.transpose()*states).transpose();
        res = 1000.0*(inf - state);

        if(jacobians != NULL)
        {
            for(int b = 0; b < nb_states_; ++b)
            {
                if(jacobians[b] != NULL)
                {
                    Eigen::Map<MatX> j_s(&(jacobians[b][0]),size_, size_);
                    j_s.setZero();
                    for(int r = 0; r < size_; ++r) j_s(r,r) = jacobian_(b);
                    j_s *= 1000.0;
                }
            }
        }
        return true;
    }





















    
    GpStateManager::GpStateManager( const GpStateType type, const double state_frequency ):
        type_(type)
        , state_period_(1.0/state_frequency)
    {
        hyper_.l2 = std::pow(kLengthScaleFactor*state_period_, 2);
        hyper_.sf2 = 1.0;
        hyper_.sz2 = 0.0;

        bias_hyper_.l2 = std::pow(kBiasLengthScaleFactor*state_period_, 2);
        bias_hyper_.sf2 = 1.0;
        bias_hyper_.sz2 = 0.0;

        cutoff_ = kLengthScaleFactor*kCutoffFactor*state_period_;
        cutoff_state_ = kLengthScaleFactor*state_period_;

        for(int i = 2; i <= 2*kLengthScaleFactor*kCutoffFactor+1; ++i)
        {
            VecX x_temp = VecX::LinSpaced(i, 0, (i-1)*state_period_);

            K_inv_.emplace(i, (seKernel(x_temp, x_temp, hyper_.l2, hyper_.sf2)+(0.00001*MatX::Identity(x_temp.size(),x_temp.size()))).inverse());
            //K_inv_.emplace(i, (seKernel(x_temp, x_temp, hyper_.l2, hyper_.sf2)).inverse());
        }
    }


    void GpStateManager::setStartingTime(const double time)
    {
        
//        for( int i = -(kLengthScaleFactor); i <= (kLengthScaleFactor); ++i)
        for( int i = 0; i <= 1; ++i)
        {
            time_buffer_.emplace(time  + (i*state_period_), time_.size());
            time_.push_back(time + (i*state_period_));
            if(type_ == GpStateType::kSe2)
            {
                rot_.push_back(SubStatePtr(new SubState(1,0)));
                pos_.push_back(SubStatePtr(new SubState(2,0)));
            }
            else if(type_ == GpStateType::kHomography)
            {
                homography_.push_back(SubStatePtr(new SubState(8,0)));
                homography_.back()->at(0) = 1;
                homography_.back()->at(4) = 1;
            }
            used_.push_back(false);
        }
    }


    void GpStateManager::updateTime( const std::vector<double>& times)
    {
//        while( time_.back() < (times.back() + cutoff_state_) )
        while( time_.back() < (times.back()) )
        {
            time_buffer_.emplace(time_.back() + state_period_, time_.size());
            time_.push_back(time_.back() + state_period_);

            if(type_ == GpStateType::kHomography)
            {
                std::vector<double> temp_homography = *(homography_.back());
                homography_.push_back(SubStatePtr(new SubState));
                *(homography_.back()) = temp_homography;
            }
            else
            {
                std::vector<double> temp_rot = *(rot_.back());
                std::vector<double> temp_pos = *(pos_.back());
                rot_.push_back(SubStatePtr(new SubState));
                pos_.push_back(SubStatePtr(new SubState));
                *rot_.back() = temp_rot;
                *pos_.back() = temp_pos;
            }
            used_.push_back(false);
        }
    }






    
    std::tuple<
        std::vector<double*>,   // Pointers to the states
        std::vector<int>,       // State sizes
        std::vector<IndexVec>,  // Rot ID mapping of the timestamps to the ptrs
        std::vector<IndexVec>,  // Pos ID mapping of the timestamps to the ptrs
        std::vector<VecX>       // ksKinv
        > GpStateManager::getStatePosePtrs(const std::vector<double>& time, bool query_only)
    {
        if( type_ == GpStateType::kHomography) throw std::range_error("GpStateManager::getPosRotVec: Cannot retrieve pose pointers for homography state");



        if(!query_only)  updateTime(time);



        std::vector<double*> ptrs;
        std::vector<int> sizes;
        std::vector<IndexVec> rot_id(time.size());
        std::vector<IndexVec> pos_id(time.size());
        std::vector<VecX> ks_K_inv(time.size());



        std::vector<std::vector<double> > state_times(time.size());

        double min_state_time = ( *std::min_element(time.begin(), time.end()) ) - cutoff_;
        double max_state_time = ( *std::max_element(time.begin(), time.end()) ) + cutoff_;

        // Select the states that correspond to the queried timestamps
        int i = time_buffer_.selectLastBeforeOrFirst(min_state_time);
        while( (i < time_.size()) && (time_.at(i) <= max_state_time))
        {
            bool add_ptr = false;
            for(int j = 0; j < time.size(); ++j)
            {
                if( (time_.at(i) <= (time.at(j) + cutoff_) ) 
                    && (time_.at(i) >= (time.at(j) - cutoff_) ) )
                {
                    add_ptr = true;
                    rot_id.at(j).push_back(ptrs.size());
                    pos_id.at(j).push_back(ptrs.size()+1);

                    state_times.at(j).push_back(time_.at(i));

                }
            }

            if(add_ptr)
            {
                if(!query_only) used_.at(i) = true;
                ptrs.push_back(rot_.at(i)->data());
                ptrs.push_back(pos_.at(i)->data());
                if(type_ == GpStateType::kSe2)
                {
                    sizes.push_back(1);
                    sizes.push_back(2);
                }
                else
                {
                    sizes.push_back(3);
                    sizes.push_back(3);
                }
            }
            i++;
        }
        

        // Compute the ks_K_inv vectors for the GP inferrence
        for(int j = 0; j < time.size(); ++j)
        {
            if(K_inv_.count(state_times.at(j).size()) )
            {
                VecX eig_time(1);
                eig_time(0) = time.at(j);
                Eigen::Map<VecX> eig_state_time(state_times.at(j).data(),state_times.at(j).size());
                MatX ks = seKernel(eig_time, eig_state_time, hyper_.l2, hyper_.sf2);
                ks_K_inv.at(j) = (ks * K_inv_.at(state_times.at(j).size())).transpose();
            }
            else
            {
                throw std::range_error("GpStateManager::getStatePosePtrs: Seems that the number of state times does not corresponds to one of the two values");
            }
        }

        return {ptrs, sizes, rot_id, pos_id, ks_K_inv};
    }







    std::tuple<
        std::vector<double*>,   // Pointers to the states
        std::vector<IndexVec>,  // Id of the states for each of the timestamps
        std::vector<VecX>       // ksKinv
        > GpStateManager::getStateHomographyPtrs(const std::vector<double>& time, bool query_only)
    {
        if( type_ != GpStateType::kHomography) throw std::range_error("GpStateManager::getPosRotVec: Cannot retrieve homography pointers for non homography state");



        if(!query_only)  updateTime(time);



        std::vector<double*> ptrs;
        std::vector<IndexVec> ids(time.size());
        std::vector<std::vector<double> > ptrs_time(time.size());
        std::vector<VecX> ks_K_inv(time.size());

        double min_state_time = ( *std::min_element(time.begin(), time.end()) ) - cutoff_;
        double max_state_time = ( *std::max_element(time.begin(), time.end()) ) + cutoff_;
        int i = time_buffer_.selectLastBeforeOrFirst(min_state_time);

        while( (i < time_.size()) && (time_.at(i) <= max_state_time))
        {
            bool add_ptr = false;
            for(int j = 0; j < time.size(); ++j)
            {
                if( (time_.at(i) <= (time.at(j) + cutoff_) ) 
                    && (time_.at(i) >= (time.at(j) - cutoff_) ) )
                {
                    add_ptr = true;
                    ids.at(j).push_back(ptrs.size());
                    ptrs_time.at(j).push_back(time_.at(i));
                }
            }
            
            if(add_ptr)
            {
                if(!query_only) used_.at(i) = true;
                ptrs.push_back(homography_.at(i)->data());
            }
            i++;

        }
        // Compute the ks_K_inv vectors for the GP inferrence
        for(int j = 0; j < time.size(); ++j)
        {
            if(K_inv_.count(ptrs_time.at(j).size()) )
            {
                VecX eig_time(1);
                eig_time(0) = time.at(j);
                Eigen::Map<VecX> eig_state_time(ptrs_time.at(j).data(),ptrs_time.at(j).size());
                MatX ks = seKernel(eig_time, eig_state_time, hyper_.l2, hyper_.sf2);
                ks_K_inv.at(j) = (ks * K_inv_.at(ptrs_time.at(j).size())).transpose();
            }
            else
            {
                throw std::range_error("GpStateManager::getStateHomographyPtrs: Seems that the number of state times does not corresponds to one of the two values");
            }
        }
        return {ptrs, ids, ks_K_inv};
    }









    std::vector<double*> GpStateManager::getClosestStatePtrs(const double time)
    {
        std::vector<double*> output;
        int id = time_buffer_.selectIterator(time, state_period_)->second;
        if(type_ == GpStateType::kHomography)
        {
            output.push_back(homography_.at(id)->data());
        }
        else
        {
            output.push_back(rot_.at(id)->data());
            output.push_back(pos_.at(id)->data());
        }
        return output;
    }
    
    std::vector<double*> GpStateManager::getStatePtrsUpTo(const double time)
    {
        std::vector<double*> output;
        int id = time_buffer_.selectIterator(time, state_period_)->second;
        bool loop = true;
        while(loop)
        {
            if(type_ == GpStateType::kHomography)
            {
                output.push_back(homography_.at(id)->data());
            }
            else
            {
                output.push_back(rot_.at(id)->data());
                output.push_back(pos_.at(id)->data());
            }
            if(id == 0)
            {
                loop = false;
            }
            else
            {
                id--;
                if(time_.at(id) < (time - cutoff_) )
                {
                    loop = false;
                }
            }
        }
        return output;
    }


    std::pair<VecX, VecX> GpStateManager::getPoseRotVec(const double time)
    {
        if( type_ == GpStateType::kHomography) throw std::range_error("GpStateManager::getPosRotVec: Cannot retrieve pose for homography state");        


        auto [states, sizes, rot_ids, pos_ids, ksKinv] = GpStateManager::getStatePosePtrs(std::vector<double>(1, time), true);

        int rot_size = sizes.at(rot_ids.at(0).at(0));
        int pos_size = sizes.at(pos_ids.at(0).at(0));

        // First get the rot and pos of the slice based on GP trajectories
        MatX rot_state(rot_ids.at(0).size(), rot_size);
        MatX pos_state(rot_ids.at(0).size(), pos_size);

        for(int j = 0; j < rot_ids.at(0).size(); ++j)
        {
            for(int r = 0; r < rot_size; ++r)
            {
                rot_state(j,r) = states[rot_ids.at(0).at(j)][r];
            }
            for(int p = 0; p < pos_size; ++p)
            {
                pos_state(j,p) = states[pos_ids.at(0).at(j)][p];
            }
        }


        VecX pos = (ksKinv.at(0).transpose()*pos_state).transpose();
        VecX rot = (ksKinv.at(0).transpose()*rot_state).transpose();

        return {rot, pos};
    }


    Mat3 GpStateManager::getHomography(const double time)
    {
        VecX h = getHomographyVec(time); 
        
        Mat3 output;
        output << h(0), h(3), h(6),
                  h(1), h(4), h(7),
                  h(2), h(5), h(8);

        return output;
    }

    Vec9 GpStateManager::getHomographyVec(const double time)
    {
        if( type_ != GpStateType::kHomography) throw std::range_error("GpStateManager::getPosRotVec: Cannot retrieve pose for homography state");        

        auto [states, ids, ksKinv] = GpStateManager::getStateHomographyPtrs(std::vector<double>(1, time), true);

        MatX state(ksKinv.at(0).size(), 8);

        for(int j = 0; j < ksKinv.at(0).size(); ++j)
        {
            Eigen::Map<Vec8> temp_state(states.at(j));
            state.row(j) = temp_state.transpose();
        }
        Vec9 h;
        h.segment<8>(0) = (ksKinv.at(0).transpose()*state).transpose();
        h(8) = 1.0;

        return h;
    }

    void GpStateManager::plotPose(const std::string window_prefix)
    {
        if( type_ == GpStateType::kHomography) throw std::range_error("GpStateManager::plotPose: Cannot plot pose for homography state");

        double resolution = (time_.back() - time_.at(0))/(15*time_.size());
        std::vector<Plotter> plots(rot_.at(0)->size()+pos_.at(0)->size());
        std::vector<std::vector<double> > y_inf(plots.size());
        std::vector<double> x_inf;
        std::vector<std::vector<double> > y_val(plots.size());
        for(double t = time_.at(0); t <= time_.back(); t+=resolution)
        {
            auto [rot, pos] = getPoseRotVec(t);

            int counter = 0;
            for(int r = 0; r < rot_.at(0)->size(); ++r)
            {
                y_inf.at(counter).push_back(rot(r));
                counter++;
            }
            for(int p = 0; p < pos_.at(0)->size(); ++p)
            {
                y_inf.at(counter).push_back(pos(p));
                counter++;
            }
            x_inf.push_back(t);
        }
        for(int i = 0; i < time_.size(); ++i)
        {
            int counter = 0;
            for(int r = 0; r < rot_.at(0)->size(); ++r)
            {
                y_val.at(counter).push_back(rot_.at(i)->at(r));
                counter++;
            }
            for(int p = 0; p < pos_.at(0)->size(); ++p)
            {
                y_val.at(counter).push_back(pos_.at(i)->at(p));
                counter++;
            }
        }
        int counter = 0;
        for(int r = 0; r < rot_.at(0)->size(); ++r)
        {
            plots.at(counter).setTitles(window_prefix + "Rot axis " + std::to_string(r) + " till " + std::to_string(time_.back()), "Time", "Rad");
            counter++;
        }
        for(int p = 0; p < pos_.at(0)->size(); ++p)
        {
            plots.at(counter).setTitles(window_prefix + "Pos axis " + std::to_string(p) + " till " + std::to_string(time_.back()), "Time", "pix");
            counter++;
        }

        for(int i = 0; i < plots.size(); ++i)
        {
            plots.at(i).addPlot(x_inf, y_inf.at(i), "Inference", celib::Color::BLACK, celib::PointStyle::NONE, celib::LineStyle::FULL);
            plots.at(i).addPlot(time_, y_val.at(i), "Inference", celib::Color::RED, celib::PointStyle::CROSS, celib::LineStyle::NONE);
            plots.at(i).plot();
        }

    }

    void GpStateManager::plotHomography(const std::string window_prefix)
    {
        if( type_ != GpStateType::kHomography) throw std::range_error("GpStateManager::plotHomography: Cannot plot homography for non homography state");

        double resolution = (time_.back() - time_.at(0))/(15*time_.size());
        Plotter plot(window_prefix + " Homography plot at " + std::to_string(time_.back()), "Time index", "Homography value");
        std::vector<std::vector<double> > y_inf(8);
        std::vector<double> x_inf;
        std::vector<std::vector<double> > y_val(8);
        for(double t = time_.at(0); t <= time_.back(); t+=resolution)
        {
            Vec9 h = getHomographyVec(t);
            for(int i = 0; i < 8; ++i)
            {
                y_inf.at(i).push_back(h(i));
            }
            x_inf.push_back(t);
        }
        for(int i = 0; i < time_.size(); ++i)
        {
            for(int j = 0; j < 8; ++j)
            {
                y_val.at(j).push_back(homography_.at(i)->at(j));
            }
        }


        plot.addPlot(x_inf, y_inf.at(0), "Inference h 0", celib::Color::BLACK, celib::PointStyle::NONE, celib::LineStyle::FULL);
        plot.addPlot(time_, y_val.at(0), "Value h 0"    , celib::Color::BLACK, celib::PointStyle::CROSS, celib::LineStyle::NONE);
        
        plot.addPlot(x_inf, y_inf.at(1), "Inference h 1", celib::Color::BLUE, celib::PointStyle::NONE, celib::LineStyle::FULL);
        plot.addPlot(time_, y_val.at(1), "Value h 1"    , celib::Color::BLUE, celib::PointStyle::CIRCLE, celib::LineStyle::NONE);

        plot.addPlot(x_inf, y_inf.at(2), "Inference h 2", celib::Color::ORANGE, celib::PointStyle::NONE, celib::LineStyle::DOTS);
        plot.addPlot(time_, y_val.at(2), "Value h 2"    , celib::Color::ORANGE, celib::PointStyle::CROSS, celib::LineStyle::NONE);

        plot.addPlot(x_inf, y_inf.at(3), "Inference h 3", celib::Color::PURPLE, celib::PointStyle::NONE, celib::LineStyle::DOTS);
        plot.addPlot(time_, y_val.at(3), "Value h 3"    , celib::Color::PURPLE, celib::PointStyle::CROSS, celib::LineStyle::NONE);

        plot.addPlot(x_inf, y_inf.at(4), "Inference h 4", celib::Color::CYAN, celib::PointStyle::NONE, celib::LineStyle::DOTS);
        plot.addPlot(time_, y_val.at(4), "Value h 4"    , celib::Color::CYAN, celib::PointStyle::CIRCLE, celib::LineStyle::NONE);

        plot.addPlot(x_inf, y_inf.at(5), "Inference h 5", celib::Color::RED, celib::PointStyle::NONE, celib::LineStyle::FULL);
        plot.addPlot(time_, y_val.at(5), "Value h 5"    , celib::Color::RED, celib::PointStyle::CROSS, celib::LineStyle::NONE);
        
        plot.addPlot(x_inf, y_inf.at(6), "Inference h 6", celib::Color::DARK_GREEN, celib::PointStyle::NONE, celib::LineStyle::FULL);
        plot.addPlot(time_, y_val.at(6), "Value h 6"    , celib::Color::DARK_GREEN, celib::PointStyle::CIRCLE, celib::LineStyle::NONE);

        plot.addPlot(x_inf, y_inf.at(7), "Inference h 7", celib::Color::GREY, celib::PointStyle::NONE, celib::LineStyle::DOTS);
        plot.addPlot(time_, y_val.at(7), "Value h 7"    , celib::Color::GREY, celib::PointStyle::CROSS, celib::LineStyle::NONE);

        plot.plot();

    }


    void GpStateManager::addNormalisersToProblem(ceres::Problem& problem) 
    {
        std::vector<std::vector<double*> > normaliser_state;
        std::vector<std::vector<int> > normaliser_size;
        std::vector<VecX> normaliser_ksKinv;
        std::vector<double*> state;
        std::vector<bool> temp_used = used_;

        for(int j = 0; j < used_.size(); ++j)
        {
            if(temp_used.at(j))
            {
                std::vector<std::vector<double*> > normaliser_state;
                std::vector<std::vector<int> > normaliser_size;
                std::vector<VecX> normaliser_ksKinv;
                std::vector<double*> state;
                std::vector<double> times;

                if(type_ == kHomography)
                {
                    state.push_back(homography_.at(j)->data());
                    normaliser_state.resize(1);
                }
                else
                {
                    state.push_back(rot_.at(j)->data());
                    state.push_back(pos_.at(j)->data());
                    normaliser_state.resize(2);
                }
                normaliser_size.resize(normaliser_state.size());
                normaliser_ksKinv.resize(normaliser_state.size());

                int i = time_buffer_.selectLastBeforeOrFirst(time_.at(j) - cutoff_);
                double max_time = time_.at(j) + cutoff_;
                while( (i < time_.size()) && (time_.at(i) <= max_time))
                {
                    if( (time_.at(i) <= (time_.at(j) + cutoff_) ) 
                        && (time_.at(i) >= (time_.at(j) - cutoff_) ) )
                    {
                        times.push_back(time_.at(i));
                        used_.at(i) = true;

                        if(type_ == kHomography)
                        {
                            normaliser_state.at(0).push_back(homography_.at(i)->data());
                            normaliser_size.at(0).push_back(8);
                        }
                        else
                        {
                            normaliser_state.at(0).push_back(rot_.at(i)->data());
                            normaliser_state.at(1).push_back(pos_.at(i)->data());
                            if(type_ == kSe2)
                            {
                                normaliser_size.at(0).push_back(1);
                                normaliser_size.at(1).push_back(2);
                            }
                            else
                            {
                                normaliser_size.at(0).push_back(3);
                                normaliser_size.at(1).push_back(3);
                            }
                        }
                    }
                    ++i;
                }
                if(K_inv_.count(times.size()) )
                {
                    VecX eig_ptr_time(1);
                    eig_ptr_time(0) = time_.at(j);
                    Eigen::Map<VecX> eig_times(times.data(), times.size());
                    for(int c = 0; c < normaliser_state.size(); ++c)
                    {
                        MatX ks;
                        ks = seKernel(eig_ptr_time, eig_times, hyper_.l2, hyper_.sf2);
                        
                        normaliser_ksKinv.at(c) = (ks * K_inv_.at(times.size())).transpose();
                    }
                }
                else
                {
                    throw std::range_error("GpStateManager::addNormalisersToProblem: Seems that the number of state times for normalisation does not corresponds to one of the two values");
                }

                for(int r = 0; r < normaliser_state.size(); ++r)
                {
                    CostFunctionGPNormaliser* gp_norm_cost_function = new CostFunctionGPNormaliser(state.at(r), normaliser_state.at(r), normaliser_size.at(r), normaliser_ksKinv.at(r));
                    problem.AddResidualBlock(gp_norm_cost_function, NULL, normaliser_state.at(r));
                }
            }

        }

    }

    void GpStateManager::resetUsedFlags()
    {
        for(int i = 0; i < used_.size(); ++i) used_.at(i) = false;
    }



    // Set constant the state variables that have a timestamp inferrior to "time"
    void GpStateManager::setConstantTill(ceres::Problem& problem, const double time)
    {
        // Get the variables used with the appropriate timestamp
        std::vector<double*> state_to_make_constant;
        for(int i = 0; i < used_.size(); ++i)
        {
            if(time_.at(i) > time) break;

            if(type_ == GpStateType::kHomography)
            {
                state_to_make_constant.push_back(homography_.at(i)->data());
            }
            else
            {
                state_to_make_constant.push_back(rot_.at(i)->data());
                state_to_make_constant.push_back(pos_.at(i)->data());
            }
        }


        // Set to constant the state variables
        for(int i = 0; i < state_to_make_constant.size(); ++i)
        {
            problem.SetParameterBlockConstant(state_to_make_constant.at(i));
        }
    }



    // Set constant the state variables that have the closest timestamp to "time"
    void GpStateManager::setConstantCloseTo(ceres::Problem& problem, const double time) const
    {
        // Get the timestamps of state variables that have been used
        std::vector<double> used_time;
        std::vector<double> index_mapping;
        for(int i = 0; i < used_.size(); ++i)
        {
            if(used_.at(i))
            {
                used_time.push_back(time_.at(i));
                index_mapping.push_back(i);
            }
        }

        // Find the used timestamp the closest to the provided time
        int id_closest = 0;
        int i = 0;
        bool loop = true;
        while(loop)
        {
            if( std::abs(time - used_time.at(i)) > std::abs(time - used_time.at(i+1)) )
            {
                id_closest = i+1;
            }
            else
            {
                loop = false;
            }

            ++i;
            if(i == (used_time.size()-1)) loop = false;
        }
        id_closest = index_mapping.at(id_closest);


        // Get the corresponding state variables
        std::vector<double*> state_to_make_constant;
        if(type_ == GpStateType::kHomography)
        {
            state_to_make_constant.push_back(homography_.at(id_closest)->data());
        }
        else
        {
            state_to_make_constant.push_back(rot_.at(id_closest)->data());
            state_to_make_constant.push_back(pos_.at(id_closest)->data());
        }


        // Set to constant the state variables
        for(int i = 0; i < state_to_make_constant.size(); ++i)
        {
            problem.SetParameterBlockConstant(state_to_make_constant.at(i));
        }
    }


} // namespace celib
