/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#include "event_gp/cost_function.h"
#include "event_gp/se2_frontend.h"
#include "cilantro/core/kd_tree.hpp"
#include "common/random.h"



namespace celib
{



// TODO To make more efficient:
// . Investigate the use of KISS GP
    CostFunctionSe2FrontEnd::CostFunctionSe2FrontEnd(
            const std::vector<EventPtr>& events_in,
            const std::vector<int>& event_ids_in,
            const std::vector<EventPtr>& fixed_events,
            GpStateManager& state_manager,
            const HyperParam& hyper,
            std::vector<double*>& state_ptrs,
            const bool filter,
            const int downsample
            )
                : xy_hyper_(hyper)
    {


        std::vector<EventPtr> events;
        std::vector<int> event_ids;

        events = events_in;
        event_ids = event_ids_in;

        if(downsample > 0)
        {
            // Downsample the events to match 'downsample' but keep the first and last 20 events
            std::vector<EventPtr> events_ds;
            std::vector<int> event_ids_ds;
            events_ds.reserve(downsample);
            for(int i=0; i<20;++i)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }

            std::vector<int> downsampled_ids = generateRandomIndexes(20, events.size()-20, downsample-40);
            for(const auto& i: downsampled_ids)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }

            for(int i=events.size()-20; i<events.size();++i)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }
            events = events_ds;
            event_ids = event_ids_ds;
        }


        // Cast the indexes to match the state_manager type
        std::vector<double> event_time(event_ids.begin(), event_ids.end());

        std::tie(state_ptrs, state_size_, rot_ids_, pos_ids_, ks_K_inv_) = state_manager.getStatePosePtrs(event_time);


        // Fill the information needed by Ceres
        state_to_event_.resize(state_ptrs.size());
        set_num_residuals(1);
        std::vector<int>* block_sizes = mutable_parameter_block_sizes();
        for(const auto& s: state_size_) block_sizes->push_back(s);


        // Store the event locations                
        events_xy_.resize(2,events.size());
        fixed_events_.resize(2,fixed_events.size());
        for(int i = 0; i < fixed_events.size(); ++i)
        {
            fixed_events_(0,i) = fixed_events.at(i)->x;
            fixed_events_(1,i) = fixed_events.at(i)->y;
        }
        for(int i = 0; i < events.size(); ++i)
        {
            events_xy_(0,i) = events.at(i)->x;
            events_xy_(1,i) = events.at(i)->y;

            for(int j = 0; j < rot_ids_.at(i).size(); ++j)
            {
                state_to_event_.at(rot_ids_.at(i).at(j)).push_back(i);
                state_to_event_.at(pos_ids_.at(i).at(j)).push_back(i);
            }
        }


    }

    bool CostFunctionSe2FrontEnd::Evaluate(double const* const* parameters,
                                double* residuals,
                                double** jacobians) const
    {
        // Store the jacobians of the reprojected events 
        // with std::map< std::pair<event id, paramter block id>, d_xy_d_param_block
        std::vector<std::vector<MatX> > jacobian_xy_storage(events_xy_.cols() + fixed_events_.size(), std::vector<MatX>(state_size_.size()));

        // Store the reprojected events
        MatX xy(2, events_xy_.cols() + fixed_events_.cols());
        if(fixed_events_.cols() != 0) xy.block(0, events_xy_.cols(), 2, fixed_events_.cols()) = fixed_events_;

        // Reproject events in first time stamp.
        #pragma omp parallel for
        for(int i = 0; i < events_xy_.cols(); ++i)
        {
            Vec2 pos = ceresParamToEigen(parameters, pos_ids_.at(i), ks_K_inv_.at(i), 2);
            double rot = ceresParamToEigen(parameters, rot_ids_.at(i), ks_K_inv_.at(i), 1)(0);


            if(jacobians != NULL)
            {
                Mat2 d_xy_d_pos;
                Vec2 d_xy_d_rot;
                Vec2 temp_xy;
                std::tie(temp_xy, d_xy_d_pos, d_xy_d_rot) = se2TransJacobian(events_xy_.col(i), pos, rot);
                xy.col(i) = temp_xy;
                for(int j = 0; j < rot_ids_.at(i).size(); ++j)
                {
                    if(jacobians[rot_ids_.at(i).at(j)] != NULL)
                    {

                        jacobian_xy_storage[i][rot_ids_.at(i).at(j)] = d_xy_d_rot*ks_K_inv_.at(i)(j);
                    }
                    if(jacobians[pos_ids_.at(i).at(j)] != NULL)
                    {
                        jacobian_xy_storage[i][pos_ids_.at(i).at(j)] = d_xy_d_pos*ks_K_inv_.at(i)(j);
                    }
                }
            }
            else
            {
                xy.col(i) = se2Trans(events_xy_.col(i), pos, rot);
            }
        }


        VecX Y(xy.cols());
        Y.setOnes();

        MatX& X = xy;

        // Compute the covariance matrix
        MatX K;
        std::vector<std::vector<Row4, Eigen::aligned_allocator<Row4> > > jacobian_K_storage;
        if(jacobians != NULL)
        {
            jacobian_K_storage = std::vector<std::vector<Row4, Eigen::aligned_allocator<Row4> > >(xy.cols(), std::vector<Row4, Eigen::aligned_allocator<Row4> >(xy.cols()));
            #pragma omp parallel for
            for(int r = 0; r < xy.cols(); ++r)
            {
                for(int c = r; c < xy.cols(); ++c)
                {
                    jacobian_K_storage[r][c] = seKernelJacobians(xy.col(r), xy.col(c), xy_hyper_.l2, xy_hyper_.sf2);
                    jacobian_K_storage[c][r].segment<2>(0) = jacobian_K_storage[r][c].segment<2>(2);
                    jacobian_K_storage[c][r].segment<2>(2) = jacobian_K_storage[r][c].segment<2>(0);
                }
            }

            K = seKernelCov(X, xy_hyper_.l2, xy_hyper_.sf2) + 0.001*MatX::Identity(X.cols(), X.cols());
        }
        else
        {
            K = seKernelCov(X, xy_hyper_.l2, xy_hyper_.sf2) + 0.001*MatX::Identity(X.cols(), X.cols());
        }

        VecX alpha = solveKinvY(K,Y);

        // Compute the cost function
        residuals[0] = alpha.sum();


        if(jacobians != NULL)
        {
            // Loop through all the parameter blocks
            #pragma omp parallel for
            for(int b = 0; b < state_to_event_.size(); ++b)
            {
                if(jacobians[b] != NULL)
                {
                    Eigen::Map<RowX> j_r(&jacobians[b][0], 1, state_size_.at(b));

                    // Loop through all the components of the parameter block
                    for(int j = 0; j < state_size_.at(b); ++j)
                    {
                        MatX d_K_d_theta_j(xy.cols(), xy.cols());
                        std::vector<std::vector<bool> > filled(xy.cols(),std::vector<bool>(xy.cols(), false));
                        d_K_d_theta_j.setZero();

                        // Loop through all the events
                        for(int i = 0; i < state_to_event_.at(b).size(); ++i)
                        {
                            int r = state_to_event_.at(b).at(i);

                            for(int c = 0; c < xy.cols(); ++c)
                            {
                                if(!filled[r][c])
                                {
                                    Vec4 d_xr_xc_d_theta_j;
                                    if(jacobian_xy_storage[r][b].cols()!=0)
                                    {
                                        d_xr_xc_d_theta_j.segment<2>(0) = jacobian_xy_storage[r][b].col(j);
                                    }
                                    else
                                    {
                                        d_xr_xc_d_theta_j.segment<2>(0) = Vec2::Zero();
                                    }

                                    if(jacobian_xy_storage[c][b].cols()!=0)
                                    {
                                        d_xr_xc_d_theta_j.segment<2>(2) = jacobian_xy_storage[c][b].col(j);
                                    }
                                    else
                                    {
                                        d_xr_xc_d_theta_j.segment<2>(2) = Vec2::Zero();
                                    }
                                    double temp = jacobian_K_storage[r][c] * d_xr_xc_d_theta_j;

                                    d_K_d_theta_j(r,c) = temp;
                                    d_K_d_theta_j(c,r) = temp;
                                    filled[r][c] = true;
                                    filled[c][r] = true;
                                }
                            }




                        }
                        j_r(j) = -alpha.segment(0,xy.cols()).transpose()*d_K_d_theta_j*alpha.segment(0,xy.cols());

                    }
                }
            }
        }

        return true;
    }















    CostFunctionSe2ApproxFrontEnd::CostFunctionSe2ApproxFrontEnd(
            const std::vector<EventPtr>& events_in,
            const std::vector<int>& event_ids_in,
            GpStateManager& state_manager,
            const HyperParam& hyper,
            std::vector<double*>& state_ptrs,
            const bool filter,
            const int downsample
            )
                : xy_hyper_(hyper)
    {


        std::vector<EventPtr> events;
        std::vector<int> event_ids;

        events = events_in;
        event_ids = event_ids_in;

        if(downsample > 0)
        {
            // Downsample the events to match 'downsample' but keep the first and last 20 events
            std::vector<EventPtr> events_ds;
            std::vector<int> event_ids_ds;
            events_ds.reserve(downsample);
            for(int i=0; i<20;++i)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }

            std::vector<int> downsampled_ids = generateRandomIndexes(20, events.size()-20, downsample-40);
            for(const auto& i: downsampled_ids)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }

            for(int i=events.size()-20; i<events.size();++i)
            {
                events_ds.push_back(events.at(i));
                event_ids_ds.push_back(event_ids.at(i));
            }
            events = events_ds;
            event_ids = event_ids_ds;
        }


        // Cast the indexes to match the state_manager type
        std::vector<double> event_time(event_ids.begin(), event_ids.end());

        std::tie(state_ptrs, state_size_, rot_ids_, pos_ids_, ks_K_inv_) = state_manager.getStatePosePtrs(event_time);


        // Fill the information needed by Ceres
        state_to_event_.resize(state_ptrs.size());
        set_num_residuals(event_time.size());
        std::vector<int>* block_sizes = mutable_parameter_block_sizes();
        for(const auto& s: state_size_) block_sizes->push_back(s);


        // Store the event locations                
        events_xy_.resize(2,events.size());
        for(int i = 0; i < events.size(); ++i)
        {
            events_xy_(0,i) = events.at(i)->x;
            events_xy_(1,i) = events.at(i)->y;

            for(int j = 0; j < rot_ids_.at(i).size(); ++j)
            {
                state_to_event_.at(rot_ids_.at(i).at(j)).push_back(i);
                state_to_event_.at(pos_ids_.at(i).at(j)).push_back(i);
            }
        }


    }


    bool CostFunctionSe2ApproxFrontEnd::Evaluate(double const* const* parameters,
                                double* residuals,
                                double** jacobians) const
    {
        // Store the jacobians of the reprojected events 
        // with std::map< std::pair<event id, paramter block id>, d_xy_d_param_block
        std::vector<std::vector<MatX> > jacobian_xy_storage(events_xy_.cols(), std::vector<MatX>(state_size_.size()));

        // Store the reprojected events
        MatX xy(2, events_xy_.cols());

        // Reproject events in first time stamp.
        //#pragma omp parallel for
        for(int i = 0; i < events_xy_.cols(); ++i)
        {
            Vec2 pos = ceresParamToEigen(parameters, pos_ids_.at(i), ks_K_inv_.at(i), 2);
            double rot = ceresParamToEigen(parameters, rot_ids_.at(i), ks_K_inv_.at(i), 1)(0);


            if(jacobians != NULL)
            {
                Mat2 d_xy_d_pos;
                Vec2 d_xy_d_rot;
                Vec2 temp_xy;
                std::tie(temp_xy, d_xy_d_pos, d_xy_d_rot) = se2TransJacobian(events_xy_.col(i), pos, rot);
                xy.col(i) = temp_xy;
                for(int j = 0; j < rot_ids_.at(i).size(); ++j)
                {
                    if(jacobians[rot_ids_.at(i).at(j)] != NULL)
                    {

                        jacobian_xy_storage[i][rot_ids_.at(i).at(j)] = d_xy_d_rot*ks_K_inv_.at(i)(j);
                    }
                    if(jacobians[pos_ids_.at(i).at(j)] != NULL)
                    {
                        jacobian_xy_storage[i][pos_ids_.at(i).at(j)] = d_xy_d_pos*ks_K_inv_.at(i)(j);
                    }
                }
            }
            else
            {
                xy.col(i) = se2Trans(events_xy_.col(i), pos, rot);
            }
        }


        VecX Y(xy.cols());
        Y.setOnes();

        MatX& X = xy;

        // Create a kd-tree in 2D for the events in xy using cilantro
        std::vector<Eigen::Vector2f> xy_vec;
        xy_vec.reserve(xy.cols());
        for(int i = 0; i < xy.cols(); ++i)
        {
            xy_vec.push_back(xy.col(i).cast<float>());
        }
        cilantro::KDTree2f<> kd_tree(xy_vec);

        // For each event, find the closest events in a radius of 5 times the lenghtscale
        std::vector<std::vector<int> > closest_events(xy.cols());
        float radius = 20.0*xy_hyper_.l2*xy_hyper_.l2;
        for(int i = 0; i < xy.cols(); ++i)
        {
            cilantro::NeighborSet<float> nn = kd_tree.radiusSearch(xy_vec[i], radius);
            for (int j = 0; j < nn.size(); j++) {
                closest_events[i].push_back(nn[j].index);
            }
        }

        // Compute the covariance vector for each event and its neighbours
        for(int i = 0; i < xy.cols(); ++i)
        {
            MatX X_local(2, closest_events[i].size());
            for(int j = 0; j < closest_events[i].size(); ++j)
            {
                X_local.col(j) = xy.col(closest_events[i].at(j));
            }
            MatX temp_xy(2,1);
            temp_xy = xy.col(i);
            MatX k = seKernelXD(X_local, temp_xy, xy_hyper_.l2, xy_hyper_.sf2);
            double k_sum = k.sum();
            residuals[i] = 1.0/k_sum;

            if(jacobians != NULL)
            {
                double factor = -1.0/(k_sum*k_sum);
                // Loop through all the entries in the covariance vector
                for(int h = 0; h < closest_events[i].size(); ++h)
                {
                    Row4 d_k_d_xi = seKernelJacobians(X_local.col(h), xy.col(i), xy_hyper_.l2, xy_hyper_.sf2);
                    
                    int c = closest_events[i].at(h);
                    // Loop through all the parameter blocks
                    for(int j = 0; j < rot_ids_.at(i).size(); ++j)
                    {
                        if(jacobians[rot_ids_.at(i).at(j)] != NULL)
                        {
                            if(h==0)
                            {
                                jacobians[rot_ids_.at(i).at(j)][i] = 0.0;
                            }
                            else
                            {
                                MatX temp = (d_k_d_xi.segment<2>(0)*jacobian_xy_storage[c][rot_ids_.at(c).at(j)] + d_k_d_xi.segment<2>(2)*jacobian_xy_storage[i][rot_ids_.at(i).at(j)]);
                                jacobians[rot_ids_.at(i).at(j)][i] += factor*temp(0,0);
                            }

                        }
                        if(jacobians[pos_ids_.at(i).at(j)] != NULL)
                        {
                            if(h==0)
                            {
                                jacobians[pos_ids_.at(i).at(j)][2*i] = 0.0;
                                jacobians[pos_ids_.at(i).at(j)][2*i+1] = 0.0;
                            }
                            else
                            {
                                MatX temp = (d_k_d_xi.segment<2>(0)*jacobian_xy_storage[c][pos_ids_.at(c).at(j)] + d_k_d_xi.segment<2>(2)*jacobian_xy_storage[i][pos_ids_.at(i).at(j)]);
                                jacobians[pos_ids_.at(i).at(j)][2*i] += factor*temp(0,0);
                                jacobians[pos_ids_.at(i).at(j)][2*i+1] += factor*temp(0,1);
                            }

                        }
                    }

                }

            }
        }

        return true;
    }

    bool CostFunctionSe2ApproxFrontEnd::testJacobian()
    {
        std::vector<int> block_sizes = parameter_block_sizes();

        RandomGenerator r_gen;

        std::vector<double*> parameters(block_sizes.size());
        for(int i = 0; i < block_sizes.size(); ++i)
        {
            std::cout << "Block " << i << " size " << block_sizes[i] << std::endl;
            parameters[i] = new double[block_sizes[i]];
            for(int j = 0; j < block_sizes[i]; ++j)
            {
                parameters[i][j] = r_gen.randGauss(0,1);
            }
        }

        double** jacobians = new double*[block_sizes.size()];
        for(int i = 0; i < block_sizes.size(); ++i)
        {
            jacobians[i] = new double[block_sizes[i]*num_residuals()];
        }

        double* residuals = new double[num_residuals()];

        Evaluate(parameters.data(), residuals, jacobians);

        double* residuals_pert = new double[num_residuals()];

        double pert = 1e-6;

        for(int i = 0; i < block_sizes.size(); ++i)
        {
            for(int j = 0; j < block_sizes[i]; ++j)
            {
                parameters[i][j] += pert;
                Evaluate(parameters.data(), residuals_pert, NULL);
                parameters[i][j] -= pert;

                std::cout << "\n\n\n\n\n\n\nBlock " << i << " parameter " << j << std::endl;

                for(int k = 0; k < num_residuals(); ++k)
                {
                    double diff = (residuals_pert[k] - residuals[k])/pert;
                    std::cout << "num diff: " << diff << " \t jacobians[i][j]: " << jacobians[i][k*block_sizes[i] +j] << " \t error " << diff - jacobians[i][k*block_sizes[i]+j] << std::endl;
                }
            }
        }

        // Free the memory
        for(int i = 0; i < block_sizes.size(); ++i)
        {
            delete[] parameters[i];
            delete[] jacobians[i];
        }
        delete[] jacobians;
        delete[] residuals;
        delete[] residuals_pert;
    

        return true;
        
    }















    CostFunctionDiscreteInvHomographyFrontEnd::CostFunctionDiscreteInvHomographyFrontEnd(
                    const EventPtr& event,
                    const Mat2XPtr& template_events,
                    const VecXPtr& alpha,
                    const HyperParam& hyper,
                    const double weight
                    )
                    : template_(template_events)
                    , alpha_(alpha)
                    , event_xy_(event->x, event->y)
                    , xy_hyper_(hyper)
                    , weight_(weight)
    {
        set_num_residuals(1);

        std::vector<int>* block_sizes = mutable_parameter_block_sizes();
        block_sizes->push_back(8);
    }



    bool CostFunctionDiscreteInvHomographyFrontEnd::Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const
    {

        const Eigen::Map<const Vec8> vec_homography(parameters[0]);

        Vec2 ev;
        Mat2_8 d_ev_d_h;
        if(jacobians != NULL)
        {
            std::tie(ev, d_ev_d_h) = homographyInvProjectionJacobian(event_xy_, vec_homography);
        }
        else
        {
            ev = homographyInvProjection(event_xy_, vec_homography);
        }


        MatX ks = seKernelND(ev, *template_, xy_hyper_.l2, xy_hyper_.sf2);

        double inferrence = (ks*(*alpha_))(0,0);

        residuals[0] = (inferrence <= 0) ? 1000.0 : (-std::log(inferrence));
        residuals[0] *= weight_;

        if(jacobians != NULL)
        {
            Vec2 d_res_d_ev = (ksJacobian(ks.row(0), ev, *template_, xy_hyper_.l2) * (*alpha_)) / (-inferrence);

            if(jacobians[0] != NULL)
            {
                Eigen::Map<Vec8> d_res_d_s(jacobians[0]);
                d_res_d_s.setZero();
                if(inferrence != 0)
                {
                    d_res_d_s = (d_res_d_ev.transpose() * d_ev_d_h).transpose();
                    d_res_d_s *= weight_;
                }
            }
        }
        return true;
    }



















    CostFunctionDiscreteHomographyFrontEnd::CostFunctionDiscreteHomographyFrontEnd(
                    const EventPtr& event,
                    const Mat2XPtr& template_events,
                    const VecXPtr& alpha,
                    const HyperParam& hyper,
                    const double weight
                    )
                    : template_(template_events)
                    , alpha_(alpha)
                    , event_xy_(event->x, event->y)
                    , xy_hyper_(hyper)
                    , weight_(weight)
    {
        set_num_residuals(1);

        std::vector<int>* block_sizes = mutable_parameter_block_sizes();
        block_sizes->push_back(8);
    }



    bool CostFunctionDiscreteHomographyFrontEnd::Evaluate(double const* const* parameters,
                                        double* residuals,
                                        double** jacobians) const
    {

        const Eigen::Map<const Vec8> vec_homography(parameters[0]);

        Vec2 ev;
        Mat2_8 d_ev_d_h;
        if(jacobians != NULL)
        {
            std::tie(ev, d_ev_d_h) = homographyProjectionJacobian(event_xy_, vec_homography);
        }
        else
        {
            ev = homographyProjection(event_xy_, vec_homography);
        }

        MatX ks = seKernelND(ev, *template_, xy_hyper_.l2, xy_hyper_.sf2);

        double inferrence = (ks*(*alpha_))(0,0);

        residuals[0] = (inferrence <= 0) ? 1000.0 : (-std::log(inferrence));
        residuals[0] *= weight_;


        if(jacobians != NULL)
        {
            Vec2 d_res_d_ev = (ksJacobian(ks.row(0), ev, *template_, xy_hyper_.l2) * (*alpha_)) / (-inferrence);

            if(jacobians[0] != NULL)
            {
                Eigen::Map<Vec8> d_res_d_s(jacobians[0]);
                d_res_d_s.setZero();
                if(inferrence != 0)
                {
                    d_res_d_s = (d_res_d_ev.transpose() * d_ev_d_h).transpose();
                    d_res_d_s *= weight_;
                }
            }
        }

        return true;
    }

} //namespace