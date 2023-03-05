/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#include "event_gp/se2_frontend.h"
#include "common/internal_utils.h"
#include <iostream>
#include <fstream>
#include <future>

#include "cilantro/utilities/point_cloud.hpp"
#include "cilantro/visualization/visualizer.hpp"
#include "cilantro/visualization/common_renderables.hpp"
#include <omp.h>

namespace celib
{
    
    Se2GpFrontEnd::Se2GpFrontEnd(YAML::Node& config)
    {

        image_width_ = readRequiredField<double>(config, "cam_width");
        image_height_ = readRequiredField<double>(config, "cam_height");
        track_width_ = readField<double>(config, "track_width", 20);
        recenter_seeds_ = readField<bool>(config, "recenter_seeds", false);
        downsample_ = readField<int>(config, "downsample", 0);


        std::string feature_type = readRequiredField<std::string>(config, "feature_format");
        if(feature_type == "tracks")
        {
            use_tracks_after_ = readField<double>(config, "use_track_after", std::numeric_limits<double>::min());
            use_tracks_before_ = readField<double>(config, "use_track_before", std::numeric_limits<double>::max());
            loadTrackFile(readRequiredField<std::string>(config, "feature_path"), readField<std::string>(config, "time_offset", ""));
        }
        else if(feature_type == "debug")
        {
            feature_px_list_.push_back(Vec2(readField<double>(config, "debug_x_start", image_width_/2.0), readField<double>(config, "debug_y_start", image_height_/2.0)));
            feature_t_list_.push_back(readField<double>(config, "debug_t_start", 0));
            feature_ids_.push_back(0);
        }
        else
        {
            printErrorValueOption("feature_format");
        }

        std::string tracker_type_str = readField<std::string>(config, "tracker_type", "homography_template");
        if(tracker_type_str == "se2_only") tracker_type_ = TrackerType::kSE2;
        else if(tracker_type_str == "homography_frame") tracker_type_ = TrackerType::kHF2F;
        else tracker_type_ = TrackerType::kHF2T;


        max_track_length_ = readField<double>(config, "limit_track_length", std::numeric_limits<double>::max());

        filter_ = readField<bool>(config, "filter_events", true);
        if(filter_) filter_thr_ = readField<double>(config, "filter_events_thr", 1e-5);

        state_every_ = readField<int>(config, "state_every_n_events", 60);
        nb_state_opt_ = readField<int>(config, "nb_state_opt", 5);
        lengthscale_factor_ = readField<int>(config, "lengthscale_factor", 3);
        nb_lengthscale_cutoff_ = readField<int>(config, "nb_lengthscale_cutoff", 3);
        event_noise_ = readField<double>(config, "add_noise", 0.0);
        result_path_ = readRequiredField<std::string>(config, "result_path");
        write_events_ = readField<bool>(config, "write_events", false);
        only_one_ = readField<bool>(config, "only_one", false);

        visualise_ = readField<bool>(config, "visualisation", false);

        sae_.resize(image_height_, image_width_);
        sae_.setOnes();
        sae_ = -sae_;



        //viz_ = createVisualizer("SAE");

        hyper_.l2 = std::pow(lengthscale_factor_*state_every_,2);
        hyper_.sf2 = 100;
        hyper_.sz2 = 0.01;

        K_inv_ = std::shared_ptr<std::map<int, MatX> >(new std::map<int, MatX>);
        for(int i = 2*(lengthscale_factor_*nb_lengthscale_cutoff_+1); i <= (2*lengthscale_factor_*nb_lengthscale_cutoff_+1+nb_state_opt_); ++i)
        {
            VecX temp_x = VecX::LinSpaced(i, 0.0, (i-1)*state_every_);
            MatX K_inv = seKernel(temp_x, temp_x, hyper_.l2, hyper_.sf2).inverse();
            K_inv_->emplace(i,K_inv);
        }
    }

    bool Se2GpFrontEnd::isEventInImage(const EventPtr& event)
    {
        return ( (event->x >= 0) && (event->x <= (image_width_-1)) 
                && (event->y >= 0) && (event->y <= (image_height_-1)) );
    }


    void Se2GpFrontEnd::startTrack(const double t, const double x, const double y, const int id)
    {
        EventTrackPtr new_track(new EventTrack(
                id,
                x,y,
                t,
                track_width_,
                image_width_, image_height_,
                nb_state_opt_,
                state_every_,
                downsample_,
                visualise_,
                tracker_type_,
                result_path_,
                max_track_length_,
                recenter_seeds_,
                only_one_,
                write_events_
                ));
        std::cout << "Start track " << id << std::endl;
        active_tracks_.emplace(new_track);
        track_buffer_.emplace(t, new_track);
    }

    void Se2GpFrontEnd::addEvent(const EventPtr& event)
    {
        if(isEventInImage(event))
        {
            
            EventPtr temp_event(new Event);
            temp_event->x = event->x + 0.5;
            temp_event->y = event->y + 0.5;
            temp_event->t = event->t;
            temp_event->pol = event->pol;

            bool new_event = false;
            if(filter_)
            {
                int row = std::max(std::min(std::round(event->y), image_height_-1),0.0);
                int col = std::max(std::min(std::round(event->x), image_width_-1),0.0);
                if( (event->t - sae_(row, col)) > filter_thr_ )
                {
                    new_event = true;
                }
                sae_(row,col) = event->t;
            }
            else
            {   
                int row = std::max(std::min(std::round(event->y), image_height_-1),0.0);
                int col = std::max(std::min(std::round(event->x), image_width_-1),0.0);
                sae_(row,col) = event->t;
                new_event = true;
            }

            if(new_event)
            {
                event_buffer_.emplace(event->t, temp_event);

                if( (feature_list_ptrs_ < feature_px_list_.size()) && (event->t >= feature_t_list_.at(feature_list_ptrs_)) )
                {
                    Vec2& temp_px = feature_px_list_.at(feature_list_ptrs_);
                    startTrack(feature_t_list_.at(feature_list_ptrs_), temp_px(0), temp_px(1), feature_ids_.at(feature_list_ptrs_));
                    feature_list_ptrs_++;
                }

                std::vector<EventTrackPtr> inactive_tracks;
                std::vector<EventTrackPtr> temp_tracks(active_tracks_.begin(), active_tracks_.end());
                std::vector<std::future<bool> > futures;
                for(int i = 0; i < temp_tracks.size(); ++i)
                {   
                    if(!temp_tracks[i]->addEvent(temp_event, event_noise_))
                    {   
                        std::cout << "Tracks " << temp_tracks[i]->getId() << " terminated" << std::endl;
                        inactive_tracks.push_back(temp_tracks[i]);
                    }
                }
                for(auto& t: inactive_tracks)
                {   
                    active_tracks_.erase(t);
                }

            }
        }
    }

    void Se2GpFrontEnd::displaySae()
    {
        auto temp = saeToEvents(sae_);
        display(viz_, temp, true, true);
    }


    std::vector<EventPtr> saeToEvents(const MatX& sae)
    {
        // Initialise output structure
        std::vector<EventPtr> output;

        // Fill output structure
        for(int c = 0; c < sae.cols(); ++c)
        {
            for(int r = 0; r < sae.rows(); ++r)
            {
                if(sae(r,c) > -0.5)
                {
                    EventPtr temp(new Event);
                    temp->x = c;
                    temp->y = r;
                    temp->t = sae(r,c);
                    output.push_back(temp);
                }
            }
        }
        return output;
    }


    void Se2GpFrontEnd::loadTrackFile(std::string path, std::string time_offset)
    {
        double nsec = 0;
        int sec = 0;
        if(time_offset != "")
        {
            std::tie(sec, nsec) = stringToSec(time_offset);
        }
        std::ifstream file;
        file.open(path, std::ios::in);
        std::string line;
        std::set<int> feature_id_set;
        while(getline(file, line))
        {
            std::vector<double> line_values = stringToNumVector<double>(line, ' ');
            
            if(feature_id_set.count(line_values.at(0)) == 0)
            {
                feature_id_set.insert(line_values.at(0));
                std::vector<std::string> line_strings = stringToStringVector(line, ' ');
                auto [temp_sec, temp_nsec] = stringToSec(line_strings.at(1));

                double t = (temp_sec - sec);
                t = t + (temp_nsec - nsec);

                if( (t>= use_tracks_after_) && (t < use_tracks_before_ ) )
                {
                    feature_ids_.push_back(line_values.at(0));
                    feature_px_list_.push_back(Vec2(line_values.at(2),line_values.at(3)));
                    feature_t_list_.push_back(t);

                    std::cout.precision( std::numeric_limits<double>::digits10 + 1);
                    std::cout << "Track " << feature_ids_.back() << " added at " << feature_px_list_.back().transpose() << " from time " << feature_t_list_.back() << std::endl;
                }
            }
        }

        std::cout << feature_px_list_.size() << " track seeds have been added from time " << feature_t_list_.at(0) << " to " << feature_t_list_.back() << std::endl;
    }
    






    EventTrack::EventTrack(
                const int id,
                const double x,
                const double y,
                const double start_t,
                const double track_width,
                const double image_width,
                const double image_height,
                const int nb_state_opt,
                const int state_every,
                const int downsample,
                const bool visualise,
                const TrackerType type,
                const std::string write_path,
                const double limit_track_length,
                const bool recenter_seeds,
                const bool only_one,
                const bool write_events
                ):
                    id_(id)
                    ,nb_event_opt_((nb_state_opt)*state_every)
                    , track_half_width_(track_width/2.0)
                    , image_height_(image_height)
                    , image_width_(image_width)
                    , image_half_height_(image_height/2.0)
                    , image_half_width_(image_width/2.0)
                    , state_every_(state_every)
                    , nb_state_opt_(nb_state_opt)
                    , downsample_(downsample)
                    , cutoff_(state_every*kLengthScaleFactor*kCutoffFactor)
                    , cutoff_state_(kLengthScaleFactor*kCutoffFactor)
                    , track_width_(track_width)
                    , first_template_xy_(new Mat2X())
                    , first_template_alpha_(new VecX())
                    , last_template_xy_(new Mat2X())
                    , last_template_alpha_(new VecX())
                    , sae_(-MatX::Ones(image_height, image_width))
                    , connectivity_(image_height, std::vector<int>(image_width, 0))
                    , visualise_(visualise)
                    , write_path_(write_path)
                    , start_time_(start_t)
                    , d_template_(track_width, kGridFactor)
                    , only_one_(only_one)
                    , write_events_(write_events)
                    , max_length_(limit_track_length)
                    , disp_str_("Track " + std::to_string(id))
                    , type_(type)
                    , recenter_(recenter_seeds)
    {
        origin_[0] = x;
        origin_[1] = y;

        // Set hyperparameters
        hyper_.l2 = 0.25;
        hyper_.sf2 = 1;
        hyper_.sz2 = 0.001;

        centroid_dist_ = 1.2;        

        field_hyper_.l2 = 0.25; // 0.25
        field_hyper_.sf2 = 1;
        field_hyper_.sz2 = 0.1;
        

        // Init the homography to identity
        current_homography_.setZero();
        current_homography_(0) = 1.0;
        current_homography_(4) = 1.0;

        homographies_.push_back(current_homography_);


        // Set the optimisation options
        se2_opts_.minimizer_progress_to_stdout = false;
        se2_opts_.minimizer_type = ceres::LINE_SEARCH;
        se2_opts_.line_search_direction_type = ceres::BFGS;
        se2_opts_.max_num_iterations = 200;
        se2_opts_.function_tolerance = 1e-10;
        se2_opts_.num_threads = 1;

        h_opts_.minimizer_progress_to_stdout = false;
        h_opts_.num_threads = 14;

        std::cout << "Created new track with id " << id_ << std::endl;
    }



    bool EventTrack::isEventInTrack(const EventPtr& event)
    {
        int row = std::floor(event->y);
        int col = std::floor(event->x);
        bool keep = false;

        // Test if event projection is in original pattern area
        Vec2 temp_event = homographyInvProjection(Vec2::Zero(), current_homography_) + origin_;
        //if( first_ && ((event->pixVec() - temp_event).norm() < track_half_width_)) keep = true;
        if( ((event->pixVec() - temp_event).norm() < track_half_width_)) keep = true;

        // Test for connectivity
        if(!keep)
        {
            for(int dr = -1; dr <= 1; ++dr)
            {
                int row_nn = row + dr;
                if( (row_nn >= 0) && (row_nn < image_height_) )
                {
                    for(int dc = -1; dc <= 1; ++dc)
                    {
                        int col_nn = col + dc;
                        if( (col_nn >= 0) && (col_nn < image_width_) )
                        {
                            if(connectivity_[row_nn][col_nn] 
                                && ( (event->t - sae_(row_nn,col_nn)) < kAssociationTimeThr) )
                            {
                                keep = true;
                                connectivity_[row_nn][col_nn] -= 1;
                            }
                        }
                    }
                }
            }
        }

        // If keep the event, update the connectivity
        if(keep)
        {
            connectivity_[row][col] = kMaxConnectivity;
            sae_(row, col) = event->t;
        }
        return keep;


    }




    bool EventTrack::isTrackInImage()
    {
        Vec2 temp_head = homographyInvProjection(Vec2(0,0), current_homography_) + origin_;
        return (   (temp_head[0] >= 1.5*track_half_width_)
                && (temp_head[0] <= (image_width_ - 1.0 - (1.5*track_half_width_))) 
                && (temp_head[1] >= 1.5*track_half_width_)
                && (temp_head[1] <= (image_height_ - 1.0 - (1.5*track_half_width_))) );
    }


    bool EventTrack::addEvent(const EventPtr event, const double add_noise)
    {
        bool output = true;
        if(isEventInTrack(event))
        {
            events_.push_back(EventPtr(new Event));
            *(events_.back()) = *event;
            events_.back()->x -= origin_(0);
            events_.back()->y -= origin_(1);
            if(add_noise != 0)
            {
                events_.back()->x += rand_gen_.randGauss(0, add_noise);
                events_.back()->y += rand_gen_.randGauss(0, add_noise);
            }

            if( events_.size() >= nb_event_opt_)
            {
                if( ((events_.size()-nb_event_opt_)%(nb_state_opt_*state_every_)) == 0) output = optimise();
            }
        }
        return output;
    }



    // Helper functions for visualisation and debug
    void EventTrack::displayRawProjection(std::string window_name)
    {
        display(vizs_[1], events_, true);
    }

    void EventTrack::displayProjection(std::string window_name)
    {
        display(vizs_[2], undistorted_events_, true);
    }






    bool EventTrack::optimise()
    {
        if(first_ && recenter_)
        {
            Vec2 centroid = Vec2::Zero();
            for(const auto& e : events_) centroid += e->pixVec();
            centroid /= events_.size();
            for(int i = 0; i < events_.size(); ++i)
            {
                events_.at(i)->x -= centroid(0);
                events_.at(i)->y -= centroid(1);
            }
            origin_ += centroid;
        }


        disp_str_ = "Track " + std::to_string(id_) + " at nb event " + std::to_string(events_.size());


        // Perform the SE2 motion compensation and update the homography prior
        auto [event_batch, event_ids] = getLastEventBatch();
        

        auto [event_corrected, pose, valid, event_poses, ev_offset] = se2MotionCorrection(event_batch, event_ids, first_, only_one_);


        if(!valid) return false;





        
        std::cout << "**************************************************" << std::endl;
        std::cout << disp_str_ <<" length = " << (events_.back()->t - events_[0]->t) << "s" << std::endl;
        std::cout << "**************************************************" << std::endl;

        if(first_)
        {
            current_homography_ = homographyMatToVec(homographyVecToMat(current_homography_)*pose);
            // Create the template
            //createTemplate(event_corrected);

            // Add the SE2 corrected events to the list of undistorted events
            for(const auto& e: event_corrected) undistorted_events_.push_back(e);

            // Store the first event's timestamp to later recover the track
            homographies_t_.push_back(event_corrected.at(0)->t);
            first_pose_ = pose;
            poses_temp_ = event_poses;
            first_offset_ = ev_offset;

            // Do some filtering and add the events to the dynamic template
            auto[temp_events, temp_id] = filterGaussianDensity(event_corrected,event_ids,1);
            d_template_.addFrame(sampleCentroid(temp_events, 1.2), frame_id_);

            // Create the first GP template from the dynamic template
            createTemplateD();


            //auto [track_pt, track_valid] = d_template_.getMostViewed();
            //track_seed_ = track_pt;
            track_seed_ = Vec2::Zero();
            
        }
        else if(type_ == TrackerType::kSE2)
        {
            std::cout << "SE2 ONLY" << std::endl;
            current_homography_ = homographyMatToVec(homographyVecToMat(current_homography_)*pose);
            for(const auto& e: event_corrected) undistorted_events_.push_back(e);
            homographies_.push_back(current_homography_);
            homographies_t_.push_back(event_corrected.back()->t);
            interpolateTrack(event_corrected, event_poses, ev_offset);
            first_homography_ = false;
        }
        else
        {


            // Compute the local homography prior
            Vec8 frame_frame_prior = first_homography_ ? homographyMatToVec( first_pose_ * pose ): homographyMatToVec(pose);
            auto[delta_homography, valid_homography] = frameToFrameHomography(previous_se2_, event_corrected, field_hyper_, frame_frame_prior);
            current_homography_ = first_homography_ ? homographyMatToVec(delta_homography) :  homographyMatToVec(homographyVecToMat(current_homography_) * delta_homography);

            if(type_ == TrackerType::kHF2F)
            {
                std::cout << "HOMOGRAPHY FRAME TO FRAME" << std::endl;
                if(!se2HomoCoherence(pose, ev_offset)) return false;
                homographies_.push_back(current_homography_);
                homographies_t_.push_back(event_corrected.back()->t);
                interpolateTrack(event_corrected, event_poses, ev_offset);
                for(const auto& e: event_corrected) undistorted_events_.push_back(EventPtr(new Event(homographyProjection(e->pixVec(), current_homography_), e->t)));

            }
            else if(type_ == TrackerType::kHF2T)
            {
                std::cout << "HOMOGRAPHY FRAME TO TEMPLATE" << std::endl;

                // Perform the homography registration and update the template
                auto[event_homographed, homography] = discreteHomographyRegistration(event_corrected);

                if(!se2HomoCoherence(pose, ev_offset)) return false;

                homographies_.push_back(current_homography_);
                homographies_t_.push_back(event_corrected.back()->t);

                // Add the homography corrected events to the list of undistorted events
                for(const auto& e: event_homographed) undistorted_events_.push_back(e);

                interpolateTrack(event_corrected, event_poses, ev_offset);


                // Filter and add events to dynamic template object
                auto[temp_events, temp_id] = filterGaussianDensity(event_homographed,event_ids,1);
                d_template_.addFrame(sampleCentroid(temp_events, 1.2), frame_id_);

                // Update the dynamic template GP
                updateTemplateD();
            }
            first_homography_ = false;
        }
        if(visualise_) d_template_.displayFrameOccupancy(disp_str_ + " d template");

        // Write to file
        if(write_events_) writeEventsToFile(write_path_);
        if(only_one_)
        {
            for(int i = 0; i < event_corrected.size(); ++i)
                track_.push_back(EventPtr (new Event(se2InvTrans(-ev_offset, homographyVecToMat(event_poses.at(i))), event_corrected.at(i)->t)));
        }
        writeTrack();
        if(only_one_) return false;


        // Update info for next iteration and check if track is leaving the image
        previous_se2_ = event_corrected;
        cleanConnectiveSAE();
        index_to_optimise_ = events_.size();
        first_ = false;
        frame_id_++;
        bool output = isTrackInImage();
        if((events_.back()->t - events_[0]->t)> max_length_) output = false;
        if(!output) std::cout << "Finish track with id " << id_ << std::endl;
        return output;
    }





    // Return the last batch of events as well as the corresponding indexes
    std::tuple<std::vector<EventPtr>, std::vector<int> > EventTrack::getLastEventBatch()
    {
        std::vector<EventPtr> events;
        std::vector<int> event_ids;
        for(int i = index_to_optimise_; i < events_.size(); ++i)
        {
            events.push_back(events_.at(i));
            event_ids.push_back(i);
        }
        return {events, event_ids};
    }



    std::tuple<std::vector<EventPtr>, Mat3, bool, std::vector<Vec8>, Vec2 > EventTrack::se2MotionCorrection(
            const std::vector<EventPtr>& events, const std::vector<int>& event_ids, bool project_start, bool project_all)
    {


        std::vector<EventPtr> events_centered(events.size());
        events_centered.reserve(events.size());
        Vec2 sum_coor = Vec2::Zero();
        for(int i = 0; i < events.size(); ++i) {sum_coor(0) += events.at(i)->x; sum_coor(1) += events.at(i)->y;}
        sum_coor = sum_coor / events.size();


        for(int i = 0; i < events.size(); ++i) events_centered.at(i) = EventPtr(new Event(events.at(i)->pixVec() - sum_coor, events.at(i)->t));

        // Create the cost function
        // TODO can be optimised by reusing the inverted matrices
        GpStateManager state_manager(GpStateType::kSe2, 1.0/state_every_);
        state_manager.setStartingTime(event_ids.at(0));

        auto [pose, valid] = se2MotionOptimisation(events_centered, event_ids, state_manager);
        if(!valid)
        {
            RandomGenerator rg;
            state_manager = GpStateManager(GpStateType::kSe2, 1.0/state_every_);
            state_manager.setStartingTime(event_ids.at(0));
            std::vector<EventPtr> events_temp;
            events_temp.reserve(events_centered.size());
            for(const auto& e: events_centered) events_temp.push_back(EventPtr(new Event(e->pixVec()+Vec2(rg.randGauss(0,kEventStd),rg.randGauss(0,kEventStd)), e->t) ));
            std::tie(pose, valid) = se2MotionOptimisation(events_temp, event_ids, state_manager, 0.7);
        }

        // Reproject and filter events
        std::vector<EventPtr> temp_corrected;
        temp_corrected.reserve(events_centered.size());
        auto [rot_temp, pos_temp] = state_manager.getPoseRotVec(event_ids.at(0));
        Mat3 inv_pose_temp = invHomogeneous(rotPosToHomogeneous(angleToRotMat(rot_temp[0]), pos_temp));
        for(int i = 0; i < events_centered.size(); ++i)
        {
            auto [rot, pos] = state_manager.getPoseRotVec(event_ids.at(i));
            temp_corrected.push_back(se2Trans(events_centered.at(i), inv_pose_temp*rotPosToHomogeneous(angleToRotMat(rot[0]), pos)));
        }
        std::vector<double> score = getGaussianDensity(temp_corrected, 0.5);
        std::vector<EventPtr> temp_filtered;
        std::vector<int> temp_id;
        for(int i = 0; i < events_centered.size(); ++i)
        {
            if(score.at(i) > kDensityThr)
            {
                temp_filtered.push_back(events_centered.at(i));
                temp_id.push_back(event_ids.at(i));
            }
        }

        std::vector<Vec8> output_pose;
        if(!valid)
        {
            return {temp_corrected, pose, false, output_pose, sum_coor};
        }


        if( ((double)(temp_filtered.size()) / events_centered.size()) < 0.80)
        {
            
            std::tie(pose, valid) = se2MotionOptimisation(temp_filtered, temp_id, state_manager);
        }


        if(project_all)
        {
            temp_filtered = events_centered;
            temp_id = event_ids;
        }


        // Reproject the events
        auto [rot_0, pos_0] = state_manager.getPoseRotVec(temp_id.at(0));
        auto [rot_1, pos_1] = state_manager.getPoseRotVec(temp_id.back());
        Mat3 pose_0 = rotPosToHomogeneous(angleToRotMat(rot_0[0]), pos_0);
        Mat3 pose_1 = rotPosToHomogeneous(angleToRotMat(rot_1[0]), pos_1);
        Mat3 inv_pose_0 = invHomogeneous(pose_0);
        std::vector<EventPtr> event_reproj;
        event_reproj.reserve(temp_filtered.size());
        output_pose.reserve(temp_filtered.size());
        for(int i = 0; i < temp_filtered.size(); ++i)
        {
            auto [rot, pos] = state_manager.getPoseRotVec(temp_id.at(i));
            Mat3 temp_pose = rotPosToHomogeneous(angleToRotMat(rot[0]), pos);
            output_pose.push_back(homographyMatToVec(inv_pose_0*temp_pose));
            temp_pose = project_start ? inv_pose_0*temp_pose: invHomogeneous(invHomogeneous(temp_pose)*pose_1);
            event_reproj.push_back(se2Trans(temp_filtered.at(i), temp_pose));
            event_reproj.back()->x += sum_coor(0);
            event_reproj.back()->y += sum_coor(1);
        }

        return {event_reproj, pose, true, output_pose, sum_coor};
    }



    
    // Returns the undistorted events, the pose difference between the first and last event
    std::tuple<Mat3, bool> EventTrack::se2MotionOptimisation(
            const std::vector<EventPtr>& events, const std::vector<int>& event_ids, GpStateManager& state_manager, const double thr)
    {
        std::vector<EventPtr> fixed_events;
        std::vector<double*> state;
        std::vector<int> state_size;
        CostFunctionSe2FrontEnd* cost_function = new CostFunctionSe2FrontEnd(events, event_ids, fixed_events, state_manager, hyper_, state, false, downsample_);


        ceres::Problem problem;
        problem.AddResidualBlock(cost_function, NULL, state);

        // Add normalisers to the problem (in the case of )
        // and constant the previous states
        state_manager.addNormalisersToProblem(problem);
        state_manager.setConstantCloseTo(problem, event_ids.at(0));
        state_manager.resetUsedFlags();


        // Optimisation
        ceres::Solver::Summary summary;
        ceres::Solve(se2_opts_, &problem, &summary);
        //std::cout << summary.FullReport() << std::endl;

        bool valid = (summary.final_cost / summary.initial_cost) < thr;

        auto [rot_0, pos_0] = state_manager.getPoseRotVec(event_ids.at(0));
        auto [rot_1, pos_1] = state_manager.getPoseRotVec(event_ids.back());
        Mat3 pose_0 = rotPosToHomogeneous(angleToRotMat(rot_0[0]), pos_0);
        Mat3 pose_1 = rotPosToHomogeneous(angleToRotMat(rot_1[0]), pos_1);
        Mat3 inv_pose_0 = invHomogeneous(pose_0);
        Mat3 delta_pose = inv_pose_0 * pose_1;

        return {delta_pose, valid};
    }





    // Perfom homography registration of the given events (previously undistorted) with respect to the template (stored as private attribute)
    std::tuple<std::vector<EventPtr>, Mat3> EventTrack::discreteHomographyRegistration(
            const std::vector<EventPtr>& se2_corrected)
    {
        //ceres::ArctanLoss* loss_function = NULL;//new ceres::ArctanLoss(1);
        ceres::CauchyLoss* loss_function = new ceres::CauchyLoss(1);
        ceres::Problem problem;
        problem.AddParameterBlock(current_homography_.data(), 8);

        auto[event_filtered, temp_id] = filterGaussianDensity(se2_corrected, std::vector<int>(se2_corrected.size()), 1);

        for(int i = 0; i < event_filtered.size(); ++i)
        {
            CostFunctionDiscreteHomographyFrontEnd* cost_function = new CostFunctionDiscreteHomographyFrontEnd(event_filtered.at(i), first_template_xy_, first_template_alpha_, field_hyper_);
            problem.AddResidualBlock(cost_function, loss_function, current_homography_.data());

            CostFunctionDiscreteHomographyFrontEnd* cost_function_last = new CostFunctionDiscreteHomographyFrontEnd(event_filtered.at(i), last_template_xy_, last_template_alpha_, field_hyper_);
            problem.AddResidualBlock(cost_function_last, loss_function, current_homography_.data());
        }


        auto [ev_alpha, ev_xy] = fromEventsToGP(event_filtered, field_hyper_);
        for(int i = 0; i < first_template_xy_->cols(); ++i)
        {

            EventPtr temp_event(new Event(first_template_xy_->col(i), 0));
            CostFunctionDiscreteInvHomographyFrontEnd* cost_function = new CostFunctionDiscreteInvHomographyFrontEnd(temp_event, ev_xy, ev_alpha, field_hyper_);
            problem.AddResidualBlock(cost_function, loss_function, current_homography_.data());
        }
        for(int i = 0; i < last_template_xy_->cols(); ++i)
        {

            EventPtr temp_event(new Event(last_template_xy_->col(i), 0));
            CostFunctionDiscreteInvHomographyFrontEnd* cost_function = new CostFunctionDiscreteInvHomographyFrontEnd(temp_event, ev_xy, ev_alpha, field_hyper_);
            problem.AddResidualBlock(cost_function, loss_function, current_homography_.data());

        }


        // Optimisation
        ceres::Solver::Summary summary;
        ceres::Solve(h_opts_, &problem, &summary);
        //std::cout << summary.FullReport() << std::endl;

        std::vector<EventPtr> events_homographed = homographyProjection(se2_corrected, current_homography_);
        std::cout << "Homography\n" << homographyVecToMat(current_homography_) << std::endl << std::endl;

        return {events_homographed, homographyVecToMat(current_homography_)};


    }

    bool EventTrack::se2HomoCoherence(const Mat3& pose, const Vec2& offset)
    {
        std::vector<double> delta = {-0.75*track_half_width_, 0.75*track_half_width_};

        Vec2 center = homographyInvProjection(track_seed_, homographies_.back());

        bool output = true;
        int counter = 0;
        for(int i = 0; i < delta.size(); ++i)
        {
            for(int j = 0; j < delta.size(); ++j)
            {
                Vec2 temp = center + Vec2(delta[i], delta[j]);
                Vec2 temp_se2;
                if(first_homography_)
                {   
                    temp_se2 = se2InvTrans(temp - first_offset_, first_pose_) + first_offset_;
                    temp_se2 = se2InvTrans(temp_se2 - offset, pose) + offset;
                }
                else
                {
                    temp_se2 = se2InvTrans(temp - offset, pose) + offset;
                }
                Vec2 temp_homo = homographyInvProjection( homographyProjection(temp, homographies_.back()), current_homography_);
                if((temp_homo-temp_se2).norm() > kCoherenceThr) counter++;

                std::cout << "Coherence error " << (temp_homo-temp_se2).norm() << std::endl;
            }
        }
        std::cout << " Coherence test : " << (counter < 2) << std::endl;
        return (counter < 2);
    }

    // Create the template from undistorted events
    void EventTrack::createTemplateD()
    {
        // Create the vectors and matrices needed for the GP inferrence
        (*first_template_xy_) = d_template_.getSkeletonPts(1.0);
        int nb_pts = first_template_xy_->cols();
        MatX K = seKernelCov(*first_template_xy_, field_hyper_.l2, field_hyper_.sf2) + (field_hyper_.sz2*MatX::Identity(nb_pts, nb_pts));
        *first_template_alpha_ = solveKinvY(K, VecX::Ones(nb_pts));

        // Create a KNN search tree for later use (only picking event that are close to the template when creating homography residuals)
        first_template_points_.reserve(nb_pts);
        for(int i = 0; i < first_template_xy_->cols(); ++i) first_template_points_.push_back(first_template_xy_->col(i));
        first_template_tree_ = std::shared_ptr<cilantro::KDTree2d<> >(new cilantro::KDTree2d<>(first_template_points_));
    }






    // Create the template from undistorted events
    void EventTrack::createTemplate(const std::vector<EventPtr>& undistorted_events)
    {
        // Filter the events
        std::vector<EventPtr> temp_events = sampleCentroid(undistorted_events, centroid_dist_);

        // Create the vectors and matrices needed for the GP inferrence
        (*first_template_xy_) = eventVecToMat(temp_events);
        MatX K = seKernelCov(*first_template_xy_, field_hyper_.l2, field_hyper_.sf2) + (field_hyper_.sz2*MatX::Identity(temp_events.size(), temp_events.size()));
        *first_template_alpha_ = solveKinvY(K, VecX::Ones(temp_events.size()));

        // Create a KNN search tree for later use (only picking event that are close to the template when creating homography residuals)
        first_template_points_.reserve(undistorted_events.size());
        for(const auto& e: undistorted_events) first_template_points_.push_back(e->pixVec());
        first_template_tree_ = std::shared_ptr<cilantro::KDTree2d<> >(new cilantro::KDTree2d<>(first_template_points_));
    }


    MatX EventTrack::inferTemplate(const Mat2X& XY, const VecX& alpha)
    {
        MatX inf((int) (2*kGridFactor*track_width_), (int) (2*kGridFactor*track_width_) );

        for(int r = 0; r < 2*kGridFactor*track_width_; ++r)
        {
            for(int c = 0; c < 2*kGridFactor*track_width_; ++c)
            {
                Vec2 xy(
                    - 2*track_half_width_ + ((double)(c)/(double)(kGridFactor)),
                    - 2*track_half_width_ + ((double)(r)/(double)(kGridFactor)));
                
                MatX ks = seKernelND(xy, XY, field_hyper_.l2, field_hyper_.sf2);
                double temp = (ks*(alpha))(0,0);
                inf(r,c) = temp;
            }
        }
        return inf;
    }

    void EventTrack::displayLastTemplate()
    {
        display("Template at " + std::to_string(events_.size()), inferTemplate(*last_template_xy_, *last_template_alpha_), true);
    }
    void EventTrack::displayTemplate()
    {
        display("Template first at " + std::to_string(events_.size()), inferTemplate(*first_template_xy_, *first_template_alpha_), true);
    }

    void EventTrack::updateTemplate(const std::vector<EventPtr>& undistorted_events)
    {

        // Filter the events
        std::vector<EventPtr> temp_events = sampleCentroid(undistorted_events, centroid_dist_);

        // Create the vectors and matrices needed for the GP inferrence
        (*last_template_xy_) = eventVecToMat(temp_events);
        MatX K = seKernelCov(*last_template_xy_, field_hyper_.l2, field_hyper_.sf2) + (field_hyper_.sz2*MatX::Identity(temp_events.size(), temp_events.size()));
        *last_template_alpha_ = solveKinvY(K, VecX::Ones(temp_events.size()));

        // Create a KNN search tree for later use (only picking event that are close to the template when creating homography residuals)
        last_template_points_.clear();
        last_template_points_.reserve(undistorted_events.size());
        for(const auto& e: undistorted_events) last_template_points_.push_back(e->pixVec());
        last_template_tree_ = std::shared_ptr<cilantro::KDTree2d<> >(new cilantro::KDTree2d<>(last_template_points_));
    }

    void EventTrack::updateTemplateD()
    {
        (*last_template_xy_) = d_template_.getSkeletonPts(1.0);
        int nb_pts = last_template_xy_->cols();
        MatX K = seKernelCov(*last_template_xy_, field_hyper_.l2, field_hyper_.sf2) + (field_hyper_.sz2*MatX::Identity(nb_pts, nb_pts));
        *last_template_alpha_ = solveKinvY(K, VecX::Ones(nb_pts));

        // Create a KNN search tree for later use (only picking event that are close to the template when creating homography residuals)
        last_template_points_.clear();
        last_template_points_.reserve(nb_pts);
        for(int i = 0; i < last_template_xy_->cols(); ++i) last_template_points_.push_back(last_template_xy_->col(i));
        last_template_tree_ = std::shared_ptr<cilantro::KDTree2d<> >(new cilantro::KDTree2d<>(last_template_points_));
    }


    void EventTrack::cleanConnectiveSAE()
    {
        Vec2 temp_origin = homographyInvProjection(Vec2::Zero(), current_homography_) + origin_;
        for(int r = 0; r < sae_.rows(); ++r)
        {
            for(int c = 0; c < sae_.cols(); ++c)
            {
                Vec2 cr(c,r);
                if((cr-temp_origin).norm() > track_half_width_)
                {
                    connectivity_[r][c] = 0;
                    sae_(r,c) = -1;
                }
            }
        }
    }


    void EventTrack::displayHomographyDiscreteResults(const std::vector<EventPtr>& undistorted_events)
    {
        std::vector<EventPtr> display_events;
        for(int i = 0; i < first_template_xy_->cols(); ++i) display_events.push_back(EventPtr(new Event(first_template_xy_->col(i), -35)));
        for(int i = 0; i < last_template_xy_->cols(); ++i) display_events.push_back(EventPtr(new Event(last_template_xy_->col(i), -25)));
        for(const auto& e: undistorted_events) display_events.push_back(EventPtr(new Event(e->pixVec(), -1)));
        display_events.at(0)->t = -45;
        display(disp_str_ + " Homography results at " + std::to_string(events_.size()), display_events ,true, false);
    }


    void EventTrack::interpolateTrack(const std::vector<EventPtr>& ev, const std::vector<Vec8>& poses, const Vec2& offset)
    {
        
        double t0 = homographies_t_.at(homographies_.size()-2);
        Vec2 seed_t0 = homographyInvProjection(track_seed_, homographies_.at(homographies_.size()-2));
        Vec2 seed_t1 = homographyInvProjection(track_seed_, homographies_.back());
        Vec2 seed_prediction = seed_t0;
        if(first_homography_)
        {
            seed_prediction = se2InvTrans(seed_t0 - first_offset_, first_pose_) + first_offset_;
        }
        seed_prediction = se2InvTrans(seed_prediction - offset, homographyVecToMat(poses.back())) + offset;
        Vec2 alpha_seed = (seed_t1 - seed_prediction) / (ev.back()->t - t0);


        Vec2 seed_temp = seed_t0;
        if(first_homography_)
        {
            for(int i = 0; i < poses_temp_.size(); ++i)
            {
                Vec2 temp_seed = se2InvTrans(seed_t0 - first_offset_, homographyVecToMat(poses_temp_.at(i))) + first_offset_ + (alpha_seed*(undistorted_events_.at(i)->t - t0));
                track_.push_back(EventPtr (new Event(temp_seed, undistorted_events_.at(i)->t)));
            }
            poses_temp_.clear();

            seed_temp = se2InvTrans(seed_t0 - first_offset_, first_pose_) + first_offset_;
        }
        for(int i = 0; i < ev.size(); ++i)
        {
            Vec2 temp_seed = se2InvTrans(seed_temp - offset, homographyVecToMat(poses.at(i))) + offset + (alpha_seed*(ev.at(i)->t - t0));
            track_.push_back(EventPtr (new Event(temp_seed, ev.at(i)->t)));
        }

    }





    std::vector<EventPtr> EventTrack::filter(const std::vector<EventPtr>& events)
    {
        std::vector<EventPtr> output;
        output.push_back(events.at(0));

        std::vector<Vec2> points;
        points.reserve(events.size());
        for(const auto& e: events) points.push_back(e->pixVec());
        cilantro::KDTree2d<> tree(points);

        for(int i = 1; i < events.size(); ++i)
        {
            cilantro::NeighborSet<double> nn = tree.kNNInRadiusSearch(events.at(i)->pixVec(), 2, 1.3);
            if(nn.size() > 1)
            {
                bool keep = true;
                for(int j = 0; j < output.size(); ++j)
                {
                    if( (output.at(j)->pixVec() - points.at(nn.at(1).index)).norm() < 2.0)
                    {
                        keep = false;
                        break;
                    }
                }
                if(keep) output.push_back(events.at(nn.at(1).index));
            }
        }
        return output;
    }

    std::vector<EventPtr> EventTrack::sampleGMM(const std::vector<EventPtr>& events, const double sigma, const double nb_samples)
    {
        RandomGenerator rand_gen;
        
        std::vector<EventPtr> output;
        output.reserve(nb_samples);
        for(int i = 0; i < nb_samples; ++i)
        {
            int id = std::rand() % (events.size());
            double x_mean = events.at(id)->x;
            double y_mean = events.at(id)->y;
            output.push_back(EventPtr(new Event()));
            output.back()->x = rand_gen.randGauss(x_mean, sigma);
            output.back()->y = rand_gen.randGauss(y_mean, sigma);
            output.back()->t = 0.0;
        }
        return output;
    }

    std::vector<EventPtr> sampleCentroid(const std::vector<EventPtr>& events, const double sigma)
    {
        std::vector<EventPtr> output;

        std::vector<Vec2> points;
        points.reserve(events.size());
        for(const auto& e: events) points.push_back(e->pixVec());
        cilantro::KDTree2d<> tree(points);

        for(int i = 0; i < events.size(); ++i)
        {
            cilantro::NeighborSet<double> nn = tree.kNNInRadiusSearch(events.at(i)->pixVec(), 10, sigma);
            if(nn.size() > 3)
            {   
                Vec2 centroid(0,0);
                for(int j = 1; j < nn.size(); ++j)
                {
                    centroid += points.at(nn.at(j).index);
                }
                centroid /= (nn.size()-1);
                EventPtr temp_ev(new Event(centroid, 0.0));
                output.push_back(temp_ev);
            }
        }
        return output;
    }


    std::vector<EventPtr> EventTrack::sampleGrid(const std::vector<EventPtr>& events, const double grid_size)
    {
        std::vector<EventPtr> output;

        cilantro::PointCloud2d temp_pc;
        temp_pc.points.resize(2, events.size());
        for(int i = 0; i < events.size(); ++i) temp_pc.points.col(i) = events.at(i)->pixVec();
        cilantro::PointCloud2d down_pc(temp_pc.gridDownsample(grid_size, 3));
        output.reserve(down_pc.size());
        for(int i = 0; i < down_pc.size(); ++i)
        {
            EventPtr temp_ev(new Event());
            temp_ev->x = down_pc.points(0,i);
            temp_ev->y = down_pc.points(1,i);
            temp_ev->t = 0.0;
            output.push_back(temp_ev);
        }
        return output;
    }


    std::vector<double> getGaussianDensity(const std::vector<EventPtr>& events, const double sigma)
    {
        std::vector<double> output;
        const double neg_inv2s2 = - 1.0 / (2.0*sigma*sigma);
        const double max_radius = 4*sigma;

        std::vector<Vec2> points;
        points.reserve(events.size());
        for(const auto& e: events) points.push_back(e->pixVec());
        cilantro::KDTree2d<> tree(points);

        for(int i = 0; i < events.size(); ++i)
        {
            cilantro::NeighborSet<double> nn = tree.kNNInRadiusSearch(events.at(i)->pixVec(), 100, max_radius);
            double sum = 0;
            if(nn.size() > 0)
            {
                for(int n=1; n<nn.size(); ++n)
                {
                    double dist = (events.at(i)->pixVec() - points.at(nn.at(n).index)).norm();
                    sum += std::exp(dist*dist*neg_inv2s2);
                }
            }
            output.push_back(sum);
        }
        return output;
    }



    std::vector<EventPtr> EventTrack::gaussianColoring(const std::vector<EventPtr>& events, const double sigma)
    {
        std::vector<double> score = getGaussianDensity(events, sigma);

        std::vector<EventPtr> output;
        for(int i = 0; i < events.size(); ++i)
        {
            double temp = (score.at(i) < kDensityThr) ? -10 : score.at(i);
            output.push_back(EventPtr(new Event(events.at(i)->pixVec(), temp)));
        }
        return output;
    }


    std::tuple<std::vector<EventPtr>, std::vector<int> > filterGaussianDensity(const std::vector<EventPtr>& events, const std::vector<int>& event_ids, const double sigma)
    {
        std::vector<double> score = getGaussianDensity(events, sigma);

        std::vector<EventPtr> output;
        std::vector<int> ids;
        for(int i = 0; i < events.size(); ++i)
        {
            if(score.at(i) > kDensityThr)
            {
                output.push_back(EventPtr(new Event(events.at(i)->pixVec3())));
                ids.push_back(event_ids.at(i));
            }
        }
        return {output,ids};
    }




    // Functions to write the tracks to file
    void EventTrack::writeEventsToFile(const std::string& path_base)
    {
        eventsToCsv(path_base + "track_" + std::to_string(id_) + "_raw_events.txt", events_, origin_);
        eventsToCsv(path_base + "track_" + std::to_string(id_) + "_undistorted_events.txt", undistorted_events_, origin_);
    }

    void EventTrack::writeTrack()
    {
        if(write_path_ != "") eventsToCsv(write_path_ + "track_" + std::to_string(id_) + ".txt", track_, origin_);
    }

    void EventTrack::eventsToCsv(const std::string& file_path, const std::vector<EventPtr>& events, Vec2 offset)
    {
        std::ofstream out_stream;
        out_stream.open(file_path);

        for(int i = 0; i < events.size(); ++i)
        {
            if(i!=0) out_stream << std::endl;
            out_stream << std::to_string(events.at(i)->t) << " " << std::to_string(events.at(i)->x + offset(0)-0.5) << " " << std::to_string(events.at(i)->y + offset(1)-0.5) << " " << std::to_string(events.at(i)->pol);
        }
        out_stream.close();
    }




    void display(VisualizerPtr& viz, const std::vector<EventPtr>& events, bool flat, bool block, bool colorise)
    {
        cilantro::PointCloud3f pc;
        double normalise_to = 200;

        pc.points.resize(3, events.size());

        double t_1 = events.front()->t;
        double t_2 = events.front()->t;
        for(const auto e:events)
        {
            if(e->t < t_1) t_1 = e->t;
            if(e->t > t_2) t_2 = e->t;
        }
        double t_range = t_2 - t_1;

        std::vector<float> scalars(events.size());

        for(int i = 0; i < events.size(); ++i)
        {
            pc.points(0,i) = (events.at(i)->x) / normalise_to;
            pc.points(1,i) = (events.at(i)->y) / normalise_to;
            double temp = (events.at(i)->t - t_1)/t_range;
            if(!flat)
            {
                pc.points(2,i) = temp;
            }
            else
            {
                pc.points(2,i) = 0;
            }
            
            scalars.at(i) = temp;
        }

    }

    void clearDisplay(VisualizerPtr& viz)
    {
        viz->clearRenderArea();
        viz->addObject<cilantro::CoordinateFrameRenderable>("axis", Eigen::Matrix4f::Identity(), 0.4f, cilantro::RenderingProperties().setLineWidth(5.0f));
    }

    void display(std::string window_prefix, const std::vector<EventPtr>& events, bool flat, bool block, bool colorise)
    {
        const std::string window_name = window_prefix;
        pangolin::CreateWindowAndBind(window_name, 640, 480);
        VisualizerPtr viz(new cilantro::Visualizer(window_name, window_name));
        clearDisplay(viz);
        display(viz, events, flat, block, colorise);
    }

    void waitOnVisualizers(std::vector<VisualizerPtr>& vizs)
    {
        bool loop = true;
        while (loop)
        {
            for(auto& v : vizs)
            {
                v->spinOnce();
                if(v->wasStopped()) loop = false;
            }
        }
    }

    void display(std::string window_prefix, const MatX& mat, bool flat, bool block, bool colorise)
    {
        cilantro::PointCloud3f pc;
        double normalise_to = mat.cols()*4;

        pc.points.resize(3, mat.cols()*mat.rows());

        std::vector<float> scalars(mat.cols()*mat.rows());

        double min_val = mat.minCoeff();
        double range_val = mat.maxCoeff() - min_val;

        for(int i = 0; i < mat.cols(); ++i)
        {
            for(int j = 0; j < mat.rows(); ++j)
            {
                int index = i*mat.rows() + j;
                pc.points(0,index) = (i / normalise_to) - 0.125;
                pc.points(1,index) = (j / normalise_to) - 0.125;
                if(!flat)
                {
                    pc.points(2,index) = (mat(j,i) - min_val)/range_val;
                }
                else
                {
                    pc.points(2,index) = 0;
                }
                
                scalars.at(index) = (mat(j,i) - min_val)/range_val;
            }
        }

        const std::string window_name = window_prefix;
        pangolin::CreateWindowAndBind(window_name, 640, 480);
        cilantro::Visualizer viz(window_name, window_name);
        viz.addObject<cilantro::CoordinateFrameRenderable>("axis", Eigen::Matrix4f::Identity(), 0.1f, cilantro::RenderingProperties().setLineWidth(2.0f));
        if(colorise)
        {
            viz.addObject<cilantro::PointCloudRenderable>("pcd", pc, cilantro::RenderingProperties().setColormapType(cilantro::ColormapType::JET))->setPointValues(scalars);
        }
        else
        {
            viz.addObject<cilantro::PointCloudRenderable>("pcd", pc, cilantro::RenderingProperties().setPointColor(0,0,0));
        }
        viz.setCameraPose(0,0,-0.4, 0,0,0, 0,-1,0);
        if(block)
        {
            while (!viz.wasStopped())
            {
                viz.spinOnce();
            }
        }
        else
        {
            viz.spinOnce();
        }
    }



    std::tuple<VecXPtr, Mat2XPtr> fromEventsToGP(const std::vector<EventPtr>& events, const HyperParam& hyper)
    {
        int nb_ev = events.size();
        Mat2XPtr ev_xy(new Mat2X(2,nb_ev));
        for(int i = 0; i < nb_ev; ++i) ev_xy->col(i) = events.at(i)->pixVec();

        MatX K = seKernelCov(*ev_xy, hyper.l2, hyper.sf2) + (hyper.sz2*MatX::Identity(nb_ev, nb_ev));
        VecXPtr alpha(new VecX);
        *alpha = solveKinvY(K, VecX::Ones(nb_ev));

        return {alpha, ev_xy};
    }

    std::tuple<Mat3, bool> frameToFrameHomography(const std::vector<EventPtr>& target_ev, const std::vector<EventPtr>& source_ev, const HyperParam& hyper, const Vec8& homo_prior)
    {
        // Create the loss function and the state to estimate
        ceres::CauchyLoss* loss_function = new ceres::CauchyLoss(2);
        Vec8 local_homo = homo_prior;
        ceres::Problem pb;

        // Filter the events
        auto[target_ev_temp, target_id] = filterGaussianDensity(target_ev, std::vector<int>(target_ev.size()), 1);
        target_ev_temp = sampleCentroid(target_ev_temp, 1.2);

        auto[source_ev_temp, source_id] = filterGaussianDensity(source_ev, std::vector<int>(source_ev.size()), 1);
        source_ev_temp = sampleCentroid(source_ev_temp, 1.2);


        // Get the GP out of the events
        auto[target_alpha, target_xy] = fromEventsToGP(target_ev_temp, hyper);
        auto[source_alpha, source_xy] = fromEventsToGP(source_ev_temp, hyper);

        // Add the residuals of the source
        for(int i = 0; i < source_ev_temp.size(); ++i)
        {
            CostFunctionDiscreteHomographyFrontEnd* cost_function = new CostFunctionDiscreteHomographyFrontEnd(source_ev_temp.at(i), target_xy, target_alpha, hyper);
            pb.AddResidualBlock(cost_function, loss_function, local_homo.data());
        }

        // Add the residuals of the target 
        for(int i = 0; i < target_ev_temp.size(); ++i)
        {
            CostFunctionDiscreteInvHomographyFrontEnd* cost_function = new CostFunctionDiscreteInvHomographyFrontEnd(target_ev_temp.at(i), source_xy, source_alpha, hyper);
            pb.AddResidualBlock(cost_function, loss_function, local_homo.data());
        }


        // Optimisation
        ceres::Solver::Summary summary;
        ceres::Solver::Options pb_opts;
        pb_opts.minimizer_progress_to_stdout = false;
        pb_opts.num_threads = 14;
        ceres::Solve(pb_opts, &pb, &summary);
        //std::cout << summary.FullReport() << std::endl;


        // Check if estimation is valid
        bool valid = true;
        double ratio_source = (seKernelND(homographyMatProjection(*source_xy, local_homo), *target_xy, hyper.l2, hyper.sf2)*(*target_alpha)).sum() / source_ev_temp.size();
        if(ratio_source < 0.65) valid = false;
        double ratio_target = (seKernelND(homographyInvMatProjection(*target_xy, local_homo), *source_xy, hyper.l2, hyper.sf2)*(*source_alpha)).sum() / target_ev_temp.size();
        if(ratio_target < 0.65) valid = false;


        return {homographyVecToMat(local_homo), valid};
    }



} // namespace celib
