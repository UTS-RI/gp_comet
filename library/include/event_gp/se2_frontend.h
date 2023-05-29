/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef SE2_FRONTEND_H
#define SE2_FRONTEND_H

#include <yaml-cpp/yaml.h>
#include <unordered_set>

#include "common/types.h"
#include "common/math_utils.h"
#include "common/gp_state_manager.h"
#include "common/random.h"
#include "cost_function.h"

#include <opencv2/opencv.hpp>

// FOR DEBUG
#include "common/utils.h"
#include "cilantro/utilities/point_cloud.hpp"
#include "cilantro/visualization/visualizer.hpp"
#include "cilantro/visualization/common_renderables.hpp"


namespace celib
{

    const int kGridFactor = 4;
    const double kAssociationTimeThr = 0.05;
    const int kMaxConnectivity = 2;
    const double kDensityThr = 1;
    const double kSe2CostValidThr = 0.90;
    const double kEventStd = 0.25;
    const double kCoherenceThr = 10;


    enum TrackerType{kSE2, kHF2F, kHF2T};

    std::vector<EventPtr> saeToEvents(const MatX& sae);

    typedef std::shared_ptr<cilantro::Visualizer> VisualizerPtr;
    // Helper functions for displaying events 
    void display(std::string window_prefix, const MatX& mat, bool flat = false, bool block = false, bool colorise = true);

    void display(std::string window_prefix, const std::vector<EventPtr>& events, bool flat = false, bool block = false, bool colorise = true);
    void display(VisualizerPtr& viz, const std::vector<EventPtr>& events, bool flat = false, bool block = false, bool colorise = true);
    void clearDisplay(VisualizerPtr& viz);
    void waitOnVisualizers(std::vector<VisualizerPtr>& vizs);




    class DynamicTemplate
    {
        private:
            const int kNbSigma = 3;
            const double kMinFrame = 5.0;
            const double kMinFrameRatio = 0.7;
            int nb_frames_ = 0;
            int last_frame_ = -1;
            // Size of the matrix/arrays to maintain
            int size_;
            double width_;
            double ratio_;
            double ratio_inv_;
            double offset_;
            double half_res_;
            // Parameter of the gaussian for the occupancy smoothing
            double inv_sigma2_;
            // For each "pixel" stores the list of frame ids
            std::vector<std::vector< std::vector<int > > > ev_collec_;
            // Stores the occucancy of the events
            MatX occupancy_;
            MatX frame_occupancy_;
            MatX skeleton_;


            bool isPixIn(int r, int c)
            {
                return ( (r >= 0)&&(r < size_)&&(c >= 0)&&(c < size_) );
            }

            std::tuple<int, int> coor2Pix(const Vec2& v)
            {
                int r = (v(1) * ratio_) + offset_;
                int c = (v(0) * ratio_) + offset_;
                return {r,c};
            }
             
            Vec2 pix2Coor(int r, int c)
            {
                return Vec2( (ratio_inv_*(c - offset_)) + half_res_, (ratio_inv_*(r - offset_)) + half_res_);
            }

            MatX skeletonise(MatX& input, double thr)
            {
                MatX output(input.rows(), input.cols());
                cv::Mat img(input.rows(), input.cols(), CV_8UC1, cv::Scalar(0));
                for(int r = 0; r < input.rows(); ++r)
                {
                    for(int c = 0; c < input.cols(); ++c)
                    {
                        if(input(r,c) >= thr)
                        {
                            img.at<uchar>(r,c) = 255;
                        }
                    }
                }
                cv::Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
                cv::Mat temp;
                cv::Mat eroded;
                cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
                bool done;		
                do
                {
                    cv::erode(img, eroded, element);
                    cv::dilate(eroded, temp, element); // temp = open(img)
                    cv::subtract(img, temp, temp);
                    cv::bitwise_or(skel, temp, skel);
                    eroded.copyTo(img);
                    done = (cv::countNonZero(img) == 0);
                } while (!done);
                
                for(int r = 0; r < input.rows(); ++r)
                {
                    for(int c = 0; c < input.cols(); ++c)
                    {
                        output(r,c) = skel.at<uchar>(r,c) > 127 ? 1.0 : 0.0;
                    }
                }
                return output;
            }

            // Grid size zero for no downsampling
            Mat2X getSkeletonPts(const MatX& input, const double downsample_grid_size = 0)
            {
                int nb_pts = input.size() - std::count(input.data(), input.data()+input.size(), 0.0);
                cilantro::PointCloud2d temp_pc;
                temp_pc.points.resize(2, nb_pts);
                int counter = 0;
                for(int r = 0; r < input.rows(); ++r)
                {
                    for(int c = 0; c < input.cols(); ++c)
                    {
                        if(input(r,c)>0.0)
                        {
                            temp_pc.points.col(counter) = pix2Coor(r,c);
                            counter++;
                        }
                    }
                }
                if(downsample_grid_size > 0)
                {
                    temp_pc = temp_pc.gridDownsample(downsample_grid_size, 1);
                }
                return temp_pc.points;
            }

        public:
            DynamicTemplate(const int width, const int res_factor)
                : size_(2*width*res_factor)
                , width_(2.0*width)
                , ratio_(size_/width_)
                , ratio_inv_(1.0/ratio_)
                , offset_(((double)size_)/2.0)
                , half_res_(1.0/(2.0*res_factor))
                , inv_sigma2_(1.0/std::pow(2.0 / (double)res_factor,2))
                , ev_collec_(size_, std::vector<std::vector<int> >(size_))
                , occupancy_(size_, size_)
                , frame_occupancy_(size_, size_)
                , skeleton_(size_, size_)
            {
                occupancy_.setZero();
                frame_occupancy_.setZero();
            }

            void addEvent(const EventPtr& ev, const int frame_id)
            {
                auto [r,c] = coor2Pix(ev->pixVec());
                if(isPixIn(r,c))
                {   
                    bool new_frame_o = false;
                    if(last_frame_ != frame_id)
                    {
                        nb_frames_++;
                        last_frame_ = frame_id;
                    }
                    if((ev_collec_[r][c].empty())||(frame_id != ev_collec_[r][c].back()))
                    {
                        ev_collec_[r][c].push_back(frame_id);
                        new_frame_o = true;
                    }
                    for(int i = -kNbSigma; i <= kNbSigma; ++i)
                    {
                        for(int j = -kNbSigma; j <= kNbSigma; ++j)
                        {
                            int temp_r = r+i;
                            int temp_c = c+j;
                            if(isPixIn(temp_r, temp_c))
                            {
                                Vec2 temp_v = pix2Coor(temp_r, temp_c);
                                double d2 = (temp_v - ev->pixVec()).squaredNorm();
                                double temp_incr = std::exp(-d2*inv_sigma2_);
                                occupancy_(temp_r, temp_c) += temp_incr;
                                skeleton_(temp_r, temp_c) += temp_incr;
                                if(new_frame_o) frame_occupancy_(temp_r, temp_c) += temp_incr;
                            }
                        }
                    }
                }
            }
            void addFrame(const std::vector<EventPtr>& evs, const int frame_id)
            {
                bool copy = (nb_frames_ == 0);
                for(const auto& e: evs) addEvent(e, frame_id);
                if(copy) skeleton_ = occupancy_;
            }

            std::tuple<Vec2, bool> getMostViewed()
            {
                int max_r = 0;
                int max_c = 0;
                double max_occup = occupancy_(max_r, max_c);
                //int nb_observed = 0;
                for(int r = 0; r < size_; ++r)
                {
                    for(int c = 0; c < size_; ++c)
                    {
                        if(occupancy_(r, c) >= max_occup) 
                        {
                            max_r = r;
                            max_c = c;
                            max_occup = occupancy_(r,c);
                        }
                    }
                }
                //return {pix2Coor(max_r, max_c), (nb_observed != 0)};
                return {pix2Coor(max_r, max_c), true};
            }

            void display(std::string str = "", bool block = true)
            {
                celib::display(str + "Dynamic template at frame " + std::to_string(nb_frames_), occupancy_, true, block);
            }

            void displayFrameOccupancy(std::string str = "", bool block = true)
            {
                celib::display(str + "Dynamic frame at " + std::to_string(nb_frames_), frame_occupancy_, true, block);
            }
            void displayFrameOccupancySkeleton(std::string str = "", bool block = true)
            {
                celib::display(str + "Dynamic frame at " + std::to_string(nb_frames_), skeletonise(frame_occupancy_, std::min(kMinFrameRatio*nb_frames_, kMinFrame)), true, block);
            }


            // Grid size zero for no downsampling
            Mat2X getSkeletonPts(const double downsample_grid_size = 0)
            {
                return getSkeletonPts(skeletonise(frame_occupancy_, std::min(kMinFrameRatio*nb_frames_, kMinFrame)), downsample_grid_size);
            }
            
    };





    class EventTrack
    {
        public:
            EventTrack(
                const int id,
                const double x,
                const double y,
                const double start_t,
                const double track_width,
                const double image_width,
                const double image_height,
                const int nb_state_opt = 10,
                const int state_every = 30,
                const int downsample = 0,
                const bool visualise = true,
                const TrackerType type = TrackerType::kHF2T,
                const std::string write_path = "",
                const double limit_track_length = std::numeric_limits<double>::max(),
                const bool recenter_seeds = false,
                const bool only_one = false,
                const bool write_events = false,
                const bool approx_only = false
                );

            bool isTrackInImage();
            void writeEventsToFile(const std::string& path_base);

            bool addEvent(const EventPtr event, const double add_noise = 0);

            // Function to visualise the event with and without trajectory correction
            void displayRawProjection(std::string window_name);
            void displayProjection(std::string window_name);

            int getId(){return id_;};

        private:
            bool isEventInTrack(const EventPtr& event);

            bool optimise();
            void eventsToCsv(const std::string& file_path, const std::vector<EventPtr>& events, const Vec2 offest = Vec2::Zero());

            void createTemplate(const std::vector<EventPtr>& undistorted_events);
            void createTemplateD();
            MatX inferTemplate(const Mat2X& XY, const VecX& alpha);
            void displayLastTemplate();
            void displayTemplate();
            void updateTemplate(const std::vector<EventPtr>& undistorted_events);
            void updateTemplateD();

            void addSe2Costs(ceres::Problem& problem);

            void writeTrack();

            std::vector<EventPtr> filter(const std::vector<EventPtr>& events);
            std::vector<EventPtr> sampleGMM(const std::vector<EventPtr>& events, const double sigma, const double nb_samples);
            std::vector<EventPtr> sampleGrid(const std::vector<EventPtr>& events, const double grid_size);
            std::vector<EventPtr> gaussianColoring(const std::vector<EventPtr>& events, const double sigma);

            void displayHomographyDiscreteResults(const std::vector<EventPtr>& undistorted_events);

            void cleanConnectiveSAE();

            // Return the last batch of events as well as the corresponding indexes
            std::tuple<std::vector<EventPtr>, std::vector<int> > getLastEventBatch();

            // Returns the undistorted events, the pose difference between the first and last event
            // Needs the events, the ids, whether it porjects event at the start of the end, and is it uses a field prior (from template in private attributes)
            std::tuple<std::vector<EventPtr>, Mat3, bool, std::vector<Vec8>, Vec2 > se2MotionCorrection(
                    const std::vector<EventPtr>& events, const std::vector<int>& event_ids, bool project_start = false, bool project_all = false);

            // Low level function for se2MotionCorrection
            std::tuple<Mat3, bool> se2MotionOptimisation(
                    const std::vector<EventPtr>& events, const std::vector<int>& event_ids, GpStateManager& state_manager, const double thr = kSe2CostValidThr);

            // Perfom homography registration of the given events (previously undistorted) with respect to the template (stored as private attribute)
            std::tuple<std::vector<EventPtr>, Mat3> discreteHomographyRegistration(
                    const std::vector<EventPtr>& events);

            void interpolateTrack(const std::vector<EventPtr>& ev, const std::vector<Vec8>& poses, const Vec2& offset);


            bool se2HomoCoherence(const Mat3& pose, const Vec2& offset);



            const int id_;

            RandomGenerator rand_gen_;
            const double track_half_width_;
            const double track_width_;
            const double image_width_;
            const double image_height_;
            const int image_half_width_;
            const int image_half_height_;
            const int nb_event_opt_;
            const int state_every_;
            const bool only_one_;
            const bool write_events_;
            const double max_length_;
            const int downsample_;
            const bool approx_only_;

            // Hyper parameters
            HyperParam hyper_;
            HyperParam field_hyper_;
            double centroid_dist_;


            bool first_ = true;
            TrackerType type_;
            

            Vec2 origin_;
            double start_time_;
            std::vector<EventPtr> events_;
            std::vector<EventPtr> undistorted_events_;
            std::vector<EventPtr> track_;
            bool visualise_;
            std::string write_path_;
            bool recenter_;


            // Template for homography
            Mat2XPtr first_template_xy_;
            VecXPtr first_template_alpha_;
            Mat2XPtr last_template_xy_;
            VecXPtr last_template_alpha_;
            std::shared_ptr<cilantro::KDTree2d<> > first_template_tree_;
            std::vector<Vec2> first_template_points_;
            std::shared_ptr<cilantro::KDTree2d<> > last_template_tree_;
            std::vector<Vec2> last_template_points_;


            // Optimisation related stuff
            ceres::Solver::Options se2_opts_;
            ceres::Solver::Options h_opts_;
            Vec8 current_homography_;
            int index_to_optimise_ = 0;

            // Surface of active events and connectivity for event "selection/clustering"
            MatX sae_;
            std::vector<std::vector<int> > connectivity_;


            // Accumulation of undistorted events_;
            DynamicTemplate d_template_;

            // Stores the homography estimates
            std::vector<Vec8> homographies_;
            std::vector<double> homographies_t_;
            Mat3 first_pose_;
            std::vector<Vec8> poses_temp_;
            Vec2 first_offset_;
            int frame_id_ = 0;
            Vec2 track_seed_;

            std::vector<EventPtr> previous_se2_;
            bool first_homography_ = true;


            std::string disp_str_;

            std::vector<VisualizerPtr> vizs_;
            
            


            double cutoff_;
            int cutoff_state_;
            int nb_state_opt_;
    };

    typedef std::shared_ptr<EventTrack> EventTrackPtr;


    class Se2GpFrontEnd
    {

        public:
            Se2GpFrontEnd(YAML::Node& config);

            void addEvent(const EventPtr& event);

        private:
            
            void displaySae();

            void startTrack(const double t, const double x, const double y, const int id);
            bool isEventInImage(const EventPtr& event);

            void processTracks();

            // Load a track file (output of other tracking methods) to start tracks at same locations
            void loadTrackFile(std::string path, std::string time_offeset);

            EventTrack getTracks();


            double track_width_;
            double image_width_;
            double image_height_;
            double nb_event_for_opt_;
            bool filter_;
            double filter_thr_;
            double max_track_length_;

            int nb_state_opt_;
            int state_every_;
            int lengthscale_factor_;
            int nb_lengthscale_cutoff_;
            double event_noise_;
            int downsample_;

            bool visualise_;
            bool recenter_seeds_;

            Buffer<EventPtr> event_buffer_;
            Buffer<EventTrackPtr> track_buffer_;

            std::set<EventTrackPtr> active_tracks_;

            MatX sae_;
            VisualizerPtr viz_;


            std::shared_ptr<std::map<int, MatX> > K_inv_;

            HyperParam hyper_;

            std::string result_path_;
            bool only_one_;
            bool write_events_;
            bool approx_only_;

            // Preloaded features (start of tracks)
            std::vector<Vec2> feature_px_list_;
            std::vector<int> feature_ids_;
            std::vector<double> feature_t_list_;
            double feature_list_ptrs_ = 0;

            TrackerType tracker_type_;

            double use_tracks_after_;
            double use_tracks_before_;

    };




    std::tuple<VecXPtr, Mat2XPtr> fromEventsToGP(const std::vector<EventPtr>& events, const HyperParam& hyper);

    std::tuple<Mat3, bool> frameToFrameHomography(const std::vector<EventPtr>& target_ev, const std::vector<EventPtr>& source_ev, const HyperParam& hyper, const Vec8& homo_prior);


    std::tuple<std::vector<EventPtr>, std::vector<int> > filterGaussianDensity(const std::vector<EventPtr>& events, const std::vector<int>& event_ids, const double sigma);

    std::vector<EventPtr> sampleCentroid(const std::vector<EventPtr>& events, const double sigma);

    std::vector<double> getGaussianDensity(const std::vector<EventPtr>& events, const double sigma);







    inline std::shared_ptr<cilantro::Visualizer> createVisualizer(std::string window_name)
    {
        pangolin::CreateWindowAndBind(window_name, 640, 480);
        std::shared_ptr<cilantro::Visualizer> viz(new cilantro::Visualizer(window_name, "disp"+window_name));
        viz->spinOnce();
        return viz;
    }



    // Square exponential kernel N dimension input (self cov: aka same x on both axis)
    // X is N*M with M the number of samples
    inline MatX seKernelND(const MatX& X1, const MatX& X2, const double l2, const double sf2)
    {
        int nb_2 = X2.cols();
        int nb_1 = X1.cols();
        MatX K(nb_1, nb_2);
        for(int c = 0; c < nb_2; ++c)
        {
            for(int r = 0; r < nb_1; ++r)
            {
                VecX temp = X1.col(r) - X2.col(c);
                K(r,c) = temp.squaredNorm();
            }
        }
        return ((K*(-0.5/l2)).array().exp() *sf2).matrix();
    }

    // Square exponential kernel N dimension input (self cov: aka same x on both axis)
    // X is N*M with M the number of samples
    inline MatX seKernelCov(const MatX& X, const double l2, const double sf2)
    {
        int nb_data = X.cols();
        MatX K(nb_data, nb_data);
        #pragma omp parallel for
        for(int c = 0; c < nb_data; ++c)
        {
            for(int r = c; r < nb_data; ++r)
            {
                VecX temp = X.col(c) - X.col(r);
                K(r,c) = temp.squaredNorm();
                if(r!=c) K(c,r) = K(r,c);
            }
        }
        return ((K*(-0.5/l2)).array().exp() *sf2).matrix();
    }

    inline VecX seKernelNDJacobians(const double k, const VecX& x_r, const VecX& x_c ,const double l2)
    {
        int nb_dim = x_r.size();
        VecX output(2*nb_dim);
        VecX x_diff = x_r - x_c;
        double k_temp = k / l2;
        for(int i = 0; i < x_r.size(); ++i)
        {
            double temp  = k_temp*x_diff(i);
            output(i) = -temp;
            output(nb_dim + i) = temp;
        }

        return output;
    }

    inline Row4 seKernelJacobians(const Vec2& x_r, const Vec2& x_c ,const double l2, const double sf2)
    {
        Row4 output;
        Vec2 x_diff = x_c - x_r;
        double k = sf2*std::exp( - (x_diff.squaredNorm())/(2.0*l2));
        output <<
             2.0*k*x_diff(0)/(2.0*l2),
             2.0*k*x_diff(1)/(2.0*l2),
            -2.0*k*x_diff(0)/(2.0*l2),
            -2.0*k*x_diff(1)/(2.0*l2);

        return output;
    }
    inline Mat2X ksJacobian(const VecX& k, const Vec2& x_r, const Mat2X& x_c ,const double l2)
    {
        Mat2X output(2,x_c.cols());
        for(int i = 0; i < x_c.cols(); ++i)
        {
            Vec2 x_diff = x_r - x_c.col(i);
            output(0, i) = -k(i)*x_diff(0)/l2;
            output(1, i) = -k(i)*x_diff(1)/l2;
        }
        return output;
    }

    inline void testKsJacobian()
    {
        double quantum = 0.0001;
        Mat2X pattern(2,10);
        pattern.setRandom();
        Vec2 xy = Vec2::Random();
        double l2 = std::abs(xy(0));
        double sf2 = std::abs(xy(1));
        xy = Vec2::Random();


        MatX ks = seKernelND(xy, pattern, l2, sf2);

        Mat2X jacobian = ksJacobian(ks.row(0), xy, pattern, l2);
        Mat2X num_jacobian(2,10);
        for(int i = 0; i < 2; ++i)
        {
            Vec2 xy_shift = xy;
            xy_shift(i) += quantum;

            MatX ks_shift = seKernelND(xy_shift, pattern, l2, sf2);
            
            num_jacobian.row(i) = (ks_shift - ks)/quantum;
        }
        std::cout << "========= TestKsJacobian" << std::endl;
        std::cout << "Analytical\n" << jacobian << std::endl;
        std::cout << "Numerical\n" << num_jacobian << std::endl;
    }



    inline void testSeKernelJacobian()
    {
        Vec2 x_1 = Vec2::Random();
        Vec2 x_2 = Vec2::Random();
        double sf2 = std::pow(x_2(0),2);
        double l2 = std::pow(x_2(1),2);
        x_2 = Vec2::Random();
        double quantum = 0.001;

        double k = sf2*std::exp( - (x_1-x_2).squaredNorm()/(2*l2));
        Row4 d_k = seKernelJacobians(x_1, x_2, l2, sf2);

        Row4 d_k_num;

        for(int i = 0; i < 2; ++i)
        {
            Vec2 x_1_diff = x_1;
            x_1_diff(i) += quantum;
            double k_1_diff = sf2*std::exp( - (x_1_diff-x_2).squaredNorm()/(2*l2));
            d_k_num(i) = (k_1_diff - k) / quantum;

            Vec2 x_2_diff = x_2;
            x_2_diff(i) += quantum;
            double k_2_diff = sf2*std::exp( - (x_1-x_2_diff).squaredNorm()/(2*l2));
            d_k_num(2+i) = (k_2_diff - k) / quantum;
        }

        std::cout << "========= TestSeKernelJacobian" << std::endl;
        std::cout << "Analytical " << d_k << std::endl;
        std::cout << "Numerical  " << d_k_num << std::endl;
    }

    inline Vec2 se2InvTrans(const Vec2& event, const Vec2& pos, const double rot)
    { 
        Vec2 output = angleToRotMat(rot).transpose() * (event - pos);
        return output;
    }
    inline EventPtr se2InvTrans(const EventPtr& event, const Vec2& pos, const double rot)
    { 
        EventPtr output(new Event());
        Vec2 temp = se2InvTrans(event->pixVec(), pos, rot);
        output->x = temp(0);
        output->y = temp(1);
        output->t = event->t;
        output->pol = event->pol;
        return output;
    }



    inline Vec2 se2Trans(const Vec2& event, const Vec2& pos, const double rot)
    { 
        Mat2 rot_mat = angleToRotMat(rot);
        return (rot_mat*event) + pos;
    }
    inline EventPtr se2Trans(const EventPtr& event, const Vec2& pos, const double rot)
    { 
        EventPtr output(new Event());
        Vec2 temp = se2Trans(event->pixVec(), pos, rot);
        output->x = temp(0);
        output->y = temp(1);
        output->t = event->t;
        output->pol = event->pol;
        return output;
    }


    inline Vec2 se2Trans(const Vec2& event, const Mat3& pose)
    {
        Vec3 temp_e;
        temp_e << event(0), event(1), 1.0;
        temp_e = pose*temp_e;
        return temp_e.segment<2>(0);
    }
    inline EventPtr se2Trans(const EventPtr& event, const Mat3& pose)
    {
        Vec2 temp_xy = se2Trans(event->pixVec(), pose);
        EventPtr output(new Event);
        output->x = temp_xy(0);
        output->y = temp_xy(1);
        output->t = event->t;
        output->pol = event->pol;
        return output;
    }


    inline Vec2 se2InvTrans(const Vec2& event, const Mat3& pose)
    {
        return se2Trans(event, invHomogeneous(pose));
    }
    inline EventPtr se2InvTrans(const EventPtr& event, const Mat3& pose)
    {
        return se2Trans(event, invHomogeneous(pose));
    }



    inline std::tuple<Vec2, Mat2, Vec2> se2TransJacobian(const Vec2& event, const Vec2& pos, const double rot)
    { 
        Mat2 rot_mat = angleToRotMat(rot);
        Vec2 temp = rot_mat*event;
        return {temp + pos, Mat2::Identity(), Vec2(-temp(1), temp(0))};
    }

    inline void testSe2TransJacobian()
    {
        Vec2 event = Vec2::Random();
        Vec2 pos = Vec2::Random();
        double rot = pos(0);
        pos = Vec2::Random();

        auto [xy, pos_J, rot_J] = se2TransJacobian(event, pos, rot);

        double quantum = 0.001;
        Vec2 xy_shift = se2Trans(event, pos, rot+quantum);
        

        std::cout << "========= TestSe2TransJacobian" << std::endl;
        std::cout << "Rot analytical" << std::endl << rot_J << std::endl;
        std::cout << "Rot numerical" << std::endl << (xy_shift-xy)/quantum << std::endl;

        Mat2 pos_J_num;
        Vec2 pos_shift = pos;
        pos_shift(0) += quantum;
        xy_shift = se2Trans(event, pos_shift, rot);
        pos_J_num.col(0) = (xy_shift - xy)/quantum;
        pos_shift = pos;
        pos_shift(1) += quantum;
        xy_shift = se2Trans(event, pos_shift, rot);
        pos_J_num.col(1) = (xy_shift - xy)/quantum;

        std::cout << "Pos analytical" << std::endl << pos_J << std::endl;
        std::cout << "Pos numerical" << std::endl << pos_J_num << std::endl;
    }

    inline Mat2X se2GridInvTrans(const Mat2X& grid, const Vec2& pos, const double rot)
    {
        Mat2X output = grid;
        output.colwise() -= pos;
        output = angleToRotMat(rot).transpose()*output;

        return output;
    }

    // Returns projected grid, jacobian w.r. pos, jacobian w.r. rot
    inline std::tuple<Mat2X, Mat2, Mat2X > se2GridInvTransJacobian(const Mat2X grid, const Vec2& pos, const double rot)
    {
        Mat2X d_g_d_r(2,grid.cols());
        Mat2X output = grid;
        output.colwise() -= pos;

        Mat2 rot_mat = angleToRotMat(rot);
        Mat2 d_g_d_p = - rot_mat.transpose();
        
        output = angleToRotMat(rot).transpose()*output;

        d_g_d_r.row(0) = output.row(1);
        d_g_d_r.row(1) = -output.row(0);

        return {output, d_g_d_p, d_g_d_r};
    }


    inline void testSe2GridTransJacobian()
    {
        std::cout << "========= TestSe2GridTransJacobian" << std::endl;
        int grid_size = 7;
        double quantum = 0.001;
        Mat2X grid = MatX::Random(2,7);
        Vec2 pos = Vec2::Random();
        double rot = pos(0);
        pos = Vec2::Random();

        auto [grid_proj, d_g_d_p, d_g_d_r] = se2GridInvTransJacobian(grid, pos, rot);
        
        Mat2 d_g_d_p_num;
        Vec2 pos_shift = pos;
        pos_shift(0) += quantum;
        Mat2X grid_shift = se2GridInvTrans(grid, pos_shift, rot);
        std::cout << "Numerical jacobian pos(0)\n" << (grid_shift-grid_proj)/quantum << std::endl;
        d_g_d_p_num.col(0) = (grid_shift.col(0)-grid_proj.col(0))/quantum;
        pos_shift = pos;
        pos_shift(1) += quantum;
        grid_shift = se2GridInvTrans(grid, pos_shift, rot);
        std::cout << "Numerical jacobian pos(1)\n" << (grid_shift-grid_proj)/quantum << std::endl;
        d_g_d_p_num.col(1) = (grid_shift.col(0)-grid_proj.col(0))/quantum;
        std::cout << "Numerical pos \n" << d_g_d_p_num << std::endl; 
        std::cout << "Analytical pos \n" << d_g_d_p << std::endl; 
        double rot_shift = rot + quantum;
        grid_shift = se2GridInvTrans(grid, pos, rot_shift);
        std::cout << "Numerical rot\n" << (grid_shift-grid_proj)/quantum << std::endl;
        std::cout << "Analytical rot\n" << d_g_d_r << std::endl;

    }

    inline std::tuple<Vec2, Mat2, Vec2> se2InvTransJacobian(const Vec2& event, const Vec2& pos, const double rot)
    {
        Vec2 temp = event - pos;
        Vec2 output;
        double cos_r = std::cos(rot);
        double sin_r = std::sin(rot);
        output(0) = (temp(0)*cos_r) - (temp(1)*sin_r);
        output(1) = (temp(0)*sin_r) + (temp(1)*cos_r);
        Mat2 pos_jacobian;
        Vec2 rot_jacobian;

        rot_jacobian(0) = -(temp(0)*sin_r) - (temp(1)*cos_r);
        rot_jacobian(1) = output(0);

        pos_jacobian(0,0) = - cos_r;
        pos_jacobian(0,1) = sin_r;
        pos_jacobian(1,0) = -sin_r;
        pos_jacobian(1,1) = - cos_r;


        return {output, pos_jacobian, rot_jacobian};
    }


    inline void testSe2TransInvJacobian()
    {
        Vec2 event = Vec2::Random();
        Vec2 pos = Vec2::Random();
        double rot = pos(0);
        pos = Vec2::Random();

        auto [xy, pos_J, rot_J] = se2InvTransJacobian(event, pos, rot);

        double quantum = 0.001;
        Vec2 xy_shift = se2InvTrans(event, pos, rot+quantum);
        

        std::cout << "========= TestSe2TransInvJacobian" << std::endl;
        std::cout << "Rot analytical" << std::endl << rot_J << std::endl;
        std::cout << "Rot numerical" << std::endl << (xy_shift-xy)/quantum << std::endl;

        Mat2 pos_J_num;
        Vec2 pos_shift = pos;
        pos_shift(0) += quantum;
        xy_shift = se2InvTrans(event, pos_shift, rot);
        pos_J_num.col(0) = (xy_shift - xy)/quantum;
        pos_shift = pos;
        pos_shift(1) += quantum;
        xy_shift = se2InvTrans(event, pos_shift, rot);
        pos_J_num.col(1) = (xy_shift - xy)/quantum;

        std::cout << "Pos analytical" << std::endl << pos_J << std::endl;
        std::cout << "Pos numerical" << std::endl << pos_J_num << std::endl;
    }





    inline Mat3 homographyVecToMat(const Vec8& h)
    {
        Mat3 H;
        H << h(0), h(3), h(6),
             h(1), h(4), h(7),
             h(2), h(5), 1.0;
        return H;
    }

    inline Vec8 homographyMatToVec(const Mat3& H)
    {
        Vec8 h;
        h << H(0,0), H(1,0), H(2,0), H(0,1), H(1,1), H(2,1), H(0,2), H(1,2);
        h /= H(2,2);
        return h;
    }


    inline Vec2 homographyProjection(const Vec2& xy, const Vec8& h)
    {
        Mat3 H = homographyVecToMat(h);
        Vec3 xy1(xy(0), xy(1), 1.0);
        Vec3 temp_ev = H * xy1;
        Vec2 ev = temp_ev.segment<2>(0) / temp_ev(2);
        return ev;
    }

    inline Mat2X homographyMatProjection(const Mat2X& xy, const Vec8& h)
    {
        Mat2X output(2, xy.cols());
        for(int i = 0; i < xy.cols(); ++i) output.col(i) = homographyProjection(xy.col(i),h);
        return output;
    }

    inline EventPtr homographyProjection(const EventPtr& ev, const Vec8& h)
    {
        Vec2 xy(ev->x, ev->y);
        Vec2 proj = homographyProjection(xy, h);
        EventPtr ev_out(new Event);
        ev_out->x = proj(0);
        ev_out->y = proj(1);
        ev_out->t = ev->t;
        ev_out->pol = ev->pol;
        return ev_out;
    }

    inline std::vector<EventPtr> homographyProjection(const std::vector<EventPtr>& evs, const Vec8& h)
    {
        std::vector<EventPtr> output;
        output.reserve(evs.size());
        for(const auto& e:evs) output.push_back( homographyProjection(e,h));
        return output;
    }



    inline Vec2 homographyInvProjection(const Vec2& xy, const Vec8& h)
    {
        Mat3 H = homographyVecToMat(h);
        Vec3 xy1(xy(0), xy(1), 1.0);
        Vec3 temp_ev = H.inverse() * xy1;
        Vec2 ev = temp_ev.segment<2>(0) / temp_ev(2);
        return ev;
    }

    inline Mat2X homographyInvMatProjection(const Mat2X& xy, const Vec8& h)
    {
        Mat2X output(2, xy.cols());
        for(int i = 0; i < xy.cols(); ++i) output.col(i) = homographyInvProjection(xy.col(i),h);
        return output;
    }



    inline std::tuple<Vec2, Mat2_8> homographyInvProjectionJacobian(const Vec2& xy, const Vec8& h)
    {
        Mat3 H = homographyVecToMat(h);

        Vec3 xy1(xy(0), xy(1), 1.0);

        Vec3 temp_ev = H.inverse() * xy1;

        Vec2 ev = temp_ev.segment<2>(0) / temp_ev(2);

        Mat9_8 d_H_inv_d_h;
        // (h1*h5 - h2*h4 - h1*h6*h8 + h2*h6*h7 + h3*h4*h8 - h3*h5*h7)^2
        double temp_a = std::pow(h(0)*h(4) - h(1)*h(3) - h(0)*h(5)*h(7) + h(1)*h(5)*h(6) + h(2)*h(3)*h(7) - h(2)*h(4)*h(6),2.0);
        // (h5 - h6*h8)
        double temp_b = (h(4) - h(5)*h(7));
        // (h4 - h6*h7)
        double temp_c = (h(3) - h(5)*h(6));
        // (h4*h8 - h5*h7)
        double temp_d = (h(3)*h(7) - h(4)*h(6));
        // (h2 - h3*h8)
        double temp_e = (h(1) - h(2)*h(7));
        // (h2*h6 - h3*h5)
        double temp_f = (h(1)*h(5) - h(2)*h(4));
        // (h1 - h3*h7)
        double temp_g = (h(0) - h(2)*h(6));
        // (h1*h6 - h3*h4)
        double temp_h = (h(0)*h(5) - h(2)*h(3));
        // (h1*h5 - h2*h4)
        double temp_i = (h(0)*h(4) - h(1)*h(3));
        // (h1*h8 - h2*h7)
        double temp_j = (h(0)*h(7) - h(1)*h(6));

        d_H_inv_d_h << 
            -std::pow(temp_b,2.0),(temp_c*temp_b),-(temp_d*temp_b),(temp_e*temp_b),-(temp_e*temp_c),(temp_d*temp_e),-(temp_f*temp_b),(temp_f*temp_c),(temp_e*temp_b),-(temp_g*temp_b),(temp_j*temp_b),-std::pow(temp_e,2.0),(temp_g*temp_e),-(temp_j*temp_e),(temp_f*temp_e),-(temp_f*temp_g),-(temp_f*temp_b),(temp_h*temp_b),-(temp_i*temp_b),(temp_f*temp_e),-(temp_h*temp_e),(temp_i*temp_e),-std::pow(temp_f,2.0),(temp_h*temp_f),(temp_c*temp_b),-std::pow(temp_c,2.0),(temp_d*temp_c),-(temp_g*temp_b),(temp_g*temp_c),-(temp_d*temp_g),(temp_h*temp_b),-(temp_h*temp_c),-(temp_e*temp_c),(temp_g*temp_c),-(temp_j*temp_c),(temp_g*temp_e),-std::pow(temp_g,2.0),(temp_j*temp_g),-(temp_h*temp_e),(temp_h*temp_g),(temp_f*temp_c),-(temp_h*temp_c),(temp_i*temp_c),-(temp_f*temp_g),(temp_h*temp_g),-(temp_i*temp_g),(temp_h*temp_f),-std::pow(temp_h,2.0),-(temp_d*temp_b),(temp_d*temp_c),-std::pow(temp_d,2.0),(temp_j*temp_b),-(temp_j*temp_c),(temp_j*temp_d),-(temp_i*temp_b),(temp_i*temp_c),(temp_d*temp_e),-(temp_d*temp_g),(temp_j*temp_d),-(temp_j*temp_e),(temp_j*temp_g),-std::pow(temp_j,2.0),(temp_i*temp_e),-(temp_i*temp_g),-(temp_f*temp_d),(temp_h*temp_d),-(temp_i*temp_d),(temp_f*temp_j),-(temp_h*temp_j),(temp_i*temp_j),-(temp_i*temp_f),(temp_i*temp_h);

        Mat3_9 d_vec3_d_H_inv;
        d_vec3_d_H_inv.block<3,3>(0,0) = xy(0)*Mat3::Identity();
        d_vec3_d_H_inv.block<3,3>(0,3) = xy(1)*Mat3::Identity();
        d_vec3_d_H_inv.block<3,3>(0,6) = Mat3::Identity();

        Mat2_3 d_vec2_d_vec3;
        double inv_e2 =  1.0 / temp_ev(2);
        double sq_inv_e2 = std::pow(inv_e2, 2);
        d_vec2_d_vec3 <<
            inv_e2, 0, -temp_ev(0)*sq_inv_e2,
            0, inv_e2, -temp_ev(1)*sq_inv_e2;

        return {ev, (d_vec2_d_vec3*d_vec3_d_H_inv*d_H_inv_d_h)/temp_a};
    }



    inline std::tuple<Vec2, Mat2_8> homographyProjectionJacobian(const Vec2& xy, const Vec8& h)
    {
        Mat3 H = homographyVecToMat(h);
        Vec3 xy1(xy(0), xy(1), 1.0);

        Vec3 temp_ev = H * xy1;

        Vec2 ev = temp_ev.segment<2>(0) / temp_ev(2);

        double inv_e2 =  1.0 / temp_ev(2);
        double x_inv_e2 = xy(0)*inv_e2;
        double y_inv_e2 = xy(1)*inv_e2;
        double sq_inv_e2 = std::pow(inv_e2, 2);

        Mat2_8 d_ev_d_h;
        d_ev_d_h << x_inv_e2,   0,   -xy(0)*temp_ev(0)*sq_inv_e2,  y_inv_e2,  0,   -xy(1)*temp_ev(0)*sq_inv_e2,  inv_e2, 0,
                    0,   x_inv_e2,   -xy(0)*temp_ev(1)*sq_inv_e2,  0,  y_inv_e2,   -xy(1)*temp_ev(1)*sq_inv_e2,  0, inv_e2;

        return {ev, d_ev_d_h};
    }

    inline void testHomographyInvProjectionJacobian()
    {
        Vec2 xy = Vec2::Random();

        Vec8 h = Vec8::Random();

        auto [proj, jacobian] = homographyInvProjectionJacobian(xy, h);

        Mat2_8 num_jacobian;

        double quantum = 0.001;
        for(int i = 0; i < 8; ++i)
        {
            Vec8 h_shift = h;
            h_shift(i) += quantum;
            Vec2 proj_shift = homographyInvProjection(xy, h_shift);
            num_jacobian.col(i) = (proj_shift - proj)/quantum;
        }


        std::cout << "========= TestHomographyProjectionJacobian" << std::endl;
        std::cout << "Analytical" << std::endl << jacobian << std::endl;
        std::cout << "Numerical" << std::endl << num_jacobian << std::endl;

    }



    inline void testHomographyProjectionJacobian()
    {
        Vec2 xy = Vec2::Random();

        Vec8 h = Vec8::Random();

        auto [proj, jacobian] = homographyProjectionJacobian(xy, h);

        Mat2_8 num_jacobian;

        double quantum = 0.001;
        for(int i = 0; i < 8; ++i)
        {
            Vec8 h_shift = h;
            h_shift(i) += quantum;
            Vec2 proj_shift = homographyProjection(xy, h_shift);
            num_jacobian.col(i) = (proj_shift - proj)/quantum;
        }


        std::cout << "========= TestHomographyProjectionJacobian" << std::endl;
        std::cout << "Analytical" << std::endl << jacobian << std::endl;
        std::cout << "Numerical" << std::endl << num_jacobian << std::endl;

    }



    inline Mat2X eventVecToMat(const std::vector<EventPtr>& events)
    {
        Mat2X mat(2, events.size());
        for(int i = 0; i < events.size(); ++i)
        {
            mat.col(i) = events.at(i)->pixVec();
        }
        return mat;
    }

} // namespace name









#endif