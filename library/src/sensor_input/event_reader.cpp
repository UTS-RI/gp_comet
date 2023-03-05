/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/


#include "sensor_input/event_reader.h"
#include "common/internal_utils.h"
#include <iostream>

#include <stdexcept>

namespace celib
{


    EventReader::EventReader(YAML::Node& config)
    {
        std::string reader_type = readRequiredField<std::string>(config, "event_data_format");


        std::cout << "----- ImuReader: Loading the IMU data" << std::endl;


        if(reader_type == "csv")
        {

            // Open the rosbag
            event_file_.open(readRequiredField<std::string>(config, "event_path"), std::ios::in);
            if(!event_file_.is_open())
            {
                printErrorValueOption("event_path");
            }

            std::ifstream calib_file;
            calib_file.open(readRequiredField<std::string>(config, "calib_path"), std::ios::in);
            std::string line;
            if(getline(calib_file, line))
            {
                std::vector<double> values = stringToNumVector<double>(line, ' ');
                if(values.size() != 9) throw std::runtime_error("Calibration file is not with the right format");
                cam_mat_ = Mat3::Zero();
                cam_mat_(0,0) = values.at(0);
                cam_mat_(1,1) = values.at(1);
                cam_mat_(0,2) = values.at(2);
                cam_mat_(1,2) = values.at(3);
                k1_ = values.at(4);
                k2_ = values.at(5);
                p1_ = values.at(6);
                p2_ = values.at(7);
                k3_ = values.at(8);
            }
            else
            {
                throw std::runtime_error("Calibration file seems to be empty");
            }

        }
        else
        {
            printErrorValueOption("event_data_format");
        }

    }

    EventPtr EventReader::getEvent()
    {
        std::string line;
        if(getline(event_file_, line))
        {
            const double& fx = cam_mat_(0,0);
            const double& fy = cam_mat_(1,1);
            const double& cx = cam_mat_(0,2);
            const double& cy = cam_mat_(1,2);

            EventPtr event(new Event());
            
            std::vector<double> values = stringToNumVector<double>(line, ' ');
            event->t = values[0];
            double x = (values[1] - cx)/fx;
            double y = (values[2] - cy)/fy;
            event->pol = (values[3] > 0.5);

            double x0 = x;
            double y0 = y;

            for(int i = 0; i < 3; ++i)
            {
                double r2 = (x*x) + (y*y);
                double k_inv = 1.0 / (1 + (k1_*r2) + (k2_*r2*r2) + (k3_*r2*r2*r2));
                double delta_x = (2*p1_*x*y) + (p2_*(r2+(2*x*x)));
                double delta_y = (p1_*(r2+(2*y*y))) + (2*p2_*x*y);
                x = (x0 - delta_x) * k_inv;
                y = (y0 - delta_y) * k_inv;
            }

            event->x = x*fx + cx;
            event->y = y*fy + cy;
            return event;
        }
        else
        {
            std::cout << "Event Reader reached end of files" << std::endl;
            return nullptr;
        }
    }


} // End namespace