/**
 *  Author: Cedric LE GENTIL
 *
 *  Copyright 2021 Cedric LE GENTIL
 *
 *  This is a simple example of how to generate the UGPMs
 *  
 *  Disclaimer:
 *  This code is not optimised neither for performance neither for maintainability.
 *  There might still be errors in the code, logic or procedure. Every feedback is welcomed.
 *
 *  For any further question or if you see any problem in that code
 *  le.gentil.cedric@gmail.com
 **/


#include <iostream>
#include <string>
#include <yaml-cpp/yaml.h>

#include <boost/program_options.hpp>

#include "sensor_input/event_reader.h"
#include "event_gp/se2_frontend.h"


int main(int argc, char* argv[]){

    // Program options
    boost::program_options::options_description opt_description("Allowed options");
    opt_description.add_options()
        ("help,h", "Produce help message")
        ("config,c", boost::program_options::value< std::string >(), "Absolute path to the config file");
        ;

    boost::program_options::variables_map var_map;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, opt_description), var_map);
    boost::program_options::notify(var_map);    

    if(var_map.count("help")) {
        std::cout << opt_description << std::endl;
        return 1;
    }

    YAML::Node config;
    if(var_map.count("config"))
    {
        config = YAML::LoadFile(var_map["config"].as<std::string>());
    }
    else
    {
        std::cout << "Error, no config file have been provided. Terminating now." << std::endl;
        return 0;
    }


    celib::EventReader event_reader(config);

    celib::Se2GpFrontEnd gp_frontend(config);

    bool loop = true;
    while(loop)
    {
        celib::EventPtr event_temp = event_reader.getEvent();
        if(event_temp)
        {
            gp_frontend.addEvent(event_temp);
        }
        else
        {
            std::cout << "No more event in the data" << std::endl;
            return 0;
        }
    }

//    // As IN2LAAMA stores time in `double`, the reader
//    // probably need to offset the data timestamps to not lose precision
//    // (these data readers and data structures are built accrodingly,
//    // but that's something to be aware of if using the In2laama object
//    // in a different context)
//    celib::ImuReader imu_reader(config);
//    celib::LidarReader lidar_reader(config, imu_reader.getOffset()); 
//
//    double start_time = std::max(imu_reader.firstTimestamp(), lidar_reader.firstTimestamp());
//    
//
//    //celib::LidarDataPtr lidar = lidar_reader.get(start_time+start_offset);
//    celib::In2laama in2laama(config, start_time);
//
//
//    bool loop = true;
//
//    while(loop)
//    {
//        double imu_start, imu_end, lidar_end;
//        loop = in2laama.getDataWindows(imu_start, imu_end, lidar_end);
//
//        if(loop)
//        {
//            auto imu_data = imu_reader.get(imu_start, imu_end);
//            auto lidar_data = lidar_reader.get(lidar_end);
//
//            if(lidar_data && imu_data)
//            {
//                in2laama.processFrame(imu_data, lidar_data);
//            }
//            else
//            {
//                loop = false;
//                in2laama.finalise();
//            }
//        }
//    }




    return 0;
}
