/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2022 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef EVENT_READER_H
#define EVENT_READER_H


#include <yaml-cpp/yaml.h>
#include <stdexcept>
#include <fstream>

#include "common/types.h"

namespace celib
{



    class EventReader
    {
        public:
            EventReader(YAML::Node& config);

            ~EventReader()
            {
                if(event_file_.is_open()) event_file_.close();
            }

            EventPtr getEvent();

        private:
            
            bool first_ = true;
            std::ifstream event_file_;
            Event first_event_;
            Mat3 cam_mat_;
            double k1_, k2_, k3_;
            double p1_, p2_;

    };






} // End namespace celib

#endif