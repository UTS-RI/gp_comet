/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2021 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef COMMON_INTERNAL_UTILS_H
#define COMMON_INTERNAL_UTILS_H

#include "common/math_utils.h"
#include "yaml-cpp/yaml.h"
#include <utility>
#include <numeric>

namespace celib
{
    

    // Function to sort indexes
    template <typename T>
    std::vector<int> sortIndexes(const std::vector<T> &v)
    {

        // initialize original index locations
        std::vector<int> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0); 

        // sort indexes based on comparing values in v
        std::sort(idx.begin(), idx.end(),
            [&v](int i1, int i2) {return v[i1] < v[i2];});

        return idx;
    }

    template <typename T>
    std::vector<T> organise(const std::vector<int>& indexes,const std::vector<T> &v)
    {
        if(indexes.size() != v.size()) throw std::range_error("Trying to orgainse a vector given a vector of indexes with a different size");
        std::vector<T> output(v.size());
        for(int i = 0; i < indexes.size(); ++i)
        {
            output[i] = v[indexes[i]];
        }
        return output;
    }

    template <typename T>
    inline void concatenateVector(std::vector<T>& v_base, const std::vector<T>& v_add)
    {
        v_base.insert( v_base.end(), v_add.begin(), v_add.end() );
    }



    template <class T>
    inline void printOption(bool user_defined, std::string field, T value){
        if(user_defined){
            std::cout << "[Config file] User defined value for " << field << " : " << value << std::endl;
        }else{
            std::cout << "[Config file] Default value for " << field << " : " << value << std::endl;
        }
    }

    inline void printErrorValueOption(std::string field){
        std::cout << "[Config file] ERROR: User defined value for " << field << " is not valid" << std::endl;
        throw std::invalid_argument("Invalid configuration value");
    }

    template <class T>
    inline T readField(YAML::Node& config, std::string field, T default_value){
        T output;
        bool user_defined;
        if(config[field]){
            output = config[field].as<T>();
            user_defined = true;
        }else{
            output = default_value;
            user_defined = false;
        }
        printOption<T>(user_defined, field, output);
        return output;
    }

    template <class T>
    inline T readRequiredField(YAML::Node& config, std::string field){
        T output;
        if(config[field]){
            output = config[field].as<T>();
        }else{
            std::cout << "[Config file] ERROR: It seems that the configuration file does not contain the required field " << field << std::endl; 
            throw std::invalid_argument("Information missing in congig file");
        }
        printOption<T>(true, field, output);
        return output;
    }


    inline std::vector<std::string> stringToStringVector(std::string s, char separator){
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream token_stream(s);
        while (std::getline(token_stream, token, separator)) 
        {
            tokens.push_back(token);
        }
        return tokens;
    }


    template <class T>
    inline std::vector<T> stringToNumVector(std::string s, char separator){
        std::vector<std::string> tokens = stringToStringVector(s,separator);
        std::vector<T> values;
        for(int i = 0; i<tokens.size(); ++i){
            if(tokens[i] != ""){
                values.push_back(std::stof(tokens[i]));
            }
        }
        return values;
    }

    inline std::tuple<int, double> stringToSec(std::string s)
    {
        std::vector<std::string> time_strings = stringToStringVector(s, '.');
        if(time_strings.size() != 2) throw std::invalid_argument("stringToSec: The string does not have the right format (\"sec.nsec\")");
        double nsec = std::stod("0."+time_strings.at(1));
        int sec = std::stoi(time_strings.at(0));
        return{sec, nsec};
    }


}


#endif
