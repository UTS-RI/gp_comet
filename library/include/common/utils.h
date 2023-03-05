/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2021 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <ctime>
#include <chrono>

namespace celib
{
    class StopWatch
    {
        public:
            StopWatch()
            {
                duration_ = std::chrono::high_resolution_clock::duration::zero();
            }

            //Start the counting (also used to get time while running);
            double start()
            {
                double output = 0.0;
                auto temp = std::chrono::high_resolution_clock::now();
                if(!stopped_)
                {
                    duration_ += (temp - last_);
                    output = std::chrono::duration_cast<std::chrono::microseconds>(duration_).count();
                }
                last_ = temp;
                stopped_ = false;
                return output/1000.0;
            }

            // Since the first "start after the last reset"
            double getTotal()
            {
                auto temp = std::chrono::high_resolution_clock::now();
                return std::chrono::duration_cast<std::chrono::microseconds>(duration_ + (temp - last_)).count()/1000.0;
            }

            // Since the 
            double getLast()
            {
                auto temp = std::chrono::high_resolution_clock::now();
                return std::chrono::duration_cast<std::chrono::microseconds>(temp - last_).count()/1000.0;

            }

            double stop()
            {
                if(!stopped_)
                {
                    auto temp = std::chrono::high_resolution_clock::now();
                    stopped_ = true;
                    duration_ += (temp - last_);
                    return std::chrono::duration_cast<std::chrono::microseconds>(duration_).count()/1000.0;
                }
                else
                {
                    std::cout << "WARNING: Stopping a StopWatch that is not running, returning negative time" << std::endl;
                    return -1;
                }
            }

            void reset()
            {
                stopped_ = true;
                duration_ = std::chrono::high_resolution_clock::duration::zero();
            }


            void print()
            {
                print("Time elapsed:");
            }

            void print(std::string str)
            {
                std::cout << str << " " << std::chrono::duration_cast<std::chrono::microseconds>(duration_).count()/1000.0 << " ms" << std::endl;
            }

        private:
            double counter_ = 0;
            bool stopped_ = true;
            std::chrono::high_resolution_clock::duration duration_;
            std::chrono::high_resolution_clock::time_point last_;
    };


    inline void pause()
    {
        do 
        {
            std::cout << '\n' << "Press a key to continue...";
        } while (std::cin.get() != '\n');
    }


    // Gerenate random indexes between start and end (excluded) without repetition
    inline std::vector<int> generateRandomIndexes(int start, int end, int nb)
    {
        std::vector<int> output;
        if(nb > end - start)
        {
            std::cout << "ERROR: generateRandomIndexes: nb > end - start" << std::endl;
            return output;
        }
        std::vector<int> indexes;
        for(int i = start; i < end; i++)
        {
            indexes.push_back(i);
        }
        std::random_shuffle(indexes.begin(), indexes.end());
        for(int i = 0; i < nb; i++)
        {
            output.push_back(indexes[i]);
        }
        return output;
    }

}



#endif