/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2021 Cedric LE GENTIL
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include "Eigen/Dense"
#include "Eigen/Sparse"





namespace celib{


    typedef Eigen::Matrix<double, 2, 2> Mat2;
    typedef Eigen::Matrix<double, 3, 3> Mat3;
    typedef Eigen::Matrix<double, 3, 1> Vec3;
    typedef Eigen::Matrix<double, 4, 1> Vec4;
    typedef Eigen::Matrix<double, 1, 3> Row3;
    typedef Eigen::Matrix<double, 1, 4> Row4;
    typedef Eigen::Matrix<double, 1, 2> Row2;
    typedef Eigen::Matrix<double, 6, 6> Mat6;
    typedef Eigen::Matrix<double, 4, 4> Mat4;
    typedef Eigen::Matrix<double, 9, 9> Mat9;
    typedef Eigen::Matrix<double, 12, 12> Mat12;
    typedef Eigen::Matrix<double, 2, 1> Vec2;
    typedef Eigen::Matrix<double, 6, 1> Vec6;
    typedef Eigen::Matrix<double, 8, 1> Vec8;
    typedef Eigen::Matrix<double, 1, 8> Row8;
    typedef Eigen::Matrix<double, 9, 1> Vec9;
    typedef Eigen::Matrix<double, 1, 9> Row9;
    typedef Eigen::Matrix<double, 1, 12> Row12;
    typedef Eigen::Matrix<double, 3, 6> Mat3_6;
    typedef Eigen::Matrix<double, 3, 2> Mat3_2;
    typedef Eigen::Matrix<double, 3, 4> Mat3_4;
    typedef Eigen::Matrix<double, 2, 6> Mat2_6;
    typedef Eigen::Matrix<double, 9, 6> Mat9_6;
    typedef Eigen::Matrix<double, 3, 9> Mat3_9;
    typedef Eigen::Matrix<double, 2, 8> Mat2_8;
    typedef Eigen::Matrix<double, 2, 9> Mat2_9;
    typedef Eigen::Matrix<double, 2, 3> Mat2_3;
    typedef Eigen::Matrix<double, 9, 3> Mat9_3;
    typedef Eigen::Matrix<double, 9, 8> Mat9_8;
    typedef Eigen::Matrix<double, 3, 12> Mat3_12;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecX;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowX;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatX;
    typedef Eigen::SparseMatrix<double> sMatX;
    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Mat3X;
    typedef Eigen::Matrix<double, 2, Eigen::Dynamic> Mat2X;


    // Data structure to store individual event measurements
    struct Event
    {
        double t;
        double x;
        double y;
        bool pol;

        Vec2 pixVec()
        {
            Vec2 output(x, y);
            return output;
        }
        Vec3 pixVec3()
        {
            Vec3 output(x, y, 1);
            return output;
        }
        Vec3 pt3D()
        {
            Vec3 output(x, y, t);
            return output;
        }

        void print()
        {
            std::cout << "Event a t = " << t << "    : x = " << x << "    y : " << y << "    pol: " << pol << std::endl;
        }

        Event(){}

        Event(Vec2 xy)
        {
            x = xy(0);
            y = xy(1);
            t = 0;
            pol = false;
        }
        Event(Vec2 xy, double _t)
        {
            x = xy(0);
            y = xy(1);
            t = _t;
            pol = false;
        }

        Event(Vec3 xyt)
        {
            x = xyt(0);
            y = xyt(1);
            t = xyt(2);
            pol = false;
        }
        Event(double _x, double _y, double _t)
        {
            x = _x;
            y = _y;
            t = _t;
            pol = false;
        }
    };
    typedef std::shared_ptr<Event> EventPtr;





    inline Vec3 stdToCelibVec3(const std::vector<double>& v)
    {
        Vec3 output;
        output << v[0], v[1], v[2];
        return output;
    }


    inline std::array<double, 3> vec3ToArray(const Vec3& vec)
    {
        return {vec(0), vec(1), vec(2)};
    }
    inline Vec3 arrayToVec3(const std::array<double,3>& array)
    {
        return Vec3{array[0], array[1], array[2]};
    }




    // Buffer structure from WOLF
    template <typename T>
    class Buffer
    {
    public:

        typedef typename std::map<double,T>::iterator Iterator; // buffer iterator

        Buffer(){};
        ~Buffer(void){};

        /**\brief Select an element from the buffer
        *
        *  Select from the buffer the closest element (w.r.t. time stamp),
        * respecting a defined time tolerances
        */
        T select(const double& _time_stamp, const double& _time_tolerance);
        
        /**\brief Select an element iterator from the buffer
        *
        *  Select from the buffer the iterator pointing to the closest element (w.r.t. time stamp),
        * respecting a defined time tolerances
        */
        Iterator selectIterator(const double& _time_stamp, const double& _time_tolerance);

        T selectLastBefore(const double& _time_stamp);
        T selectLastBeforeOrFirst(const double& _time_stamp);

        T selectFirstBefore(const double& _time_stamp, const double& _time_tolerance);
        
        T selectLastAfter(const double& _time_stamp, const double& _time_tolerance);

        T selectFirst();

        T selectLast();

        /**\brief Buffer size
        *
        */
        int size(void);

        /**\brief Add a element to the buffer
        *
        */
        void emplace(const double& _time_stamp, const T& _element);

        /** \brief returns the container with elements of the buffer
        *
        * elements are ordered from most recent to oldest
        */
        std::map<double,T>& getContainer();

        /**\brief Remove all elements in the buffer with a time stamp older than the specified
        *
        */
        void removeUpTo(const double& _time_stamp);

        /**\brief Remove all elements in the buffer with a time stamp older than the specified
        *
        */
        void removeUpToLower(const double& _time_stamp);

        /**\brief Clear the buffer
        *
        */
        void clear();

        /**\brief is the buffer empty ?
        *
        */
        bool empty();

    protected:
        /**\brief Check time tolerance
        *
        * Check if the time distance between two time stamps is smaller than
        * the time tolerance.
        */
        static bool checkTimeTolerance(const double& _time_stamp1, const double& _time_stamp2, const double& _time_tolerance);

        /**\brief Check time tolerance
        *
        * Check if the time distance between two time stamps is smaller than
        * the minimum time tolerance of the two frames.
        */
        static bool doubleCheckTimeTolerance(const double& _time_stamp1, const double& _time_tolerance1, const double& _time_stamp2, const double& _time_tolerance2);

    protected:

        std::map<double,T> container_; // Main buffer container
    };



    template <typename T>
    typename Buffer<T>::Iterator Buffer<T>::selectIterator(const double& _time_stamp, const double& _time_tolerance)
    {
        Buffer<T>::Iterator post = container_.upper_bound(_time_stamp);

        bool prev_exists = (post != container_.begin());
        bool post_exists = (post != container_.end());

        bool post_ok = post_exists && checkTimeTolerance(post->first, _time_stamp, _time_tolerance);

        if (prev_exists)
        {
            Buffer<T>::Iterator prev = std::prev(post);

            bool prev_ok = checkTimeTolerance(prev->first, _time_stamp, _time_tolerance);

            if (prev_ok && !post_ok)
                return prev;

            else if (!prev_ok && post_ok)
                return post;

            else if (prev_ok && post_ok)
            {
                if (std::fabs(post->first - _time_stamp) < std::fabs(prev->first - _time_stamp))
                    return post;
                else
                    return prev;
            }
        }
        else if (post_ok)
            return post;

        return container_.end();
    }

    template <typename T>
    T Buffer<T>::select(const double& _time_stamp, const double& _time_tolerance)
    {
        if (container_.empty())
            return nullptr;

        Buffer<T>::Iterator it = selectIterator(_time_stamp, _time_tolerance);

        // end is returned from selectIterator if an element of the buffer complying with the time stamp
        // and time tolerance has not been found
        if (it != container_.end()){
            return it->second;
        }
        
        return nullptr;
    }

    template <typename T>
    T Buffer<T>::selectLastBefore(const double& _time_stamp)
    {
        // There is no element
        if (container_.empty())
            throw std::range_error("Buffer::selectLastBefore: Querying empty buffer");

        // Checking on begin() since elements are ordered in time
        // Return first element if is older than time stamp
        if (container_.begin()->first > _time_stamp)
            throw std::range_error("Buffer::selectLastBefore: Querying before first entry");

        return container_.lower_bound(_time_stamp)->second;
    }

    template <typename T>
    T Buffer<T>::selectLastBeforeOrFirst(const double& _time_stamp)
    {
        // There is no element
        if (container_.empty())
            throw std::range_error("Buffer::selectLastBefore: Querying empty buffer");

        // Checking on begin() since elements are ordered in time
        // Return first element if is older than time stamp
        if (container_.begin()->first > _time_stamp)
            container_.begin()->second;

        return container_.lower_bound(_time_stamp)->second;
    }


    template <typename T>
    T Buffer<T>::selectFirstBefore(const double& _time_stamp, const double& _time_tolerance)
    {
        // There is no element
        if (container_.empty())
            return nullptr;

        // Checking on begin() since elements are ordered in time
        // Return first element if is older than time stamp
        if (container_.begin()->first < _time_stamp)
            return container_.begin()->second;

        // Return first element if despite being newer, it is within the time tolerance
        if (checkTimeTolerance(container_.begin()->first, _time_stamp, _time_tolerance))
            return container_.begin()->second;

        // otherwise return nullptr (no element before the provided ts or within the tolerance was found)
        return nullptr;
    }


    template <typename T>
    T Buffer<T>::selectLastAfter(const double& _time_stamp, const double& _time_tolerance)
    {
        // There is no element
        if (container_.empty())
            return nullptr;

        // Checking on rbegin() since elements are ordered in time
        // Return last element if is newer than time stamp
        if (container_.rbegin()->first > _time_stamp)
            return container_.rbegin()->second;

        // Return last element if despite being older, it is within the time tolerance
        if (checkTimeTolerance(container_.rbegin()->first, _time_stamp, _time_tolerance))
            return container_.rbegin()->second;

        // otherwise return nullptr (no element after the provided ts or within the tolerance was found)
        return nullptr;
    }

    template <typename T>
    T Buffer<T>::selectFirst()
    {
        // There is no element
        if (container_.empty())
            return nullptr;

        // Returning first map element
        return container_.begin()->second;
    }

    template <typename T>
    T Buffer<T>::selectLast()
    {
        // There is no element
        if (container_.empty())
            return nullptr;

        // Returning last map element
        return container_.rbegin()->second;
    }

    template <typename T>
    void Buffer<T>::emplace(const double& _time_stamp, const T& _element)
    {
        container_.emplace(_time_stamp, _element);
    }

    template <typename T>
    std::map<double,T>& Buffer<T>::getContainer()
    {
        return container_;
    }

    template <typename T>
    inline void Buffer<T>::clear()
    {
        container_.clear();
    }

    template <typename T>
    inline bool Buffer<T>::empty()
    {
        return container_.empty();
    }

    template <typename T>
    inline int Buffer<T>::size(void)
    {
        return container_.size();
    }

    template <typename T>
    inline void Buffer<T>::removeUpTo(const double& _time_stamp)
    {
        Buffer::Iterator post = container_.upper_bound(_time_stamp);
        container_.erase(container_.begin(), post); // erasing by range
    }

    template <typename T>
    inline void Buffer<T>::removeUpToLower(const double& _time_stamp)
    {
        Buffer::Iterator post = container_.lower_bound(_time_stamp);
        container_.erase(container_.begin(), post); // erasing by range
    }

    template <typename T>
    inline bool Buffer<T>::doubleCheckTimeTolerance(const double& _time_stamp1,
                                                    const double& _time_tolerance1,
                                                    const double& _time_stamp2,
                                                    const double& _time_tolerance2)
    {
        double time_diff = std::fabs(_time_stamp1 - _time_stamp2);
        double time_tol  = std::min(_time_tolerance1, _time_tolerance2);
        bool pass = time_diff <= time_tol;
        return pass;
    }

    template <typename T>
    inline bool Buffer<T>::checkTimeTolerance(const double& _time_stamp1,
                                                    const double& _time_stamp2,
                                                    const double& _time_tolerance)
    {
        double time_diff = std::fabs(_time_stamp1 - _time_stamp2);
        bool pass = time_diff <= _time_tolerance;
        return pass;
    }


    struct HyperParam
    {
        double sf2;
        double l2;
        double sz2;
    };

}

#endif
