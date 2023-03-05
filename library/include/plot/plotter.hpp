/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2017 Cedric LE GENTIL
 *
 *  This code is the declaration of a simple wrapper for gnuplot-iostream.
 *  This code is initially designed to get hands-on interaction with gnuplot
 *  and gnuplot-iostream.
 *
 *  To work, this wrapper needs gnuplot-iostream (header-only library) and
 *  gnuplot to be installed on your computer
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/

#ifndef PLOTTER_H
#define PLOTTER_H

#include "gnuplot-iostream/gnuplot-iostream.h"
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

namespace celib {

// Default values for the figure titles
const std::string kDefaultTitle = "Title";
const std::string kDefaultXTitle = "X Title";
const std::string kDefaultYTitle = "Y Title";

// Color choices
enum class Color{
    RED, GREEN, BLUE, ORANGE, CYAN, PURPLE, DARK_RED, DARK_GREEN, DARK_BLUE, BLACK, GREY, DEFAULT
};

// Point style choices (do not change the numerical values, directly maps to gnuplot parameters)
enum class PointStyle{
    NONE = 0, PLUS = 1, CROSS = 2, SQUARE = 4, SQUARE_FULL = 5, CIRCLE = 6, CIRCLE_FULL = 7
};

// Line style choices (do not change the numerical values, directly maps to gnuplot parameters)
enum class LineStyle{
    DOTS = 0, FULL = 1, NONE
};

// Class Plotter
// Provides a "high-level" API to interact with gnuplot through gnuplot-iostream
// An instance of that class correspond to a window/figure.
// Curves are added to that window either from the constructor or from the method addPlot(...).
// Once all the desired curves are added to the Plotter instance, call the methode plot() to
// actually display the curves. Calling plot() empties the data buffers. Therfore to plot new
// data in the same windows, you need to addPlot(...) and use plot() again (that will erase the
// previously plotted data).
class Plotter{

    public :
        // Constructor setting the figure title and the axis titles
        Plotter(std::string title = kDefaultTitle,
                std::string x_title = kDefaultXTitle,
                std::string y_title = kDefaultYTitle);

        // Constructor setting the figure title, the axis titles,
        // the first set of data and the associated style (optional).
        Plotter(std::string title,
                std::string x_title,
                std::string y_title,
                std::vector<double> x_data,
                std::vector<double> y_data,
                std::string legend,
                Color color = Color::DEFAULT,
                PointStyle point_style = PointStyle::NONE,
                LineStyle line_style = LineStyle::FULL,
                int point_size = 1,
                int line_width = 1);

        // Setter for the figure title
        void setTitle(std::string title);
       
        // Setter for the figure title and axis titles
        void setTitles(std::string title, std::string x_title, std::string y_title);

        // Add a new curve to be plotted
        // Arguments provide the data to plot, the legend and the optional style of the plot
        void addPlot(
                std::vector<double> x_data,
                std::vector<double> y_data,
                std::string legend = "legend",
                Color color = Color::DEFAULT,
                PointStyle point_style = PointStyle::NONE,
                LineStyle line_style = LineStyle::FULL,
                int point_size = 1,
                int line_width = 1);

        // Add a new curve to be plotted
        // Arguments provide the data to plot, the legend and the optional style of the plot
        void addPlot(
                Eigen::VectorXd x_data,
                Eigen::VectorXd y_data,
                std::string legend = "legend",
                Color color = Color::DEFAULT,
                PointStyle point_style = PointStyle::NONE,
                LineStyle line_style = LineStyle::FULL,
                int point_size = 1,
                int line_width = 1);

        // Plot the curves stored and clean the storage
        void plot();

    private :
        // Gnuplot-iostream interface to send commands to gnuplot
        gnuplotio::Gnuplot gplot_;
      
        // Storage of the figure titles
        std::string title_;
        std::string x_title_;
        std::string y_title_;
        
        // Static map to the Color enum values to the gnuplot color strings
        static const std::map<Color, std::string> kColorMap;
        
        // Data structure of one buffered curve
        struct Plot{
            std::pair<std::vector<double>, std::vector<double> > data;
            std::string legend;
            std::string style;
        };

        // Buffer of curves
        std::vector< Plot > plots_;
};



} // namespace

#endif // PLOTTER_H
