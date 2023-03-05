/**
 *  Author: Cedric LE GENTIL 
 *
 *  Copyright 2017 Cedric LE GENTIL
 *
 *  This code is the implementation of a simple wrapper for gnuplot-iostream.
 *  This code is initially designed to get hands-on interaction with gnuplot
 *  and gnuplot-iostream.
 *
 *  To work, this wrapper needs gnuplot-iostream (header-only library) and
 *  gnuplot to be installed on your computer
 *
 *  For any further question, recommendation or contribution
 *  le.gentil.cedric@gmail.com
 **/
#include <iostream>

#include "plot/plotter.hpp"

namespace celib {

    // Map the Color enum values to the gnuplot color strings
    const std::map<Color, std::string> Plotter::kColorMap = { 
        {Color::RED, "red"},
        {Color::GREEN, "green"},
        {Color::BLUE, "blue"},
        {Color::ORANGE, "orange"},
        {Color::CYAN, "cyan"},
        {Color::PURPLE, "violet"},
        {Color::DARK_RED, "dark-red"},
        {Color::DARK_GREEN, "dark-green"},
        {Color::DARK_BLUE, "dark-blue"},
        {Color::BLACK, "black"},
        {Color::GREY, "grey"}
        };

    // Constructor setting the figure title and the axis titles
    Plotter::Plotter(std::string title, std::string x_title, std::string y_title){
        setTitles(title, x_title, y_title);
    }

    // Constructor setting the figure title, the axis titles,
    // the first set of data and the associated style.
    Plotter::Plotter(
            std::string title,
            std::string x_title,
            std::string y_title,
            std::vector<double> x_data,
            std::vector<double> y_data,
            std::string legend,
            Color color,
            PointStyle point_style,
            LineStyle line_style,
            int point_size,
            int line_width){

        setTitles(title, x_title, y_title);
        addPlot(x_data, y_data, legend, color, point_style, line_style, point_size, line_width);
    }



    // Setter for the figure title
    void Plotter::setTitle(std::string title){
        title_ = title;
    }

    // Setter for the figure title and axis titles
    void Plotter::setTitles(std::string title, std::string x_title, std::string y_title){
        title_ = title;
        x_title_ = x_title;
        y_title_ = y_title;
    }

    // Add a new curve to be plotted
    // Arguments provide the data to plot, the legend and the style of the plot
    void Plotter::addPlot(
            std::vector<double> x_data,
            std::vector<double> y_data,
            std::string legend,
            Color color,
            PointStyle point_style,
            LineStyle line_style,
            int point_size,
            int line_width){

        // Print an error message if the x and y data don't have the same length
        if(x_data.size() != y_data.size()){
            std::cerr << "ADDING plot to the plotter, x and y data vector "
                << "don't have the same lenght" << std::endl;
        }else{
            // Populate the vector of plots
            Plot temp_plot;
            temp_plot.data = std::make_pair(x_data, y_data);
            temp_plot.legend = legend;

            // Generate the style string corresping to the arguments
            temp_plot.style = "";
            if(line_style != LineStyle::NONE){
                temp_plot.style += "linespoints lt '";
                temp_plot.style += std::to_string(static_cast<int>(line_style));
                temp_plot.style += "' lw ";
                temp_plot.style += std::to_string(line_width);
            }else{
                temp_plot.style += "points";
            }
            temp_plot.style += " pt ";
            temp_plot.style += std::to_string(static_cast<int>(point_style));
            temp_plot.style += " ps ";
            temp_plot.style += std::to_string(point_size);
            temp_plot.style += " lc ";
            if(color != Color::DEFAULT){
                temp_plot.style += "rgb '";
                temp_plot.style += kColorMap.at(color);
                temp_plot.style += "'";
            }else{
                temp_plot.style += std::to_string(plots_.size() + 1);
            }


            plots_.push_back(temp_plot);
        }
    }



    void Plotter::addPlot(
            Eigen::VectorXd x_data,
            Eigen::VectorXd y_data,
            std::string legend,
            Color color,
            PointStyle point_style,
            LineStyle line_style,
            int point_size,
            int line_width){

        if(x_data.size() != y_data.size())
        {
            throw std::range_error("Plotter: The X and Y axis data do not have the same lenght");
        }
        int nb_data = x_data.size();
        std::vector<double> x;
        std::vector<double> y;
        x.reserve(nb_data);
        y.reserve(nb_data);
        for(int i = 0; i < nb_data; ++i)
        {
            x.push_back(x_data[i]);
            y.push_back(y_data[i]);
        }
        addPlot(x,y,legend,color,point_style,line_style,point_size,line_width);
    }

    // Plot the curves stored and clean the storage
    void Plotter::plot(){
        // If there is data to plot send the gnuplot commands
        if(plots_.size() > 0){

            // Send gnuplot commands
            gplot_<< "set title '" << title_ <<"'" << std::endl;
            gplot_<< "set xlabel '" << x_title_ <<"'" << std::endl;
            gplot_<< "set ylabel '" << y_title_ <<"'" << std::endl;
            gplot_<< "plot '-' with " << plots_[0].style << " title '"
                << plots_[0].legend << "' ";
            for(int i = 1; i < plots_.size(); i++){
                gplot_<< ",'-' with " << plots_[i].style << " title '" << plots_[i].legend << "' ";
            }
            gplot_<< std::endl;

            gplot_.send1d(plots_[0].data);
            for(int i = 1; i < plots_.size(); i++){
                gplot_.send1d(plots_[i].data);
            }
            
        }else{
            std::cerr << "PLOTTING, No data to plot, or number of legends different "
                << "from the number of curves" << std::endl;
        }



    }

} // namespace
