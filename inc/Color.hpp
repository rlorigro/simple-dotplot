#ifndef SIMPLE_DOTPLOT_COLOR_HPP
#define SIMPLE_DOTPLOT_COLOR_HPP


#include <vector>
#include <array>
#include <string>

using std::vector;
using std::array;
using std::string;


class ColorMap {
    /// Methods ///
    virtual array<double,3> get_rgb(double x)=0;
};


class MplRainbow: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    MplRainbow()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class MplGnuplot: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    MplGnuplot()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class Viridis: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    Viridis()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};






#endif //SIMPLE_DOTPLOT_COLOR_HPP
