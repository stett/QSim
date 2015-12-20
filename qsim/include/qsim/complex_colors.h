#ifndef COMPLEX_COLORS_H
#define COMPLEX_COLORS_H

#include <SFML/Graphics.hpp>
#include <gsl/gsl_complex.h>
//#include <cmath>

namespace qsim {

    sf::Color color_complex(double imag) {
        double arg = imag + M_PI;
        double hue = arg / (2.0 * M_PI);
        double r, g, b;
        double one_sixth = 1.0 / 6.0;
        
        if (hue > 5.0 * one_sixth) {
            r = 1; g = 0; b = 6.0 * (1.0 - hue);
        } else if (hue > 4.0 * one_sixth) {
            g = 0; b = 1; r = 6.0 * (hue - 4.0 * one_sixth);
        } else if (hue > 3.0 * one_sixth) {
            b = 1; r = 0; g = 6.0 * (4.0 * one_sixth - hue);
        } else if (hue > 2.0 * one_sixth) {
            r = 0; g = 1; b = 6.0 * (hue - 2.0 * one_sixth);
        } else if (hue > one_sixth) {
            g = 1; b = 0; r = 6.0 * (2.0 * one_sixth - hue);
        } else {
            b = 0; r = 1; g = 6.0 * hue;
        }

        return sf::Color(
            (unsigned int)(255.0 * r),
            (unsigned int)(255.0 * g),
            (unsigned int)(255.0 * b),
            255);
    }

    sf::Color color_complex(const gsl_complex& z) {
        return color_complex(gsl_complex_arg(z));
    }

    sf::Color color_complex_circle(const gsl_complex& z) {

        // Define color base vectors (splitting unit circle into 3rds)
        double rx = 0;
        double ry = 1;
        double gx = cos(M_PI * 7.0 / 6.0);
        double gy = sin(M_PI * 7.0 / 6.0);
        double bx = cos(M_PI * 11.0 / 6.0);
        double by = sin(M_PI * 11.0 / 6.0);

        // Get position on unit circle from complex argument
        double arg = gsl_complex_arg(z);
        double vx = cos(arg);
        double vy = sin(arg);

        // Project position onto color base vectors
        double r = std::max(vx * rx + vy * ry, 0.0);
        double g = std::max(vx * gx + vy * gy, 0.0);
        double b = std::max(vx * bx + vy * by, 0.0);

        // 
        return sf::Color(
            (unsigned int)(255.0 * r),
            (unsigned int)(255.0 * g),
            (unsigned int)(255.0 * b),
            255);
    }
};

#endif