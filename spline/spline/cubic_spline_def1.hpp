//
//  cubic_spline_def1.hpp
//  spline
//
//  Created by Александр Никифоров on 7/18/16.
//  Copyright © 2016 thealeksandr. All rights reserved.
//

#ifndef cubic_spline_def1_hpp
#define cubic_spline_def1_hpp

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <limits>

class cubic_spline {
    
private:
    // Cupic spline structure.
    struct spline_tuple
    {
        double a, b, c, d, x;
    };
    
    spline_tuple *splines; // Spline.
    std::size_t n; // Number of nodes.
    
    void free_mem(); // Free memory.
    
public:
    cubic_spline(); //Constructor.
    ~cubic_spline(); //Deconstructor.
    
    // Build spline
    // x - nodes order by ascending
    // y - function values in nodes
    // n - number of nodes.
    void build_spline(const double *x, const double *y, std::size_t n);
    
    // Calculate spline value in x
    double f(double x) const;
};

#endif /* cubic_spline_def1_hpp */
