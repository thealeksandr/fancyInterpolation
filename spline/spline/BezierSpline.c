//
//  BezierSpline.c
//  spline
//
//  Created by Aleksandr Nikiforov on 7/22/16.
//  Copyright © 2016 Aleksandr Nikiforov. All rights reserved.
//

#include "BezierSpline.h"

//https://www.particleincell.com/2012/bezier-splines/
//https://habrahabr.ru/post/247235/
void computeControlPoints(float *points, int n) {
    
    float *p1;
    float *p2;
    
    /*rhs vector*/
    float *a;
    float *b;
    float *c;
    float *r;
    //TODO: alloc arrays. What size?
    
    /*left most segment*/
    a[0] = 0;
    b[0] = 2;
    c[0] = 1;
    r[0] = points[0] + 2 * points[1];
    
    /*internal segments*/
    for (int i = 1; i < n - 1; i++) {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        r[i] = 4 * points[i] + 2 * points[i+1];
    }
    
    /*right segment*/
    a[n-1] = 2;
    b[n-1] = 7;
    c[n-1] = 0;
    r[n-1] = 8 * points[n-1]+points[n];
    
    /*solves Ax=b with the Thomas algorithm (from Wikipedia)*/
    for (int i = 1; i < n; i++) {
        int m = a[i] / b[i-1];
        b[i] = b[i] - m * c[i - 1];
        r[i] = r[i] - m * r[i-1];
    }
    
    p1[n-1] = r[n-1]/b[n-1];
    for (int i = n - 2; i >= 0; --i)
        p1[i] = (r[i] - c[i] * p1[i+1]) / b[i];
    
    /*we have p1, now compute p2*/
    for (int i = 0; i < n - 1; i++)
        p2[i] = 2 * points[i+1] - p1[i+1];
    
    p2[n-1]=0.5 * (points[n] + p1[n-1]);
    
    return {p1:p1, p2:p2};
}