//
//  norm_rad.h
//  UKFProject
//
//  Created by 刘威 on 2017/3/23.
//  Copyright © 2017年 wsix. All rights reserved.
//

#ifndef norm_rad_h
#define norm_rad_h

#include <cmath>

/*
 * Bring the 'difference' between two angles into [-pi; pi].
 * See http://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi
 */
template <typename T>
T normalizeRadiansPiToMinusPi(T rad) {
  // Copy the sign of the value in radians to the value of pi.
  T signed_pi = std::copysign(M_PI,rad);
  // Set the value of difference to the appropriate signed value between pi and -pi.
  rad = std::fmod(rad + signed_pi,(2 * M_PI)) - signed_pi;
  return rad;
}

#endif /* norm_rad_h */
