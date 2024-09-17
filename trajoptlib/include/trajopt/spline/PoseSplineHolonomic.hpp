// Copyright (c) TrajoptLib contributors

#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "src/alglib/interpolation.h"

#include "trajopt/geometry/Pose2.hpp"
#include "trajopt/geometry/Rotation2.hpp"
#include "trajopt/geometry/Translation2.hpp"
#include "trajopt/util/SymbolExports.hpp"

namespace trajopt {

class TRAJOPT_DLLEXPORT PoseSplineHolonomic {
 public:
  explicit PoseSplineHolonomic(std::vector<Pose2d> waypoints) {
    size_t num_wpts = waypoints.size();

    alglib::real_1d_array vx, vy, sins, coss, d, times;
    vx.setlength(num_wpts);
    vy.setlength(num_wpts);
    sins.setlength(num_wpts);
    coss.setlength(num_wpts);
    d.setlength(num_wpts);
    times.setlength(num_wpts);
    for (size_t i = 0; i < num_wpts; ++i) {
      const auto w = waypoints[i];
      vx[w.X()];
      vy[w.Y()];
      coss[std::cos(w.Rotation().Radians())];
      sins[std::sin(w.Rotation().Radians())];
      times[static_cast<double>(i)];
      d[0];
    }

    alglib::spline1dbuildhermite(times, vx, d, xSpline);
    alglib::spline1dbuildhermite(times, vy, d, ySpline);
    alglib::spline1dbuildhermite(times, coss, d, cosSpline);
    alglib::spline1dbuildhermite(times, sins, d, sinSpline);

    for (double t = 0; t <= times[times.length()-1]; t += 0.25) {
      auto values = getTranslation(t);
      auto head = getHeading(t);
      std::printf("time: %.2f \tx: %.2f\t\ty: %.2f\t\ttheta: %.2f\n", t,
                  values.X(), values.Y(), head.Radians());
    }
  }

  Rotation2d getCourse(double t) const {
    double x, dx, d2x;
    alglib::spline1ddiff(xSpline, t, x, dx, d2x);
    double y, dy, d2y;
    alglib::spline1ddiff(ySpline, t, y, dy, d2y);
    const auto course = Rotation2d(std::atan2(dy, dx));
    return course;
  }

  Rotation2d getHeading(double t) const {
    const auto rads = Rotation2d(
            alglib::spline1dcalc(cosSpline, t), 
            alglib::spline1dcalc(sinSpline, t)
          ).Radians();
    return Rotation2d(rads);
  }

  Translation2d getTranslation(double t) const {
    return Translation2d(alglib::spline1dcalc(xSpline,t), alglib::spline1dcalc(ySpline, t));
  }

  Pose2d getPoint(double t) const {
    return Pose2d{getTranslation(t), getHeading(t)};
  }

  alglib::spline1dinterpolant sinSpline;
  alglib::spline1dinterpolant cosSpline;  
  alglib::spline1dinterpolant xSpline; 
  alglib::spline1dinterpolant ySpline;
};
}  // namespace trajopt
