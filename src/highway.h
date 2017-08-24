#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"
#include "planner.h"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

double euclidean_distance(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x1-x2, 2) + pow(y2-y1, 2));
}


double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

class Car {
 public:
  double x, y, s, d, yaw, speed;
  const double conv = 0.44703;

  Car(double x, double y, double s, double d, double yaw, double speed) {
    this->x = x;
    this->y = y;
    this->s = s;
    this->d = d;
    this->yaw = yaw;
    // mph to m/s
    this->speed = speed*conv;
  }


};

class Trajectory {
 public:

};

class JMT;

class Highway {

 public:

  vector<double> xs, ys, ss, dxs, dys;
  int size;

  Highway(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys) {
    this->xs = xs;
    this->ys = ys;
    this->ss = ss;
    this->dxs = dxs;
    this->dys = dys;
    this->size = xs.size();
  }

  vector<double> slice(vector<double> vec, int start, int c) {
    if (start + c < vec.size())
      return vector<double>(vec.begin() + start, vec.begin() + start + c);
    int extra = (start + c) - vec.size();
    vector<double> res = vector<double>(vec.begin() + start, vec.end());
    res.insert(res.end(), vec.begin(), vec.begin() + extra);
    return res;
  }

  int closest_waypoint(Car car) {
    vector<double> distances(this->size);
    transform(this->xs.begin(), this->xs.end(), this->ys.begin(), distances.begin(),
              [car](const double& x, const double& y) { return euclidean_distance(x, y, car.x, car.y); });
    return std::distance(distances.begin(), min_element(distances.begin(), distances.end()));
  }

  int next_waypoint(Car car) {
    int next = closest_waypoint(car);
    double wx = this->xs[next], wy = this->ys[next];
    double angle = atan2(wy-car.y, wx-car.x);
    if (fabs(angle-car.yaw) > M_PI/4) {
      next ++;
    }
    return next;
  }

  std::pair<vector<double>, vector<double>> next_traj(Car car) {
    vector<double> xpts, ypts;
    double inc = 0.5;
    for (int i = 0; i < 50; ++i) {
      double s = car.s + (i+1)*inc;
      auto xy = this->getXY(s, car.d);
      xpts.push_back(xy[0]);
      ypts.push_back(xy[1]);
    }
    return std::make_pair(xpts, ypts);
  }

  JMT next_path(Car car, double speed_limit);

  vector<double> getXY(double s, double d) {
    return ::getXY(s, d, this->ss, this->xs, this->ys);
  }

  void plot_highway(Car car, int section = 10000, int width = 50) {
    plt::plot({car.x}, {car.y}, "r*");
    this->plot_highway(section, width);
  }

  void plot_highway(int section = 10000, int width = 50) {
    if (section > dxs.size()) {
      plt::plot(xs, ys);
      section = dxs.size();
    }
    for (int i = 0; i < section; ++i) {
      plt::plot({xs[i], xs[i]+width*dxs[i]}, {ys[i], ys[i]+width*dys[i]}, "g");
      plt::plot({xs[i]+width*dxs[i]}, {ys[i]+width*dys[i]}, "g.");
    }
  }

};

class Poly {
 private:
  vector<double> coords;

 public:
  Poly(vector<double> coords) {
    this->coords = coords;
  }

  double value_at(double t) {
    double val = 0.0;
    for (int j=0; j<this->coords.size(); j++) val += this->coords[j]*pow(t,j);
    return val;
  }

  void plot(double start, double end, double inc) {
    double t = start;
    vector<double> c = this->coords;
    vector<double> p1, p2, p3;
    while(t<end){
      double t2=t*t, t3=t2*t, t4=t3*t, t5=t4*t;
      double f = c[0] + c[1]*t + c[2]*t2 + c[3]*t3 + c[4]*t4 + c[5]*t5;
      double f_dot = c[1] + 2*c[2]*t + 3*c[3]*t2 + 4*c[4]*t3 + 5*c[5]*t4;
      double f_ddot = 2*c[2] + 6*c[3]*t + 12*c[4]*t2 + 20*c[5]*t3;
      p1.push_back(f);
      p2.push_back(f_dot);
      p3.push_back(f_ddot);
      t += inc;
    }
    plt::plot(p1, "r");
    plt::plot(p2, "g");
    plt::plot(p3, "b");
  }
};

class JMT {
 private:
  vector<double> s_vec, d_vec;

 public:
  JMT(vector<double> s_vec, vector<double> d_vec) {
    this->s_vec = s_vec;
    this->d_vec = d_vec;
  }

  std::pair<double, double> length() {
    double slen = this->s_vec[3] - this->s_vec[0];
    double dlen = this->d_vec[3] - this->d_vec[0];
    return std::make_pair(slen, dlen);
  }

  std::pair<Poly, Poly> interpolate(double tf);

  vector<vector<double>> generate_points(double tf, Highway& highway);
};


std::pair<Poly, Poly> JMT::interpolate(double tf) {
  double t = tf, t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t;
  MatrixXd Tmat(3, 3);
  Tmat << t3, t4, t5,
          3*t2, 4*t3, 5*t4,
          6*t, 12*t2, 20*t3;
  MatrixXd Tinv = Tmat.inverse();
  VectorXd Smat(3);
  Smat << this->s_vec[3],
          this->s_vec[4],
          this->s_vec[5];
  VectorXd Dmat(3);
  Dmat << this->d_vec[3],
          this->d_vec[4],
          this->d_vec[5];
  VectorXd Svec = Tinv*Smat;
  VectorXd Dvec = Tinv*Dmat;
  vector<double> s_poly = { this->s_vec[0], this->s_vec[1], this->s_vec[2]/2.0, Svec(0), Svec(1), Svec(2) };
  vector<double> d_poly = { this->d_vec[0], this->d_vec[1], this->d_vec[2]/2.0, Dvec(0), Dvec(1), Dvec(2) };
  cout << "spoly - ";
  for (int i = 0; i < s_poly.size(); ++i) {
    cout << s_poly[i] << ", ";
  }
  cout << endl;
  return std::make_pair(Poly(s_poly), Poly(d_poly)) ;
}

vector<vector<double>> JMT::generate_points(double tf, Highway& highway) {
  auto polys = interpolate(tf);
  vector<double> xpts, ypts, ss, ds, ts;
  double inc = (tf/50.0);
  cout << "inc - "<< inc << endl;
  for (int i=0; i<50; i++) {
    double t=inc*(i+1), s=0, d=0;
    s += polys.first.value_at(t);
    d += polys.second.value_at(t);
    ss.push_back(s);
    ds.push_back(d);
    ts.push_back(t);
    auto xy = highway.getXY(s, d);
    xpts.push_back(xy[0]);
    ypts.push_back(xy[1]);
  }
  return { xpts, ypts, ss, ds, ts};
}

JMT Highway::next_path(Car car, double max_speed) {
  int wp = next_waypoint(car);
  double s_delta = this->ss[wp] - car.s;
  vector<double> s_vec = {car.s, car.speed, 0, this->ss[wp], 8, 8.0};
  vector<double> d_vec = {car.d, 0, 0, car.d, 0, 0};
  return JMT(s_vec, d_vec);
}
