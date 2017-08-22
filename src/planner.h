#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"


#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

double euclidean_distance(double x1, double y1, double x2, double y2) {
  return sqrt(pow(x1-x2, 2) + pow(y2-y1, 2));
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-M_PI/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

class Car {
 public:
  double x, y, s, d, yaw, speed;

  Car(double x, double y, double s, double d, double yaw, double speed) {
    this->x = x;
    this->y = y;
    this->s = s;
    this->d = d;
    this->yaw = yaw;
    this->speed = speed;
  }
};

class Path;

class Trajectory {

 private:
  vector<double> xs, ys, ss, dxs, dys;
  int size;

 public:

  Trajectory(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys);

  vector<double> slice(vector<double> vec, int start, int c);

  int closest_waypoint(Car car);

  int next_waypoint(Car car);

  Path next_path(Car car, double speed_limit);

  vector<double> getXY(double s, double d);

  void plot_highway(Car car, int section, int width);

  void plot_highway(int section, int width);

};

class Path {
 private:
  vector<double> s_vec, d_vec;

 public:
  Path(vector<double> s_vec, vector<double> d_vec) {
    this->s_vec = s_vec;
    this->d_vec = d_vec;
  }

  std::pair<double, double> length() {
    double slen = this->s_vec[3] - this->s_vec[0];
    double dlen = this->d_vec[3] - this->d_vec[0];
    return std::make_pair(slen, dlen);
  }

  vector<vector<double>> interpolate(double tf);

  std::pair<vector<double>, vector<double>> generate_points(double tf, Trajectory& traj);
};


vector<vector<double>> Path::interpolate(double tf) {
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
  return { s_poly, d_poly };
}

std::pair<vector<double>, vector<double>> Path::generate_points(double tf, Trajectory& traj) {
  auto poly = interpolate(tf);
  vector<double> xpts, ypts;
  double inc = tf/50.0;
  cout << "inc - "<< inc << endl;
  for (int i=0; i<50; i++) {
    double t=inc*(i+1), s=0, d=0;
    for (int j=0; j<poly[0].size(); j++) s += poly[0][j]*pow(t,j);
    for (int j=0; j<poly[1].size(); j++) d += poly[1][j]*pow(t,j);
    auto xy = traj.getXY(s, d);
    xpts.push_back(xy[0]);
    ypts.push_back(xy[1]);
  }
  return std::make_pair(xpts, ypts);
}



Trajectory::Trajectory(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys) {
  this->xs = xs;
  this->ys = ys;
  this->ss = ss;
  this->dxs = dxs;
  this->dys = dys;
  this->size = xs.size();
}

vector<double> Trajectory::slice(vector<double> vec, int start, int c) {
  if (start + c < vec.size())
    return vector<double>(vec.begin() + start, vec.begin() + start + c);
  int extra = (start + c) - vec.size();
  vector<double> res = vector<double>(vec.begin() + start, vec.end());
  res.insert(res.end(), vec.begin(), vec.begin() + extra);
  return res;
}

int Trajectory::closest_waypoint(Car car) {
  vector<double> distances(this->size);
  transform(this->xs.begin(), this->xs.end(), this->ys.begin(), distances.begin(),
            [car](const double& x, const double& y) { return euclidean_distance(x, y, car.x, car.y); });
  return std::distance(distances.begin(), min_element(distances.begin(), distances.end()));
}

int Trajectory::next_waypoint(Car car) {
  int next = closest_waypoint(car);
  double wx = this->xs[next], wy = this->ys[next];
  double angle = atan2(wy-car.y, wx-car.x);
  if (fabs(angle-car.yaw) > M_PI/4) {
    next ++;
  }
  return next;
}

Path Trajectory::next_path(Car car, double speed_limit) {
  int wp = next_waypoint(car);
  vector<double> s_vec = {car.s, car.speed, 0, this->ss[wp], speed_limit, 8.0};
  vector<double> d_vec = {car.d, 0, 0, car.d, 0, 0};
  return Path(s_vec, d_vec);
}

vector<double> Trajectory::getXY(double s, double d) {
  return ::getXY(s, d, this->ss, this->xs, this->ys);
}

void Trajectory::plot_highway(Car car, int section = 10000, int width = 50) {
  plt::plot({car.x}, {car.y}, "r*");
  this->plot_highway(section, width);
}

void Trajectory::plot_highway(int section = 10000, int width = 50) {
  if (section > dxs.size()) {
    plt::plot(xs, ys);
    section = dxs.size();
  }
  for (int i = 0; i < section; ++i) {
    plt::plot({xs[i], xs[i]+width*dxs[i]}, {ys[i], ys[i]+width*dys[i]}, "g");
    plt::plot({xs[i]+width*dxs[i]}, {ys[i]+width*dys[i]}, "g.");
  }
}
