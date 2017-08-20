#include <iostream>
#include <fstream>
#include <sstream>


#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;


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


class Trajectory {

 private:
  vector<double> xs, ys, ss, dxs, dys;
  int size;

 public:

  Trajectory(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys) {
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

  vector<vector<double>> next_path(Car car) {
    int i=0, j=1;
    while(true) {
      if (this->ss[i]<= car.s && car.s < this->ss[j]) {
        vector<double> xs = slice(this->xs, j, 50);
        vector<double> ys = slice(this->ys, j, 50);
        return { xs, ys };
      }
      if (this->ss[j]<= car.s && car.s < this->ss[i]) {
        vector<double> xs = slice(this->xs, i, 50);
        vector<double> ys = slice(this->ys, i, 50);
        return { xs, ys };
      }
      i = (i+1) % this->size;
      j = (j+1) % this->size;
    }

  }

  void plot_highway() {
    int width = 50;
    plt::plot(xs, ys);
    for (int i = 0; i < dxs.size(); ++i) {
      plt::plot({xs[i], xs[i]+width*dxs[i]}, {ys[i], ys[i]+width*dys[i]}, "g");
      plt::plot({xs[i]+width*dxs[i]}, {ys[i]+width*dys[i]}, "g.");
    }
  }

  void plot_car_in_highway(Car car) {
    plt::plot({car.x}, {car.y}, "r*");
  }

};

