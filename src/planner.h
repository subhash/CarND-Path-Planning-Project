#include "spline.h"
#include <map>

/*
 * planner.h
 *
 *  Created on: 23-Aug-2017
 *      Author: subhash
 */

#ifndef PLANNER_H_
#define PLANNER_H_

using namespace std;

const double speed_conv = 0.44703, lane_width = 4.0;


// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


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


double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2( (map_y-y),(map_x-x) );

  double angle = abs(theta-heading);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}


class Planner {

 private:
  vector<double> xs, ys, ss, dxs, dys;
  int size;
  tk::spline WP_spline_x, WP_spline_y, WP_spline_dx, WP_spline_dy;

 public:

  Planner(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys) {
    this->xs = xs;
    this->ys = ys;
    this->ss = ss;
    this->dxs = dxs;
    this->dys = dys;
    this->size = xs.size();

    WP_spline_x.set_points(ss, xs);
    WP_spline_y.set_points(ss, ys);
    WP_spline_dx.set_points(ss, dxs);
    WP_spline_dy.set_points(ss, dys);
  }


  vector<double> getXY(double s, double d) {
    double x = WP_spline_x(s), y = WP_spline_y(s);
    double dx = WP_spline_dx(s) * d, dy = WP_spline_dy(s) * d;
    return { x+dx, y+dy };
  }


};

class Vehicle {
 public:
  double x, y, s, d, yaw, v, a, j;
  vector<double> xs, ys, ss;
  vector<double> vs, as, js;
  vector<double> verr_xs, verr_ys, verr, aerr_xs, aerr_ys, aerr;
  int lane;
  int iter = 0;
  bool initialized = false;
  Planner& planner;

  Vehicle(Planner& p): planner(p) {
    this->initialized = false;
  }

  void init(double x, double y, double s, double d, double yaw, double speed) {
    this->initialized = true;
    this->x = x;
    this->y = y;
    this->s = s;
    this->d = d;
    this->yaw = yaw;
    this->v = speed;
    this->a = 0;
    this->j = 0;

    this-> lane = int(d/lane_width);

    this->xs.push_back(x);
    this->ys.push_back(y);
    this->ss.push_back(s);
  }

  void trim(int n) {
    if (n <= this->ss.size()) {
     this->ss.erase(this->ss.begin(), this->ss.begin()+n);
     this->xs.erase(this->xs.begin(), this->xs.begin()+n);
     this->ys.erase(this->ys.begin(), this->ys.begin()+n);
    }
  }

  double inc_velocity(double time, double speed_limit, double acc_limit, double jerk_limit) {
//    double max_j = jerk_limit * time;
//
//    double a = this->a;
//    double max_a_inc = max_j * time;
//    double a_inc = acc_limit - a;
//    if (a_inc > max_a_inc) a_inc = max_a_inc;
//    if (a_inc > 0) a += a_inc;
//
//    double v = this->v;
//    double max_v_inc = a * time;
//    double v_inc = speed_limit - v;
//    if (v_inc > max_v_inc) v_inc = max_v_inc;
//    if (v_inc > 0) v += v_inc;
//    return v;

    double a = this->a;
    double a_diff = jerk_limit*time*time;
    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
    double v = this->v;
    double v_diff = a*time;
    if (v + v_diff < speed_limit) v += v_diff;

    cout << "inc v - "<< this->iter<< "  " << v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
    return v;
  }

  double dec_velocity(double time, double speed_limit, double acc_limit, double jerk_limit) {
//    double max_j = jerk_limit * time;
//
//    double a = this->a;
//    double max_a_inc = max_j * time;
//    double a_inc = acc_limit - a;
//    if (a_inc > max_a_inc) a_inc = max_a_inc;
//    if (a_inc > 0) a += a_inc;
//
//    double v = this->v;
//    double max_v_inc = a * time;
//    double v_inc = speed_limit - v;
//    if (v_inc > max_v_inc) v_inc = max_v_inc;
//    if (v_inc > 0) v += v_inc;
//    return v;

    acc_limit *= -1;
    jerk_limit *= -1;

    double a = this->a;
    double a_diff = jerk_limit*time*time;
    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
    double v = this->v;
    double v_diff = a*time;
    if (v + v_diff > speed_limit) v += v_diff;

    cout << "dec v - "<< v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
    return v;
  }


  double velocity(double x, double y, double time) {
    double dist = sqrt(pow(x-this->x, 2) + pow(y-this->y, 2));
    double v = dist/time;
    return v;
  }


  void move(int steps, double dest_d, double time, double speed_limit, double acc_limit, double jerk_limit) {
    this->trim(steps);
    for (int i = 0; i < steps; ++i) {
      double v;
      if (this->v < speed_limit)
        v = this->inc_velocity(time, speed_limit, acc_limit, jerk_limit);
      else
        v = this->dec_velocity(time, speed_limit, acc_limit, jerk_limit);
      this->step(time, dest_d, v, speed_limit, acc_limit, this->planner);
    }
  }

  void step(double time, double d, double v, double speed_limit, double acc_limit, Planner planner) {
    this->iter++;
    double s_diff = v * time;
    double d_diff = (d - this->d) * time * 0.5;
    vector<double> xy = planner.getXY(this->s + s_diff, this->d + d_diff);
    double corr_v = velocity(xy[0], xy[1], time);

    if (corr_v > speed_limit) {
      this->verr_xs.push_back(xy[0]);
      this->verr_ys.push_back(xy[1]);
      this->verr.push_back(corr_v-speed_limit);
    }

    double a = (v - this->v)/time;

    double j = (a - this->a)/time;
    double x = xy[0], y = xy[1];

    this->s += s_diff;
    this->d += d_diff;
    this->x = x;
    this->y = y;
    this->v = v;
    if(this->iter > 1) this->a = a;
    if(this->iter > 2) this->j = j;

    this->xs.push_back(x);
    this->ys.push_back(y);
    this->ss.push_back(s);
    this->vs.push_back(v);
    if(this->iter > 1) this->as.push_back(a);
    if(this->iter > 2) this->js.push_back(j);
  }

  void print_stats() {
    int min_acc = min_element(this->as.begin()+2, this->as.end()) - this->as.begin();
    int max_acc = max_element(this->as.begin()+2, this->as.end()) - this->as.begin();
    int min_jerk = min_element(this->js.begin()+3, this->js.end()) - this->js.begin();
    int max_jerk = max_element(this->js.begin()+3, this->js.end()) - this->js.begin();
    cout << "Acc: " << this->as[min_acc] << "["<< min_acc << "], " << this->as[max_acc] << "["<< max_acc << "]" << endl;
    cout << "Jrk: " << this->js[min_jerk] << "["<< min_jerk << "], " << this->js[max_jerk] << "["<< max_jerk << "]" << endl;
  }
};



class Obstacle {
 public:
  int id;
  double x, y, vx, vy, s, d;

  void update(int id, double x, double y, double vx, double vy, double s, double d) {
    this->id = id;
    this->x = x;
    this->y = y;
    this->vx = vx;
    this->vy = vy;
    this->s = s;
    this->d = d;
  }

  double speed() {
    return sqrt(pow(this->vx, 2) + pow(this->vy, 2));
  }

};


class Environment {

 public:

  map<int, Obstacle> obstacles;

  void update(int id, double x, double y, double vx, double vy, double s, double d) {
    obstacles[id].update(id, x, y, vx, vy, s, d);
  }

  int leading_obstacle(int lane, double s) {
    double min_dist = std::numeric_limits<double>::max();
    int min_id = -1;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      double dist = o.s - s;
      if (dist > 0 && dist < min_dist && o.d > lane * lane_width && o.d < (lane+1)*lane_width) {
        min_dist = dist;
        min_id = el.first;
      }
    }
    return min_id;
  }

  double leading_distance(int lane, double s) {
    int id = leading_obstacle(lane, s);
    if (id == -1) {
      return std::numeric_limits<double>::max();
    } else {
      return obstacles[id].s - s;
    }
  }

  double lane_speed(int lane, double s) {
    int id = leading_obstacle(lane, s);
    if (id == -1) {
      return 0;
    } else {
      return obstacles[id].speed();
    }
  }

};



class Behaviour {

 protected:
  const double time = 0.02, speed_limit = 49.0 * speed_conv, acc_limit = 10, jerk_limit = 10;

 public:
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) {
    double dest_d = vehicle.d;
//    const double lane_speed = env.lane_speed(vehicle.lane, vehicle.s);
//    if (lane_speed < speed_limit) {
//      cout << "Limiting speed to "<< min(speed_limit, lane_speed) << endl;
//    }
//    vehicle.move(nsteps, dest_d, time, min(speed_limit, lane_speed), acc_limit, jerk_limit);
    vehicle.move(nsteps, dest_d, time, speed_limit, acc_limit, jerk_limit);
  }

  virtual ~Behaviour() {}

};

class LeftChangeBehaviour: public Behaviour {
 public:
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override {
    double dest_d = (vehicle.lane - 1)*lane_width + 2;
    vehicle.move(nsteps, dest_d, time, speed_limit, acc_limit, jerk_limit);
    if (fabs(vehicle.d - dest_d) < 0.0001){
      vehicle.lane -= 1;
    }
  }
};

class SlowingDownBehaviour: public Behaviour {
 public:
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override {
    vehicle.move(nsteps, vehicle.d, time, 5, acc_limit, jerk_limit);
  }
};


class BehaviourPlanner {

 private:
  Behaviour keepLane;
  LeftChangeBehaviour changeLane;
  SlowingDownBehaviour slowDown;

 public:
  Behaviour& behaviour(Vehicle& vehicle) {
    if (vehicle.iter > 1000 && vehicle.lane != 0) {
      cout << "Slow down at "<< vehicle.v << " and " << vehicle.a <<endl;
      return slowDown;
    }
    return keepLane;
  }
};



#endif /* PLANNER_H_ */
