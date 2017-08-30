#include "spline.h"
#include <map>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"


#ifndef PLANNER_H_
#define PLANNER_H_

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  int lane() {
    return int(d/lane_width);
  }

  double inc_velocity(double time, double init_v, double init_a, double v_desired, double acc_limit, double jerk_limit) {
    double a = init_a;
    double a_diff = jerk_limit*time;
    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
    double v = init_v;
    double v_diff = a*time;
    if (v + v_diff < v_desired) v += v_diff;

    //cout << "inc v - "<< this->iter<< "  " << v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
    return v;
  }

  double dec_velocity(double time, double init_v, double init_a, double v_desired, double acc_limit, double jerk_limit) {
    acc_limit *= -1;
    jerk_limit *= -1;

    double a = init_a;
    double a_diff = jerk_limit*time;
    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
    double v = init_v;
    double v_diff = a*time;
    //v_diff = -0.1;
    if (v + v_diff > v_desired) v += v_diff;

    //cout << "dec v - "<< v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
    return v;
  }


  double velocity(double x, double y, double time) {
    double dist = sqrt(pow(x-this->x, 2) + pow(y-this->y, 2));
    double v = dist/time;
    return v;
  }

  void move(int steps, double dest_d, double time, double speed_limit, double acc_limit, double jerk_limit, vector<double> d_poly = {}, int nsteps = 0) {
    this->trim(steps);
    for (int i = 0; i < steps; ++i) {
      double v;
      if (this->v < speed_limit)
        v = this->inc_velocity(time, this->v, this->a, speed_limit, acc_limit, jerk_limit);
      else
        v = this->dec_velocity(time, this->v, this->a, speed_limit, acc_limit, jerk_limit);
      double d_diff = 0;
      if (!d_poly.empty()) {
        double tf = time*(i+nsteps);
        for (int c = 0; c < d_poly.size(); ++c) {
          d_diff += d_poly[c]*pow(tf, c);
        }
      }
      this->step(time, dest_d, v, d_diff, this->planner);
    }
  }

  void step(double time, double d, double v, double d_diff, Planner planner) {
    this->iter++;

    double s_diff = v * time;
    // Just to support turns
    if (d_diff > 0) this->d = d_diff;
    vector<double> xy = planner.getXY(this->s + s_diff, this->d);

//    double corr_v = velocity(xy[0], xy[1], time);
//    if (corr_v > speed_limit) {
//      this->verr_xs.push_back(xy[0]);
//      this->verr_ys.push_back(xy[1]);
//      this->verr.push_back(corr_v-speed_limit);
//    }

    double a = (v - this->v)/time;
    double j = (a - this->a)/time;
    double x = xy[0], y = xy[1];

    this->s += s_diff;
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

  int lane() {
    return int(d/lane_width);
  }

  double speed() {
    return sqrt(pow(this->vx, 2) + pow(this->vy, 2));
  }

};


class Environment {

 private:
  const double time = 0.02, speed_limit = 49.0 * speed_conv, acc_limit = 2, jerk_limit = 2;

 public:

  map<int, Obstacle> obstacles;

  void update(int id, double x, double y, double vx, double vy, double s, double d) {
    obstacles[id].update(id, x, y, vx, vy, s, d);
  }

  int closest_obstacle(int lane, double s) {
    double min_dist = std::numeric_limits<double>::max();
    int min_id = -1;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      double dist = fabs(o.s - s);
      if (o.lane() == lane && dist < min_dist) {
        min_dist = dist;
        min_id = el.first;
      }
    }
    return min_id;
  }

  bool will_vehicle_collide(int lane, Vehicle& vehicle) {
    double band = vehicle.v * 2, min_s = vehicle.s - band, max_s = vehicle.s + band;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      if (o.lane() == lane && o.s > min_s and o.s < max_s) {
        return true;
      }
    }
    return false;
  }

  bool will_obstacle_collide(int lane, Vehicle& vehicle) {
    double band = vehicle.v * 2, min_s = vehicle.s, max_s = vehicle.s + band;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      if (o.lane() == lane && o.s > min_s and o.s < max_s and fabs(o.d - vehicle.d) < 0.9*lane_width) {
        return true;
      }
    }
    return false;
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

  double lane_speed(Vehicle vehicle, int lane) {
    int id = leading_obstacle(lane, vehicle.s);
    if (id != -1) {
      double gap = obstacles[id].s - vehicle.s;
      if (gap < vehicle.v*3) {
        return obstacles[id].speed();
      }
    }
    return speed_limit;
  }

};

class BehaviourPlanner;

class Behaviour {

 protected:
  double time = 0.02, speed_limit = 49.0 * speed_conv, acc_limit = 2, jerk_limit = 2;
  BehaviourPlanner& behaviour_planner;

  double decide_speed(Vehicle& vehicle, Environment& env, double speed_limit) {
    return env.lane_speed(vehicle, vehicle.lane());
  }

 public:
  string name;
  Behaviour(string name, BehaviourPlanner& bp): name(name), behaviour_planner(bp) {}
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) { }
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env) { return { *this }; }
  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) = 0;
  virtual ~Behaviour() {}
};


class KeepLane: public Behaviour {
 public:
  KeepLane(BehaviourPlanner& bp, string name): Behaviour(name, bp) {};
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env) override;

  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
    if(debug) cout << this->name << ": ";
    int lane = vehicle.lane(), llane = lane-1, rlane = lane+1;
    double s = vehicle.s, d = vehicle.d;
    double speed_cost = env.lane_speed(vehicle, vehicle.lane()) < speed_limit;
    double safety_cost = (llane >= 0) && env.will_obstacle_collide(llane, vehicle) +
        (rlane <= 2) && env.will_obstacle_collide(rlane, vehicle);
    if (debug) cout << "speed_cost=" << speed_cost << ", safety_cost=" << safety_cost;
    if(debug) cout << ", llane=" << env.will_obstacle_collide(llane, vehicle) << ", rlane=" << env.will_obstacle_collide(rlane, vehicle);
    return 2*speed_cost + 10*safety_cost;
  }
};

class LaneChange: public Behaviour {
  int dest_lane = 0;
  vector<double> d_poly;
  int steps = 0;
 public:
  LaneChange(BehaviourPlanner& bp, string name): Behaviour(name, bp) {}
  void set_lane(Vehicle& vehicle, int dest_lane) {
    this->dest_lane = dest_lane;
    double dest_d = dest_lane*lane_width + 2;
    d_poly = make_poly(3, vehicle.d, dest_d);
    steps = 0;
  }
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env) override;

  vector<double> make_poly(double tf, double d_init, double d_final) {
    double t = tf, t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t;
    MatrixXd Tmat(3, 3);
    Tmat << t3, t4, t5,
            3*t2, 4*t3, 5*t4,
            6*t, 12*t2, 20*t3;
    MatrixXd Tinv = Tmat.inverse();
    VectorXd Dmat(3);
    Dmat << d_final - d_init, 0, 0;
    VectorXd Dvec = Tinv*Dmat;
    vector<double> d_poly = { d_init, 0, 0, Dvec(0), Dvec(1), Dvec(2) };
    return d_poly;
  }

  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
    if(debug) cout << this->name << ": ";
    double s = vehicle.s, d = vehicle.d;
    double speed_cost = env.lane_speed(vehicle, dest_lane) < speed_limit;
    double safety_cost = env.will_vehicle_collide(dest_lane, vehicle);
    if (debug) cout << "speed_cost=" << speed_cost << ", safety_cost=" << safety_cost;
    return 1+2*speed_cost + 10*safety_cost;
  }
};

class SlowDown: public Behaviour {
 public:
  SlowDown(BehaviourPlanner& bp, string name): Behaviour(name, bp) {};
  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
    double speed_cost = 1;
    if (debug) cout << "speed_cost=" << speed_cost;
    return speed_cost;
  }
};


class BehaviourPlanner {

 private:
  KeepLane kl;
  SlowDown sd;
  LaneChange lc;
  LaneChange rc;
  Behaviour* last_action = &kl;

 public:
  BehaviourPlanner(): kl(*this, "KL"), sd(*this, "SD"), lc(*this, "LC"), rc(*this, "RC") {
  }

  Behaviour& keepLane() {
    return this->kl;
  }

  Behaviour& leftChange(Vehicle& vehicle, int lane) {
    this->lc.set_lane(vehicle, lane);
    return this->lc;
  }

  Behaviour& rightChange(Vehicle& vehicle, int lane) {
    this->rc.set_lane(vehicle, lane);
    return this->rc;
  }

  Behaviour& slowDown() {
    return this->sd;
  }

  Behaviour& behaviour(Vehicle& vehicle, Environment& env) {
    auto actions = last_action->next_actions(vehicle, env);
    auto min_action = *min_element(actions.begin(), actions.end(),
         [&vehicle, &env](const std::reference_wrapper<Behaviour> a, const std::reference_wrapper<Behaviour> b) {
            return a.get().cost(vehicle, env) < b.get().cost(vehicle, env);
          } );
    auto next_action = &min_action.get();
    bool changed = last_action!=next_action;
    if (changed) {
      cout << last_action->name << endl;
      for (auto a: actions) {
        Behaviour& b = a.get();
        b.cost(vehicle, env, changed);
        cout << endl;
      }
    }
    last_action = next_action;
    return *last_action;
  }
};

void KeepLane::effect(Vehicle& vehicle, Environment& env, int nsteps) {
  double dest_d = vehicle.lane()*lane_width + 2;
  vehicle.move(nsteps, dest_d, time, decide_speed(vehicle, env, speed_limit), acc_limit, jerk_limit);
}

void LaneChange::effect(Vehicle& vehicle, Environment& env, int nsteps) {
  double dest_d = dest_lane*lane_width + 2;
  vehicle.move(nsteps, vehicle.d, time, decide_speed(vehicle, env, speed_limit), acc_limit, jerk_limit, d_poly, this->steps);
  this->steps += nsteps;
}

void SlowDown::effect(Vehicle& vehicle, Environment& env, int nsteps) {
  vehicle.move(nsteps, vehicle.d, time, 5, acc_limit, jerk_limit);
}

vector<std::reference_wrapper<Behaviour>> KeepLane::next_actions(Vehicle& vehicle, Environment& env) {
  vector<reference_wrapper<Behaviour>> actions = { behaviour_planner.keepLane() };
//  if (vehicle.iter > 500) {
//    return { behaviour_planner.leftChange(vehicle, 0) };
//  }
  if (vehicle.lane() > 0) {
    Behaviour& lc = behaviour_planner.leftChange(vehicle, vehicle.lane() - 1);
    actions.push_back(lc);
  }
  if (vehicle.lane() < 2) {
    Behaviour& lc = behaviour_planner.rightChange(vehicle, vehicle.lane() + 1);
    actions.push_back(lc);
  }
  return actions;
}

vector<std::reference_wrapper<Behaviour>> LaneChange::next_actions(Vehicle& vehicle, Environment& env) {
  double dest_d = dest_lane*lane_width + 2;
  if (this->steps < 3*50){
    return { *this };
  }
  return { behaviour_planner.keepLane() };
}


#endif /* PLANNER_H_ */
