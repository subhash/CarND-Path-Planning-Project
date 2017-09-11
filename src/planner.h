#include "spline.h"
#include <map>
#include <queue>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;


#ifndef PLANNER_H_
#define PLANNER_H_

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const double speed_conv = 0.44703, lane_width = 4.0, speed_limit = 49.0 * speed_conv, timestep = 0.02, max_s = 6945.554;
//  double time = 0.02, speed_limit = 49.0 * speed_conv, acc_limit = 2, jerk_limit = 2;



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

 public:

  vector<double> xs, ys, ss, dxs, dys, thetas;
  int size;
  tk::spline WP_spline_x, WP_spline_y, WP_spline_dx, WP_spline_dy, WP_spline_theta;
  std::vector<tk::spline> lane_spline_s, lane_spline_rev_s;

  Planner(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys):
    lane_spline_s(3), lane_spline_rev_s(3) {

//    xs.erase(xs.begin());
//    ys.erase(ys.begin());
//    ss.erase(ss.begin());
//    dxs.erase(dxs.begin());
//    dys.erase(dys.begin());

    this->xs = xs;
    this->ys = ys;
    this->ss = ss;
    this->dxs = dxs;
    this->dys = dys;
    this->size = xs.size();

    this->xs.push_back(this->xs[0]);
    this->ys.push_back(this->ys[0]);
    this->ss.push_back(max_s);
    this->dxs.push_back(this->dxs[0]);
    this->dys.push_back(this->dys[0]);
    this->size += 1;

    WP_spline_x.set_points(this->ss, this->xs);
    WP_spline_y.set_points(this->ss, this->ys);
    WP_spline_dx.set_points(this->ss, this->dxs);
    WP_spline_dy.set_points(this->ss, this->dys);


    for (int i = 0; i < this->size; ++i) {
      double theta = atan2(this->dys[i], this->dxs[i]);
      this->thetas.push_back(theta);
    }

    WP_spline_theta.set_points(this->ss, this->thetas);

    for (int i=0; i < 3; i++) {
      this->setup_lane_splines(i);
    }

  }


  void setup_lane_splines(int lane) {
    vector<double> ss1, xs1, ys1, dxs1, dys1;
    double lane_offset = lane*lane_width + lane_width/2.0;
    double err = 0.0;

    ss1.push_back(this->ss[0]);
    for (int i = 1; i < this->size; ++i) {
      int prev_pt = i-1;
      if (prev_pt < 0) {
        prev_pt = this->size-1;
      }
      double prev_s = this->ss[prev_pt], s = this->ss[i];
      double theta1 = WP_spline_theta(prev_s), theta2 = WP_spline_theta(s);
      err += lane_offset*tan(theta2 - theta1);
      double proj_s = s + err;
      if (proj_s < 0) cout << "s - " << s << "proj_s - " << proj_s << ", " << i << endl;
      ss1.push_back(proj_s);
    }

//    vector<double> new_ss(this->ss);
//    {
//      new_ss.push_back(max_s);
//      cout << new_ss.size() << ", "<< this->ss.size() << endl;
//      double prev_s = this->ss[this->size-1], s = max_s;
//      double theta1 = WP_spline_theta(prev_s), theta2 = WP_spline_theta(s);
//      err += lane_offset*tan(theta2 - theta1);
//      double proj_s = s + err;
//      ss1.push_back(proj_s);
//    }


    lane_spline_s[lane].set_points(ss1, this->ss);
    lane_spline_rev_s[lane].set_points(this->ss, ss1);
  }

  vector<double> getXY(double s, double d) {
    s = fmod(s, max_s);
    double x = WP_spline_x(s), y = WP_spline_y(s);
    double dx = WP_spline_dx(s) * d, dy = WP_spline_dy(s) * d;
    return { x+dx, y+dy };
  }

  double waypoint_s(double ln_s, int ref_lane) {
    double lane_max_s = lane_s(max_s, ref_lane);
    ln_s = fmod(ln_s, lane_max_s);
    return lane_spline_s[ref_lane](ln_s);
  }

  double lane_s(double wp_s, int ref_lane) {
    return lane_spline_rev_s[ref_lane](wp_s);
  }

  void plot(int from = 0, int to = 10000, int width = 50) {
    if (to > dxs.size()) {
      //plt::plot(xs, ys);
      to = dxs.size();
    }
    for (int i = from; i < to; ++i) {
      plt::plot({xs[i], xs[i]+width*dxs[i]}, {ys[i], ys[i]+width*dys[i]}, "g");
      plt::plot({xs[i]+width*dxs[i]}, {ys[i]+width*dys[i]}, "g.");
    }
    auto start_xy = this->getXY(0, 0);
    plt::plot({start_xy[0]}, {start_xy[1]}, "rs");
//    double start_s = from*30, end_s = to*30;
//    for (int d = 0; d < 16; d+=4) {
//      vector<double> track_x, track_y;
//      for (double s = start_s; s < end_s; s+=10.0) {
//        vector<double> xy = getXY(s, d);
//        track_x.push_back(xy[0]);
//        track_y.push_back(xy[1]);
//      }
//      plt::plot(track_x, track_y, "b-");
//    }
  }

  void plot_splines(double from_s, double to_s) {
    vector<double> xpts1, ypts1, xpts2, ypts2, xpts3, ypts3;
    for (double s=from_s; s<to_s; s+=10){
      //double throttled_s = fmod(s, max_s);
      auto xy1 = this->getXY(s, 0);
      xpts1.push_back(xy1[0]);
      ypts1.push_back(xy1[1]);

      double s2 = this->waypoint_s(s, 1);
      auto xy2 = this->getXY(s2, 6);
      xpts2.push_back(xy2[0]);
      ypts2.push_back(xy2[1]);

      double s3 = this->waypoint_s(s, 2);
      auto xy3 = this->getXY(s3, 10);
      xpts3.push_back(xy3[0]);
      ypts3.push_back(xy3[1]);

    }
    plt::plot(xpts1, ypts1, "g-");
    plt::plot(xpts2, ypts2, "r-");
    plt::plot(xpts3, ypts3, "b-");

  }

};

class Traj {

 private:

  double poly_eval(vector<double> poly, double dt) {
    double v = 0.0;
    for (int i = 0; i < poly.size(); ++i) {
      v += poly[i] * pow(dt, i);
    }
    return v;
  }

  double v_eval(vector<double> poly, double dt) {
    double v = 0.0;
    for (int i = 1; i < poly.size(); ++i) {
      v += i * poly[i] * pow(dt, i-1);
    }
    return v;
  }

  double a_eval(vector<double> poly, double dt) {
    double a = 0.0;
    for (int i = 2; i < poly.size(); ++i) {
      a += i * poly[i] * (i-1) * pow(dt, i-2);
    }
    return a;
  }


  std::pair<vector<double>, vector<double>> make_poly(
      double tf, vector<double> s_vec, vector<double> d_vec) {
    double t = tf, t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t;
    double si = s_vec[0], si_dot = s_vec[1], si_ddot = s_vec[2];
    double sf = s_vec[3], sf_dot = s_vec[4], sf_ddot = s_vec[5];
    double di = d_vec[0], di_dot = d_vec[1], di_ddot = d_vec[2];
    double df = d_vec[3], df_dot = d_vec[4], df_ddot = d_vec[5];
    MatrixXd Tmat(3, 3);
    Tmat << t3, t4, t5,
            3*t2, 4*t3, 5*t4,
            6*t, 12*t2, 20*t3;
    MatrixXd Tinv = Tmat.inverse();
    VectorXd Smat(3);
    Smat << sf - (si + si_dot * t + (si_ddot/2) * t2),
            sf_dot - (si_dot + si_ddot * t),
            sf_ddot - (si_ddot);
    VectorXd Dmat(3);
    Dmat << df - (di + di_dot * t + (di_ddot/2) * t2),
            df_dot - (di_dot + di_ddot * t),
            df_ddot - (di_ddot);
    VectorXd Svec = Tinv*Smat;
    VectorXd Dvec = Tinv*Dmat;
    vector<double> s_poly = { si, si_dot, si_ddot/2.0, Svec(0), Svec(1), Svec(2) };
    vector<double> d_poly = { di, di_dot, di_ddot/2.0, Dvec(0), Dvec(1), Dvec(2) };
    return std::make_pair(s_poly, d_poly);
  }

  vector<double> spoly, dpoly;


public:

  double duration;
  int ref_lane=-1;

  vector<double> s_vec, d_vec;

  Traj(double duration, vector<double> s_vec, vector<double> d_vec, int ref_lane, Planner& p): ref_lane(ref_lane) {
    this->s_vec = s_vec;
    this->d_vec = d_vec;
    this->duration = duration;
    auto polys = make_poly(duration, s_vec, d_vec);
    this->spoly = polys.first;
    this->dpoly = polys.second;
  }

  vector<pair<double, double>> points(double dt, Planner& p) {
    int steps = this->duration/dt;
    vector<pair<double, double>> pts;
    for (int i = 0; i < steps; ++i) {
      double s = poly_eval(this->spoly, (i+1) * dt);
      double d = poly_eval(this->dpoly, (i+1) * dt);
      double proj_s = p.waypoint_s(s, ref_lane);
      pts.push_back(make_pair(proj_s, d));
    }
    return pts;
  }

  vector<pair<double, double>> samples(int n, Planner& p) {
    int total_pts = this->duration/timestep;
    int dt = int(total_pts/n) * timestep;
    return this->points(dt, p);
  }

  vector<vector<double>> motion_vector(double dt, Planner& p) {
    int steps = this->duration/dt;
    vector<vector<double>> pts;
    for (int i = 0; i < steps; ++i) {
      double s = poly_eval(this->spoly, (i+1) * dt);
      double d = poly_eval(this->dpoly, (i+1) * dt);
      double proj_s = p.waypoint_s(s, ref_lane);
      double v = v_eval(this->spoly, (i+1)*dt);
      double a = a_eval(this->spoly, (i+1)*dt);
      pts.push_back( { proj_s, d, v, a });
    }
    return pts;
  }


  void plot(double interval, Planner& p) {
    auto pts = this->points(interval, p);
    vector<double> xpts, ypts;
    for (auto pt: pts) {
      auto xy = p.getXY(pt.first, pt.second);
      xpts.push_back(xy[0]);
      ypts.push_back(xy[1]);
    }
    plt::plot(xpts, ypts, "b+");

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

  pair<double, double> velocity(double x, double y, double time) {
    double dist = sqrt(pow(x-this->x, 2) + pow(y-this->y, 2));
    double v = dist/time;
    double a = (v-this->v)/time;
    return make_pair(v, a);
  }

  void move_to(double proj_s, double d, double v, double a, double dt, Planner planner) {
    this->iter++;

    vector<double> xy = planner.getXY(proj_s, d);
    double s = proj_s; //planner.lane_spline_rev_s[this->lane()](proj_s);
    double x = xy[0], y = xy[1];

    this->s = s;
    this->d = d;
    this->x = x;
    this->y = y;
    this->v = v;
    this->a = a;

    this->xs.push_back(x);
    this->ys.push_back(y);
    this->ss.push_back(s);
    this->vs.push_back(v);
  }

  void plot() {
    plt::plot({this->x}, {this->y}, "g+");
  }

};


class Obstacle {
 public:
  int id;
  double x, y, vx, vy, s, d;
  vector<double> xs, ys, ss, ds, svs, dvs, speeds, tss;
  double sv, sa, dv, da;

  void update(int id, double x, double y, double vx, double vy, double s, double d, double ts) {
    this->tss.push_back(ts);

    // has history?
    if (this->ss.size() > 0) {
      double sv = (s - this->s)/ts;
      double dv = (d - this->d)/ts;
      if (this->ss.size() > 1) {
        double sa = (sv - this->sv)/ts;
        double da = (dv - this->dv)/ts;
        this->sa = sa;
        this->da = da;
      } else {
        this->sa = 0;
        this->da = 0;
      }
      this->sv = sv;
      this->dv = dv;
    } else {
      this->sv = 0;
      this->dv = 0;
    }


    this->id = id;
    this->x = x;
    this->y = y;
    this->vx = vx;
    this->vy = vy;
    this->s = s;
    this->d = d;

    this->xs.push_back(this->x);
    this->ys.push_back(this->y);
    this->ss.push_back(this->s);
    this->ds.push_back(this->d);
    this->svs.push_back(this->sv);
    this->dvs.push_back(this->dv);

    this->speeds.push_back(this->speed());
  }

  int lane() const {
    return int(d/lane_width);
  }

  double speed() const {
    return sqrt(pow(this->vx, 2) + pow(this->vy, 2));
  }

};


class Prediction {

 private:
  const double s_tol = 20, d_tol = 1.5;
  map<int, vector<Obstacle>> predictions;
  map<int, vector<Obstacle>> states;

 public:
  void update(int nt, Obstacle o) {
    predictions[nt].push_back(o);
    states[o.id].push_back(o);
  }

  double dist(Obstacle o, double s, double d) {
    return sqrt(pow(o.s-s,2) + pow(o.d-d,2));
  }

  bool too_close(int nt, double s, double d) {
    nt += 50;
    vector<Obstacle> obstacles = predictions[nt];
    for (auto const & o: obstacles) {
      if (dist(o, s, d) < 6)
        return true;
    }
    return false;
  }

  bool tail_ends(int nt, double s, double d) {
    nt += 50;
    vector<Obstacle> obstacles = predictions[nt];
    for (auto const & o: obstacles) {
      // same lane, in front
      if (o.s > s && (o.lane() == int(d/lane_width))) {
        if (dist(o, s, d) < 5)
          return true;
      }
    }
    return false;
  }


  int leading_obstacle(int nt, int lane, double s) {
    double min_dist = std::numeric_limits<double>::max();
    int min_id = -1;
    for (auto const& o : predictions[nt]) {
      double dist = o.s - s;
      if (dist > 0 && dist < min_dist && o.lane() == lane) {
        min_dist = dist;
        min_id = o.id;
      }
    }
    return min_id;
  }

  double lane_speed(int nt, Vehicle vehicle, int lane) {
    nt += 50;
    int id = leading_obstacle(nt, lane, vehicle.s);
    if (id != -1) {
      for (Obstacle& o: states[id]) {
        // 2 second rule
        if (dist(o, vehicle.s, vehicle.d) < vehicle.v * 2) {
          return min(o.speed(), speed_limit);
        }
      }
    }
    return speed_limit;
  }

  void plot(Planner& p, int nt = 100000) {
    for (auto const& el: states) {
      auto obs = el.second;
      if (nt > obs.size()) nt = obs.size();
      for (int i=0; i<nt; i++) {
        Obstacle& o = obs[i];
        auto xy = p.getXY(o.s, o.d);
        plt::plot({xy[0]}, {xy[1]}, "r.");
      }
    }
  }

};


class Environment {

 private:
  const double time = 0.02, acc_limit = 2, jerk_limit = 2;

 public:
  map<int, Obstacle> obstacles;

  double latency;

  void update(int id, double x, double y, double vx, double vy, double s, double d, double ts) {
    obstacles[id].update(id, x, y, vx, vy, s, d, ts);
    this->latency = ts;
  }

  Environment predict(double dt) {
    Environment predicted;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      double d_diff = o.dv*dt;
      predicted.update(o.id, o.x+o.vx*dt, o.y+o.vy*dt, o.vx, o.vy, o.s+o.speed()*dt, o.d+d_diff, dt);
      //predicted.update(o.id, o.x+o.vx*dt, o.y+o.vy*dt, o.vx, o.vy, o.s+o.speed()*dt, o.d, dt);
    }
    return predicted;
  }

  Prediction prediction(double n_timesteps) {
    Prediction pred;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      for (int nt = 0; nt < n_timesteps; ++nt) {
        double dt = (nt+1) * timestep;
        double s_diff = o.speed() * dt;
        double d_diff = o.dv * dt;
        Obstacle po;
        po.update(o.id, o.x+o.vx*dt, o.y+o.vy*dt, o.vx, o.vy, o.s+s_diff, o.d+d_diff, dt);
        pred.update(nt, po);
      }
    }
    return pred;
  }


  void plot() {
    cout << "Plotting "<< obstacles.size() << " obstacles " << endl;
    int size = obstacles.size();
    int rows = 4;

    Environment predicted = this->predict(3.0);

    int i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.ss);
      plt::plot({200 + 3/0.02}, {predicted.obstacles[o.id].s}, "r.");
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.svs);
      plt::plot({200 + 3/0.02}, {predicted.obstacles[o.id].d}, "r.");
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.speeds);
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.ds);
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.dvs);
    }
    plt::show();


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

  double lane_speed(Vehicle vehicle, int lane) {
    int id = leading_obstacle(lane, vehicle.s);
    if (id != -1) {
      double gap = obstacles[id].s - vehicle.s;
      if (gap < 50) {
        return min(obstacles[id].speed(), speed_limit);
      }
    }
    return speed_limit;
  }

  bool too_close(double s, double d, double s_tol, double d_tol) {
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      bool s_clear = ((o.s < s) && (s - o.s < 2*s_tol)) || ((o.s > s) && (o.s - s < s_tol));
      bool d_clear = fabs(o.d - d) < d_tol;
      if (s_clear && d_clear)
        return true;
    }
    return false;
  }

  bool collision_space(double from_s, double to_s, int from_d, double to_d) {
    int from_lane = int(from_d/lane_width), to_lane = int(to_d/lane_width);
    if (from_lane > to_lane) {
      int tmp = from_lane;
      from_lane = to_lane;
      to_lane = from_lane;
    }
    double min_d = from_lane*lane_width, max_d = (to_lane+1)*lane_width;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      if (o.s >= from_s && o.s <= to_s && o.d >= min_d && o.d <= max_d)
        return true;
    }
    return false;
  }

};

class BehaviourPlanner;

class Behaviour {

 protected:

  vector<double> long_vector(double s_init, double v_init, double a_init, double v_final, double &tf) {
    double si = s_init, si_dot = v_init, si_ddot = a_init;
    double sf_dot = v_final;
    if (si_dot < sf_dot/2.0) {
      tf += 2;
      sf_dot *= 0.75;
    }
    double v_ave = 0.5*(si_dot + sf_dot);
    double sf = si + v_ave*tf;
    return { si, si_dot, si_ddot, sf, sf_dot, 0 };
  }

 public:
  string name;
  Behaviour(string name): name(name) {}
  virtual Traj convolute(Vehicle& vehicle, Prediction& pred, Planner& p) = 0;
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) {
    return {};
  }

//  virtual bool will_collide(Prediction& pred, int nt, double s, double d) {
//    return pred.too_close(nt, s, d);
//  }

  map<int, pair<double, double>> path_in_time(Traj traj, Planner& p) {
    int n_pts = 4;
    int total_pts = traj.duration/timestep;
    int interval_steps = int(total_pts/n_pts);
    double interval =  interval_steps * timestep;
    vector<pair<double, double>> pts = traj.points(interval, p);
    map<int, pair<double, double>> path;
    for (int i = 0; i < pts.size(); ++i) {
      path[interval_steps*(i+1)] = pts[i];
    }
    return path;
  }

  virtual bool collides_with(Prediction& pred, int nt, double s, double d) {
    return pred.too_close(nt, s, d);
  }

  virtual double safety_cost(Vehicle& vehicle, Prediction& pred, Traj& traj, Planner& p) {
    map<int, pair<double, double>> path = path_in_time(traj, p);
    double cost = 0;
    for (auto el: path) {
      auto pt = el.second;
      double traj_s = pt.first, traj_d = pt.second;
      if (this->collides_with(pred, el.first, traj_s, traj_d)) {
        cost += 1.0;
      }
    }

//    if (cost > 2 && (this->name == "LC" || this->name == "RC")) {
//      int total_pts = traj.duration/timestep;
//      int interval_steps = int(total_pts/4);
//      double interval =  interval_steps * timestep;
//      //p.plot(0, 10);
//
//      pred.plot(p, 50);
//      traj.plot(interval, p);
//      plt::show();
//      exit(0);
//    }


    return cost;
  }

  virtual double efficiency_cost(Vehicle& vehicle, Prediction& pred, Traj& traj) {
    int dest_lane = int(traj.d_vec[3]/lane_width);
    double dest_v = traj.s_vec[4];
    int last_timestep = int(traj.duration/timestep);
    //double lane_speed = pred.lane_speed(last_timestep, vehicle, dest_lane);
    return (dest_v < 0.9 * vehicle.v) + (dest_v < 0.95*speed_limit);
  }

  virtual double cost(Vehicle& vehicle, Prediction& pred, Planner& p) {
    Traj traj = this->convolute(vehicle, pred, p);
    return 10*this->safety_cost(vehicle, pred, traj, p) + this->efficiency_cost(vehicle, pred, traj);
  }

  void plot(Vehicle& vehicle, Environment& env, Traj& traj, Planner& p) {
    plt::plot(vehicle.xs, vehicle.ys, "r");
    plt::plot({vehicle.x}, {vehicle.y}, "r+");
    Obstacle& o = env.obstacles[1];
    Obstacle& po = env.predict(3.0).obstacles[1];
    vector<double> opt = p.getXY(po.s, po.d);
    double s1 = traj.s_vec[0], s2 = traj.s_vec[3], d1 = traj.d_vec[0], d2= traj.d_vec[3];
    cout << "s diff " << s2-s1 << ", d diff - "<< d2 - d2 << endl;
    vector<double> pt1 = p.getXY(s1, d1), pt2 = p.getXY(s2, d1), pt3= p.getXY(s1, d2), pt4 = p.getXY(s2, d2) ;
    plt::plot({pt1[0], pt2[0], pt3[0], pt4[0]}, {pt1[1], pt2[1], pt3[1], pt4[1]}, "g.");
    plt::show();
    exit(0);
  }

  virtual ~Behaviour() {}
};


class KeepLane: public Behaviour {

 public:

  KeepLane(): Behaviour("KL") {};

  virtual Traj convolute(Vehicle& vehicle, Prediction& pred, Planner& p) override {
    double tf = 3.0; // sec
    double df = vehicle.lane() * lane_width + lane_width/2.0;
    int nt = tf/timestep;
    double sf_dot = pred.lane_speed(0, vehicle, vehicle.lane());

    double si = vehicle.s;
    si = p.lane_s(si, vehicle.lane());

    vector<double> s_vec = this->long_vector(si, vehicle.v, vehicle.a, sf_dot, tf);
    vector<double> d_vec = { vehicle.d, 0, 0, df, 0, 0, 0};
    Traj traj(tf, s_vec, d_vec, vehicle.lane(), p);
    return traj;
  }

  virtual bool collides_with(Prediction& pred, int nt, double s, double d) override {
    return pred.tail_ends(nt, s, d);
  }


  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) override;
};


class LaneChange: public Behaviour {
  int dest_lane = 0;
  vector<double> d_poly;
  int steps = 0;

 public:

  LaneChange(string name, int lane): Behaviour(name), dest_lane(lane) {}

  virtual Traj convolute(Vehicle& vehicle, Prediction& pred, Planner& p) override {
    double tf = 3.0; // sec
    int nt = tf/timestep;
    double df = dest_lane * lane_width + lane_width/2.0;
    double sf_dot = pred.lane_speed(0, vehicle, dest_lane);

    double si = vehicle.s;
    si = p.lane_s(si, vehicle.lane());

    vector<double> s_vec = this->long_vector(si, vehicle.v, vehicle.a, sf_dot, tf);
    vector<double> d_vec = { vehicle.d, 0, 0, df, 0, 0, 0};
    Traj traj(tf, s_vec, d_vec, vehicle.lane(), p);
    return traj;
  }

  void set_lane(Vehicle& vehicle, int dest_lane) {
    this->dest_lane = dest_lane;
  }

  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) override;
};


class SlowDown: public Behaviour {

  double ref_speed = 0.0;

 public:
  SlowDown(): Behaviour("SD") {};
  virtual Traj convolute(Vehicle& vehicle, Prediction& pred, Planner& p) override {
    double tf = 2.0; // sec
    double sf_dot = ref_speed;

    double si = vehicle.s;
    si = p.lane_s(si, vehicle.lane());

    vector<double> s_vec = this->long_vector(si, vehicle.v, vehicle.a, sf_dot, tf);
    double df = vehicle.lane() * lane_width + lane_width/2.0;
    vector<double> d_vec = { vehicle.d, 0, 0, df, 0, 0, 0};
    Traj traj(tf, s_vec, d_vec, vehicle.lane(), p);
    return traj;
  }

  virtual bool collides_with(Prediction& pred, int nt, double s, double d) override {
    return pred.tail_ends(nt, s, d);
  }

  void set_speed(double speed) {
    this->ref_speed = speed;
  }

  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) override;
};


class BehaviourPlanner {

 private:
  KeepLane kl;
  SlowDown sd;
  LaneChange lc;
  LaneChange rc;
  Behaviour* last_action = &kl;

 public:

  BehaviourPlanner(): lc("LC", 0), rc("RC", 0) {}

  Traj convolute(Vehicle& vehicle, Prediction& pred, Planner& p) {
    vector<reference_wrapper<Behaviour>> actions = last_action->next_actions(vehicle, pred, *this);
    cout << "Costs: ";
    for(reference_wrapper<Behaviour> a: actions) {
      Behaviour& b = a.get();
      cout << "[" << b.name << "-" << b.cost(vehicle, pred, p) << "]";
    }
    auto min_action = *min_element(actions.begin(), actions.end(),
         [&vehicle, &pred, &p](const std::reference_wrapper<Behaviour> a, const std::reference_wrapper<Behaviour> b) {
            return a.get().cost(vehicle, pred, p) < b.get().cost(vehicle, pred, p);
          } );
    last_action = &min_action.get();
    cout << " ["<< last_action->name << "]";
    Traj traj =  last_action->convolute(vehicle, pred, p);
    return traj;
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

  Behaviour& slowDown(double speed) {
    this->sd.set_speed(speed);
    return this->sd;
  }

};

vector<std::reference_wrapper<Behaviour>> KeepLane::next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) {
  vector<std::reference_wrapper<Behaviour>> behaviours = { bp.keepLane() };
  int lane = vehicle.lane();
  if (lane > 0) {
    behaviours.push_back(bp.leftChange(vehicle, lane - 1));
  }
  if (lane < 2) {
    behaviours.push_back(bp.rightChange(vehicle, lane + 1));
  }
  if (vehicle.v > speed_limit/2.0) {
    behaviours.push_back(bp.slowDown(vehicle.v * 0.9));
  }
  return behaviours;
}

vector<std::reference_wrapper<Behaviour>> LaneChange::next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) {
  return { bp.keepLane(), bp.slowDown(vehicle.v * 0.9) };
}

vector<std::reference_wrapper<Behaviour>> SlowDown::next_actions(Vehicle& vehicle, Prediction& pred, BehaviourPlanner& bp) {
  return { bp.keepLane(), bp.slowDown(vehicle.v * 0.9) };
}


class TrajectoryGenerator {
 private:
  std::deque<vector<double>> points;
  BehaviourPlanner& bp;
  double step_duration;

 public:
  TrajectoryGenerator(BehaviourPlanner& bp, double step_duration = 0.02): bp(bp), step_duration(step_duration) {
  }

  void refresh_trajectory(Vehicle& vehicle, Prediction& pred, Planner& p) {
    cout << "Refreshing "<< vehicle.s << " at "<< vehicle.iter ;
    Traj traj = bp.convolute(vehicle, pred, p);
    cout<<  " with ";
    for(auto s: traj.s_vec) cout << s << ", ";
    cout<<  " and ";
    for(auto d: traj.d_vec) cout << d << ", ";
    cout << endl;
    auto pts = traj.motion_vector(this->step_duration, p);
    for (auto p: pts) this->points.push_back(p);
  }

  void effect(Vehicle& vehicle, int nsteps, Prediction& pred, Planner& p) {
    for (int i = 0; i < nsteps; ++i) {
      if (points.size() <= 0) {
        this->refresh_trajectory(vehicle, pred, p);
//        if (vehicle.iter > 800) {
//          p.plot(5, 25, 2);
//          vehicle.plot();
//          pred.plot(p, 100);
//          this->plot(p, 100);
//          plt::show();
//          exit(0);
//        }
      }
      auto pt = points.front();
      points.pop_front();
      vehicle.move_to(pt[0], pt[1], pt[2], pt[3], 0.02, p);
    }
  }


  void plot(Planner& p, int s = 10000) {
    cout << "Points in gen - "<< points.size() << endl;
    if (s>this->points.size()) s = points.size();
    for (int i = 0; i < s; ++i) {
      auto pt = this->points[i];
      auto xy = p.getXY(pt[0], pt[1]);
      plt::plot({xy[0]}, {xy[1]}, "g.");
    }
  }


};


#endif /* PLANNER_H_ */
