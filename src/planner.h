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

const double speed_conv = 0.44703, lane_width = 4.0, speed_limit = 49.0 * speed_conv;


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
  std::vector<tk::spline> lane_spline_x, lane_spline_y, lane_spline_dx, lane_spline_dy;
  std::vector<tk::spline> lane_spline_s, lane_spline_rev_s;


  Planner(vector<double> xs, vector<double> ys, vector<double> ss, vector<double> dxs, vector<double> dys):
    lane_spline_x(3), lane_spline_y(3), lane_spline_dx(3), lane_spline_dy(3), lane_spline_s(3), lane_spline_rev_s(3) {
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

    for (int i = 0; i < this->size; ++i) {
      this->thetas.push_back(atan2(this->dys[i], this->dxs[i]));
    }
    WP_spline_theta.set_points(ss, this->thetas);

    for (int i=0; i < 3; i++) {
      this->setup_lane_splines(i);
    }
  }


  void setup_lane_splines(int lane) {
    vector<double> ss1, xs1, ys1, dxs1, dys1;
    double lane_offset = lane*lane_width + lane_width/2.0;
    double err = 0.0;

//    ss1.push_back(this->ss[0]);
//    xs1.push_back(WP_spline_x(this->ss[0]));
//    ys1.push_back(WP_spline_y(this->ss[0]));
//    dxs1.push_back(WP_spline_dx(this->ss[0]));
//    dys1.push_back(WP_spline_dy(this->ss[0]));

    for (int i = 0; i < this->size; ++i) {
      int prev_pt = i-1;
      if (prev_pt < 0) {
        prev_pt = this->size-1;
      }
      double prev_s = this->ss[prev_pt], s = this->ss[i];
      double theta1 = WP_spline_theta(prev_s), theta2 = WP_spline_theta(s);
      err += lane_offset*tan(theta2 - theta1);
      double proj_s = s + err;
      ss1.push_back(proj_s);
      xs1.push_back(WP_spline_x(s));
      ys1.push_back(WP_spline_y(s));
      dxs1.push_back(WP_spline_dx(s));
      dys1.push_back(WP_spline_dy(s));
    }

    lane_spline_x[lane].set_points(ss1, xs1);
    lane_spline_y[lane].set_points(ss1, ys1);
    lane_spline_dx[lane].set_points(ss1, dxs1);
    lane_spline_dy[lane].set_points(ss1, dys1);

    lane_spline_s[lane].set_points(ss1, this->ss);
    lane_spline_rev_s[lane].set_points(this->ss, ss1);
  }


//    vector<double> getXY(double s, double d) {
//      int lane = int(d/lane_width);
//      double wp_s = lane_spline_s[lane](s);
//      double x = WP_spline_x(wp_s), y = WP_spline_y(wp_s);
//      double dx = WP_spline_dx(wp_s) * d, dy = WP_spline_dy(wp_s) * d;
//      return { x+dx, y+dy };
//    }

//  vector<double> getXY(double s, double d) {
//    int lane = int(d/lane_width);
//    double x = lane_spline_x[lane](s), y = lane_spline_y[lane](s);
//    double dx = lane_spline_dx[lane](s) * d, dy = lane_spline_dy[lane](s) * d;
//    return { x+dx, y+dy };
//  }

  vector<double> getXY(double s, double d) {
    double x = WP_spline_x(s), y = WP_spline_y(s);
    double dx = WP_spline_dx(s) * d, dy = WP_spline_dy(s) * d;
    return { x+dx, y+dy };
  }


//
//  double theta(double s) {
//    double dx = WP_spline_dx(s), dy = WP_spline_dy(s);
//    return atan2(dy, dx);
//  }
//
//  double curve_error(double s1, double s2, double d) {
//    double theta1 = this->theta(s1);
//    double theta2 = this->theta(s2);
//    return d*tan(theta2 - theta1);
//  }
//
//  double projected_s(double s1, double s2, double d) {
//    double theta1 = this->theta(s1);
//    double theta2 = this->theta(s2);
//    return s2+d*tan(theta2 - theta1);
//  }

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
    //s_vec[0] = p.lane_spline_rev_s[ref_lane](s_vec[0]);
    //s_vec[3] = p.lane_spline_rev_s[ref_lane](s_vec[3]);
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
      double proj_s = p.lane_spline_s[ref_lane](s);
      pts.push_back(make_pair(proj_s, d));
    }
    return pts;
  }

  vector<vector<double>> motion_vector(double dt, Planner& p) {
    int steps = this->duration/dt;
    vector<vector<double>> pts;
    for (int i = 0; i < steps; ++i) {
      double s = poly_eval(this->spoly, (i+1) * dt);
      double d = poly_eval(this->dpoly, (i+1) * dt);
      double proj_s = p.lane_spline_s[ref_lane](s);
      double v = v_eval(this->spoly, (i+1)*dt);
      double a = a_eval(this->spoly, (i+1)*dt);
      pts.push_back( { proj_s, d, v, a });
    }
    return pts;
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
  const double max_s = 6945.554;

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

    if (proj_s > max_s) proj_s -= 6945.554;
    vector<double> xy = planner.getXY(proj_s, d);
    double s = proj_s; //planner.lane_spline_rev_s[this->lane()](proj_s);

//    double corr_v, corr_a;
//    tie(corr_v, corr_a) = velocity(xy[0], xy[1], dt);
//    if (corr_v > 50*speed_conv) {
//      //this->verr_xs.push_back(xy[0]);
//      //this->verr_ys.push_back(xy[1]);
//      //this->verr.push_back(corr_v-speed_limit);
//      // cout << this->s <<": > correct velocity "<< corr_v/speed_conv << " at iter " << this->iter << endl;
//    }

//    double v = (s - this->s)/dt;
//    double a = (v - this->v)/dt;
//    double j = (a - this->a)/dt;
    double x = xy[0], y = xy[1];

//    if (corr_a > 10 && this->iter > 1) {
//      this->aerr_xs.push_back(xy[0]);
//      this->aerr_ys.push_back(xy[1]);
//      this->aerr.push_back(corr_a);
//      // cout << this->s <<": > correct acc "<< corr_a << " and "<< corr_v << " instead of " << a <<" at iter " << this->iter << endl;
//      // cout << "previous - " << this->s << ", " << this->ss[this->ss.size() -1] << ", "<< this->ss[this->ss.size() -2] <<", "<< this->ss[this->ss.size() -3] << endl;
//    }

    this->s = s;
    this->d = d;
    this->x = x;
    this->y = y;
    this->v = v;
    this->a = a;
    //if(this->iter > 1) this->a = a;
    //if(this->iter > 2) this->j = j;

    this->xs.push_back(x);
    this->ys.push_back(y);
    this->ss.push_back(s);
    this->vs.push_back(v);
    //if(this->iter > 1) this->as.push_back(a);
    //if(this->iter > 2) this->js.push_back(j);
  }

//
//
//  double inc_velocity(double time, double init_v, double init_a, double v_desired, double acc_limit, double jerk_limit) {
//    double a = init_a;
//    double a_diff = jerk_limit*time;
//    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
//    double v = init_v;
//    double v_diff = a*time;
//    if (v + v_diff < v_desired) v += v_diff;
//
//    //cout << "inc v - "<< this->iter<< "  " << v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
//    return v;
//  }
//
//  double dec_velocity(double time, double init_v, double init_a, double v_desired, double acc_limit, double jerk_limit) {
//    acc_limit *= -1;
//    jerk_limit *= -1;
//
//    double a = init_a;
//    double a_diff = jerk_limit*time;
//    if (fabs(a + a_diff) < fabs(acc_limit)) a += a_diff;
//    double v = init_v;
//    double v_diff = a*time;
//    //v_diff = -0.1;
//    if (v + v_diff > v_desired) v += v_diff;
//
//    //cout << "dec v - "<< v << " adiff - "<< a_diff << " a - "<< a << " bool - "<< (fabs(a + a_diff) < fabs(acc_limit)) << " vdiff " << v_diff << " speed_limit " << speed_limit <<  endl;
//    return v;
//  }
//
//
//
//  void move(int steps, double dest_d, double time, double speed_limit, double acc_limit, double jerk_limit, vector<double> d_poly = {}, int nsteps = 0) {
//    this->trim(steps);
//    for (int i = 0; i < steps; ++i) {
//      double v;
//      if (this->v < speed_limit)
//        v = this->inc_velocity(time, this->v, this->a, speed_limit, acc_limit, jerk_limit);
//      else
//        v = this->dec_velocity(time, this->v, this->a, speed_limit, acc_limit, jerk_limit);
//      double d_diff = 0;
//      if (!d_poly.empty()) {
//        double tf = time*(i+nsteps);
//        for (int c = 0; c < d_poly.size(); ++c) {
//          d_diff += d_poly[c]*pow(tf, c);
//        }
//      }
//      this->step(time, dest_d, v, d_diff, this->planner);
//    }
//  }
//
//  void step(double time, double d, double v, double d_diff, Planner planner) {
//    this->iter++;
//
//    double s_diff = v * time;
//    // Just to support turns
//    if (d_diff > 0) this->d = d_diff;
//    vector<double> xy = planner.getXY(this->s + s_diff, this->d);
//
////    double corr_v = velocity(xy[0], xy[1], time);
////    if (corr_v > speed_limit) {
////      this->verr_xs.push_back(xy[0]);
////      this->verr_ys.push_back(xy[1]);
////      this->verr.push_back(corr_v-speed_limit);
////    }
//
//    double a = (v - this->v)/time;
//    double j = (a - this->a)/time;
//    double x = xy[0], y = xy[1];
//
//    this->s += s_diff;
//    this->x = x;
//    this->y = y;
//    this->v = v;
//    if(this->iter > 1) this->a = a;
//    if(this->iter > 2) this->j = j;
//
//    this->xs.push_back(x);
//    this->ys.push_back(y);
//    this->ss.push_back(s);
//    this->vs.push_back(v);
//    if(this->iter > 1) this->as.push_back(a);
//    if(this->iter > 2) this->js.push_back(j);
//  }
//
//  void print_stats() {
//    int min_acc = min_element(this->as.begin()+2, this->as.end()) - this->as.begin();
//    int max_acc = max_element(this->as.begin()+2, this->as.end()) - this->as.begin();
//    int min_jerk = min_element(this->js.begin()+3, this->js.end()) - this->js.begin();
//    int max_jerk = max_element(this->js.begin()+3, this->js.end()) - this->js.begin();
//    cout << "Acc: " << this->as[min_acc] << "["<< min_acc << "], " << this->as[max_acc] << "["<< max_acc << "]" << endl;
//    cout << "Jrk: " << this->js[min_jerk] << "["<< min_jerk << "], " << this->js[max_jerk] << "["<< max_jerk << "]" << endl;
//  }
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

  int lane() {
    return int(d/lane_width);
  }

  double speed() {
    return sqrt(pow(this->vx, 2) + pow(this->vy, 2));
  }

};


class Environment {

 private:
  const double time = 0.02, acc_limit = 2, jerk_limit = 2;

 public:


//  int closest_obstacle(int lane, double s) {
//    double min_dist = std::numeric_limits<double>::max();
//    int min_id = -1;
//    for (auto const& el : obstacles) {
//      Obstacle o = el.second;
//      double dist = fabs(o.s - s);
//      if (o.lane() == lane && dist < min_dist) {
//        min_dist = dist;
//        min_id = el.first;
//      }
//    }
//    return min_id;
//  }
//
//  bool will_vehicle_collide(int lane, Vehicle& vehicle) {
//    double band = vehicle.v * 2, min_s = vehicle.s - band, max_s = vehicle.s + band;
//    for (auto const& el : obstacles) {
//      Obstacle o = el.second;
//      if (o.lane() == lane && o.s > min_s and o.s < max_s) {
//        return true;
//      }
//    }
//    return false;
//  }
//
//  bool will_obstacle_collide(int lane, Vehicle& vehicle) {
//    double band = vehicle.v * 2, min_s = vehicle.s, max_s = vehicle.s + band;
//    for (auto const& el : obstacles) {
//      Obstacle o = el.second;
//      if (o.lane() == lane && o.s > min_s and o.s < max_s and fabs(o.d - vehicle.d) < 0.9*lane_width) {
//        return true;
//      }
//    }
//    return false;
//  }
//
//  double leading_distance(int lane, double s) {
//    int id = leading_obstacle(lane, s);
//    if (id == -1) {
//      return std::numeric_limits<double>::max();
//    } else {
//      return obstacles[id].s - s;
//    }
//  }


  map<int, Obstacle> obstacles;

  void update(int id, double x, double y, double vx, double vy, double s, double d, double ts) {
    obstacles[id].update(id, x, y, vx, vy, s, d, ts);
  }

  Environment predict(double dt) {
    Environment predicted;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      predicted.update(o.id, o.x+o.vx*dt, o.y+o.vy*dt, o.vx, o.vy, o.s+o.speed()*dt, o.d+o.dv*dt, dt);
    }
    return predicted;
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
      plt::plot(o.ds);
      plt::plot({200 + 3/0.02}, {predicted.obstacles[o.id].d}, "r.");
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.svs);
    }
    plt::show();

    i = 0;
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      plt::subplot(rows,3,++i);
      plt::plot(o.tss);
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
      if (gap < 40) {
        return obstacles[id].speed();
      }
    }
    return speed_limit;
  }

  bool too_close(double s, double d, double s_tol, double d_tol) {
    for (auto const& el : obstacles) {
      Obstacle o = el.second;
      if (fabs(o.s - s) < s_tol && fabs(o.d - d) < d_tol)
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
//    cout << "Scanning from "<< from_s << ", " << min_d << " to "<< to_s << ", " << max_d << endl;
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
  double time = 0.02, speed_limit = 49.0 * speed_conv, acc_limit = 2, jerk_limit = 2;

  double decide_speed(Vehicle& vehicle, Environment& env, double speed_limit) {
    return env.lane_speed(vehicle, vehicle.lane());
  }

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
//  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) { }
//  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env) { return { *this }; }
//  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) = 0;
// virtual Behaviour& next_behaviour(Vehicle& vehicle, Environment& env, Planner& planner);
  virtual Traj convolute(Vehicle& vehicle, Environment& env, Planner& p) {
    cout << "Base convolute "<< endl;
    return Traj(0, vector<double>(), vector<double>(), 0, p);
  }
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) {
    return {};
  }

  virtual double safety_cost(Vehicle& vehicle, Environment& env, Traj& traj, Planner& p) {
    vector<pair<double, double>> pts = traj.points(traj.duration, p);
    double dt = traj.duration/pts.size();
    double cost = 0;
    for (int i = 0; i < pts.size(); ++i) {
      double traj_s = pts[i].first, traj_d = pts[i].second;
      Environment predicted = env.predict(dt*(i+1));
      if (predicted.too_close(traj_s, traj_d, 30, 2.5))
        cost += 1.0;
    }
    return cost;
    //return predicted.collision_space(traj.s_vec[0], traj.s_vec[3], traj.d_vec[0], traj.d_vec[3]);
  }

  virtual double efficiency_cost(Vehicle& vehicle, Environment& env, Traj& traj) {
    int dest_lane = int(traj.d_vec[3]/lane_width);
    double lane_speed = env.lane_speed(vehicle, dest_lane);
    return lane_speed < speed_limit;
  }

  virtual double cost(Vehicle& vehicle, Environment& env, Planner& p) {
    Traj traj = this->convolute(vehicle, env, p);
    return 10*this->safety_cost(vehicle, env, traj, p) + this->efficiency_cost(vehicle, env, traj);
  }

  void plot(Vehicle& vehicle, Environment& env, Traj& traj, Planner& p) {
    plt::plot(vehicle.xs, vehicle.ys, "r");
    plt::plot({vehicle.x}, {vehicle.y}, "r+");
    Obstacle& o = env.obstacles[1];
    //plt::plot(o.xs, o.ys, "b.");
    Obstacle& po = env.predict(3.0).obstacles[1];
    vector<double> opt = p.getXY(po.s, po.d);
    //plt::plot({opt[0]}, {opt[1]}, "b+");
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
//  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
//
//
//  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
//    if(debug) cout << this->name << ": ";
//    int lane = vehicle.lane(), llane = lane-1, rlane = lane+1;
//    double s = vehicle.s, d = vehicle.d;
//    double speed_cost = env.lane_speed(vehicle, vehicle.lane()) < speed_limit;
//    double safety_cost = (llane >= 0) && env.will_obstacle_collide(llane, vehicle) +
//        (rlane <= 2) && env.will_obstacle_collide(rlane, vehicle);
//    if (debug) cout << "speed_cost=" << speed_cost << ", safety_cost=" << safety_cost;
//    if(debug) cout << ", llane=" << env.will_obstacle_collide(llane, vehicle) << ", rlane=" << env.will_obstacle_collide(rlane, vehicle);
//    return 2*speed_cost + 10*safety_cost;
//  }


  virtual Traj convolute(Vehicle& vehicle, Environment& env, Planner& p) override {
    double tf = 3.0; // sec
    double sf_dot = env.lane_speed(vehicle, vehicle.lane());

    double si = vehicle.s;
    si = p.lane_spline_rev_s[vehicle.lane()](si);

    vector<double> s_vec = this->long_vector(si, vehicle.v, vehicle.a, sf_dot, tf);
    vector<double> d_vec = { vehicle.d, 0, 0, vehicle.d, 0, 0, 0};
    Traj traj(tf, s_vec, d_vec, vehicle.lane(), p);
    return traj;
  }

  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) override;
  //virtual Behaviour& next_behaviour(Vehicle& vehicle, Environment& env, Planner& planner) override;
};

class LaneChange: public Behaviour {
  int dest_lane = 0;
  vector<double> d_poly;
  int steps = 0;
 public:
  LaneChange(string name, int lane): Behaviour(name), dest_lane(lane) {}

  virtual Traj convolute(Vehicle& vehicle, Environment& env, Planner& p) override {
    double tf = 7.0; // sec
    double sf_dot = env.lane_speed(vehicle, dest_lane);

    double si = vehicle.s;
    si = p.lane_spline_rev_s[vehicle.lane()](si);

    vector<double> s_vec = this->long_vector(si, vehicle.v, vehicle.a, sf_dot, tf);
    double df = dest_lane * lane_width + lane_width/2.0;
    vector<double> d_vec = { vehicle.d, 0, 0, df, 0, 0, 0};
    Traj traj(tf, s_vec, d_vec, vehicle.lane(), p);
    return traj;
  }

  void set_lane(Vehicle& vehicle, int dest_lane) {
    this->dest_lane = dest_lane;
  }

  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) override;

//  virtual Behaviour& next_behaviour(Vehicle& vehicle, Environment& env, Planner& planner) override;
//  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
//
//
//  vector<double> make_poly(double tf, double d_init, double d_final) {
//    double t = tf, t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t;
//    MatrixXd Tmat(3, 3);
//    Tmat << t3, t4, t5,
//            3*t2, 4*t3, 5*t4,
//            6*t, 12*t2, 20*t3;
//    MatrixXd Tinv = Tmat.inverse();
//    VectorXd Dmat(3);
//    Dmat << d_final - d_init, 0, 0;
//    VectorXd Dvec = Tinv*Dmat;
//    vector<double> d_poly = { d_init, 0, 0, Dvec(0), Dvec(1), Dvec(2) };
//    return d_poly;
//  }
//
//  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
//    if(debug) cout << this->name << ": ";
//    double s = vehicle.s, d = vehicle.d;
//    double speed_cost = env.lane_speed(vehicle, dest_lane) < speed_limit;
//    double safety_cost = env.will_vehicle_collide(dest_lane, vehicle);
//    if (debug) cout << "speed_cost=" << speed_cost << ", safety_cost=" << safety_cost;
//    return 1+2*speed_cost + 10*safety_cost;
//  }
};

class SlowDown: public Behaviour {
 public:
  SlowDown(): Behaviour("SD") {};
  virtual vector<std::reference_wrapper<Behaviour>> next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) override;
//  virtual void effect(Vehicle& vehicle, Environment& env, int nsteps) override;
//  virtual double cost(Vehicle& vehicle, Environment& env, bool debug = false) override {
//    double speed_cost = 1;
//    if (debug) cout << "speed_cost=" << speed_cost;
//    return speed_cost;
//  }
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

  Traj convolute(Vehicle& vehicle, Environment& env, Planner& p) {
    vector<reference_wrapper<Behaviour>> actions = last_action->next_actions(vehicle, env, *this);
    cout << "Costs: ";
    for(reference_wrapper<Behaviour> a: actions) {
      Behaviour& b = a.get();
      cout << "[" << b.name << "-" << b.cost(vehicle, env, p) << "]";
    }
    auto min_action = *min_element(actions.begin(), actions.end(),
         [&vehicle, &env, &p](const std::reference_wrapper<Behaviour> a, const std::reference_wrapper<Behaviour> b) {
            return a.get().cost(vehicle, env, p) < b.get().cost(vehicle, env, p);
          } );
    last_action = &min_action.get();
    cout << " ["<< last_action->name << "]";
    Traj traj =  last_action->convolute(vehicle, env, p);
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

  Behaviour& slowDown() {
    return this->sd;
  }
//
//  Behaviour& next_action(Vehicle& vehicle, Environment& env) {
//  }
};

vector<std::reference_wrapper<Behaviour>> KeepLane::next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) {
  vector<std::reference_wrapper<Behaviour>> behaviours = { bp.keepLane() };
  int lane = vehicle.lane();
  if (lane > 0) {
    behaviours.push_back(bp.leftChange(vehicle, lane - 1));
  }
  if (lane < 2) {
    behaviours.push_back(bp.rightChange(vehicle, lane + 1));
  }
  return behaviours;
}

vector<std::reference_wrapper<Behaviour>> LaneChange::next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) {
  return { bp.keepLane() };
}

vector<std::reference_wrapper<Behaviour>> SlowDown::next_actions(Vehicle& vehicle, Environment& env, BehaviourPlanner& bp) {
  return { bp.keepLane() };
}


//void KeepLane::effect(Vehicle& vehicle, Environment& env, int nsteps) {
//  double dest_d = vehicle.lane()*lane_width + 2;
//  vehicle.move(nsteps, dest_d, time, decide_speed(vehicle, env, speed_limit), acc_limit, jerk_limit);
//}
//
//void LaneChange::effect(Vehicle& vehicle, Environment& env, int nsteps) {
//  double dest_d = dest_lane*lane_width + 2;
//  vehicle.move(nsteps, vehicle.d, time, decide_speed(vehicle, env, speed_limit), acc_limit, jerk_limit, d_poly, this->steps);
//  this->steps += nsteps;
//}
//
//void SlowDown::effect(Vehicle& vehicle, Environment& env, int nsteps) {
//  vehicle.move(nsteps, vehicle.d, time, 5, acc_limit, jerk_limit);
//}
//
//vector<std::reference_wrapper<Behaviour>> KeepLane::next_actions(Vehicle& vehicle, Environment& env) {
//  vector<reference_wrapper<Behaviour>> actions = { behaviour_planner.keepLane() };
////  if (vehicle.iter > 500) {
////    return { behaviour_planner.leftChange(vehicle, 0) };
////  }
//  if (vehicle.lane() > 0) {
//    Behaviour& lc = behaviour_planner.leftChange(vehicle, vehicle.lane() - 1);
//    actions.push_back(lc);
//  }
//  if (vehicle.lane() < 2) {
//    Behaviour& lc = behaviour_planner.rightChange(vehicle, vehicle.lane() + 1);
//    actions.push_back(lc);
//  }
//  return actions;
//}
//
//vector<std::reference_wrapper<Behaviour>> LaneChange::next_actions(Vehicle& vehicle, Environment& env) {
//  double dest_d = dest_lane*lane_width + 2;
//  if (this->steps < 3*50){
//    return { *this };
//  }
//  return { behaviour_planner.keepLane() };
//}


class TrajectoryGenerator {
 private:
  std::deque<vector<double>> points;
  BehaviourPlanner& bp;
  double step_duration;

 public:
  TrajectoryGenerator(BehaviourPlanner& bp, double step_duration = 0.02): bp(bp), step_duration(step_duration) {
  }

  void refresh_trajectory(Vehicle& vehicle, Environment& env, Planner& p) {
    cout << "Refreshing "<< vehicle.s << " at "<< vehicle.iter ;
    Traj traj = bp.convolute(vehicle, env, p);
    cout<<  " with ";
    for(auto s: traj.s_vec) cout << s << ", ";
    cout<<  " and ";
    for(auto d: traj.d_vec) cout << d << ", ";
    cout << endl;
    auto pts = traj.motion_vector(this->step_duration, p);
    for (auto p: pts) this->points.push_back(p);
  }

  void effect(Vehicle& vehicle, int nsteps, Environment& env, Planner& p) {
    for (int i = 0; i < nsteps; ++i) {
      if (points.size() <= 0) {
        this->refresh_trajectory(vehicle, env, p);
      }
      auto pt = points.front();
//      if (env.too_close(pt[0], pt[1], 10, 1)) {
//        cout << "Cancelling traj "<< endl;
//        for (int j = 0; j < points.size(); ++j) {
//          points.pop_back();
//        }
//        this->refresh_trajectory(vehicle, env, p);
//        pt = points.front();
//      }
      points.pop_front();
      vehicle.move_to(pt[0], pt[1], pt[2], pt[3], 0.02, p);
    }
  }


};


#endif /* PLANNER_H_ */
