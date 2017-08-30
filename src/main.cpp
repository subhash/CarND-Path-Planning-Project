#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "planner.h"
#include "highway.h"
#include "json.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

// for convenience
using json = nlohmann::json;


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

void plot_debug1(Planner& planner, BehaviourPlanner& bp, TrajectoryGenerator& traj_gen, Environment& env, Highway highway) {
  Vehicle vehicle(planner);
  double vehicle_s = 3900, vehicle_d = 6;
  //double vehicle_s = 124.834, vehicle_d = 6.16483;
  vector<double> xy = planner.getXY(vehicle_s, vehicle_d);
  double vehicle_x = xy[0], vehicle_y = xy[1];
  cout << "x - "<< vehicle_x << ", y - "<< vehicle_y << endl;

  vehicle.init(vehicle_x, vehicle_y, vehicle_s, vehicle_d,  0, speed_limit);
  int nsteps = 3000;
  for(int i=0; i< nsteps; i++)
    traj_gen.effect(vehicle, 1, env, planner);

  plt::subplot(2,2,1);
  Car car(vehicle_x, vehicle_y, vehicle_s, vehicle_d,  0, speed_limit);
  highway.plot_highway(car);
  plt::plot(vehicle.xs, vehicle.ys, "r.");
  plt::subplot(2,2,2);
  //vehicle.vs.erase(vehicle.vs.begin(), vehicle.vs.begin()+200);
  plt::plot(vehicle.vs);
  plt::plot({0.0, 1000.0}, {49.3*speed_conv, 49.3*speed_conv}, "g-");
  plt::subplot(2,2,3);
  //vehicle.as.erase(vehicle.as.begin(), vehicle.as.begin()+2);
  plt::plot(vehicle.aerr);
  //plt::plot({0.0, 1000.0}, {10.0, 10.0}, "g-");
  plt::subplot(2,2,4);
  //vehicle.ss.erase(vehicle.ss.begin(), vehicle.ss.begin()+200);
  //plt::plot(vector<double>(vehicle.vs.begin()+1000, vehicle.vs.end()));
  plt::plot({0.0, 1000.0}, {speed_limit, speed_limit}, "g-");

  plt::show();
  exit(0);

}

//void plot_debug2(Planner& planner, BehaviourPlanner& bp, TrajectoryGenerator& traj_gen, Environment& env, Highway highway) {
//
//  vector<double> thetas, errs1, errs2;
//  double ds = 60, d = 6;
//  for(double s=1200; s<1400.0; s=s+60.0) {
//    thetas.push_back(planner.theta(s));
//    double err = planner.curve_error(s, s+ds, d);
//    double projected_s = planner.projected_s(s, s+ds-err, d);
//    errs1.push_back(err);
//    errs2.push_back(s+ds-projected_s);
//  }
//  plt::subplot(2,2,1);
//  plt::plot(thetas);
//  plt::subplot(2,2,2);
//  highway.plot_highway();
//  plt::subplot(2,2,3);
//  plt::plot(errs1);
//  plt::subplot(2,2,4);
//  plt::plot(errs2);
//
//  plt::show();
//  exit(0);
//
//}

void plot_debug3(Planner& planner, BehaviourPlanner& bp, TrajectoryGenerator& traj_gen, Environment& env, Highway highway) {
  vector<double> xpts, ypts, xpts1, ypts1;
  for(double s=4500; s<5000; s+=0.5) {
    auto xy = planner.getXY(s, 0);
    xpts.push_back(xy[0]);
    ypts.push_back(xy[1]);
    xy = planner.getXY(s, 6);
    xpts1.push_back(xy[0]);
    ypts1.push_back(xy[1]);
  }

  highway.plot_highway(120,140);

  plt::plot(xpts, ypts, "r");
  plt::plot(xpts1, ypts1, "g");
  plt::show();
  exit(0);
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  Planner planner(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);

  BehaviourPlanner bp;

  TrajectoryGenerator traj_gen(bp);

  Environment env;

  Highway highway = Highway(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);

  //plot_debug1(planner, bp, traj_gen, env, highway);

  Vehicle vehicle(planner);



  int step = 0;
  h.onMessage([&vehicle,&bp,&env,&traj_gen,&planner,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	for (int i = 0; i < sensor_fusion.size(); ++i) {
          	  auto data = sensor_fusion[i];
              env.update(data[0], data[1], data[2], data[3], data[4], data[5], data[6]);
            }

          	//cout << "Speed - " << env.lane_speed(0, car_s) << ", " << env.lane_speed(1, car_s) << ", "<< env.lane_speed(2, car_s) << endl;


            if (!vehicle.initialized)
              vehicle.init(car_x, car_y, car_s, car_d, car_yaw, car_speed * speed_conv);

//            if (!vehicle.initialized)
//              vehicle.init(car_x, car_y, 3900, car_d, car_yaw, car_speed * speed_conv);

            int nsteps = 50 - previous_path_x.size();
            //cout << "nsteps - "<< nsteps << endl;
            //Behaviour& b = bp.behaviour(vehicle, env);
            //b.effect(vehicle, env, nsteps);

            vehicle.trim(nsteps);
            traj_gen.effect(vehicle, nsteps, env, planner);


            for (int i = 1; i < vehicle.xs.size(); ++i) {
              next_x_vals.push_back(vehicle.xs[i]);
              next_y_vals.push_back(vehicle.ys[i]);
            }


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);


          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































