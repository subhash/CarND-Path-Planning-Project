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

void plot_debug(Planner planner, Highway highway) {
  //double car_s = 124.834, car_d = 6.16483;
  //double car_s = 390, car_d = 6.16483;
  double car_s = 124.834, car_d = 6.16483;
  Car car = Car(909.48, 1128.67, car_s, car_d,  0, 0);

  double x = 0, y = 0, prev_dist = 0;
  vector<double> ddist, dvel, infl_x, infl_y;
  vector<double> xpts1, ypts1;

  Vehicle vehicle(909.48, 1128.67, car_s, car_d,  0, 20);
  for (int i=0; i<50*240; i++) {
    vehicle.step(21, 0.02, planner);
  }
  plt::subplot(4,1,1);
  highway.plot_highway(car);
  plt::plot(vehicle.xs, vehicle.ys, "r.");
  plt::subplot(4,1,2);
  plt::plot(vehicle.vs);
  plt::subplot(4,1,3);
  vehicle.as.erase(vehicle.as.begin(), vehicle.as.begin()+2);
  plt::plot(vehicle.as);
  plt::subplot(4,1,4);
  vehicle.js.erase(vehicle.js.begin(), vehicle.js.begin()+3);
  plt::plot(vehicle.js);


    cout << max_element(vehicle.as.begin(), vehicle.as.end()) - vehicle.as.begin() << ", "<< min_element(vehicle.as.begin(), vehicle.as.end()) - vehicle.as.begin() << endl;
    cout << max_element(vehicle.js.begin(), vehicle.js.end()) - vehicle.js.begin() << ", "<< min_element(vehicle.js.begin(), vehicle.js.end()) - vehicle.js.begin() << endl;
    cout << *max_element(vehicle.js.begin(), vehicle.js.end()) << ", "<< *min_element(vehicle.js.begin(), vehicle.js.end()) << endl;

  cout << "violated speed at " << vehicle.verr_xs.size() << endl;

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

  Highway highway = Highway(map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy);

  //plot_debug(planner, highway);
  //car_debug(planner, highway);

  int step = 0;
  h.onMessage([&step,&highway,&planner,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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


//          	// Straight path
//          	double inc = 0.25;
//          	for (int i = 0; i < 50; ++i) {
//              next_x_vals.push_back(car_x + i*inc*cos(car_yaw));
//              next_y_vals.push_back(car_y + i*inc*sin(car_yaw));
//            }

//            // Circular path
//            double inc = 0.5, angle = deg2rad(car_yaw);
//            for (int i = 0; i < 50; ++i) {
//              car_x += inc*cos(angle + (i+1)*deg2rad(0.2));
//              car_y += inc*sin(angle + (i+1)*deg2rad(0.2));
//              next_x_vals.push_back(car_x);
//              next_y_vals.push_back(car_y);
//            }

//          	double speed_limit = 49;
//          	if (car.speed < speed_limit) {
//          	  car.speed += 5;
//          	}
//            if (car.speed > speed_limit) {
//              car.speed -= 5;
//            }
//

//            double time = 0.02, speed_limit = 30, acc_limit = 10;
//            Car car = Car(car_x, car_y, car_s, car_d, car_yaw, speed_limit);
//            double d = 6;
//            for (int i = 0; i < 500; ++i) {
//              //car.adjust_speed(speed_limit, acc_limit, time);
//              car.step(time);
//              vector<double> xy = planner.getXY(car.s, d);
//              next_x_vals.push_back(xy[0]);
//              next_y_vals.push_back(xy[1]);
//            }


          	Vehicle vehicle(car_x, car_y, end_path_s, end_path_d, car_yaw, 20);
          	for (int i = 0; i < 2000-previous_path_x.size(); ++i) {
              vehicle.step(20, 0.02, planner);
            }
          	//vehicle.print_stats();

          	for (int i = 0; i < previous_path_x.size(); ++i) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

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
















































































