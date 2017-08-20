#include <iostream>
#include <fstream>
#include <sstream>


#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;


vector<vector<double>> get_highway_wp(std::string filepath) {
      std::ifstream file(filepath);
      vector<vector<double>> wp(5);
      while (file) {
        std::string line;
        std::getline(file, line);
        std::stringstream line_stream(line);
        if (!line.empty()) {
          for (int i = 0; i < 5; ++i) {
            std::string data;
            std::getline(line_stream, data, ' ');
            //cout << " >>" << data << "<< ";
            wp[i].push_back(std::stod(data));
          }
        }
      }
      file.close();
      return wp;
}

void plot_highway_wp() {
  auto wp = get_highway_wp("../data/highway_map.csv");
  vector<double> xs = wp[0], ys = wp[1], dxs = wp[3], dys = wp[4];
  int width = 50;
  plt::plot(xs, ys);
  for (int i = 0; i < dxs.size(); ++i) {
    plt::plot({xs[i], xs[i]+width*dxs[i]}, {ys[i], ys[i]+width*dys[i]}, "g");
    plt::plot({xs[i]+width*dxs[i]}, {ys[i]+width*dys[i]}, "g.");
  }
  plt::show();
}
