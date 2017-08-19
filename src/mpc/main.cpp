#include <chrono>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#include <uWS/uWS.h>

#include "json/json.hpp"

#include "mpc.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// for convenience
using json = nlohmann::json;

constexpr const int kDelay = 100;

double ToSpeed(double v) { return v * 0.44703888888; }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string HasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  Mpc mpc(ToSpeed(60));

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char* data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata(data, length);

    if (!length || length <= 2 || data[0] != '4' || data[1] != '2') return;

    auto s = HasData(sdata);
    if (s.empty()) {
      string msg = "42[\"manual\",{}]";
      ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      return;
    }

    auto j = json::parse(s);
    string event = j[0].get<string>();
    if (event != "telemetry") return;

    // j[1] is the data JSON object
    vector<double> ptsx = j[1]["ptsx"];
    vector<double> ptsy = j[1]["ptsy"];
    double px = j[1]["x"];
    double py = j[1]["y"];
    double psi = j[1]["psi"];
    double v = j[1]["speed"];
    double delta = j[1]["steering_angle"];
    double a = j[1]["throttle"];

    if (ptsx.size() != ptsy.size()) return;
    for (size_t i = 0; i != ptsx.size(); ++i) {
      double x = ptsx[i] - px;
      double y = ptsy[i] - py;
      ptsx[i] = x * cos(-psi) - y * sin(-psi);
      ptsy[i] = x * sin(-psi) + y * cos(-psi);
    }

    double steer_value = delta;
    double throttle_value = a;
    vector<double> mpc_x_vals;
    vector<double> mpc_y_vals;

    auto mpc_vals = mpc.Solve(0.001 * kDelay, 0, 0, 0, ToSpeed(v), -delta, a, ptsx, ptsy);

    if (mpc_vals.size() != 0) {
      steer_value = -mpc_vals(1, 6);
      throttle_value = mpc_vals(1, 7);
      
      mpc_x_vals.resize(mpc_vals.rows());
      mpc_y_vals.resize(mpc_vals.rows());
      for (int i = 0; i != mpc_vals.rows(); ++i) {
        mpc_x_vals[i] = mpc_vals(i, 0);
        mpc_y_vals[i] = mpc_vals(i, 1);
      }
    }

    json msgJson;
    msgJson["steering_angle"] = steer_value;
    msgJson["throttle"] = throttle_value;

    msgJson["mpc_x"] = mpc_x_vals;
    msgJson["mpc_y"] = mpc_y_vals;

    msgJson["next_x"] = ptsx;
    msgJson["next_y"] = ptsy;

    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
    // Latency
    // The purpose is to mimic real driving conditions where
    // the car does actuate the commands instantly.
    //
    // Feel free to play around with this value but should be to drive
    // around the track with 100ms latency.
    //
    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
    // SUBMITTING.
    std::this_thread::sleep_for(std::chrono::milliseconds(kDelay));
    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
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

  h.onDisconnection(
      [&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message,
           size_t length) { std::cout << "Disconnected" << std::endl; });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
