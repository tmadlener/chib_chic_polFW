#ifndef UTILS_JN_H
#define UTILS_JN_H

///////////////////////////////////
// Some helper functions for C++ //
///////////////////////////////////
//
// TODO:
//   - make all functions threadsafe
//
////////////////////////////////////

#include <chrono>
#include <iostream>
#include <sstream>
#include <mutex>

#include "TObject.h"


// Single functions

bool file_exists(const std::string &fileName);



// Structs and classes

struct stopwatch {
  explicit stopwatch(std::string description, std::mutex *cout_mutex = nullptr) :
    d(description),
    mtx(cout_mutex),
    start(std::chrono::system_clock::now()) { }
  ~stopwatch() {
    std::stringstream ss;
    auto end = std::chrono::system_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //ss << duration_ms.count() << '\n';
    ss << "Time for '" << d << "': ";
    unsigned long long mins = duration_ms.count() / 1000. / 60.;
    unsigned long long secs = (duration_ms.count() - mins * 1000 * 60) / 1000;
    unsigned long long ms = (duration_ms.count() - secs * 1000 - mins * 1000 * 60);
    ss << mins << ":" << (secs < 10 ? "0" : "") << secs << ',' << ms;
    ss << " [minutes:seconds,milliseconds]";
    if (mtx) {
      std::lock_guard<std::mutex> lock(*mtx);
      std::cout << ss.rdbuf() << std::endl;
    }
    else std::cout << ss.rdbuf() << std::endl;
  }
private:
  std::string d;
  std::chrono::system_clock::time_point start;
  std::mutex *mtx;
};

struct scopeLog {
  explicit scopeLog(const std::string &description) :
    desc(description)
  {
    std::cout << "START of '" << desc << "'" << std::endl;
  }
  ~scopeLog() {
    std::cout << "END of '" << desc << "'" << std::endl;
  }
private:
  std::string desc;
};

#endif // UTILS_JN_H