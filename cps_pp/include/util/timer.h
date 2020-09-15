// vim: set ts=2 sw=2 expandtab:
#ifndef INCLUDED_TIMER_H
#define INCLUDED_TIMER_H

#include <vector>
#include <string>
#include <algorithm>

#ifdef USE_PAPI
#include <papi.h>
#include <omp.h>
#endif

#include <util/time_cps.h>
#include <util/verbose.h>

CPS_START_NAMESPACE

inline double& getStartTime() {
  static double time = dclock();
  return time;
}

inline double getTotalTime() {
  return dclock() - getStartTime();
}

inline long long getTotalFlops() {
  long long flops = 0;
#ifdef USE_PAPI
  const int n_threads = omp_get_max_threads();
  long long flopses[n_threads];
  memset(flopses, 0, n_threads * sizeof(long long));
#pragma omp parallel
  {
    float rtime, ptime, mflops;
    int i = omp_get_thread_num();
    PAPI_flops(&rtime, &ptime, &flopses[i], &mflops);
  }
  for (int i = 0; i < n_threads; i++) {
    flops += flopses[i];
  }
#endif
  return flops;
}

inline void initializePAPI() {
#ifdef USE_PAPI
  static bool initialized = false;
  if (initialized) {
    return;
  }
  VRB.Result("PAPI", "initializePAPI", "Start.");
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
  initialized = true;
  VRB.Result("PAPI", "initializePAPI", "Finish.");
#endif
}

class TimerInfo {

  public :

    std::string fname;
    double dtime;
    double accumulated_time;
    long long dflops;
    long long accumulated_flops;
    int call_times;

    TimerInfo() {
      fname = "Unknown";
      dtime = 0.0 / 0.0;
      accumulated_time = 0;
      dflops = 0;
      accumulated_flops = 0;
      call_times = 0;
    }

    void showLast(const char *info = NULL) {
      double total_time = getTotalTime();
      std::string fnameCut;
      fnameCut.assign(fname, 0, 60);
      VRB.Result("Timer", NULL == info ? "" : info,
          "%60s :%5.1f%%%9d calls. Last %.3E secs%8.3f Gflops (%.3E per call)\n",
          fnameCut.c_str(),
          accumulated_time / total_time * 100, call_times,
          dtime,
          dflops / dtime / 1.0E9,
          (double)dflops);
    }

    void showAvg(const char *info = NULL) {
      double total_time = getTotalTime();
      std::string fnameCut;
      fnameCut.assign(fname, 0, 60);
      VRB.Result("Timer", NULL == info ? "" : info,
          "%60s :%5.1f%%%9d calls. Avg %.3E secs%8.3f Gflops (%.3E per call)\n",
          fnameCut.c_str(),
          accumulated_time / total_time * 100, call_times,
          accumulated_time / call_times,
          accumulated_flops / accumulated_time / 1.0E9,
          (double)accumulated_flops / (double)call_times);
    }

    void show(const char *info = NULL) {
      showAvg(info);
    }

    void clear()
    {
	dtime = 0;
	dflops = 0;
	accumulated_time = 0;
	accumulated_flops = 0;
	call_times = 0;
    }
};

inline bool compareTimeInfoP(TimerInfo *p1, TimerInfo *p2) {
  return p1->accumulated_time < p2->accumulated_time;
}


class Timer {
 protected:
  static Float dtime_begin;
  static Float dtime_last;
  static double autodisplay_interval;
  static double show_stop;

  public :

    const char *cname;

    int index;
    TimerInfo info;

    static std::vector<TimerInfo *>& getTimerDatabase() {
      static std::vector<TimerInfo *> timerDatabase;
      return timerDatabase;
    }

    double start_time;
    double stop_time;

    long long start_flops;
    long long stop_flops;

    static double& set_autodisplay_interval(Float _time) {
      autodisplay_interval=_time;
      return autodisplay_interval;
    }

    static double& minimum_autodisplay_interval() {
      return autodisplay_interval;
    }

    static double& set_duration_for_show_stop_info( Float _time) {
      show_stop = _time ;
      return show_stop;
    }
    static double& minimum_duration_for_show_stop_info() {
      return show_stop;
    }

    static double& minimum_duration_for_show_start_info() {
//      static double time = 0.1;
//      return time;
      return show_stop;
    }

    Timer() {
      init();
    }

    Timer(const std::string& fname_str) {
      init();
      init(fname_str);
    }

    Timer(const std::string& cname_str, const std::string& fname_str) {
      init();
      init(cname_str, fname_str);
    }

    void init() {
      cname = "Timer";
      getStartTime();
      initializePAPI();
      index = getTimerDatabase().size();
      getTimerDatabase().push_back(&info);
    }

    void init(const std::string& fname_str) {
      info.fname = fname_str;
    }

    void init(const std::string& cname_str, const std::string& fname_str) {
      info.fname = cname_str;
      info.fname += "::";
      info.fname += fname_str;
    }

    void start(bool verbose = false) {
      if (verbose || info.dtime >= minimum_duration_for_show_start_info()) {
        info.showLast("start");
      }
      start_time = dclock();
      start_flops = getTotalFlops();
    }

    void stop(bool verbose = false) {
      stop_flops = getTotalFlops();
      stop_time = dclock();
      info.dtime = stop_time - start_time;
      info.dflops = stop_flops - start_flops;
      info.accumulated_time += info.dtime;
      info.accumulated_flops += info.dflops;
      info.call_times++;
      if (verbose || info.dtime >= minimum_duration_for_show_stop_info()) {
        info.showLast("stop ");
      }
      autodisplay();
    }

    static void clear_all()
    {
	std::vector<TimerInfo *> db(getTimerDatabase());
	for (int i = 0; i < db.size(); i++) {
	    db[i]->clear();
	}
	getStartTime() = dclock();
    }

    static void display(const std::string& str = "") {
      double total_time = getTotalTime();
      VRB.Result("Timer", "DisplayStart", "%s ------------ total %.4e secs -----------------------\n", str.c_str(), total_time);
      std::vector<TimerInfo *> db(getTimerDatabase());
      std::sort(db.begin(), db.end(), compareTimeInfoP);
      for (int i = 0; i < db.size(); i++) {
        db[i]->showAvg("Display");
      }
      VRB.Result("Timer", "DisplayEnd  ", "%s ------------ total %.4e secs -----------------------\n", str.c_str(), total_time);
    }

    static void autodisplay() {
      static double last_time = getTotalTime();
      double time = getTotalTime();
      if (time - last_time > minimum_autodisplay_interval()) {
        last_time = time;
        display("autodisplay");
      }
    }
  static void reset();
  //Time since last reset
  static Elapsed elapsed_time();
  //Time since last call to this function
  static Elapsed relative_time();

};


CPS_END_NAMESPACE

#endif
