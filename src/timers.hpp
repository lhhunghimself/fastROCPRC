#include <time.h>
#include <sys/time.h>
void current_utc_time(struct timespec *ts);
double get_elapsed_time(const struct timespec *start_time,const struct timespec *end_time);

double get_elapsed_time(const struct timespec *start_time, const struct timespec *end_time)
{
    int64_t sec = end_time->tv_sec - start_time->tv_sec;
    int64_t nsec;
    if (end_time->tv_nsec >= start_time->tv_nsec) {
        nsec = end_time->tv_nsec - start_time->tv_nsec;
    } else {
        nsec = 1000000000 - (start_time->tv_nsec - end_time->tv_nsec);
        sec -= 1;
    }
    return ((double) sec + (double) nsec *1e-9);  
}
	
void current_utc_time(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_MONOTONIC_RAW, ts);
#endif
}
