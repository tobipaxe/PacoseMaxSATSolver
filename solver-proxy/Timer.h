/********************************************************************************************
SolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian, Nhat Minh Hoang

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************/

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>

enum TimeUnit { SECONDS, MILLISECONDS, MICROSECONDS };

double GetTimeUnitFactor(TimeUnit tunit)
{
    switch (tunit) {
    case SECONDS:
        return 1;
    case MILLISECONDS:
        return 1e3;
    case MICROSECONDS:
        return 1e6;
    }
    return 0;    
}

class CpuTimer
{
public:
    CpuTimer() { _totalElapsedTime = 0; };

    /**
     * Start the timer at the current time of execution.
     */
    void SetTimeReference()
    {
        _running = true;
        _begin = std::chrono::steady_clock::now();
    };

    /**
     * Set the total elapsed time to 0 to start a new run.
     */
    void Reset() { _totalElapsedTime = 0; }

    /**
     * Reset the timer to begin a new run.
     */
    void Restart()
    {
        Reset();
        SetTimeReference();
    }

    /**
     * Return the total elapsed time since the previous call of SetTimeReference.
     *
     * @return The total elapsed time.
     */
    double TotalTime() { return _totalElapsedTime; };

    /**
     * Return the difference between the previous measured start time and the current time.
     *
     * @param tunit Scale the total elapsed time by the factor given by the unit.
     * @return The total elapsed time.
     */
    double TimeSinceReference(TimeUnit tunit = SECONDS)
    {
        _end = std::chrono::steady_clock::now();
        _factor = GetTimeUnitFactor(tunit);
        return std::chrono::duration_cast<std::chrono::microseconds>(_end - _begin).count()
               / (1e6 / _factor);
    };

    /**
     * Stop the timer and scale the total elapsed time by a factor given by the time unit.
     *
     * @param tunit Scale the total elapsed time by the factor given by the unit.
     */
    void Stop(TimeUnit tunit)
    {
        if (_running) {
            _totalElapsedTime += TimeSinceReference(tunit);
            _running = false;
        }
    };

private:
    std::chrono::steady_clock::time_point _begin;
    std::chrono::steady_clock::time_point _end;
    double _factor;
    double _totalElapsedTime;
    bool _running;
};

struct FunctionTimer
{
    double *timeDiff;
    double timeStart;
    double timeUnitFactor;
    bool addTimeDiff;

    FunctionTimer(double *timeVar, TimeUnit tunit, bool add = false)
        : timeDiff(timeVar)
        , timeStart(0)
        , timeUnitFactor(GetTimeUnitFactor(tunit))
        , addTimeDiff(add)
    {
        struct rusage resources;
        getrusage(RUSAGE_SELF, &resources);
        timeStart = resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
    }

    ~FunctionTimer()
    {
        struct rusage resources;
        getrusage(RUSAGE_SELF, &resources);
        if (addTimeDiff && *timeDiff < 0)
            *timeDiff += static_cast<double>(
                (resources.ru_utime.tv_sec + 1.e-6 * static_cast<double>(resources.ru_utime.tv_usec))
                - timeStart);
        else
            *timeDiff = ((resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec)
                         - timeStart)
                        * timeUnitFactor;
    }
};

#endif /* TIMER_H */
