/********************************************************************************************
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************************************/

#ifndef TIMEVARIABLES_H
#define TIMEVARIABLES_H

//#include <iomanip>
#include <iostream>

namespace Pacose {
namespace DGPW {

struct TimeVariables {
  TimeVariables()
      : solvedFirst(-1),
        fillingBuckets(-1),
        createTree(-1),
        encoding(-1),
        solving(-1),
        solvingLastBucket(-1),
        solvingTares(-1),
        total(-1) {}

  double solvedFirst;
  double fillingBuckets;
  double createTree;
  double encoding;
  double solving;
  double solvingLastBucket;
  double solvingTares;
  double total;

  void AdjustSolvingTime() {
    if (solvingLastBucket < 0 && encoding < 0 && createTree < 0) {
      solvingLastBucket -= encoding;
    }
    if (solvingLastBucket > 0 && solvingTares > 0)
      solving = solvingLastBucket + solvingTares;
  }

  void DumpVariables(unsigned iteration = 0) {
    //        std::cout << "Current cascade iteration: " << iteration <<
    //        std::endl;
    std::string iter;
    std::string point;
    if (iteration == 0) {
      iter = "";
      point = "..";

    } else {
      iter = std::to_string(iteration) + " ";
      point = "";
    }
    // std::cout << __FUNCTION__ << std::endl;
    AdjustSolvingTime();
    if (solvedFirst >= 0)
      std::cout << "c time " << iter << "approximation..." << point << ": "
                << solvedFirst << std::endl;
    if (fillingBuckets >= 0)
      std::cout << "c time " << iter << "filling buckets." << point << ": "
                << fillingBuckets << std::endl;
    if (createTree >= 0)
      std::cout << "c time " << iter << "creating tree..." << point << ": "
                << createTree << std::endl;
    if (encoding >= 0)
      std::cout << "c time " << iter << "encoding........" << point << ": "
                << encoding << std::endl;
    if (solving >= 0)
      std::cout << "c time " << iter << "solving........." << point << ": "
                << solving << std::endl;
    if (solvingLastBucket >= 0)
      std::cout << "c time " << iter << "solving last bkt" << point << ": "
                << solvingLastBucket << std::endl;
    if (solvingTares >= 0)
      std::cout << "c time " << iter << "solving tares..." << point << ": "
                << solvingTares << std::endl;
    if (total >= 0)
      std::cout << "c dgpw " << iter << "time............" << point << ": "
                << total << std::endl;
  }

  void AddTimeStruct(TimeVariables *tv) {
    if (solvedFirst < 0) {
      solvedFirst = tv->solvedFirst;
    } else if (tv->solvedFirst > 0) {
      solvedFirst += tv->solvedFirst;
    }
    if (fillingBuckets < 0) {
      fillingBuckets = tv->fillingBuckets;
    } else if (tv->fillingBuckets > 0) {
      fillingBuckets += tv->fillingBuckets;
    }
    if (encoding < 0) {
      encoding = tv->encoding;
    } else if (tv->encoding > 0) {
      encoding += tv->encoding;
    }
    if (solving < 0) {
      solving = tv->solving;
    } else if (tv->solving > 0) {
      solving += tv->solving;
    }
    if (solvingLastBucket < 0) {
      solvingLastBucket = tv->solvingLastBucket;
    } else if (tv->solvingLastBucket > 0) {
      solvingLastBucket += tv->solvingLastBucket;
    }
    if (solvingTares < 0) {
      solvingTares = tv->solvingTares;
    } else if (tv->solvingTares > 0) {
      solvingTares += tv->solvingTares;
    }
    if (total < 0) {
      total = tv->total;
    } else if (tv->total > 0) {
      total += tv->total;
    }
  }
};

}  // namespace DGPW
} // Namespace Pacose

#endif  // TIMEVARIABLES_H
