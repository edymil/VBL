/*
 *  Utilities.cpp
 *  VBL
 *
 *  Created by Edoardo Milotti on 8/22/08.
 *  Copyright 2008 I.N.F.N.-Sezione di Trieste. All rights reserved.
 *
 */
 
#ifndef UTILITY_H
#define UTILITY_H  // header guard
// numeri casuali 

// versione standalone di ran2 (versione NR modificata)

// for old compilers we need some tricks
#include<features.h> // gives various information on the build system
#if __GNUC_PREREQ(4,9)
    // means gnu compiler is higher than 4.5 and fully c++11 compatible
    // we don't need to hacking
  #include <random> // std random number generator
#else
  //old fashioned
  #define OLD_RANDOM
  #include<cstdlib>
  #include<time.h>
#endif

#include "sim.h"

#include "InputFromFile.h"
#include "CellType.h"
#include "Environment.h"
#include "EnvironmentalSignals.h"
#include "geom-2.h"
#include "CellsSystem.h"
#include "Utilities.h"


double ran2(int &idum);
double gasdev(int &idum);
double gammln(const double xx);
double factln(const int n);
double bico(const int n, const int k);
double bnldev(const double pp, const int n, int &idum);

#endif //headerguard
