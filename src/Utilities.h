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
#include <random> // std random number generator
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
