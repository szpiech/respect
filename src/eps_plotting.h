/*
    Copyright (C) 2011  Zachary A Szpiech (szpiechz@umich.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __EPS_PRIMITIVES_H__
#define __EPS_PRIMITIVES_H__
#include "eps_primitives.h"
#endif
#include <iostream>

using namespace std;

void plotAxes(ostream &o, int x, int y, int size, double xlow, double xhigh,
	      double ylow, double yhigh, string xlab, string ylab);
void scatterData(ostream &o, int x, int y, int size,
		 double **data, int dataSize, short **color, int colorSize);
