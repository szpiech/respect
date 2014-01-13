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

#include <iostream>

using namespace std;

void plotVLineCenter(ostream &o, int x, int y, int size);
void plotHLineCenter(ostream &o, int x, int y, int size);
void plotBox(ostream &o, int x, int y, int size, double color);
void plotBoxCenter(ostream &o, int x, int y, int size, double color);
void plotRect(ostream &o, int x, int y, int sizex, int sizy, double color);
void plotRectCenter(ostream &o, int x, int y, int sizex, int sizey, double color);
void plotLLTriangle(ostream &o, int x, int y, int size, double color);
void plotText(ostream &o, int x, int y, int size, double rot, string text);
void plotTextRight(ostream &o, int x, int y, int size, double rot,string text);
void plotTextCenter(ostream &o, int x, int y, int size, double rot,string text);
void setUpEPS(ostream &o,int llx, int lly, int urx, int ury,string creator, string title, string date);
void finalizeEPS(ostream &o);
