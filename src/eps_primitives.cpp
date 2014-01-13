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

void plotVLineCenter(ostream &o, int x, int y, int size);
void plotHLineCenter(ostream &o, int x, int y, int size);

void plotRectCenter(ostream &o, int x, int y, int sizex, int sizey, double color)
{
  plotRect(o,x-sizex/2,y-sizey/2,sizex,sizey,color);
  return;
}

void plotBoxCenter(ostream &o, int x, int y, int size, double color)
{
  plotRectCenter(o,x,y,size,size,color);
  return;
}

void plotRect(ostream &o, int x, int y, int sizex, int sizey, double color)
{
  /*
  o << "newpath\n"
    << "\t0 setgray\n"
    << "\t" << x << " " << y << " moveto\n"
    << "\t" << 0 << " " << sizey << " rlineto\n"
    << "\t" << sizex << " " << 0 << " rlineto\n"
    << "\t" << 0 << " " << -sizey << " rlineto\n"
    << "\tclosepath\n"
    << "\t1 setlinewidth\n"
    << "stroke\n";
  if(color >= 0)
    {
      o << "\t" << x << " " << y << " moveto\n"
	<< "\t" << 0 << " " << sizey << " rlineto\n"
	<< "\t" << sizex << " " << 0 << " rlineto\n"
	<< "\t" << 0 << " " << -sizey << " rlineto\n"
	<< "\tclosepath\n"
	<< "\t" << 1-color << " setgray\n"
	<< "\tfill\n";
    }
  o << endl;
  */
  o << sizex << " " 
    << sizey << " " 
    << x << " " 
    << y << " drawrect\n";

  if(color >= 0)
    {
      o << 1-color << " " 
	<< sizex << " " 
	<< sizey << " " 
	<< x << " " 
	<< y << " drawfillrect\n";
    }  
  return;
}

void finalizeEPS(ostream &o)
{
    o << "showpage\n";
}

void setUpEPS(ostream &o,int llx, int lly, int urx, int ury,string creator, string title, string date)
{
  o << "%!PS-Adobe-3.0 EPSF-3.0\n"
    << "%%BoundingBox: " << llx << " " << lly << " " << urx << " " << ury << "\n"
    << "%%Creator: " << creator << "\n"
    << "%%Title: " << title << "\n"
    << "%%CreationDate: " << date << "\n\n";

  o << "/rightshow\n"
    << "{dup stringwidth pop\n"
    << "0 exch sub\n"
    << "0 rmoveto\n"
    << "show } def\n\n";

  o << "/centershow\n"
    << "{dup stringwidth pop\n"
    << "2 div\n"
    << "0 exch sub\n"
    << "0 rmoveto\n"
    << "show } def\n\n";

  o << "%Usage: sizex sizey x y drawrect\n"
    << "/drawrect\n"
    << "{\n"
    << "    newpath\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup 0 exch rlineto\n"
    << "    exch 0 rlineto\n"
    << "    0 exch sub\n"
    << "    0 exch rlineto\n"
    << "    closepath\n"
    << "    1 setlinewidth\n"
    << "    stroke\n"
    << "} def\n\n";

  o << "%Usage: color sizex sizey x y drawfillrect\n"
    << "/drawfillrect\n"
    << "{\n"
    << "    newpath\n"
    << "    moveto\n"
    << "    dup 0 exch rlineto\n"
    << "    exch 0 rlineto\n"
    << "    0 exch sub\n"
    << "    0 exch rlineto\n"
    << "    closepath\n"
    << "    setgray\n"
    << "    fill\n"
    << "} def\n\n";

  o << "%Usage: text rotation x y size drawtext\n"
    << "/drawtext\n"
    << "{\n"
    << "    /Helvetica findfont exch scalefont setfont\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup rotate\n"
    << "    exch show\n"
    << "    0 exch sub rotate\n"
    << "} def\n\n";

  o << "%Usage: text rotation x y size drawrighttext\n"
    << "/drawrighttext\n"
    << "{\n"
    << "    /Helvetica findfont exch scalefont setfont\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup rotate\n"
    << "    exch rightshow\n"
    << "    0 exch sub rotate\n"
    << "} def\n\n";

  o << "%Usage: text rotation x y size drawcentertext\n"
    << "/drawcentertext\n"
    << "{\n"
    << "    /Helvetica findfont exch scalefont setfont\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup rotate\n"
    << "    exch centershow\n"
    << "    0 exch sub rotate\n"
    << "} def\n\n";

  o << "%Usage: size x y drawLLtriangle\n"
    << "/drawLLtriangle\n"
    << "{\n"
    << "   newpath\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup 0 exch rlineto\n"
    << "    dup 0 exch sub rlineto\n"
    << "    closepath\n"
    << "    1 setlinewidth\n"
    << "    stroke\n"
    << "} def\n\n";
  
  o << "%Usage: color size x y drawfillLLtriangle\n"
    << "/drawfillLLtriangle\n"
    << "{\n"
    << "    newpath\n"
    << "    moveto\n"
    << "    dup 0 exch rlineto\n"
    << "    dup 0 exch sub rlineto\n"
    << "    closepath\n"
    << "    setgray\n"
    << "    fill\n"
    << "} def\n\n";

  o << "%Usage: size x y drawVLineCenter\n"
    << "/drawVLineCenter\n"
    << "{\n"
    << "    newpath\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup 2 div\n"
    << "    0 exch rmoveto\n"
    << "    0 exch sub\n"
    << "    0 exch rlineto\n"
    << "    1 setlinewidth\n"
    << "    stroke\n"
    << "} def\n\n";
  
  o << "%Usage: size x y drawHLineCenter\n"
    << "/drawHLineCenter\n"
    << "{\n"
    << "    newpath\n"
    << "    0 setgray\n"
    << "    moveto\n"
    << "    dup 2 div\n"
    << "    0 rmoveto\n"
    << "    0 exch sub\n"
    << "    0 rlineto\n"
    << "    1 setlinewidth\n"
    << "    stroke\n"
    << "} def\n\n";
  
}


void plotText(ostream &o, int x, int y, int size, double rot,string text)
{
  o << "(" << text << ") "
    << rot << " "
    << x << " "
    << y << " "
    << size << " drawtext\n";
  /*
   o << "/Helvetica findfont\n"
    << "0 setgray\n"
    << size << " scalefont\n"
    << "setfont\n"
    << x << " " << y << " moveto\n"
    << rot << " rotate\n"
    << "(" << text << ") show\n"
    << -rot << " rotate\n";
  */

  return;
}

void plotTextRight(ostream &o, int x, int y, int size, double rot,string text)
{
  o << "(" << text << ") "
    << rot << " "
    << x << " "
    << y << " "
    << size << " drawrighttext\n";
  /*
  o << "/Helvetica findfont\n"
    << "0 setgray\n"
    << size << " scalefont\n"
    << "setfont\n"
    << x << " " << y << " moveto\n"
    << rot << " rotate\n"
    << "(" << text << ") rightshow\n"
    << -rot << " rotate\n";
  */
  return;
}

void plotTextCenter(ostream &o, int x, int y, int size, double rot,string text)
{
  o << "(" << text << ") "
    << rot << " "
    << x << " "
    << y << " "
    << size << " drawcentertext\n";
  /*
  o << "/Helvetica findfont\n"
    << "0 setgray\n"
    << size << " scalefont\n"
    << "setfont\n"
    << x << " " << y << " moveto\n"
    << rot << " rotate\n"
    << "(" << text << ") centershow\n"
    << -rot << " rotate\n";
  */
  return;
}

void plotBox(ostream &o, int x, int y, int size, double color)
{
  plotRect(o,x,y,size,size,color);
    /*
  o << "newpath\n"
    << "\t0 setgray\n"
    << "\t" << x << " " << y << " moveto\n"
    << "\t" << 0 << " " << size << " rlineto\n"
    << "\t" << size << " " << 0 << " rlineto\n"
    << "\t" << 0 << " " << -size << " rlineto\n"
    << "\tclosepath\n"
    << "\t1 setlinewidth\n"
    << "stroke\n";
  if(color >= 0)
    {
      o << "\t" << x << " " << y << " moveto\n"
	<< "\t" << 0 << " " << size << " rlineto\n"
	<< "\t" << size << " " << 0 << " rlineto\n"
	<< "\t" << 0 << " " << -size << " rlineto\n"
	<< "\tclosepath\n"
	<< "\t" << 1-color << " setgray\n"
	<< "\tfill\n";
    }
  o << endl;
    */
  return;
}

void plotLLTriangle(ostream &o, int x, int y, int size, double color)
{

  o << size << " " 
    << x << " " 
    << y << " drawLLtriangle\n";
  
  if(color >= 0)
    {
      o << 1-color << " " 
	<< size << " " 
	<< x << " " 
	<< y << " drawfillLLtriangle\n";
    }  


  /*
  o << "newpath\n"
    << "\t0 setgray\n"
    << "\t" << x << " " << y << " moveto\n"
    << "\t" << 0 << " " << size << " rlineto\n"
    << "\t" << size << " " << -size << " rlineto\n"
    << "\tclosepath\n"
    << "\t1 setlinewidth\n"
    << "stroke\n"
    << "\t" << x << " " << y << " moveto\n"
    << "\t" << 0 << " " << size << " rlineto\n"
    << "\t" << size << " " << -size << " rlineto\n"
    << "\tclosepath\n"
    << "\t" << 1-color << " setgray\n"
    << "\tfill\n\n";
  */
  return;
}
