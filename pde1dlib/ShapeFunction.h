// Copyright (C) 2016-2017 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#pragma once

class ShapeFunction
{
public:
  virtual void N(double r, double *vals) const =0;
  virtual void dNdr(double r, double *vals) const  = 0;
  virtual int numNodes() const = 0;
};

class ShapeFunction2 : public ShapeFunction
{
public:
  virtual void N(double r, double *vals) const;
  virtual void dNdr(double r, double *vals) const;
  virtual int numNodes() const { return 2; }
};

class ShapeFunction3 : public ShapeFunction
{
public:
  virtual void N(double r, double *vals) const;
  virtual void dNdr(double r, double *vals) const;
  virtual int numNodes() const { return 3; }
};

