// Copyright (C) 2016 William H. Greene
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

#include <math.h>
#include <assert.h>

#include "GausLegendreIntRule.h"

// From Felippa
// http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.AppI.d/AFEM.AppI.pdf
//

/*
c1 = sqrt(5-c0)/3;
c2 = sqrt(5+c0)/3;
w0 = 13*sqrt(70);
w1 = (322 + w0)/900;
w2 = (322 - w0)/900;
%  g5={-Sqrt[5+2*Sqrt[10/7]],-Sqrt[5-2*Sqrt[10/7]],0,
% Sqrt[5-2*Sqrt[10/7]], Sqrt[5+2*Sqrt[10/7]]}/3,
% w5={322-13*Sqrt[70],322+13*Sqrt[70],512,
% 322+13*Sqrt[70],322-13*Sqrt[70]}/900,
xi = [-c2 -c1 0 c1 c2]; w = [w2 w1 512/900 w1 w2];
*/

static const double s3o3 = sqrt(3.0) / 3;
static const double s3o5 = sqrt(3.0 / 5.);
static const double xi4p0 = 2.0*sqrt(6. / 5.);
static const double xi4p1 = sqrt((3. - xi4p0) / 7.);
static const double xi4p2 = sqrt((3. + xi4p0) / 7.);
static const double wt4p1 = (18. + sqrt(30.)) / 36.;
static const double wt4p2 = (18. - sqrt(30.)) / 36.;
// constants for 5-point rule
static const double s107 = 2 * sqrt(10. / 7.);
static const double xi5p1 = sqrt(5. + s107) / 3.;
static const double xi5p2 = sqrt(5. - s107) / 3.;
static const double wt5p0 = 13. * sqrt(70.);
static const double wt5p1 = (322. - wt5p0) / 900.;
static const double wt5p2 = (322. + wt5p0) / 900.;

static GausLegendreIntRule::IRule rule1Pt[] = { 0, 2 };
static GausLegendreIntRule::IRule rule2Pt[] = {
    { -s3o3, 1 }, { s3o3, 1 }
};

static GausLegendreIntRule::IRule rule3Pt[] = {
    { -s3o5, 5. / 9. }, { 0., 8. / 9. }, { s3o5, 5. / 9. }
};
static GausLegendreIntRule::IRule rule4Pt[] = {
    { -xi4p2, wt4p2 }, { -xi4p1, wt4p1 }, { xi4p1, wt4p1 }, { xi4p2, wt4p2 }
};
static GausLegendreIntRule::IRule rule5Pt[] = {
    { -xi5p1, wt5p1 }, { -xi5p2, wt5p2 }, { 0., 512. / 900. },
    { xi5p2, wt5p2 }, { xi5p1, wt5p1 }
};


static  GausLegendreIntRule::IRule *rules[] = {
  0, rule1Pt, rule2Pt, rule3Pt, rule4Pt, rule5Pt };

GausLegendreIntRule::GausLegendreIntRule(int npts) : numPts(npts)
{
  assert(npts>0 && npts <= (sizeof(rules) / sizeof(rules[0])));
  rule = rules[npts];
}

void GausLegendreIntRule::getPoint(int iPt, double &xi, double &weight) const
{
  xi = rule[iPt].xi;
  weight = rule[iPt].wt;
}

const GausLegendreIntRule::IRule &GausLegendreIntRule::getPoint(int iPt) const
{
  return rule[iPt];
}