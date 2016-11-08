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

#include <cmath>
#include <assert.h>
#include <climits>
#include <iostream>
using std::cout;
using std::endl;

#include <Eigen/Core>

#include "GausLegendreIntRule.h"
#include "util.h"

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

// 6pt rule
static const double xi6p1 = 0.2386191860831969;
static const double xi6p2 = 0.6612093864662645;
static const double xi6p3 = 0.9324695142031521;
static const double wt6p1 = 0.4679139345726910;
static const double wt6p2 = 0.3607615730481386;
static const double wt6p3 = 0.1713244923791704;
static GausLegendreIntRule::IRule rule6Pt[] = {
    { -xi6p1, wt6p1 }, { -xi6p2, wt6p2 }, { -xi6p3, wt6p3 },
    { xi6p1, wt6p1 }, { xi6p2, wt6p2 }, { xi6p3, wt6p3 },
};

// 7pt rule
static const double xi7p1 = 0.4058451513773972;
static const double xi7p2 = 0.7415311855993945;
static const double xi7p3 = 0.9491079123427585;
static const double w7p0 = 0.4179591836734694;
static const double w7p1 = 0.3818300505051189;
static const double w7p2 = 0.2797053914892766;
static const double w7p3 = 0.1294849661688697;
static GausLegendreIntRule::IRule rule7Pt[] = {
    { -xi7p3, w7p3 }, { -xi7p2, w7p2 }, { -xi7p1, w7p1 }, { 0, w7p0 },
    { xi7p1, w7p1 }, { xi7p2, w7p2 }, { xi7p3, w7p3 }
};

// 8pt rule
static const double xi8p1 = 0.1834346424956498049394761;
static const double w8p1 = 0.3626837833783619829651504;
static const double xi8p2 = 0.5255324099163289858177390;
static const double w8p2 = 0.3137066458778872873379622;
static const double xi8p3 = 0.7966664774136267395915539;
static const double w8p3 = 0.2223810344533744705443560;
static const double xi8p4 = 0.9602898564975362316835609;
static const double w8p4 = 0.1012285362903762591525314;
static GausLegendreIntRule::IRule rule8Pt[] = {
    { -xi8p4, w8p4 }, { -xi8p3, w8p3 }, { -xi8p2, w8p2 }, { -xi8p1, w8p1 },
    { xi8p1, w8p1 }, { xi8p2, w8p2 }, { xi8p3, w8p3 }, { xi8p4, w8p4 }
};

// 9pt rule
static const double w9p0 = 0.3302393550012597631645251;
static const double xi9p1 = 0.3242534234038089290385380;
static const double w9p1 = 0.3123470770400028400686304;
static const double xi9p2 = 0.6133714327005903973087020;
static const double w9p2 = 0.2606106964029354623187429;
static const double xi9p3 = 0.8360311073266357942994298;
static const double w9p3 = 0.1806481606948574040584720;
static const double xi9p4 = 0.9681602395076260898355762;
static const double w9p4 = 0.0812743883615744119718922;
static GausLegendreIntRule::IRule rule9Pt[] = {
    { -xi9p4, w9p4 }, { -xi9p3, w9p3 }, { -xi9p2, w9p2 }, { -xi9p1, w9p1 },
    { 0, w9p0 }, { xi9p1, w9p1 }, { xi9p2, w9p2 }, { xi9p3, w9p3 }, { xi9p4, w9p4 }
};

// 10pt rule
static const double xi10p1 = 0.1488743389816312108848260;	
static const double w10p1 = 0.2955242247147528701738930;
static const double xi10p2 = 0.4333953941292471907992659;	
static const double w10p2 = 0.2692667193099963550912269;
static const double xi10p3 = 0.6794095682990244062343274;	
static const double w10p3 = 0.2190863625159820439955349;
static const double xi10p4 = 0.8650633666889845107320967;	
static const double w10p4 = 0.1494513491505805931457763;
static const double xi10p5 = 0.9739065285171717200779640;	
static const double w10p5 = 0.0666713443086881375935688;
static GausLegendreIntRule::IRule rule10Pt[] = {
    { -xi10p5, w10p5 }, { -xi10p4, w10p4 }, { -xi10p3, w10p3 }, { -xi10p2, w10p2 }, { -xi10p1, w10p1 },
    { xi10p1, w10p1 }, { xi10p2, w10p2 }, { xi10p3, w10p3 }, { xi10p4, w10p4 },
    { xi10p5, w10p5 },
};


static  GausLegendreIntRule::IRule *rules[] = {
  0, rule1Pt, rule2Pt, rule3Pt, rule4Pt, rule5Pt, rule6Pt, rule7Pt, rule8Pt,
  rule9Pt, rule10Pt
};

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


GausLegendreIntRule::IRuleList GausLegendreIntRule::computePts()
{
  IRuleList pts(numPts);
  int N = numPts - 1;
  int N1 = N + 1, N2 = N + 2;
  const double pi = 2*asin(1.);
  const double eps = std::numeric_limits<double>::epsilon();

  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(N1, N2);
  Eigen::MatrixXd Lp = Eigen::MatrixXd::Zero(N1, 1);

  // y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
  Eigen::VectorXd xu = pi*N/N2*Eigen::VectorXd::LinSpaced(N1, -1, 1);
  Eigen::VectorXd y = 2*Eigen::VectorXd::LinSpaced(N1, 0, N).array() + 1;
  y = (y*pi / (2.*N + 2.)).array().cos() + 0.27 / N1*xu.array().sin();

  //cout << y << endl;
  Eigen::VectorXd y0(N1);
  y0.setConstant(2);
  int iter = 0;
  
  while (true) {
    double maxD = (y - y0).cwiseAbs().maxCoeff();
    //printf("iter = %d, maxD=%g\n", iter++, maxD);
    if (maxD < eps) break;
    L.col(0).setConstant(1);
    L.col(1) = y;
    for (int k = 1; k < N1; k++) {
      // L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
      L.col(k + 1) = ((2 * k + 1)*y.cwiseProduct(L.col(k)) -
        k*L.col(k - 1)) / (double) (k+1);
    }
    //cout << "L\n" << L << endl;

    // Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    Eigen::VectorXd tmp = L.col(N1 - 1) - y.cwiseProduct(L.col(N2 - 1));
    Lp = N2*tmp.cwiseQuotient((1 - y.array().square()).matrix());
    //cout << "Lp\n" << Lp << endl;

    y0 = y;
    y = y0 - L.col(N2 - 1).cwiseQuotient(Lp);
  }

  printMat(y, "y", "%21.17f");

  //w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
  double N2N1 = (double)N2 / (double)N1;
  Eigen::VectorXd wts = 2*N2N1*N2N1 / 
    (1 - y.array().square()).cwiseProduct(Lp.array().square());
  printMat(wts, "wts", "%21.17f");
  return pts;
}