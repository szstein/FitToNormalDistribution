//
//  main.cpp
//  fitDIstribution
//
//  Created by Stephen Stein on 4/2/20.
//  Copyright © 2020 Stephen Stein. All rights reserved.
//
// NY deaths starting 3/14/2020

#define VERBOSE 0
#include <iostream>
#include <fstream>      // std::ifstream
#include <cmath>

double stdNorm(double x)
{
    double y = (1.0/sqrt(2.0*M_PI))*exp(-0.5*x*x);
    return y;
}

double norm(double x,double mean,double sig)
{
    double y = (1.0/sig)*stdNorm((x-mean)/sig);
    return y;
}
std::string datestring;
int monthFromDatestring(std::string str)
{
    size_t slashPos = str.find("/");
    int retVal = std::stoi(str.substr(0,slashPos));
    return retVal;
}
int dayFromDatestring(std::string str)
{
    size_t slashPos = str.find("/");
    size_t slashPos2 = str.find("/",slashPos);
    int retVal = std::stoi(str.substr(slashPos2+1));
    return retVal;
}

int daysInMonth[] = {0,31,29,31,30,31,30,31,31,30,31,30,31};
std::string dayToDate(int day)
{
    int d = dayFromDatestring(datestring);
    int m = monthFromDatestring(datestring);
    d = d+day;
    while (d>daysInMonth[m]) {
        d -= daysInMonth[m];
        m++;
    }
    char buff[100];
    sprintf(buff,"%d/%d/2020",m,d);
    return std::string(buff);
}


// USA data 4/3
//double testSeries[] = {14, 16, 18, 22, 24, 27, 36, 39, 49, 60, 71, 90, 112, 160, 219, 272, 398,
//                        471, 675, 900, 1163, 1530, 1965, 2428, 2939, 3746, 4700, 5784, 6962};
// NY 4/4
// double testSeries[] = {3, 7, 7, 12, 12, 35, 44, 114, 114, 210, 285, 385, 519, 728, 965, 1218, 1550, 1941, 2373, 2935, 3565};

// Italy 4/3
//double testSeries[] = {0, 1, 2, 3, 7, 10, 12, 17, 21, 29, 34, 52, 79, 107, 148, 197, 233, 366, 463, 631, 827, 827, 1266, 1441, 1809, 2158, 2503, 2978, 3405, 4032, 4825, 5476, 6077, 6820, 7503, 8215, 9134, 10023, 10779, 11591, 12428, 13155, 13915, 14681};
double *testSeries;
double *diffs;
double *smoothedDiffs;
double startMean = 31;
int nPoints;
double multiplier;
int maxDiff = 0;
int maxDiffPt = 0;
#define LASTDAY 120

double scoreGuess(double mean, double dev)
{
    multiplier = maxDiff/norm(maxDiffPt+1,mean,dev);
#if VERBOSE
    std::cout << maxPoint << "\t"<< maxVal << "\t" <<"multiplier = " << multiplier << "\n";
#endif
    double totDiffTotSquared = 0.0;
    double modelTotal = 0.0;
    for (int i=1; i<LASTDAY; i++) {
#if VERBOSE
        std::cout << i << "\t" << norm(i,mean,dev) << "\t" << multiplier*norm(i,mean,dev) << "\n";
#endif
        if (i<nPoints) {
            double modelDeathsThisPeriod = multiplier*norm(i,mean,dev);
            modelTotal += modelDeathsThisPeriod;
            double diffTot = modelTotal - testSeries[i-1];
            double diffTotSquared = diffTot*diffTot;
            totDiffTotSquared += diffTotSquared;
        }
    }
    return totDiffTotSquared;
}

void getPoints(const char *filename)
{
    // the input file should have a date, then values one per line
    std::ifstream ifs (filename, std::ifstream::in);
    double data;
    nPoints = 0;
    testSeries = (double*)malloc(0);
    ifs >> datestring;
    while (!ifs.eof()) {
        ifs >> data;
        testSeries = (double*)realloc(testSeries, ++nPoints*sizeof(double));
        testSeries[nPoints-1] = data;
    }
}

int main(int argc, const char * argv[]) {
    if (argc<2) {
        std::cout << "input file missing\n";
        return -1;
    }
    getPoints(argv[1]);
    diffs = (double*)malloc(nPoints*sizeof(double));
    smoothedDiffs = (double*)malloc(nPoints*sizeof(double));
    std::cout << nPoints << "\n";
    diffs[0] = testSeries[0];
    smoothedDiffs[0] = 0.75*diffs[0] + 0.25*diffs[1];
    for (int i=1;i<nPoints-1;i++) {
        diffs[i] = testSeries[i] - testSeries[i-1];
        if (i==nPoints-1) {
            smoothedDiffs[i] = 0.25*diffs[i-1] + 0.75*diffs[i];
        }
        else {
            smoothedDiffs[i] = 0.25*diffs[i-1] + 0.5*diffs[i] + 0.25*diffs[i+1];
        }
        std::cout << i << "\t" << diffs[i] << "\n";
        if (diffs[i]>maxDiff) {
            maxDiff = diffs[i];
            maxDiffPt = i;
        }
    }
    // if the last point is a max, move the ref point back one and smooth the last 3 values
    if (maxDiffPt==nPoints-2) {
        maxDiff = 0.25*diffs[maxDiffPt-2] + 0.5*diffs[maxDiffPt-1] + 0.25*diffs[maxDiffPt];
        maxDiffPt -= 1;
    }
    else {
        // smooth the max point
        maxDiff = 0.25*diffs[maxDiffPt-1] + 0.5*diffs[maxDiffPt] + 0.25*diffs[maxDiffPt+1];
    }
    double bestScore = std::numeric_limits<double>::max();
    double bestM = 0;
    double bestR = 0;
    for (int m = nPoints+20; m>nPoints-10; m--) {
        double bestScoreForMean = std::numeric_limits<double>::max();
        double bestRforMean = 0;
        for (int r=2; r<15; r++) {
            double score = scoreGuess(m, r);
#if VERBOSE
            std::cout << "mean="  << m << ", dev=" << r << ", score=" << score << "\n";
#endif
            if (score<bestScoreForMean) {
                bestScoreForMean = score;
                bestRforMean = r;
                if (score<bestScore) {
                    bestScore = score;
                    bestM = m;
                    bestR = r;
                    std::cout << "BEST: mean=" << bestM << " dev=" << bestR << " score=" << bestScore << "\n";
                }
            }
        }
        std::cout << "BEST score for mean=" << m << " dev=" << bestRforMean << " score=" << bestScoreForMean << "\n";
    }
    std::cout << "BEST: mean=" << bestM << " dev=" << bestR << " score=" << bestScore << "\n";
    double downScore;
    double upScore;
    downScore = scoreGuess(bestM,bestR-0.1);
    if (downScore < bestScore) {
        while (downScore<bestScore) {
            bestR -= 0.1;
            bestScore = downScore;
            downScore = scoreGuess(bestM,bestR-0.1);
        }
    }
    else {
        upScore = scoreGuess(bestM,bestR+0.1);
        while (upScore<bestScore) {
            bestR += 0.1;
            bestScore = upScore;
            upScore = scoreGuess(bestM,bestR+0.1);
        }
    }
    bestScore = scoreGuess(bestM,bestR);
    std::cout << "BEST: mean=" << bestM << " dev=" << bestR << " score=" << bestScore << " multiplier=" << multiplier << "\n";
    std::cout << "day\tdate\tdeaths\tdeaths/day\tdeaths(proj)\tdeaths/day(proj)\n";
    double projectedTotal = 0;
    for (int i=1; i<LASTDAY; i++) {
        double projectedPerDay = (int)multiplier*norm(i,bestM,bestR);
        projectedTotal += projectedPerDay;
        std::cout << i << "\t" << dayToDate(i-1) << "\t";
        if (i<nPoints) {
            std::cout << testSeries[i-1] <<"\t" << diffs[i-1] << "\t";
        }
        else {
            std::cout << "\t \t ";
        }
        std::cout << (int)(projectedTotal+0.5) << "\t" << (int)(projectedPerDay+0.5) << "\n";
    }
    std::cout << "total=" << projectedTotal << "\n";
    return 0;
}
