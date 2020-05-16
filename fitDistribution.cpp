//
//  fitDistribution.cpp
//  fitDistribution
//
//  Created by Stephen Stein on 4/2/20.
//  Copyright © 2020 Stephen Stein. All rights reserved.
//
// NAME
//      fitDistribution - finds the cumulative distribution that best fits a set of points
//
// SYNOPSIS
//      fitDistribution inputFile
//
// DESCRIPTION
//      fitDistribution takes a set of points and computes the best-fitting cumulative distribution.
//
//      The input file to fitDistribution is a text file with a date on the first line followed by any number of data points, one per line
//
//      After reading the date and the data points, the program computes an array of differences between successive data points.
//
//      The program then cycles through possible means and standard deviations, computing a score for each which represents how well that
//      mean and standard deviation fit the input points. The lowest score represents the best fitting mean and standard deviation.
//
//      The score is computed as follows:
//      1. Find the multiplier.  Choose a difference value ("D") and note its index ("X") in the array of differences. Compute the value "V"
//          of the standard normal distribution at that index (NORMDIST(X,Mean,StdDev)) for the given mean and stadard deviation. Compute
//          the multiplier M = D/V.
//      2. Compute a projected value array by summing the scaled values of the normal distribution M*NORMDIST(Index,Mean,StdDev)
//          for Index ranging from 0 to the length of the input array.
//      3. Find the difference between the input array and the projected array and square it for each index in the input array. Keep a
//          running total of these squares.
//      4. The total from step 3 is the "score"
//
//      NOTE - a smoothed value "D" in step 1 is the maximum of the differences array, smoothed by taking the weighted average of the
//          preceding and succeeding values (25% each) and the max point itself (50%). If the maximum of the differences array is the
//          last value, use the next-to-last value, smoothed, as "D" and decrement X by 1.
//
// OUTPUT
//      The mean, standard deviation, score and multiplier for the best fit
//      The best fit distribution, in columns separated by tabs. Each line contains:
//          Index: the index of the point in the input array
//          Date: "Index"-1 days after the input date. The input date is day 1.
//          Actual Cumulative Total: The given input array
//          Actual difference per day: The difference between the successive elements of the input array
//          Projected Cumulative Total: The cumulative total of the scaled normal distribution on that day
//          Projected difference per day: The scaled value of the normal distribution on that day
//
// RETURN VALUE
//      0 on success, non-zero on failure
//
#define VERBOSE 0
#include <iostream>
#include <fstream>      // std::ifstream
#include <cmath>

// normal distribution utility routines
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

// date computation utility routines
int daysInMonth[] = {0,31,29,31,30,31,30,31,31,30,31,30,31};
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

std::string datestring;    // the date string on line 1 of the input
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

double *testSeries;
double *diffs;
double *smoothedDiffs;      // computed but not used (yet?)
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
    // note: smoothedDiffs is computed but not used (yet?)
    smoothedDiffs = (double*)malloc(nPoints*sizeof(double));
    std::cout << nPoints << "\n";
    diffs[0] = testSeries[0];
    smoothedDiffs[0] = 0.75*diffs[0] + 0.25*diffs[1];
    for (int i=1;i<nPoints-1;i++) {
        diffs[i] = testSeries[i] - testSeries[i-1];
        if (i==nPoints-1) {
            smoothedDiffs[i] = 0.5*diffs[i-1] + 0.5*diffs[i];
        }
        else {
            smoothedDiffs[i] = (diffs[i-1] + diffs[i] + diffs[i+1])/3.0;
        }
        std::cout << i << "\t" << diffs[i] << "\n";
        if (diffs[i]>maxDiff) {
            maxDiff = diffs[i];
            maxDiffPt = i;
        }
    }
    // if the last point is a max, move the ref point back one and smooth the last 3 values
    if (maxDiffPt==nPoints-2) {
        maxDiff = (diffs[maxDiffPt-2] + diffs[maxDiffPt-1] + diffs[maxDiffPt])/3.0;
        maxDiffPt -= 1;
    }
    else {
        // smooth the max point
        maxDiff = (diffs[maxDiffPt-1] + diffs[maxDiffPt] + diffs[maxDiffPt+1])/3.0;
    }
    double bestScore = std::numeric_limits<double>::max();
    double bestM = 0;
    double bestStdDev = 0;
    for (int m = nPoints+60; m>0; m--) {
        double bestScoreForMean = std::numeric_limits<double>::max();
        double bestStdDevForMean = 0;
        for (int r=2; r<15; r++) {
            double score = scoreGuess(m, r);
#if VERBOSE
            std::cout << "mean="  << m << ", dev=" << r << ", score=" << score << "\n";
#endif
            if (score<bestScoreForMean) {
                bestScoreForMean = score;
                bestStdDevForMean = r;
                if (score<bestScore) {
                    bestScore = score;
                    bestM = m;
                    bestStdDev = r;
                    std::cout << "BEST: mean=" << bestM << " dev=" << bestStdDev << " score=" << bestScore << "\n";
                }
            }
        }
        std::cout << "BEST score for mean=" << m << " dev=" << bestStdDevForMean << " score=" << bestScoreForMean << "\n";
    }
    std::cout << "BEST: mean=" << bestM << " dev=" << bestStdDev << " score=" << bestScore << "\n";
    double downScore;
    double upScore;
    downScore = scoreGuess(bestM,bestStdDev-0.1);
    if (downScore < bestScore) {
        while (downScore<bestScore) {
            bestStdDev -= 0.1;
            bestScore = downScore;
            downScore = scoreGuess(bestM,bestStdDev-0.1);
        }
    }
    else {
        upScore = scoreGuess(bestM,bestStdDev+0.1);
        while (upScore<bestScore) {
            bestStdDev += 0.1;
            bestScore = upScore;
            upScore = scoreGuess(bestM,bestStdDev+0.1);
        }
    }
    bestScore = scoreGuess(bestM,bestStdDev);
    std::cout << "BEST: mean=" << bestM << " dev=" << bestStdDev << " score=" << bestScore << " multiplier=" << multiplier << "\n";
    std::cout << "day\tdate\tdeaths\tdeaths/day\tdeaths(proj)\tdeaths/day(proj)\n";
    double projectedTotal = 0;
    for (int i=1; i<LASTDAY; i++) {
        double projectedPerDay = (int)multiplier*norm(i,bestM,bestStdDev);
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
