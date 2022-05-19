#ifndef GENERATIONSIGNAL_H
#define GENERATIONSIGNAL_H

#endif // GENERATIONSIGNAL_H

#include <vector>
#include <complex>
#include <fstream>
#include "Generation.h"
#include "Methods.h"
#include "mainwindow.h"
using namespace std;

class GenerationSignal
{
public:
        int NumberOfSatellites=3;
        int NumberOfSources=3;

        //частота дискретизации
        double samplingFrequency = 300e6;
        //метка времени начала
        double startTimestamp = 0;
        //продолжительность
        double Duration = 1e-3;
        //начальная фаза
        double startPhase = 0;
        double nSamples = 0;
        //скорость передачи данных
        double Bitrate = 100e5;
        //дополнительный параметр, нечто, напоминающее modFreq
        double additionalParameter = 0;

        int sput;
        double clearance;
        double fr;
        double scale;
        int KolOtch;


        vector<vector<double>> spectr_long_signal;
        vector<vector<double>> BigMassOtrisovka;
        vector<vector<double>> BigMassOtshetX;
        vector<string> Stroka1;
        vector<string> Stroka2;
        vector<vector<int>> Time;
        vector<vector<int>> Freq;
        vector<bool> out;
        vector<string> str1;
        vector<string> str2;
        vector<string> str3;
        vector<double> timett;
        vector<double>freq;
        vector<string> str;
        vector<string> strr;
        vector<vector<int>> ind;
        vector<vector<int>> del;
        vector<Signal>* SummSignal = new vector<Signal>(NumberOfSatellites);
        vector<Signal2D>* A = new vector<Signal2D>(NumberOfSatellites);
        vector<vector<double>> massdoubleSignal;
        vector<Signal>* CorrelateSignal = new vector<Signal>(NumberOfSatellites);
         QVector<QString> listt;
         QVector<QString> listtVFN;
         vector<double> SummDelays;
         vector<double> DeltaK;


        vector<int> find_max_n(int n, Signal& sig, double clearance);
        vector<double> criteria(vector<vector<int>> delays, double clear);
        vector<double> Delays(int NumberOfSources, int NumberOfSatellites);
        double DoubleRand(double _max, double _min);
        vector<double> DelaysFrequency(int NumberOfSources, int NumberOfSatellites, double max_freq);
        void cutFrequency(Signal2D& insignal/*, Signal2D& outsignal*/);
        vector<int> find_max_time(int n, vector<double>& sig, vector<double>& keysX, vector<double>& keysY,bool time_freq, double clearance);
        vector<vector<int>> find_max_TF(int n, Signal2D& sig, double clearance);
        vector<double> criteria(vector<vector<int>> delays);
        void CreateVFN();
        void GenerateTwoLongSignal();
        vector<double> GenerateDeltaK(int NumberOfSatellites, double max_deltaK);

};
