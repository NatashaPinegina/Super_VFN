#pragma once
#include <vector>
#include <complex>


using namespace std;

enum SignalType
{
	//двоичная фазовая модуляция
	BPSK,
	//четырехпозиционная фазовая я
	QPSK,
	//частотная модуляция
	MSK,
	//мультиплексирование с ортогональным частотным разделением каналов
	OFDM,
	//Псевдослучайная перестройка рабочей частоты
	FHSS,
};

struct Signal
{
	vector<complex<double>> signal;
	vector<double> keys;

	string name;
	//описание
	string description;

	double sampling;
	double timestamp;
	double duration;
};

struct Signal2D
{
	vector<vector<complex<double>>> signal;
	vector<double> keysX;
	vector<double> keysY;

	string name;
	string description;

	double sampling;
};