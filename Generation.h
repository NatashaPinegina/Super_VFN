#pragma once

#include "Signal.h"

typedef vector < complex < double > > iq_signal;

struct SignalGenerationParameters
{
	//метка времени начала
	double start_timestamp;
	//начальная фаза
	double start_phase;
	//продолжительность
	double duration;

	double n_samples;
	//частота дискретизации
	double sampling_frequency;
	//скорость передачи данных
	double bitrate;
	//дополнительный параметр
	double additional_parameter;
};

class SignalGenerator
{
public:
	void GenerateSignal(SignalType type, SignalGenerationParameters params, Signal& dst);
	static double carriers[51];
	iq_signal buffer;
private:
	void GenerateBPSKSignal(SignalGenerationParameters params, Signal& dst);
	void GenerateMSKSignal(SignalGenerationParameters params, Signal& dst);
	void GenerateFHSSSignal(SignalGenerationParameters params, Signal& dst);
	bool RandomBit(double low_chance);
	void FHSSGenerate(SignalGenerationParameters params, Signal& ret_signal);

	void MSKManipulation(SignalGenerationParameters params, const vector <bool>& data, iq_signal& buffer);
};