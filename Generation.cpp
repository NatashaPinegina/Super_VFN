#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include "Generation.h"
#include <QtGlobal>

using namespace std::chrono;

double SignalGenerator::carriers[51] = { 969.,
												972.,
												975.,
												978.,
												981.,
												984.,
												987.,
												990.,
												993.,
												996.,
												999.,
												1002.,
												1005.,
												1008.,
												1053.,
												1056.,
												1059.,
												1062.,
												1065.,
												1113.,
												1116.,
												1119.,
												1122.,
												1125.,
												1128.,
												1131.,
												1134.,
												1137.,
												1140.,
												1143.,
												1146.,
												1149.,
												1152.,
												1155.,
												1158.,
												1161.,
												1164.,
												1167.,
												1170.,
												1173.,
												1176.,
												1179.,
												1182.,
												1185.,
												1188.,
												1191.,
												1194.,
												1197.,
												1200.,
												1203.,
												1209.
};


void  SignalGenerator::GenerateSignal(SignalType type, SignalGenerationParameters params, Signal& dst)
{
	switch (type)
	{
	case BPSK:	return GenerateBPSKSignal(params, dst);
	case MSK:	return GenerateMSKSignal(params, dst);
	//case FHSS:	return GenerateFHSSSignal(params, dst);
	case FHSS: return FHSSGenerate(params, dst);
	}
}

bool SignalGenerator::RandomBit(double low_chance)
{
	return (rand() < RAND_MAX * low_chance);
}

void SignalGenerator::GenerateBPSKSignal(SignalGenerationParameters params, Signal& ret_signal)
{
	ret_signal.signal.clear();
	ret_signal.sampling = params.sampling_frequency;
	ret_signal.duration = params.duration;
	ret_signal.timestamp = params.start_timestamp;
    qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	ret_signal.name = string("BPSK") + to_string(timestamp);
	ret_signal.description = string("BPSK signal: ") +
		string("Sampling = ") + to_string(params.sampling_frequency) + string("; ") +
		string("Bitrate = ") + to_string(params.bitrate) + string("; ") +
		string("Start phase = ") + to_string(params.start_phase) + string("; ") +
		string("Timestamp = ") + to_string(params.start_timestamp) + string("; ") +
		string("Duration = ") + to_string(params.duration) + string("; ") +
		string("Add Parameter = ") + to_string(params.additional_parameter) + string(";");

	int samples_in_bit = (int)(params.sampling_frequency / params.bitrate);
	double sampling_period = 1. / params.sampling_frequency;

	double dbl_index = params.start_timestamp;
	while (dbl_index < params.start_timestamp + params.duration)
	{
		double bit = RandomBit(0.5);
		for (int n_sample = 0; n_sample < samples_in_bit; n_sample++)
		{
			ret_signal.signal.push_back(
				complex<double>(
					cos(params.start_phase + bit * M_PI),
					sin(params.start_phase + bit * M_PI)
					));
			ret_signal.keys.push_back(dbl_index);
			dbl_index += sampling_period;
		}
	}
}
void SignalGenerator::GenerateMSKSignal(SignalGenerationParameters params, Signal& ret_signal)
{
	ret_signal.sampling = params.sampling_frequency;
	ret_signal.duration = params.duration;
	ret_signal.timestamp = params.start_timestamp;
    qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	ret_signal.name = string("MSK") + to_string(timestamp);
	ret_signal.description = string("MSK signal: ") +
		string("Sampling = ") + to_string(params.sampling_frequency) + string("; ") +
		string("Bitrate = ") + to_string(params.bitrate) + string("; ") +
		string("Start phase = ") + to_string(params.start_phase) + string("; ") +
		string("Timestamp = ") + to_string(params.start_timestamp) + string("; ") +
		string("Duration = ") + to_string(params.duration) + string("; ") +
		string("Add Parameter = ") + to_string(params.additional_parameter) + string(";");

	double sampling_period = 1. / params.sampling_frequency;


	int samples_in_bit = (int)(params.sampling_frequency / params.bitrate);

	double separated_bit_summ = 0;
	double dbl_index = params.start_timestamp;
	double t = 0;
	while (dbl_index < params.start_timestamp + params.duration)
	{
		double separated_bit = RandomBit(0.5) * 2. - 1;
		for (int n_sample = 0; n_sample < samples_in_bit; n_sample++)
		{
			ret_signal.signal.push_back(
				complex<double>(
					cos(M_PI / 2. * (separated_bit_summ + params.bitrate * separated_bit * t)),
					sin(M_PI / 2. * (separated_bit_summ + params.bitrate * separated_bit * t))
					));
			ret_signal.keys.push_back(dbl_index);
			t += sampling_period;
			dbl_index += sampling_period;
		}
		separated_bit_summ += separated_bit;
	}
}
double DoubleRrand(double _max, double _min)
{
	return _min + double(rand()) / RAND_MAX * (_max - _min);
}

//void SignalGenerator::FHSSGenerate(SignalGenerationParameters params, Signal& ret_signal)
//{
//	ret_signal.sampling = params.sampling_frequency;
//	ret_signal.duration = params.duration;
//	ret_signal.timestamp = params.start_timestamp;
//
//	__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
//	ret_signal.name = string("FHSS") + to_string(timestamp);
//	ret_signal.description = string("FHSS signal: ") +
//		string("Sampling = ") + to_string(params.sampling_frequency) + string("; ") +
//		string("Bitrate = ") + to_string(params.bitrate) + string("; ") +
//		string("Start phase = ") + to_string(params.start_phase) + string("; ") +
//		string("Timestamp = ") + to_string(params.start_timestamp) + string("; ") +
//		string("Duration = ") + to_string(params.duration) + string("; ") +
//		string("Add Parameter = ") + to_string(params.additional_parameter) + string(";");
//
//	double sampling_period = 1. / params.sampling_frequency;
//
//	double carriers[51] = {
//		969e6, 972e6, 975e6, 978e6, 981e6,
//		984e6, 987e6, 990e6, 993e6, 996e6,
//		999e6, 1002e6, 1005e6, 1008e6, 1053e6,
//		1056e6, 1059e6, 1062e6, 1065e6, 1113e6,
//		1116e6, 1119e6, 1122e6, 1125e6, 1128e6,
//		1131e6, 1134e6, 1137e6, 1140e6, 1143e6,
//		1146e6, 1149e6, 1152e6, 1155e6, 1158e6,
//		1161e6, 1164e6, 1167e6, 1170e6, 1173e6,
//		1176e6, 1179e6, 1182e6, 1185e6, 1188e6,
//		1191e6, 1194e6, 1197e6, 1200e6, 1203e6,
//		1209e6 };
//
//	double I = 0.0;
//	double Q = 0.0;
//	double phase = 0.0;
//	double phase_shift = 0.0;
//	int bits_per_word = 32;
//
//	int samples_in_bit = (int)(params.sampling_frequency / params.bitrate);
//	double freq_shift = 0;
//
//	double dbl_index = params.start_timestamp;
//	double t = 0;
//
//	double separated_bit_summ = 0;
//	GenerateMSKSignal(params, ret_signal);
//
//	ret_signal.keys.clear();
//	while (dbl_index < params.start_timestamp + params.duration)
//	{
//		//Генерируем по словам, сразу смещаем сигнал к 0 частоте.
//		double carrier_index = 50.0 * (double)rand() / (double)RAND_MAX;
//		freq_shift = (carriers[(int)carrier_index] - carriers[0]);
//		phase_shift = 0.0;
//		for (int j = 0; j < bits_per_word; j++)
//		{
//			//бит
//			double separated_bit = RandomBit(0.5) * 2. - 1;
//			//отсчеты в бите
//			for (int k = 0; k < samples_in_bit; k++)
//			{
//				phase_shift += 2 * M_PI * freq_shift / params.sampling_frequency;
//				phase = k + j * samples_in_bit + dbl_index * bits_per_word * samples_in_bit;
//				//phase = separated_bit_summ + params.bitrate * separated_bit * t;
//				if (separated_bit) phase_shift += M_PI * params.bitrate / (2.0 * params.sampling_frequency);
//				else phase_shift -= M_PI * params.bitrate / (2.0 * params.sampling_frequency);
//
//				/*if (phase_shift > 2.0 * M_PI) phase_shift -= 2.0 * M_PI;
//				if (phase_shift < -2.0 * M_PI) phase_shift += 2.0 * M_PI;
//
//				if (phase > (2.0 * M_PI)) phase -= (2.0 * M_PI);
//				if (phase < (-2.0 * M_PI)) phase += (2.0 * M_PI);*/
//
//
//				I = cos(phase);
//				Q = sin(phase);
//
//				t += sampling_period;
//				dbl_index += sampling_period;
//
//				ret_signal.signal.push_back(
//					complex <double>(
//						I * cos(phase_shift) - Q * sin(phase_shift),
//						I * sin(phase_shift) + Q * cos(phase_shift)
//						)
//				);
//
//				ret_signal.signal[phase] = complex < float >((float)(
//					ret_signal.signal[phase].real() * cos(phase_shift) -
//					ret_signal.signal[phase].imag() * sin(phase_shift)),
//					(float)(
//						ret_signal.signal[phase].imag() * cos(phase_shift) +
//						ret_signal.signal[phase].real() * sin(phase_shift)));
//
//				/*buffer[phase]= complex < float >((float)(
//					ret_signal.signal[phase].real() * cos(phase_shift) -
//					ret_signal.signal[phase].imag() * sin(phase_shift)),
//					(float)(
//						ret_signal.signal[phase].imag() * cos(phase_shift) +
//						ret_signal.signal[phase].real() * sin(phase_shift)));*/
//
//				ret_signal.keys.push_back(dbl_index);
//			}
//			separated_bit_summ += separated_bit;
//		}
//	}
//}


void SignalGenerator::MSKManipulation(SignalGenerationParameters params, const vector <bool>& data, iq_signal& buffer)
{
	buffer.clear();
	int samples_per_bit = (int)(params.sampling_frequency / params.bitrate);

	double phase = 0.0;
	int bit_sum = 0;

	for (int i = 0; i < data.size(); i++)
	{
		phase = (double)bit_sum;
		double separated_bit = RandomBit(0.5) * 2. - 1;
		for (int j = 0; j < samples_per_bit; j++)
		{
			if (separated_bit) phase += params.bitrate / params.sampling_frequency;
			else phase -= params.bitrate / params.sampling_frequency;

			buffer.push_back(complex<float>(cos(M_PI * phase / 2.0), sin(M_PI * phase / 2.0)));
		}

		if (data[i]) bit_sum++;
		else bit_sum--;
	}

}

void SignalGenerator::FHSSGenerate(SignalGenerationParameters params, Signal& ret_signal)
{
	ret_signal.sampling = params.sampling_frequency;
	ret_signal.duration = params.duration;
	ret_signal.timestamp = params.start_timestamp;

    qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
	ret_signal.name = string("FHSS") + to_string(timestamp);
	ret_signal.description = string("FHSS signal: ") +
		string("Sampling = ") + to_string(params.sampling_frequency) + string("; ") +
		string("Bitrate = ") + to_string(params.bitrate) + string("; ") +
		string("Start phase = ") + to_string(params.start_phase) + string("; ") +
		string("Timestamp = ") + to_string(params.start_timestamp) + string("; ") +
		string("Duration = ") + to_string(params.duration) + string("; ") +
		string("Add Parameter = ") + to_string(params.additional_parameter) + string(";");

	double sampling_period = 1. / params.sampling_frequency;

	/*double carriers[51] = {
		969e6, 972e6, 975e6, 978e6, 981e6,
		984e6, 987e6, 990e6, 993e6, 996e6,
		999e6, 1002e6, 1005e6, 1008e6, 1053e6,
		1056e6, 1059e6, 1062e6, 1065e6, 1113e6,
		1116e6, 1119e6, 1122e6, 1125e6, 1128e6,
		1131e6, 1134e6, 1137e6, 1140e6, 1143e6,
		1146e6, 1149e6, 1152e6, 1155e6, 1158e6,
		1161e6, 1164e6, 1167e6, 1170e6, 1173e6,
		1176e6, 1179e6, 1182e6, 1185e6, 1188e6,
		1191e6, 1194e6, 1197e6, 1200e6, 1203e6,
		1209e6 };*/

	double carriers[51] = {
		969e6, 972e6, 975e6, 978e6, 981e6,
		984e6, 987e6, 990e6, 993e6, 996e6,
		999e6, 1002e6, 1005e6, 1008e6, 1053e6,
		1056e6, 1059e6, 1062e6, 1065e6, 1113e6,
		1116e6, 1119e6, 1122e6, 1125e6, 1128e6,
		1131e6, 1134e6, 1137e6, 1140e6, 1143e6,
		1146e6, 1149e6, 1152e6, 1155e6, 1158e6,
		1161e6, 1164e6, 1167e6, 1170e6, 1173e6,
		1176e6, 1179e6, 1182e6, 1185e6, 1188e6,
		1191e6, 1194e6, 1197e6, 1200e6, 1203e6,
		1209e6 };

	double I = 0.0;
	double Q = 0.0;
	double phase = 0.0;
	double phase_shift = 0.0;
	int bits_per_word = 1;

	int samples_in_bit = (int)(params.sampling_frequency / params.bitrate);
	double freq_shift = 0;

	double dbl_index = params.start_timestamp;
	double t = 0;

	double separated_bit_summ = 0;
	auto comp_j = complex<double>(0, 1);



	double countOnBit = params.sampling_frequency / params.additional_parameter;

	GenerateBPSKSignal(params, ret_signal);

	int carrier_index = 0;

	for (int i = 0; i < ret_signal.signal.size(); i++)
	{
		if (i % (bits_per_word * samples_in_bit) == 0)
		{
			carrier_index = 51.0 * (double)rand() / RAND_MAX;
			freq_shift = (carriers[carrier_index] - carriers[0]);
		}
		
		ret_signal.signal[i] = ret_signal.signal[i] * exp(comp_j * 2.0*M_PI*freq_shift * ret_signal.keys[i]);
	}


	
	//ret_signal.keys.clear();
	//while (dbl_index < params.start_timestamp + params.duration)
	//{
	///*for(int i=0;i<count_words;i++)
	//{*/
	//	//Генерируем по словам, сразу смещаем сигнал к 0 частоте.
	//	double carrier_index = 50.0 * (double)rand() / (double)RAND_MAX;
	//	freq_shift = (carriers[(int)carrier_index] - carriers[0]);
	//	phase_shift = 0.0;
	//	for (int j = 0; j < bits_per_word; j++)
	//	{
	//		//бит
	//		double separated_bit = RandomBit(0.5) * 2. - 1;
	//		//отсчеты в бите
	//		for (int k = 0; k < samples_in_bit; k++)
	//		{
	//			phase_shift += 2 * M_PI * freq_shift / params.sampling_frequency;
	//			if (separated_bit) phase_shift += M_PI * params.bitrate / (2.0 * params.sampling_frequency);
	//			else phase_shift -= M_PI * params.bitrate / (2.0 * params.sampling_frequency);

	//			if (phase_shift > 2.0 * M_PI) phase_shift -= 2.0 * M_PI;
	//			if (phase_shift < -2.0 * M_PI) phase_shift += 2.0 * M_PI;

	//			phase = k + j * samples_in_bit + dbl_index * bits_per_word * samples_in_bit;
	//			//phase = separated_bit_summ + params.bitrate * separated_bit * t;
	//			
	//			if (phase > (2.0 * M_PI)) phase -= (2.0 * M_PI);
	//			if (phase < (-2.0 * M_PI)) phase += (2.0 * M_PI);

	//			I = cos(phase);
	//			Q = sin(phase);
	//			
	//			// ret_signal.signal[phase]
	//			t += sampling_period;
	//			dbl_index += sampling_period;

	//			ret_signal.signal[phase] = complex<double>(I, Q);// *exp(comp_j * phase_shift);
	//			ret_signal.keys.push_back(dbl_index);
	//		}
	//		separated_bit_summ += separated_bit;
	//	}
	//}
}

