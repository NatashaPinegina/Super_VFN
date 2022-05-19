#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <complex>
#include "Signal.h"
#include <chrono>
#include <algorithm>
#include <numeric>
#include <QtGlobal>

using namespace std;
using namespace std::chrono;


class Methods
{
public:

	static void correlate(Signal &signal1, Signal &signal2, Signal &correlation)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		correlation.signal.clear();
		correlation.keys.clear();
		correlation.name = string("CORR") + "(" + to_string(timestamp) + "):(" + signal1.name + "|" + signal2.name + ")";
		correlation.description = "Correlation of " + signal1.name + " and " + signal2.name;
		correlation.sampling = signal1.sampling;
		correlation.timestamp = 0;
		correlation.duration = signal1.duration + signal2.duration;
		
		for (unsigned int n = 0; n < (signal1.signal.size() - signal2.signal.size()); n++)
		{
			complex<double> counter = complex<double>(0, 0);
			for (unsigned int m = 0; m < signal2.signal.size(); m++)
			{
				if ((m + n) >= 0 && (m + n) < signal1.signal.size())
				{
					counter += signal2.signal[m] * conj(signal1.signal[m + n]);
				}
				else continue;
			}
			counter /= (double)signal1.signal.size();
			correlation.signal.push_back(counter);
			correlation.keys.push_back(n);
		}
	}
	static void Newcorrelate(Signal& signal1, Signal& signal2, Signal& correlation)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		correlation.signal.clear();
		correlation.keys.clear();
		correlation.name = string("CORR") + "(" + to_string(timestamp) + "):(" + signal1.name + "|" + signal2.name + ")";
		correlation.description = "Correlation of " + signal1.name + " and " + signal2.name;
		correlation.sampling = signal1.sampling;
		correlation.timestamp = 0;
		correlation.duration = signal1.duration + signal2.duration;
		for (int n = -(int)(signal2.signal.size()); n <= (int)(signal1.signal.size()-1); n++)
		{
			complex<double> counter = complex<double>(0, 0);
			for (unsigned int m = 0; m < signal2.signal.size(); m++)
			{
				if ((m + n) >= 0 && (m + n) < signal1.signal.size())
				{
					counter += signal2.signal[m] * conj(signal1.signal[m + n]);
				}
				else continue;
			}
			counter /=  (double)signal2.signal.size();//cnt;//
			correlation.signal.push_back(counter);
			correlation.keys.push_back(n);
		}
	}
	static void correlate(vector<complex<double>> &base_signal, vector<complex<double>> &analyzed_signal, vector<complex<double>> &correlation)
	{
		correlation.clear();
		for (int n = -(int)(analyzed_signal.size()); n <= (int)(base_signal.size()); n++)
		{
			complex<double> counter = complex<double>(0, 0);
			for (unsigned int m = 0; m < analyzed_signal.size(); m++)
			{
				if ((m + n) >= 0 && (m + n) < base_signal.size())
					counter += analyzed_signal[m] * conj(base_signal[m + n]);
				else continue;
			}
			counter /= (double)analyzed_signal.size();
			correlation.push_back(counter);
		}
	}
	static void FFT(Signal &in, Signal &out, int direction)
	{
		unsigned int pts = 2;
		while (in.signal.size() > pts)
		{
			pts *= 2;
		}

		int pts_to_add = pts - in.signal.size();

		out = in;
		out.name = "FFT" + in.name;
		out.signal = in.signal;

		if (pts_to_add >= 0)
		{
			for (int ttt = 0; ttt < pts_to_add; ttt++)
			{
				out.signal.push_back(complex<double>(0, 0));
			}
		}
		else
		{
			for (int ttt = 0; ttt < -pts_to_add; ttt++)
			{
				out.signal.pop_back();
			}
		}

		int n = out.signal.size();

		int i, j, istep;
		int m, mmax;
		double r, r1, theta, w_r, w_i, temp_r, temp_i;

		r = M_PI*direction;
		j = 0;

		for (i = 0; i<n; i++)
		{
			if (i<j)
			{
				temp_r = out.signal[j].real();
				temp_i = out.signal[j].imag();
				out.signal[j] = out.signal[i];
				out.signal[i] = complex<double>(temp_r, temp_i);
			}
			m = n >> 1;
			while (j >= m)
			{
				j -= m;
				m = (m + 1) / 2;
			}
			j += m;
		}
		mmax = 1;
		while (mmax<n)
		{
			istep = mmax << 1;
			r1 = r / (double)mmax;
			for (m = 0; m<mmax; m++)
			{
				theta = r1*m;
				w_r = (double)cos((double)theta);
				w_i = (double)sin((double)theta);
				for (i = m; i<n; i += istep)
				{
					j = i + mmax;
					temp_r = w_r*out.signal[j].real() - w_i*out.signal[j].imag();
					temp_i = w_r*out.signal[j].imag() + w_i*out.signal[j].real();
					out.signal[j] = complex<double>((out.signal[i].real() - temp_r), (out.signal[i].imag() - temp_i));
					out.signal[i] += complex<double>(temp_r, temp_i);
				}
			}
			mmax = istep;
		}
		if (direction > 0)
		{
			double sqn = sqrt(n);
			for (i = 0; i<n; i++)
			{
				out.signal[i] /= sqn;
			}
		}
		out.keys.clear();
		for (unsigned int i = 0; i < out.signal.size(); i++)
		{
			out.keys.push_back(i);
		}
	}
	static void FFT(vector<complex<double>> & in, vector<complex<double>> & out, int direction)
	{
		out = in;
		unsigned int pts = 2;
		while (out.size() > pts)
		{
			pts *= 2;
		}

		int pts_to_add = pts - out.size();

		for (int ttt = 0; ttt < pts_to_add; ttt++)
		{
			out.push_back(complex<double>(0, 0));
		}
		int n = out.size();

		int i, j, istep;
		int m, mmax;
		double r, r1, theta, w_r, w_i, temp_r, temp_i;

		r = M_PI*direction;
		j = 0;

		for (i = 0; i<n; i++)
		{
			if (i<j)
			{
				temp_r = out[j].real();
				temp_i = out[j].imag();
				out[j] = out[i];
				out[i] = complex<double>(temp_r, temp_i);
			}
			m = n >> 1;
			while (j >= m)
			{
				j -= m;
				m = (m + 1) / 2;
			}
			j += m;
		}
		mmax = 1;
		while (mmax<n)
		{
			istep = mmax << 1;
			r1 = r / (double)mmax;
			for (m = 0; m<mmax; m++)
			{
				theta = r1*m;
				w_r = (double)cos((double)theta);
				w_i = (double)sin((double)theta);
				for (i = m; i<n; i += istep)
				{
					j = i + mmax;
					temp_r = w_r*out[j].real() - w_i*out[j].imag();
					temp_i = w_r*out[j].imag() + w_i*out[j].real();
					out[j] = complex<double>((out[i].real() - temp_r), (out[i].imag() - temp_i));
					out[i] += complex<double>(temp_r, temp_i);
				}
			}
			mmax = istep;
		}
		if (direction > 0)
		{
			double sqn = sqrt(n);
			for (i = 0; i<n; i++)
			{
				out[i] /= sqn;
			}
		}
	}

	static void newFFT(vector<complex<double>>& in, int direction)
	{
		//out = in;
		unsigned int pts = 2;
		while (in.size() > pts)
		{
			pts *= 2;
		}

		int pts_to_add = pts - in.size();

		for (int ttt = 0; ttt < pts_to_add; ttt++)
		{
			in.push_back(complex<double>(0, 0));
		}
		int n = in.size();

		int i, j, istep;
		int m, mmax;
		double r, r1, theta, w_r, w_i, temp_r, temp_i;

		r = M_PI * direction;
		j = 0;

		for (i = 0; i < n; i++)
		{
			if (i < j)
			{
				temp_r = in[j].real();
				temp_i = in[j].imag();
				in[j] = in[i];
				in[i] = complex<double>(temp_r, temp_i);
			}
			m = n >> 1;
			while (j >= m)
			{
				j -= m;
				m = (m + 1) / 2;
			}
			j += m;
		}
		mmax = 1;
		while (mmax < n)
		{
			istep = mmax << 1;
			r1 = r / (double)mmax;
			for (m = 0; m < mmax; m++)
			{
				theta = r1 * m;
				w_r = (double)cos((double)theta);
				w_i = (double)sin((double)theta);
				for (i = m; i < n; i += istep)
				{
					j = i + mmax;
					temp_r = w_r * in[j].real() - w_i * in[j].imag();
					temp_i = w_r * in[j].imag() + w_i * in[j].real();
					in[j] = complex<double>((in[i].real() - temp_r), (in[i].imag() - temp_i));
					in[i] += complex<double>(temp_r, temp_i);
				}
			}
			mmax = istep;
		}
		if (direction > 0)
		{
			double sqn = sqrt(n);
			for (i = 0; i < n; i++)
			{
				in[i] /= sqn;
			}
		}
		
	}

	static void crivoeFFT(vector<complex<double>>& in, int direction)
	{
		//out = in;
		unsigned int pts = 2;
		while (in.size() > pts)
		{
			pts *= 2;
		}

		int pts_to_add = pts - in.size();

		for (int ttt = 0; ttt < pts_to_add; ttt++)
		{
			in.push_back(complex<double>(0, 0));
		}
		int n = in.size();

		int i, j, istep;
		int m, mmax;
		double r, r1, theta, w_r, w_i, temp_r, temp_i;

		r = M_PI * direction;
		j = 0;

		for (i = 0; i < n; i++)
		{
			if (i < j)
			{
				temp_r = in[j].real();
				temp_i = in[j].imag();
				in[j] = in[i];
				in[i] = complex<double>(temp_r, temp_i);
			}
			m = n >> 1;
			while (j >= m)
			{
				j -= m;
				m = (m + 1) / 2;
			}
			j += m;
		}
		mmax = 1;
		double deltaK = 1e-5;// Коэф-нт масштабирования спектра. 
		double deltaT = 1;// Шаг по времени.
		while (mmax < n)
		{
			istep = mmax << 1;
			//r1 = r / (double)mmax;
			r1 = r / n;
			for (m = 0; m < mmax; m++)
			{
				theta = 2 * r1 * m * deltaT * deltaK;
				w_r = (double)cos((double)theta);
				w_i = (double)sin((double)theta);
				for (i = m; i < n; i += istep)
				{
					j = i + mmax;
					temp_r = w_r * in[j].real() - w_i * in[j].imag();
					temp_i = w_r * in[j].imag() + w_i * in[j].real();
					in[j] = complex<double>((in[i].real() - temp_r), (in[i].imag() - temp_i));
					in[i] += complex<double>(temp_r, temp_i);
				}
			}
			mmax = istep;
		}
		if (direction > 0)
		{
			double sqn = sqrt(n);
			for (i = 0; i < n; i++)
			{
				in[i] /= sqn;
			}
		}

	}

	static void transformSignal(Signal &base_signal, double delay, double duration, double fshift, double scale, double SNR, Signal &ret_sig)
	{
		ret_sig.signal.clear();
		ret_sig.keys.clear();
		ret_sig.sampling = base_signal.sampling;
		ret_sig.timestamp = 0;
		ret_sig.duration = duration;
		ret_sig.name = string("TS") + "(" + base_signal.name + ")";
		ret_sig.description = "Transformed signal:" + base_signal.name +
			"(dur):" + to_string(duration) +
			"(delay):" + to_string(delay) +
			"(fshift):" + to_string(fshift) +
			"(td):" + to_string(delay);
		double time_shift = delay*base_signal.sampling;
		int size_of_sample = (int)(duration*base_signal.sampling);
		if (time_shift < base_signal.signal.size() && (time_shift + size_of_sample) < base_signal.signal.size())
		{
			for (int i = 0; i < (size_of_sample); i++)
			{
				ret_sig.signal.push_back(
					base_signal.signal[(unsigned int)((double)i / scale + time_shift)] *
					exp(complex<double>(0, 2. * M_PI * fshift * (1. / ret_sig.sampling) * (double)i)));
				ret_sig.keys.push_back((1. / ret_sig.sampling)*(double)i);
			}
		}
		addNoise(ret_sig.signal, SNR);
	}
	static void addNoise(vector<complex<double>> &buf, double SNR)
	{
		//srand(0);
		double energy = std::accumulate(buf.begin(), buf.end(), 0.0,
			[&](double a, complex <double> b)
		{
			return a + b.real() * b.real() + b.imag() * b.imag();
		});

		vector <complex<double>> noise(buf.size(), 0.0);
		vector <double> ampl(buf.size(), 0.0);
		std::transform(ampl.begin(), ampl.end(), ampl.begin(),
			[&](double a)
		{
			int count_rep = 20;
			a = 0.0;
			for (int rep = 0; rep < count_rep; rep++)
			{
				a += ((double)rand() / RAND_MAX) * 2.0 - 1.0;
			}
			a /= count_rep;
			return a;
		});
		std::transform(noise.begin(), noise.end(), ampl.begin(), noise.begin(),
			[&](complex<double> a, double ampl)
		{
			double fi = (double)rand() / RAND_MAX * 2.0 * M_PI;
			return complex <double>(ampl * cos(fi), ampl * sin(fi));
		});
		double noise_energy = std::accumulate(noise.begin(), noise.end(), 0.0,
			[&](double a, complex <double> b)
		{
			return a + b.real() * b.real() + b.imag() * b.imag();
		});
		double norm_coef = energy * pow(10.0, -SNR / 10.0) / noise_energy;
		std::transform(buf.begin(), buf.end(), noise.begin(), buf.begin(),
			[&](complex <double> a, complex <double> b)
		{
			return complex <double>(a.real() + sqrt(norm_coef) * b.real(),
				a.imag() + sqrt(norm_coef) * b.imag());
		});
	}
	static void AmbiguityFunction(Signal &insignal1, Signal &insignal2, Signal2D &outsignal)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		outsignal.name = "Ambiguity function" + to_string(timestamp);
		outsignal.description = "Ambiguity function of " + insignal1.name + " and " + insignal2.name;
		outsignal.signal.resize(insignal1.signal.size() - insignal2.signal.size());
		for (unsigned int i = 0; i < insignal1.signal.size() - insignal2.signal.size(); i++)
		{
			vector<complex<double>> mass;
			vector<complex<double>> fmass;
			mass.resize(insignal2.signal.size());

			for (unsigned int j = 0; j < insignal2.signal.size(); j++)
			{
				mass[j] = insignal2.signal[j] * conj(insignal1.signal[j + i]);
				mass[j] /= insignal2.signal.size();
			}

			FFT(mass, fmass, -1);

			outsignal.signal[i] = fmass;
			outsignal.keysX.push_back(i);
		}
		for (unsigned int i = 0; i < insignal2.signal.size(); i++)
		{
			outsignal.keysY.push_back(i);
		}
	}

	static void NewAmbiguityFunction(Signal& insignal1, Signal& insignal2, Signal2D& outsignal, double &scale)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		outsignal.name = "Ambiguity function" + to_string(timestamp);
		outsignal.description = "Ambiguity function of " + insignal1.name + " and " + insignal2.name;
		//outsignal.signal.resize(insignal1.signal.size() - insignal2.signal.size());
		outsignal.signal.resize(insignal1.signal.size() + insignal2.signal.size()/*insignal1.duration + insignal2.duration*/);
		outsignal.sampling= insignal1.sampling;
		int iter = 0;
		//vector<complex<double>> mass;
		vector<complex<double>>* mass = new vector<complex<double>>;
		(*mass).resize(insignal2.signal.size()); //insignal1.signal.size() + 
		//vector<complex<double>> fmass;
		vector<complex<double>>* fmass = new vector<complex<double>>;
		for(int i = -(int)(insignal2.signal.size()); i <= (int)(insignal1.signal.size() - 1); i++)
		{
			(*mass).clear();
			(*mass).resize(insignal2.signal.size());
			for (unsigned int j = 0; j < insignal2.signal.size(); j++)
			{
				if ((j + i) >= 0 && (j + i) < insignal1.signal.size())
				{
					(*mass)[j] = insignal2.signal[j] * conj(insignal1.signal[j + i]);
					(*mass)[j] /= insignal2.signal.size();
				}
				else continue;
			}
			//FFT(mass, fmass, -1);
			newFFT((*mass), -1);
			outsignal.signal[iter]=(*mass);
			outsignal.keysY.push_back(i);//ключи времени
			iter++;
		}
		scale = outsignal.sampling / (*mass).size();
		for (int i = 0; i < (*mass).size(); i++)
		{
			outsignal.keysX.push_back(i*scale);
		}
		delete mass;
		delete fmass;
	}
	static double Interpol(double y1, double y2,double t, double deltaT, double koeff)
	{
		double a = (y2 - y1) / deltaT;
		double b = y1 - a * t;
		return a * koeff + b;
	}

    static void SuperVFN(Signal& insignal1, Signal& insignal2, Signal2D& outsignal, double& scale,double deltaK)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		outsignal.name = "Ambiguity function" + to_string(timestamp);
		outsignal.description = "Ambiguity function of " + insignal1.name + " and " + insignal2.name;
		//outsignal.signal.resize(insignal1.signal.size() - insignal2.signal.size());
		outsignal.signal.resize(insignal1.signal.size() + insignal2.signal.size()/*insignal1.duration + insignal2.duration*/);
		outsignal.sampling = insignal1.sampling;
		int iter = 0;
		//vector<complex<double>> mass;
		vector<complex<double>>* mass = new vector<complex<double>>;
		(*mass).resize(insignal2.signal.size()); //insignal1.signal.size() + 
		//vector<complex<double>> fmass;
		vector<complex<double>>* fmass = new vector<complex<double>>;

        //double deltaK = 1e-5;// Коэф-нт масштабирования спектра.
		double deltaT = 1;// Шаг по времени.
		double koeff = 0;// Шаг для интерполяции.

		double re = 0, im=0;

		vector<complex<double>> s1(insignal1.signal.size());
		vector<complex<double>> s2(insignal2.signal.size());
		for (int i = 0; i < insignal1.signal.size()-1; i++)
		{
			koeff = i * deltaT * deltaK;

			re = Interpol(insignal1.signal[i].real(), insignal1.signal[i+(int)deltaT].real(), i, deltaT, koeff);
			im = Interpol(insignal1.signal[i].imag(), insignal1.signal[i + (int)deltaT].imag(), i, deltaT, koeff);
			complex<double> S1(re, im);
			s1[i] = S1;

			re = Interpol(insignal2.signal[i].real(), insignal2.signal[i + (int)deltaT].real(), i, deltaT, koeff);
			im = Interpol(insignal2.signal[i].imag(), insignal2.signal[i + (int)deltaT].imag(), i, deltaT, koeff);
			complex<double> S2(re, im);
			s2[i] = S2;
		}


		for (int i = -(int)(insignal2.signal.size()); i <= (int)(insignal1.signal.size() - 1); i++)
		{
			(*mass).clear();
			(*mass).resize(insignal2.signal.size());
			for (unsigned int j = 0; j < s1.size(); j++)
			{

				if ((j + i) >= 0 && (j + i) < s1.size())
				{
					(*mass)[j] =s1[j] * conj(s2[j + i]);
					(*mass)[j] /= s1.size();
				}
				else continue;
			}
			//FFT(mass, fmass, -1);
			crivoeFFT((*mass), -1);
			outsignal.signal[iter] = (*mass);
			outsignal.keysY.push_back(i);//ключи времени
			iter++;
		}


		scale = outsignal.sampling / (*mass).size();
		for (int i = 0; i < (*mass).size(); i++)
		{
			outsignal.keysX.push_back(i * scale);
		}
		delete mass;
		delete fmass;
	}


	static void SuperNewAmbiguityFunction(Signal& insignal1, Signal& insignal2, Signal& outsignal, double& scale)
	{
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		outsignal.name = "Ambiguity function" + to_string(timestamp);
		outsignal.description = "Ambiguity function of " + insignal1.name + " and " + insignal2.name;
		//outsignal.signal.resize(insignal1.signal.size() - insignal2.signal.size());
		outsignal.signal.resize(insignal1.signal.size() + insignal2.signal.size()/*insignal1.duration + insignal2.duration*/);
		outsignal.sampling = insignal1.sampling;
		outsignal.keys.clear();
		int iter = 0;
		//vector<complex<double>> mass;
		vector<complex<double>>* mass = new vector<complex<double>>();
		(*mass).resize(insignal2.signal.size()); //insignal1.signal.size() + 
		//vector<complex<double>> fmass;
		vector<complex<double>>* fmass = new vector<complex<double>>;
		for (int i = -(int)(insignal2.signal.size()); i <= (int)(insignal1.signal.size() - 1); i++)
		{
			(*mass).clear();
			(*mass).resize(insignal2.signal.size());
			for (unsigned int j = 0; j < insignal2.signal.size(); j++)
			{
				if ((j + i) >= 0 && (j + i) < insignal1.signal.size())
				{
					(*mass)[j] = insignal2.signal[j] * conj(insignal1.signal[j + i]);
					(*mass)[j] /= insignal2.signal.size();
				}
				else continue;
			}
			//FFT(mass, fmass, -1);
			newFFT((*mass), -1);
			double max = abs((*mass)[0]);
			complex<double> otsh_kurilsika = complex<double>(0, 0);
			for (int k = 0; k < (*mass).size(); k++)
			{
				double buf = abs((*mass)[k]);
				if (buf > max)
				{
					max = buf;
					otsh_kurilsika = (*mass)[k];
				}
			}
			outsignal.signal[iter] = (otsh_kurilsika);
			outsignal.keys.push_back(i);//ключи времени
			iter++;
		}
	
		delete mass;
		delete fmass;
	}


	static void coherentSumm(const vector<Signal> &sgs, Signal &result)
	{
		if (sgs.size() == 0)
		{
			return;
		}
		else if (sgs.size() == 1)
		{
			result = sgs[0];
			return;
		}
        qint64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		result.name = "COHERENT_SUMM:" + to_string(timestamp);
		result.description = "Coherent summ of " + sgs[0].name;
		result.signal = sgs[0].signal;
		result.sampling = sgs[0].sampling;
		for (unsigned int i = 1; i < sgs.size(); i++)
		{
			result.description = result.description + " and " + sgs[i].name;
			for (unsigned int j = 0; j < sgs[i].signal.size(); j++)
			{
				result.signal[j] += sgs[i].signal[j];
			}
		}
		if (result.keys.empty())
		{
			for (unsigned int j = 0; j < sgs[0].signal.size(); j++)
			{
				result.keys.push_back(sgs[0].keys[j]);
			}
		}
	}/*
	static void incoherentSumm(vector<Signal> &sgs, Signal &result)
	{
		if (sgs.size() == 0)
		{
			return;
		}
		else if (sgs.size() == 1)
		{
			result = sgs[0];
			return;
		}
		__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		result.name = "INCOHERENT_SUMM:" + to_string(timestamp);
		result.description = "Incoherent summ of " + sgs[0].name;
		result.signal.clear();

		for (unsigned int j = 0; j < sgs[0].signal.size(); j++)
		{
			result.signal.push_back(
				complex<double>(
					sqrt(sgs[0].signal[j].real()*sgs[0].signal[j].real() + sgs[0].signal[j].imag()*sgs[0].signal[j].imag()),
					0)
			);
		}

		for (unsigned int i = 1; i < sgs.size(); i++)
		{
			result.description = result.description + " and " + sgs[i].name;
			for (unsigned int j = 0; j < sgs[i].signal.size(); j++)
			{
				result.signal[j] += complex<double>(
					sqrt(sgs[i].signal[j].real()*sgs[i].signal[j].real() + sgs[i].signal[j].imag()*sgs[i].signal[j].imag()),
					0);
			}
		}
		for (unsigned int j = 0; j < sgs[0].signal.size(); j++)
		{
			result.keys.push_back(sgs[0].keys[j]);
		}
	}
	static void incoherentSummSq(vector<Signal> &sgs, Signal &result)
	{
		if (sgs.size() == 0)
		{
			return;
		}
		else if (sgs.size() == 1)
		{
			result = sgs[0];
			return;
		}
		__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		result.name = "INCOHERENT_SUMM:" + to_string(timestamp);
		result.description = "Incoherent summ of " + sgs[0].name;
		result.signal.clear();

		for (unsigned int j = 0; j < sgs[0].signal.size(); j++)
		{
			result.signal.push_back(
				complex<double>(
				(sgs[0].signal[j].real()*sgs[0].signal[j].real() + sgs[0].signal[j].imag()*sgs[0].signal[j].imag()),
					0)
			);
		}

		for (unsigned int i = 1; i < sgs.size(); i++)
		{
			result.description = result.description + " and " + sgs[i].name;
			for (unsigned int j = 0; j < sgs[i].signal.size(); j++)
			{
				result.signal[j] += complex<double>(
					(sgs[i].signal[j].real()*sgs[i].signal[j].real() + sgs[i].signal[j].imag()*sgs[i].signal[j].imag()),
					0);
			}
		}
		for (unsigned int j = 0; j < result.signal.size(); j++)
		{
			result.signal[j] = complex<double>(
				sqrt(result.signal[j].real()),
				0);
		}
		for (unsigned int j = 0; j < sgs[0].signal.size(); j++)
		{
			result.keys.push_back(sgs[0].keys[j]);
		}
	}
	static double countCriterium(Signal &sg)
	{
		double ret = 0;
		double max = 0;
		double average = 0;
		double msq = 0;
		for (unsigned int i = 0; i < sg.signal.size(); i++)
		{
			double mean = sqrt(sg.signal[i].real()*sg.signal[i].real() + sg.signal[i].imag()*sg.signal[i].imag());
			if (max < mean)
			{
				max = mean;
			}
			average += mean;
		}
		average /= sg.signal.size();

		for (unsigned int i = 0; i < sg.signal.size(); i++)
		{
			double mean = sqrt(sg.signal[i].real()*sg.signal[i].real() + sg.signal[i].imag()*sg.signal[i].imag());
			msq += (average - mean)*(average - mean);
		}
		msq = sqrt(msq / sg.signal.size());
		ret = (max - average) / msq;
		return ret;
	}
	static void NFFT(Signal in, Signal &out, int direction)
	{
		double sqn = sqrt(in.signal.size());
		double signal_size = in.signal.size();
		out = in;
		out.name = string("DFT:") + in.name;
		out.description = string("DFT of ") + in.name;
		for (unsigned int i = 0; i < in.signal.size(); i++)
		{
			out.signal[i] = complex<double>(0, 0);
			for (unsigned int j = 0; j < in.signal.size(); j++)
			{
				out.signal[i] += in.signal[j] * exp(direction*2.*M_PI*complex<double>(0, 1) / signal_size*(double)i*(double)j);
			}
			out.signal[i] /= sqn;
		}
	}
	static void NFFT(vector<complex<double>> & in, vector<complex<double>> & out, int direction)
	{
		double sqn = sqrt(in.size());
		double signal_size = in.size();
		out = in;
		for (unsigned int i = 0; i < in.size(); i++)
		{
			out[i] = complex<double>(0, 0);
			for (unsigned int j = 0; j < in.size(); j++)
			{
				out[i] += in[j] * exp(direction*2.*M_PI*complex<double>(0, 1) / signal_size*(double)i*(double)j);
			}
			out[i] /= sqn;
		}
	}
	static void AmbiguityFunctionSections(Signal &insignal1, Signal &insignal2, Signal &afs, Signal &afns, Signal &mafs, Signal &mafns)
	{
		Signal2D a_f;
		AmbiguityFunction(insignal1, insignal2, a_f);
		int view_idx = 0;
		afs.signal.clear(); afs.keys.clear();
		afns.signal.clear(); afns.keys.clear();
		mafs.signal.clear(); mafs.keys.clear();
		mafns.signal.clear(); mafns.keys.clear();
		__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		for (unsigned int i = 1; i < a_f.signal.size(); i++)
		{
			view_idx = 0;
			double view = 0;
			for (unsigned int j = 1; j < a_f.signal[i].size(); j++)
			{
				double mean = a_f.signal[i][j].real()*a_f.signal[i][j].real() + a_f.signal[i][j].imag()*a_f.signal[i][j].imag();
				if (view < mean)
				{
					view = mean;
					view_idx = j;
				}
			}
			mafs.signal.push_back(a_f.signal[i][view_idx]);
			mafs.keys.push_back(i);
		}
		mafs.name = string("MAFS") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
		mafs.description = "Ambiguity function section of " + insignal1.name + " and " + insignal2.name;

		for (unsigned int i = 1; i < a_f.signal[0].size(); i++)
		{
			view_idx = 0;
			double view = 0;
			for (unsigned int j = 1; j < a_f.signal.size(); j++)
			{
				double mean = a_f.signal[j][i].real()*a_f.signal[j][i].real() + a_f.signal[j][i].imag()*a_f.signal[j][i].imag();
				if (view < mean)
				{
					view = mean;
					view_idx = j;
				}
			}
			mafns.signal.push_back(a_f.signal[view_idx][i]);
			mafns.keys.push_back(i);
		}
		mafns.name = string("MAFS_N") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
		mafns.description = "Ambiguity function section - normal of " + insignal1.name + " and " + insignal2.name;

		double view = 0;
		int view_idx_i = 0;
		int view_idx_j = 0;

		for (unsigned int i = 1; i < a_f.signal.size(); i++)
		{
			for (unsigned int j = 1; j < a_f.signal[i].size(); j++)
			{
				double mean = a_f.signal[i][j].real()*a_f.signal[i][j].real() + a_f.signal[i][j].imag()*a_f.signal[i][j].imag();
				if (view < mean)
				{
					view = mean;
					view_idx_i = i;
					view_idx_j = j;
				}
			}
		}

		for (unsigned int i = 1; i < a_f.signal.size(); i++)
		{
			afs.signal.push_back(a_f.signal[i][view_idx_j]);
			afs.keys.push_back(i);
		}
		afs.name = string("AFS") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
		afs.description = "Ambiguity function section of " + insignal1.name + " and " + insignal2.name;

		for (unsigned int j = 1; j < a_f.signal[0].size(); j++)
		{
			afns.signal.push_back(a_f.signal[view_idx_i][j]);
			afns.keys.push_back(j);
		}
		afns.name = string("AFS_N") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
		afns.description = "Ambiguity function section normal of " + insignal1.name + " and " + insignal2.name;
	}
	static void AmbiguityFunctionSection(Signal &insignal1, Signal &insignal2, Signal &outsignal, bool vector = true, int *max_t = NULL, int *max_f = NULL)
	{
		Signal2D a_f;
		AmbiguityFunction(insignal1, insignal2, a_f);
		int view_idx = 0;
		double view = 0;
		int view_idx_i = 0;
		int view_idx_j = 0;
		outsignal.signal.clear();
		outsignal.keys.clear();

		for (unsigned int i = 0; i < a_f.signal.size(); i++)
		{
			for (unsigned int j = 0; j < a_f.signal[i].size(); j++)
			{
				double mean = a_f.signal[i][j].real()*a_f.signal[i][j].real() + a_f.signal[i][j].imag()*a_f.signal[i][j].imag();
				if (view < mean)
				{
					view = mean;
					view_idx_i = i;
					view_idx_j = j;
				}
			}
		}

		if (max_t != NULL && max_f != NULL)
		{
			(*max_t) = view_idx_i;
			(*max_f) = view_idx_j;
		}


		__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
		if (vector)
		{
			for (unsigned int i = 0; i < a_f.signal.size(); i++)
			{
				view_idx = 0;
				double view = 0;
				for (unsigned int j = 0; j < a_f.signal[i].size(); j++)
				{
					double mean = a_f.signal[i][j].real()*a_f.signal[i][j].real() + a_f.signal[i][j].imag()*a_f.signal[i][j].imag();
					if (view < mean)
					{
						view = mean;
						view_idx = j;
					}
				}
				outsignal.signal.push_back(a_f.signal[i][view_idx]);
				outsignal.keys.push_back(i);
			}
			outsignal.name = string("MAFS") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
			outsignal.description = "Ambiguity function section of " + insignal1.name + " and " + insignal2.name;
		}
		else
		{
			for (unsigned int i = 0; i < a_f.signal[0].size(); i++)
			{
				view_idx = 0;
				double view = 0;
				for (unsigned int j = 0; j < a_f.signal.size(); j++)
				{
					double mean = a_f.signal[j][i].real()*a_f.signal[j][i].real() + a_f.signal[j][i].imag()*a_f.signal[j][i].imag();
					if (view < mean)
					{
						view = mean;
						view_idx = j;
					}
				}
				outsignal.signal.push_back(a_f.signal[view_idx][i]);
				outsignal.keys.push_back(i);
			}
			outsignal.name = string("MAFS_N") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
			outsignal.description = "Ambiguity function section - normal of " + insignal1.name + " and " + insignal2.name;
		}


	}
	static void AmbiguityFunctionNormSection(Signal &insignal1, Signal &insignal2, Signal &outsignal, bool vector = true, int *max_t = NULL, int *max_f = NULL)
	{
		Signal2D a_f;
		AmbiguityFunction(insignal1, insignal2, a_f);
		double view = 0;
		int view_idx_i = 0;
		int view_idx_j = 0;
		for (unsigned int i = 0; i < a_f.signal.size(); i++)
		{
			for (unsigned int j = 0; j < a_f.signal[i].size(); j++)
			{
				double mean = a_f.signal[i][j].real()*a_f.signal[i][j].real() + a_f.signal[i][j].imag()*a_f.signal[i][j].imag();
				if (view < mean)
				{
					view = mean;
					view_idx_i = i;
					view_idx_j = j;
				}
			}
		}

		if (max_t != NULL && max_f != NULL)
		{
			(*max_t) = view_idx_i;
			(*max_f) = view_idx_j;
		}

		__int64 timestamp = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

		outsignal.signal.clear();
		if (vector)
		{
			for (unsigned int i = 0; i < a_f.signal.size(); i++)
			{
				outsignal.signal.push_back(a_f.signal[i][view_idx_j]);
				outsignal.keys.push_back(i);
			}
			outsignal.name = string("AFS") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
			outsignal.description = "Ambiguity function section of " + insignal1.name + " and " + insignal2.name;
		}
		else
		{
			for (unsigned int j = 0; j < a_f.signal[0].size(); j++)
			{
				outsignal.signal.push_back(a_f.signal[view_idx_i][j]);
				outsignal.keys.push_back(j);
			}
			outsignal.name = string("AFS_N") + "(" + to_string(timestamp) + "):(" + insignal1.name + "|" + insignal2.name + ")";
			outsignal.description = "Ambiguity function section normal of " + insignal1.name + " and " + insignal2.name;
		}

	}
	static void AmbiguityFunctionEnhanceFrequencyResolution(Signal &insignal1, Signal &insignal2, Signal &outsignal, int max_t, int max_freq)
	{
		double freq_up = insignal2.sampling / insignal2.signal.size() * (max_freq + 1);
		double freq_down = insignal2.sampling / insignal2.signal.size() * (max_freq - 1);
		complex<double> comp_j = complex<double>(0, 1);
		double new_freq_prec = 0.5;

		outsignal.signal.clear();
		outsignal.keys.clear();
		outsignal.name = "Enhance: " + insignal1.name + ":" + insignal2.name;

		double vval = 0;
		double fmax = 0;

		for (double freq = freq_down; freq < freq_up; freq += new_freq_prec)
		{
			complex<double> value(0, 0);
			for (unsigned int j = 0; j < insignal2.signal.size(); j++)
			{
				value += insignal2.signal[j] * conj(insignal1.signal[j + max_t]) * exp(-comp_j*2.*M_PI*freq*((double)j / insignal2.sampling));
			}
			if ((value.real()*value.real() + value.imag()*value.imag()) > vval)
			{
				vval = value.real()*value.real() + value.imag()*value.imag();
				fmax = freq;
			}
		}

		outsignal.signal.resize(insignal1.signal.size() - insignal2.signal.size());

		for (unsigned int i = 0; i < insignal1.signal.size() - insignal2.signal.size(); i++)
		{
			for (unsigned int j = 0; j < insignal2.signal.size(); j++)
			{
				outsignal.signal[i] += insignal2.signal[j] * conj(insignal1.signal[j + i]) * exp(-comp_j*2.*M_PI*fmax*((double)j / insignal2.sampling));

			}
			outsignal.signal[i] /= insignal2.signal.size();
			outsignal.keys.push_back(i);
		}
	}
	static void Zeroize(Signal & sig)
	{
		double average = 0;
		for (auto it = sig.signal.begin(); it != sig.signal.end(); it++)
		{
			average += (*it).real();
		}
		average /= sig.signal.size();
		for (auto it = sig.signal.begin(); it != sig.signal.end(); it++)
		{
			(*it) -= complex<double>(average, 0);
		}
	}
	static int findMaxIndex(Signal &sig)
	{
		int itermax = 0;
		double mean = 0;
		for (unsigned int i = 0; i < sig.signal.size(); i++)
		{
			double curr = sig.signal[i].real()* sig.signal[i].real() + sig.signal[i].imag()*sig.signal[i].imag();
			if (curr > mean)
			{
				itermax = i;
				mean = curr;
			}
		}
		return itermax;
	}
	static complex<double> findMaxValue(Signal & sig)
	{
		return sig.signal[findMaxIndex(sig)];
	}
	static pair<int, int> findMaxIndex(Signal2D & sig)
	{
		int itermax_i = 0;
		int itermax_j = 0;
		double mean = 0;
		for (unsigned int i = 0; i < sig.signal.size(); i++)
		{
			for (unsigned int j = 0; j < sig.signal[i].size(); j++)
			{
				double curr = sig.signal[i][j].real()* sig.signal[i][j].real() + sig.signal[i][j].imag()*sig.signal[i][j].imag();
				if (curr > mean)
				{
					itermax_i = i;
					itermax_j = j;
					mean = curr;
				}
			}
		}
		return pair<int, int>(itermax_i, itermax_j);
	}
	static complex<double> findMaxValue(Signal2D & sig)
	{
		auto p = findMaxIndex(sig);
		return sig.signal[p.first][p.second];
	}*/
};
