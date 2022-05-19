#include "GenerationSignal.h"
#include "mainwindow.h"

vector<int> GenerationSignal::find_max_n(int n, Signal& sig, double clearance)
{
    if (sig.signal.empty()) return {};
    vector<double> abs_values(sig.signal.size());
    for (int i = 0; i < abs_values.size(); i++)
    {
        abs_values[i] = abs(sig.signal[i]);
    }

    vector <int> max_inds;
    for (int max_idx = 0; max_idx < n; max_idx++)
    {
        int max_ind = 0;
        for (int i = 1; i < sig.signal.size() - 1; i++)
        {
            if (abs_values[i - 1] < abs_values[i] && abs_values[i + 1] < abs_values[i])
            {
                int peak_index = i;
                bool unique_peak = true;
                for (int j = 0; j < max_inds.size(); j++)
                {
                    unique_peak &= abs(peak_index - max_inds[j]) > clearance;	//&& max_inds[j] != 0;
                }
                if ((abs_values[peak_index] > abs_values[max_ind]) && unique_peak)
                {
                    max_ind = i;
                }
            }
        }
        max_inds.push_back(max_ind);
    }
    sort(max_inds.begin(), max_inds.end());
    for (auto& ind : max_inds)
    {
        ind = sig.keys[ind];
    }

    return max_inds;
}

vector<double> GenerationSignal::criteria(vector<vector<int>> delays, double clear)
{
    int summ = 0;
    int s = 0;
    int buf = 0;
    vector<double> SummDelays;
    string prob = " ";
        for (int i = 0; i < delays[0].size(); i++)// ïî ýëåìåíòàì ïåðâîé ñòðîêå
    {
        summ = delays[0][i];
        s = 1;
        for (int k = 0; k < delays[0].size(); k++)//ïî ñòðîêàì
        {
            summ += delays[s][k];
            int iter = 0;
                buf = summ;
                for (int m = 0; m < delays[0].size(); m++)// ïî ñòîëáöàì
                {
                    summ += delays[s + 1][m];
                    Stroka1.push_back(to_string(summ));
                    Stroka2.push_back( to_string(i) + ' ' + to_string(k) + ' ' + to_string(m));
                    if (abs(summ) <= clear)
                    {
                        SummDelays.push_back(summ);
                        string stroka1 = to_string(i) + ' ' + to_string(k) + ' ' + to_string(m);
                        string stroka2 = to_string(Freq[0][i]) + ' ' + to_string(Freq[1][k]) + ' ' + to_string(Freq[2][m]);
                        string stroka3 = to_string(Time[0][i]) + ' ' + to_string(Time[1][k]) + ' ' + to_string(Time[2][m]);
                        str1.push_back(stroka1);
                        str2.push_back(stroka2);
                        str3.push_back(stroka3);
                    }
                    summ = buf;
                }
            summ = delays[0][i];

        }

    }
    return SummDelays;
}

vector<double> GenerationSignal::Delays(int NumberOfSources, int NumberOfSatellites)
{
    vector<double> del;
        for (int i = 0; i < NumberOfSources; i++)
    {
        for (int j = 0; j < NumberOfSatellites; j++)
        {
            del.push_back(1.+ 0.5*((double)rand() / RAND_MAX - 0.5));
            //del.push_back(DoubleRand(3, 3.1));
            //del.push_back(DoubleRand(1.5, 2));
            //del.push_back(DoubleRand(0.75, 1.25));
        }
    }

    return del;
}

double GenerationSignal::DoubleRand(double _max, double _min)
{
    return _min + double(rand()) / RAND_MAX * (_max - _min);
}

vector<double> GenerationSignal::DelaysFrequency(int NumberOfSources, int NumberOfSatellites, double max_freq)
{
    vector<double> del;
    del.clear();
    for (int i = 0; i < NumberOfSources; i++)
    {
        for (int j = 0; j < NumberOfSatellites; j++)
        {
            del.push_back(DoubleRand(-max_freq, max_freq));
        }
    }
    return del;
}

vector<double>GenerationSignal::GenerateDeltaK(int NumberOfSatellites, double max_deltaK)
{
    vector<double>deltaK(NumberOfSatellites);
    double proc=30;
    double del=max_deltaK*proc/100;
    for(int i=0;i<deltaK.size();i++)
    {
        deltaK[i]=DoubleRand(max_deltaK-del,max_deltaK+del);
    }
    return deltaK;
}

//double DoubleRand(double _max, double _min)
//
  //  return _min + double(rand()) / RAND_MAX * (_max - _min);
//}

void GenerationSignal::cutFrequency(Signal2D& insignal/*, Signal2D& outsignal*/)
{
    complex<double> buf;
    for (int i = 0; i < insignal.signal.size() - 1; i++)
    {
        for (int j = 0; j < insignal.signal[i].size()/2; j++)
        {
            buf = insignal.signal[i][j];
            insignal.signal[i][j] = insignal.signal[i][j+ insignal.signal[i].size() / 2];
            insignal.signal[i][j + insignal.signal[i].size() / 2] = buf;
        }
    }
    vector<double> bufKey;
    bufKey.resize(insignal.keysX.size());
    for (int j = 0; j < insignal.keysX.size(); j++)
    {
        if(j< insignal.keysX.size()/2)
            bufKey[j] = -insignal.keysX[-j+ insignal.keysX.size()/2];
        if (j > insignal.keysX.size() / 2)
            bufKey[j] = insignal.keysX[j - insignal.keysX.size() / 2 ];
    }
    insignal.keysX = bufKey;
    bufKey.clear();
}


vector<int> GenerationSignal::find_max_time(int n, vector<double>& sig, vector<double>& keysX, vector<double>& keysY,bool time_freq, double clearance)
{
    if (sig.empty()) return {};
    vector<double> abs_values(sig.size());
    for (int i = 0; i < abs_values.size(); i++)
    {
        abs_values[i] = abs(sig[i]);
    }

    vector <int> max_inds;
    for (int max_idx = 0; max_idx < n; max_idx++)
    {
        int max_ind = 0;
        for (int i = 1; i < sig.size() - 1; i++)
        {
            if (abs_values[i - 1] < abs_values[i] && abs_values[i + 1] < abs_values[i])
            {
                int peak_index = i;
                bool unique_peak = true;
                for (int j = 0; j < max_inds.size(); j++)
                {
                    unique_peak &= abs(peak_index - max_inds[j]) > clearance;	//&& max_inds[j] != 0;
                }
                if ((abs_values[peak_index] > abs_values[max_ind]) && unique_peak)
                {
                    max_ind = i;
                }
            }
        }
        max_inds.push_back(max_ind);
    }
    sort(max_inds.begin(), max_inds.end());
    if (time_freq)
    {
        for (auto& ind : max_inds)
        {
            ind = keysX[ind];
        }
    }
    else
    {
        for (auto& ind : max_inds)
        {
            ind = keysY[ind];
        }
    }

    return max_inds;
}

vector<vector<int>> GenerationSignal::find_max_TF(int n, Signal2D& sig, double clearance)
{
    if (sig.signal.empty()) return {};
    //vector<vector<double>> abs_values(sig.signal.size());
    vector<vector<double>>* abs_values = new vector<vector<double>>(sig.signal.size());
    for (int i = 0; i < (*abs_values).size(); i++)
    {
        (*abs_values)[i].resize(sig.signal[i].size());
        for (int j = 0; j < (*abs_values)[i].size(); j++)
        {
            (*abs_values)[i][j] = abs(sig.signal[i][j]);
        }
    }

    vector<int> max_inds_freq;
    vector<int> max_inds_time;
    //max_inds.resize(n);
    for (int max_idx = 0; max_idx < n; max_idx++)
    {
        int max_ind_time = 0, max_ind_freq=0;
        for (int i = 1; i < sig.signal.size() - 1; i++)
        {
            //max_ind_freq = 0;
            for (int j = 1; j < sig.signal[i].size() - 1; j++)
            {
                if (((*abs_values)[i - 1][j] < (*abs_values)[i][j] && (*abs_values)[i + 1][j] < (*abs_values)[i][j]) &&
                    ((*abs_values)[i][j-1] < (*abs_values)[i][j] && (*abs_values)[i][j+1] < (*abs_values)[i][j]))
                {
                    int peak_indexI = i;
                    int peak_indexJ = j;
                    bool unique_peakI = true;
                    //bool unique_peakJ = true;
                    for (int k = 0; k < max_inds_time.size(); k++)
                    {
                            unique_peakI &= abs(peak_indexI - max_inds_time[k]) > clearance;	//&& max_inds[j] != 0;
                            //unique_peakJ &= abs(peak_indexJ - max_inds_freq[k]) > fr;	//???
                    }
                    if (((*abs_values)[peak_indexI][peak_indexJ] > (*abs_values)[max_ind_time][max_ind_freq]) && /*((*abs_values)[peak_indexJ] > (*abs_values)[max_ind_freq]) &&*/
                        unique_peakI /*&& unique_peakJ*/)
                    {
                        max_ind_time = i;
                        max_ind_freq = j;
                    }
                }
            }
        }
        max_inds_time.push_back(max_ind_time);
        max_inds_freq.push_back(max_ind_freq);
    }
    //sort(max_inds.begin(), max_inds.end());

    timett.resize(n);
    freq.resize(n);

    Time.resize(NumberOfSatellites);
    Freq.resize(NumberOfSatellites);

    for (int indI = 0; indI < max_inds_time.size(); indI++)
    {
        timett[indI] = sig.keysY[max_inds_time[indI]];
        freq[indI] = sig.keysX[max_inds_freq[indI]];

        Time[sput].push_back(sig.keysY[max_inds_time[indI]]);
        Freq[sput].push_back(sig.keysX[max_inds_freq[indI]]);
    }
    vector<vector<int>> max_inds;
    max_inds.resize(max_inds_time.size());
    for (int indI = 0; indI < max_inds_time.size(); indI++)
    {
        max_inds[indI].push_back(timett[indI]);
        max_inds[indI].push_back(freq[indI]);
    }

    /*vector<int> buf;
    double min = Time[sput][0];
    for (int i = 0; i < Time[sput].size(); i++)
    {
        if (i > 0) min = buf[i-1];
        for (int j = i + 1; j < Time[sput].size(); j++)
        {
            if (min > Time[sput][j]) min = Time[sput][j];
        }
        buf.push_back(min);
    }
    Time[sput] = buf;*/
    sort(Time[sput].begin(), Time[sput].end());
    return max_inds;
}

vector<double> GenerationSignal::criteria(vector<vector<int>> delays)
{
    int summ = 0;
    int s = 0;
    int buf = 0;
    int iter = 0;
    int MinSum = 0;
    vector<double> SummDelays;
    for (int i = 0; i < delays[0].size(); i++)// ïî ýëåìåíòàì ïåðâîé ñòðîêå
    {
        summ = delays[0][i];
        s = 1;
        for (int k = 0; k < delays[0].size(); k++)//ïî ñòðîêàì
        {
            summ += delays[s][k];
            int iter = 0;
            buf = summ;
            //vector<double>MinSum;
            for (int m = 0; m < delays[0].size(); m++)// ïî ñòîëáöàì
            {
                summ += delays[s + 1][m];
                //MinSum.push_back(summ);
                if (abs(summ) <= clearance)
                {
                    SummDelays.push_back(summ);
                    string stroka = to_string(i) + ' ' + to_string(k) + ' ' + to_string(m);
                    string stroka1 = to_string(i) +  to_string(k) +  to_string(m);
                    str.push_back(stroka);
                    strr.push_back(stroka1);
                    ind[iter].push_back(i);
                    ind[iter].push_back(k);
                    ind[iter].push_back(m);
                }
                summ = buf;
            }
            /*double min = MinSum[0];
            double m = 0;
            for (int s = 0; s < MinSum.size(); s++)
            {
                if (abs(min) > abs(MinSum[s]))
                {
                    min = MinSum[s];
                    m = s;
                }
            }
            SummDelays.push_back(min);
            string stroka = to_string(i) + to_string(k) + to_string(m);
            str.push_back(stroka);
            ind[iter].push_back(i);
            ind[iter].push_back(k);
            ind[iter].push_back(m);*/


            summ = delays[0][i];

        }
        iter++;
    }
    return SummDelays;
}

void GenerationSignal::GenerateTwoLongSignal()
{
    SignalType signal_type;
    srand(time(0));
        /*if (bpsk.GetCheck())
            signal_type = BPSK;
        if (msk.GetCheck())
            signal_type = MSK;*/

        signal_type = FHSS;
        //signal_type = MSK;
        //signal_type = BPSK;

        vector<Signal>* massLongSignal = new vector<Signal>(NumberOfSources);

        SignalGenerationParameters parametrs = {startTimestamp,startPhase,Duration,nSamples,
            samplingFrequency,Bitrate,additionalParameter };

        for (int i = 0; i < NumberOfSources; i++)
        {
            SignalGenerator ss;
            ss.GenerateSignal(signal_type, parametrs, (*massLongSignal)[i]);
        }

        spectr_long_signal.resize(NumberOfSources);



        for (int k = 0; k < NumberOfSources; k++)
            {
                double size = (*massLongSignal)[k].signal.size();


                vector<complex<double>>* sp = new vector<complex<double>>;
                vector<complex<double>>* outsp = new vector<complex<double>>;
                (*sp).resize(size);
                for (int i = 0; i < size; i++)
                {
                    (*sp)[i] = (*massLongSignal)[k].signal[i];
                }

                Methods::newFFT((*sp), -1);

                for (int i = 0; i < (*sp).size(); i++)
                {
                    spectr_long_signal[k].push_back(abs((*sp)[i]));
                }
                double sampling = (*massLongSignal)[k].sampling;
                double sscale = sampling / spectr_long_signal[k].size();
                vector<double> Fr;
                for (int i = 0; i < spectr_long_signal[k].size(); i++)
                {
                    Fr.push_back(i * sscale);
                }


                //UpdateData(false);
                BigMassOtrisovka.push_back(spectr_long_signal[k]);
                BigMassOtshetX.push_back(Fr);
                listt.push_back("Длиннющий спектр "+ QString::number(k+1));
                out.push_back(false);
            }

        vector<vector<Signal>>* massSignal = new vector<vector<Signal>>(NumberOfSatellites, vector<Signal>(NumberOfSources));



            vector<double> sdvigtime;
            sdvigtime.clear();
            sdvigtime = Delays(NumberOfSources, NumberOfSatellites);

            vector<double> sdvigFrequency;
            sdvigFrequency.clear();

            sdvigFrequency = DelaysFrequency(NumberOfSources, NumberOfSatellites, 20000);

            DeltaK.clear();
            double max_deltaK=0.0002;
            DeltaK= GenerateDeltaK(NumberOfSatellites,max_deltaK);

            int k = 0;
            double d = Duration / 90;



            double shum=-3;

            for (int i = 0; i < NumberOfSatellites; i++)
            {
                for (int j = 0; j < NumberOfSources; j++)
                {
                    Methods::transformSignal((*massLongSignal)[j], sdvigtime[i + NumberOfSatellites * j] * d, d, sdvigFrequency[i + NumberOfSatellites * j]*DeltaK[i], 1+DeltaK[i], shum, (*massSignal)[i][j]);
                }
                if (NumberOfSources == 1)
                {
                    (*SummSignal)[i] = (*massSignal)[i][0];
                }
                for (int j = 1; j < NumberOfSources; j++)
                {
                    if (j == 1)	Methods::coherentSumm({ (*massSignal)[i][0] , (*massSignal)[i][j] }, (*SummSignal)[i]);
                    else		Methods::coherentSumm({ (*SummSignal)[i] , (*massSignal)[i][j] }, (*SummSignal)[i]);
                }

            }

            massdoubleSignal.resize(NumberOfSatellites);
                for (int i = 0; i < NumberOfSatellites; i++)
                {
                    if (i + 1 == NumberOfSatellites) Methods::Newcorrelate((*SummSignal)[i], (*SummSignal)[0], (*CorrelateSignal)[i]);
                    else Methods::Newcorrelate((*SummSignal)[i], (*SummSignal)[i + 1], (*CorrelateSignal)[i]);
                    vector<double> Fr;
                    for (int j = 0; j < (*CorrelateSignal)[i].signal.size(); j++)
                    {
                        massdoubleSignal[i].push_back(abs((*CorrelateSignal)[i].signal[j]));
                        Fr.push_back((*CorrelateSignal)[i].keys[j]);
                    }
                    BigMassOtrisovka.push_back(massdoubleSignal[i]);
                    BigMassOtshetX.push_back(Fr);
                    listt.push_back("Корреляция просуммированных сигналов "+ QString::number(i));
                    out.push_back(false);
                }
                delete massLongSignal;
                delete massSignal;
                massdoubleSignal.clear();

}

void GenerationSignal::CreateVFN()
{
    sput = 0;

        for (int i = 0; i < NumberOfSatellites; i++)
        {
            if (i + 1 == NumberOfSatellites) Methods::SuperVFN((*SummSignal)[i], (*SummSignal)[0], (*A)[i], scale,1+DeltaK[i]);
            else Methods::SuperVFN((*SummSignal)[i], (*SummSignal)[i + 1], (*A)[i], scale,1+DeltaK[i]);
        }

        sput = 0;
        int N=NumberOfSources;
        for (int i = 0; i < NumberOfSatellites; i++)
        {
            cutFrequency((*A)[i]);
            clearance = 1.0 / Bitrate * samplingFrequency;
            fr = 2 * scale;
            del = find_max_TF(N/*NumberOfSources*//* * NumberOfSatellites*/, (*A)[i], clearance);
            sput++;
        }


        vector<vector<vector<double>>> srez;
        srez.resize(Time.size());

        for (int i = 0; i < Time.size(); i++)
        {
            srez[i].resize(N);
            for (int j = 0; j < N; j++)
            {
                int fr = 0;

                for (int k = 0; k < (*A)[i].keysX.size(); k++)
                {

                    if ((int)(*A)[i].keysX[k] == Freq[i][j])
                    {
                        fr = k;
                        Signal sig;
                        for (int m = 0; m < (*A)[i].keysY.size(); m++)
                        {
                            sig.signal.push_back((*A)[i].signal[m][fr]);
                        }
                        //Methods::newFFT(sig.signal, -1);
                        srez[i][j].resize(sig.signal.size());
                        vector<double> ttime;
                        ttime.resize(sig.signal.size());
                        for (int s = 0; s < sig.signal.size(); s++)
                        {
                            srez[i][j][s] = abs(sig.signal[s]);
                            ttime[s] = (*A)[i].keysY[s];
                        }

                        BigMassOtrisovka.push_back(srez[i][j]);
                        BigMassOtshetX.push_back(ttime);
                        listt.push_back("Срез ВФН "+ QString::number(j+1) +QString::number(i+1));

                        out.push_back(true);
                        //LBox.AddString(sstr);
                        break;
                    }

                }
            }
        }



        //delete SummSignal;
        //InitiateOPGL();
        //vector<vector<double>> modP;
        ////тут надо вывести картинку


        //firsttimedraw = 1;
        //delete newA;


        string res = "";

        res = res + "\r\n" + "Результат работы критерия по функциям неопределенности:\r\n";


        SummDelays.clear();
        str1.clear();
        str2.clear();
        str3.clear();
        SummDelays = criteria(Time, clearance);
        KolOtch = SummDelays.size();

        ofstream tout("Куча времени.txt");
        tout << "Номера задержек" << "\t" << "Сумма задержек" << "\n";
        for (int i = 0; i < SummDelays.size(); i++)
        {
            tout << str1[i] << "\t\t" << SummDelays[i] << "\n";
            tout << str3[i] << "\n\n";

            res = res + str1[i] + "\t\t" + to_string((int)SummDelays[i]) + "\r\n";
        }

        res = res + "\r\n" + "Группы задержек и их сумма:\r\n";

        for (int i = 0; i < Stroka1.size(); i++)
        {

            tout << Stroka1[i] << "\t\t" << Stroka2[i] << "\n";
            res = res + Stroka2[i] + "\t\t" + Stroka1[i] + "\r\n";
        }
        tout.close(); // закрываем файл
        //Info = res.c_str();
        /*str1.clear();
        vector<double> SummFreq;
        SummFreq.clear();
        //scale = (*outsignal)[0].sampling / (*mass).size();
        SummFreq = criteria(Freq, 2 * scale);

        ofstream fout("Куча частот.txt");
        fout << "Номера задержек" << "\t" << "Сумма частот" << "\n";
        for (int i = 0; i < SummFreq.size(); i++)
        {
            fout << str1[i] << "\t\t" << SummFreq[i] << "\n";
            fout << str2[i] << "\n\n";
        }
        fout.close(); // закрываем файл*/
}
