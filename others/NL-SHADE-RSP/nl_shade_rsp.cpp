#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <random>

#include <iomanip>

using namespace std;
unsigned seed1 = 0;
std::mt19937 generator_uni_i(seed1);
std::mt19937 generator_uni_r(seed1 + 100);
std::mt19937 generator_norm(seed1 + 200);
std::mt19937 generator_cachy(seed1 + 300);
std::mt19937 generator_uni_i_2(seed1 + 400);
std::uniform_int_distribution<int> uni_int(0, 32768);
std::uniform_real_distribution<double> uni_real(0.0, 1.0);
std::normal_distribution<double> norm_dist(0.0, 1.0);
std::cauchy_distribution<double> cachy_dist(0.0, 1.0);

int IntRandom(int target)
{
    if (target == 0)
        return 0;
    return uni_int(generator_uni_i) % target;
}
double Random(double minimal, double maximal) { return uni_real(generator_uni_r) * (maximal - minimal) + minimal; }
double NormRand(double mu, double sigma) { return norm_dist(generator_norm) * sigma + mu; }
double CachyRand(double mu, double sigma) { return cachy_dist(generator_cachy) * sigma + mu; }

void qSort1(double *Mass, int low, int high)
{
    int i = low;
    int j = high;
    double x = Mass[(low + high) >> 1];
    do
    {
        while (Mass[i] < x)
            ++i;
        while (Mass[j] > x)
            --j;
        if (i <= j)
        {
            double temp = Mass[i];
            Mass[i] = Mass[j];
            Mass[j] = temp;
            i++;
            j--;
        }
    } while (i <= j);
    if (low < j)
        qSort1(Mass, low, j);
    if (i < high)
        qSort1(Mass, i, high);
}

void qSort2int(double *Mass, int *Mass2, int low, int high)
{
    int i = low;
    int j = high;
    double x = Mass[(low + high) >> 1];
    do
    {
        while (Mass[i] < x)
            ++i;
        while (Mass[j] > x)
            --j;
        if (i <= j)
        {
            double temp = Mass[i];
            Mass[i] = Mass[j];
            Mass[j] = temp;
            int temp2 = Mass2[i];
            Mass2[i] = Mass2[j];
            Mass2[j] = temp2;
            i++;
            j--;
        }
    } while (i <= j);
    if (low < j)
        qSort2int(Mass, Mass2, low, j);
    if (i < high)
        qSort2int(Mass, Mass2, i, high);
}

void cec22_test_func(double *x, double *f, int nx, int mx, int func_num);
double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;
int GNVars;
double stepsFEval[16];
double ResultsArray[12][30][16];
int LastFEcount;
int NFEval = 0;
int MaxFEval = 0;
double tempF[1];
double fopt;
char buffer[200];
double globalbest;
bool globalbestinit;
bool initfinished;
vector<double> FitTemp3;

void GenerateNextRandUnif(const int num, const int Range, int *Rands, const int Prohib)
{
    for (int j = 0; j != 25; j++)
    {
        bool generateagain = false;
        Rands[num] = IntRandom(Range);
        for (int i = 0; i != num; i++)
            if (Rands[i] == Rands[num])
                generateagain = true;
        if (!generateagain)
            break;
    }
}
void GenerateNextRandUnifOnlyArch(const int num, const int Range, const int Range2, int *Rands, const int Prohib)
{
    for (int j = 0; j != 25; j++)
    {
        bool generateagain = false;
        Rands[num] = IntRandom(Range2) + Range;
        for (int i = 0; i != num; i++)
            if (Rands[i] == Rands[num])
                generateagain = true;
        if (!generateagain)
            break;
    }
}
bool CheckGenerated(const int num, int *Rands, const int Prohib)
{
    if (Rands[num] == Prohib)
        return false;
    for (int j = 0; j != num; j++)
        if (Rands[j] == Rands[num])
            return false;
    return true;
}
void SaveBestValues(int funcN, int RunN, double newbestfit)
{
    for (int stepFEcount = LastFEcount; stepFEcount < 16; stepFEcount++)
    {
        if (NFEval == int(stepsFEval[stepFEcount] * MaxFEval))
        {
            double temp = globalbest - fopt;
            if (temp <= 10E-8)
                temp = 0;
            ResultsArray[funcN - 1][RunN][stepFEcount] = temp;
            LastFEcount = stepFEcount;
        }
    }
}
void FindLimits(double *Ind, double *Parent, int CurNVars, double CurLeft, double CurRight)
{
    for (int j = 0; j < CurNVars; j++)
    {
        if (Ind[j] < CurLeft)
            Ind[j] = Random(CurLeft, CurRight);
        if (Ind[j] > CurRight)
            Ind[j] = Random(CurLeft, CurRight);
    }
}
class Optimizer
{
public:
    bool FitNotCalculated;
    int Int_ArchiveSizeParam;
    int MemorySize;
    int MemoryIter;
    int SuccessFilled;
    int MemoryCurrentIndex;
    int NVars;
    int NInds;
    int NIndsMax;
    int NIndsMin;
    int besti;
    int func_num;
    int bench_num;
    int RunN;
    int Generation;
    int ArchiveSize;
    int CurrentArchiveSize;
    double F;
    double Cr;
    double bestfit;
    double ArchiveSizeParam;
    double Right;
    double Left;
    int *Rands;
    int *Indexes;
    int *BackIndexes;
    double *Weights;
    double *Donor;
    double *Trial;
    double *FitMass;
    double *FitMassTemp;
    double *FitMassCopy;
    double *BestInd;
    double *tempSuccessCr;
    double *tempSuccessF;
    double *FGenerated;
    double *CrGenerated;
    double *MemoryCr;
    double *MemoryF;
    double *FitDelta;
    double *ArchUsages;
    double **Popul;
    double **PopulTemp;
    double **Archive;

    void Initialize(int newNInds, int newNVars, int new_func_num, int newRunN, int new_bench_num, int NewMemSize, double NewArchSizeParam);
    void Clean();
    void MainCycle();
    void FindNSaveBest(bool init, int ChosenOne);
    inline double GetValue(const int index, const int NInds, const int j);
    void CopyToArchive(double *RefusedParent, double RefusedFitness);
    void SaveSuccessCrF(double Cr, double F, double FitD);
    void UpdateMemoryCrF();
    double MeanWL_general(double *Vector, double *TempWeights, int Size, double g_p, double g_m);
    void RemoveWorst(int NInds, int NewNInds);
};

double cec_22_(double *HostVector, int func_num, int bench_num)
{
    /*double* HVT = new double[GNVars];
    for(int i=0;i!=GNVars;i++)
    {
        HVT[i] = HostVector[i] - 10000.0;
    }
    switch(bench_num)
    {           //Bias,Shift,Rotation
        case 0:{cec21_basic_func(HVT,  tempF, GNVars, 1, func_num);          break;}  //000
        case 1:{cec21_bias_func(HVT,  tempF, GNVars, 1, func_num);           break;}  //100
        case 2:{cec21_shift_func(HVT,  tempF, GNVars, 1, func_num);          break;}  //010
        case 3:{cec21_rot_func(HVT,  tempF, GNVars, 1, func_num);            break;}  //001
        case 4:{cec21_bias_shift_func(HVT,  tempF, GNVars, 1, func_num);     break;}  //110
        case 5:{cec21_bias_rot_func(HVT,  tempF, GNVars, 1, func_num);       break;}  //101
        case 6:{cec21_shift_rot_func(HVT,  tempF, GNVars, 1, func_num);      break;}  //011
        case 7:{cec21_bias_shift_rot_func(HVT,  tempF, GNVars, 1, func_num); break;}  //111
    }
    delete HVT;*/

    cec22_test_func(HostVector, tempF, GNVars, 1, func_num);
    NFEval++;
    return tempF[0];
}
void Optimizer::Initialize(int newNInds, int newNVars, int new_func_num, int newRunN,
                           int new_bench_num, int NewMemSize, double NewArchSizeParam)
{
    FitNotCalculated = true;
    NInds = newNInds;
    NIndsMax = NInds;
    NIndsMin = 4;
    NVars = newNVars;
    RunN = newRunN;
    Left = -100; //+10000;
    Right = 100; //+10000;
    Cr = 0.2;
    F = 0.2;
    besti = 0;
    Generation = 0;
    CurrentArchiveSize = 0;
    ArchiveSizeParam = NewArchSizeParam;
    Int_ArchiveSizeParam = ceil(ArchiveSizeParam);
    ArchiveSize = NIndsMax * ArchiveSizeParam;
    func_num = new_func_num;
    bench_num = new_bench_num;
    for (int steps_k = 0; steps_k != 15; steps_k++)
        stepsFEval[steps_k] = pow(double(GNVars), double(steps_k) / 5.0 - 3.0);
    stepsFEval[15] = 1.0;
    Popul = new double *[NIndsMax];
    for (int i = 0; i != NIndsMax; i++)
        Popul[i] = new double[NVars];
    PopulTemp = new double *[NIndsMax];
    for (int i = 0; i != NIndsMax; i++)
        PopulTemp[i] = new double[NVars];
    Archive = new double *[NIndsMax * Int_ArchiveSizeParam];
    for (int i = 0; i != NIndsMax * Int_ArchiveSizeParam; i++)
        Archive[i] = new double[NVars];
    FitMass = new double[NIndsMax];
    FitMassTemp = new double[NIndsMax];
    FitMassCopy = new double[NIndsMax];
    Indexes = new int[NIndsMax];
    BackIndexes = new int[NIndsMax];
    BestInd = new double[NVars];
    for (int i = 0; i < NIndsMax; i++)
        for (int j = 0; j < NVars; j++)
            Popul[i][j] = Random(Left, Right);
    Donor = new double[NVars];
    Trial = new double[NVars];
    Rands = new int[NIndsMax];
    tempSuccessCr = new double[NIndsMax];
    tempSuccessF = new double[NIndsMax];
    FitDelta = new double[NIndsMax];
    FGenerated = new double[NIndsMax];
    CrGenerated = new double[NIndsMax];
    for (int i = 0; i != NIndsMax; i++)
    {
        tempSuccessCr[i] = 0;
        tempSuccessF[i] = 0;
    }
    MemorySize = NewMemSize;
    MemoryIter = 0;
    SuccessFilled = 0;
    ArchUsages = new double[NIndsMax];
    Weights = new double[NIndsMax];
    MemoryCr = new double[MemorySize];
    MemoryF = new double[MemorySize];
    for (int i = 0; i != MemorySize; i++)
    {
        MemoryCr[i] = 0.2;
        MemoryF[i] = 0.2;
    }
}
void Optimizer::SaveSuccessCrF(double Cr, double F, double FitD)
{
    tempSuccessCr[SuccessFilled] = Cr;
    tempSuccessF[SuccessFilled] = F;
    FitDelta[SuccessFilled] = FitD;
    SuccessFilled++;
}
void Optimizer::UpdateMemoryCrF()
{
    if (SuccessFilled != 0)
    {
        MemoryCr[MemoryIter] = MeanWL_general(tempSuccessCr, FitDelta, SuccessFilled, 2, 1);
        MemoryF[MemoryIter] = MeanWL_general(tempSuccessF, FitDelta, SuccessFilled, 2, 1);
        MemoryIter++;
        if (MemoryIter >= MemorySize)
            MemoryIter = 0;
    }
    else
    {
        MemoryF[MemoryIter] = 0.5;
        MemoryCr[MemoryIter] = 0.5;
    }
}
double Optimizer::MeanWL_general(double *Vector, double *TempWeights, int Size, double g_p, double g_m)
{
    double SumWeight = 0;
    double SumSquare = 0;
    double Sum = 0;
    for (int i = 0; i != SuccessFilled; i++)
        SumWeight += TempWeights[i];
    for (int i = 0; i != SuccessFilled; i++)
        Weights[i] = TempWeights[i] / SumWeight;
    for (int i = 0; i != SuccessFilled; i++)
        SumSquare += Weights[i] * pow(Vector[i], g_p);
    for (int i = 0; i != SuccessFilled; i++)
        Sum += Weights[i] * pow(Vector[i], g_p - g_m);
    if (fabs(Sum) > 0.000001)
        return SumSquare / Sum;
    else
        return 0.5;
}
void Optimizer::CopyToArchive(double *RefusedParent, double RefusedFitness)
{
    if (CurrentArchiveSize < ArchiveSize)
    {
        for (int i = 0; i != NVars; i++)
            Archive[CurrentArchiveSize][i] = RefusedParent[i];
        CurrentArchiveSize++;
    }
    else if (ArchiveSize > 0)
    {
        int RandomNum = IntRandom(ArchiveSize);
        for (int i = 0; i != NVars; i++)
            Archive[RandomNum][i] = RefusedParent[i];
    }
}
void Optimizer::FindNSaveBest(bool init, int ChosenOne)
{
    if (FitMass[ChosenOne] <= bestfit || init)
    {
        bestfit = FitMass[ChosenOne];
        besti = ChosenOne;
        for (int j = 0; j != NVars; j++)
            BestInd[j] = Popul[besti][j];
    }
    if (bestfit < globalbest)
        globalbest = bestfit;
}
void Optimizer::RemoveWorst(int NInds, int NewNInds)
{
    int PointsToRemove = NInds - NewNInds;
    for (int L = 0; L != PointsToRemove; L++)
    {
        double WorstFit = FitMass[0];
        int WorstNum = 0;
        for (int i = 1; i != NInds; i++)
        {
            if (FitMass[i] > WorstFit)
            {
                WorstFit = FitMass[i];
                WorstNum = i;
            }
        }
        for (int i = WorstNum; i != NInds - 1; i++)
        {
            for (int j = 0; j != NVars; j++)
                Popul[i][j] = Popul[i + 1][j];
            FitMass[i] = FitMass[i + 1];
        }
    }
}
inline double Optimizer::GetValue(const int index, const int NInds, const int j)
{
    if (index < NInds)
        return Popul[index][j];
    return Archive[index - NInds][j];
}
void Optimizer::MainCycle()
{
    double ArchSuccess;
    double NoArchSuccess;
    double NArchUsages;
    double ArchProbs = 0.5;
    for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
    {
        FitMass[TheChosenOne] = cec_22_(Popul[TheChosenOne], func_num, bench_num);
        FindNSaveBest(TheChosenOne == 0, TheChosenOne);
        if (!globalbestinit || bestfit < globalbest)
        {
            globalbest = bestfit;
            globalbestinit = true;
        }
        SaveBestValues(func_num, RunN, bestfit);
    }
    do
    {
        double minfit = FitMass[0];
        double maxfit = FitMass[0];
        for (int i = 0; i != NInds; i++)
        {
            FitMassCopy[i] = FitMass[i];
            Indexes[i] = i;
            if (FitMass[i] >= maxfit)
                maxfit = FitMass[i];
            if (FitMass[i] <= minfit)
                minfit = FitMass[i];
        }
        if (minfit != maxfit)
            qSort2int(FitMassCopy, Indexes, 0, NInds - 1);
        for (int i = 0; i != NInds; i++)
            for (int j = 0; j != NInds; j++)
                if (i == Indexes[j])
                {
                    BackIndexes[i] = j;
                    break;
                }
        FitTemp3.resize(NInds);
        for (int i = 0; i != NInds; i++)
            FitTemp3[i] = exp(-double(i) / (double)NInds);
        std::discrete_distribution<int> ComponentSelector3(FitTemp3.begin(), FitTemp3.end());
        int psizeval = max(2.0, NInds * (0.2 / (double)MaxFEval * (double)NFEval + 0.2));
        int CrossExponential = 0;
        if (Random(0, 1) < 0.5)
            CrossExponential = 1;
        for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
        {
            MemoryCurrentIndex = IntRandom(MemorySize);
            Cr = min(1.0, max(0.0, NormRand(MemoryCr[MemoryCurrentIndex], 0.1)));
            do
            {
                F = CachyRand(MemoryF[MemoryCurrentIndex], 0.1);
            } while (F <= 0);
            FGenerated[TheChosenOne] = min(F, 1.0);
            CrGenerated[TheChosenOne] = Cr;
        }
        qSort1(CrGenerated, 0, NInds - 1);
        for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
        {
            Rands[0] = Indexes[IntRandom(psizeval)];
            for (int i = 0; i != 25 && !CheckGenerated(0, Rands, TheChosenOne); i++)
                Rands[0] = Indexes[IntRandom(psizeval)];
            GenerateNextRandUnif(1, NInds, Rands, TheChosenOne);
            if (Random(0, 1) > ArchProbs || CurrentArchiveSize == 0)
            {
                Rands[2] = Indexes[ComponentSelector3(generator_uni_i_2)];
                for (int i = 0; i != 25 && !CheckGenerated(2, Rands, TheChosenOne); i++)
                    Rands[2] = Indexes[ComponentSelector3(generator_uni_i_2)];
                ArchUsages[TheChosenOne] = 0;
            }
            else
            {
                GenerateNextRandUnifOnlyArch(2, NInds, CurrentArchiveSize, Rands, TheChosenOne);
                ArchUsages[TheChosenOne] = 1;
            }
            for (int j = 0; j != NVars; j++)
                Donor[j] = Popul[TheChosenOne][j] +
                           FGenerated[TheChosenOne] * (GetValue(Rands[0], NInds, j) - Popul[TheChosenOne][j]) +
                           FGenerated[TheChosenOne] * (GetValue(Rands[1], NInds, j) - GetValue(Rands[2], NInds, j));

            int WillCrossover = IntRandom(NVars);
            Cr = CrGenerated[BackIndexes[TheChosenOne]];
            double CrToUse = 0;
            if (NFEval > 0.5 * MaxFEval)
                CrToUse = (double(NFEval) / double(MaxFEval) - 0.5) * 2;
            if (CrossExponential == 0)
            {
                for (int j = 0; j != NVars; j++)
                {
                    if (Random(0, 1) < CrToUse || WillCrossover == j)
                        PopulTemp[TheChosenOne][j] = Donor[j];
                    else
                        PopulTemp[TheChosenOne][j] = Popul[TheChosenOne][j];
                }
            }
            else
            {
                int StartLoc = IntRandom(NVars);
                int L = StartLoc + 1;
                while (Random(0, 1) < Cr && L < NVars)
                    L++;
                for (int j = 0; j != NVars; j++)
                    PopulTemp[TheChosenOne][j] = Popul[TheChosenOne][j];
                for (int j = StartLoc; j != L; j++)
                    PopulTemp[TheChosenOne][j] = Donor[j];
            }
            FindLimits(PopulTemp[TheChosenOne], Popul[TheChosenOne], NVars, Left, Right);
            FitMassTemp[TheChosenOne] = cec_22_(PopulTemp[TheChosenOne], func_num, bench_num);
            if (FitMassTemp[TheChosenOne] <= globalbest)
                globalbest = FitMassTemp[TheChosenOne];

            if (FitMassTemp[TheChosenOne] < FitMass[TheChosenOne])
                SaveSuccessCrF(Cr, F, fabs(FitMass[TheChosenOne] - FitMassTemp[TheChosenOne]));
            FindNSaveBest(false, TheChosenOne);
            SaveBestValues(func_num, RunN, bestfit);
        }
        ArchSuccess = 0;
        NoArchSuccess = 0;
        NArchUsages = 0;
        for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
        {
            if (FitMassTemp[TheChosenOne] <= FitMass[TheChosenOne])
            {
                if (ArchUsages[TheChosenOne] == 1)
                {
                    ArchSuccess += (FitMass[TheChosenOne] - FitMassTemp[TheChosenOne]) / FitMass[TheChosenOne];
                    NArchUsages += 1;
                }
                else
                    NoArchSuccess += (FitMass[TheChosenOne] - FitMassTemp[TheChosenOne]) / FitMass[TheChosenOne];
                CopyToArchive(Popul[TheChosenOne], FitMass[TheChosenOne]);
                for (int j = 0; j != NVars; j++)
                    Popul[TheChosenOne][j] = PopulTemp[TheChosenOne][j];
                FitMass[TheChosenOne] = FitMassTemp[TheChosenOne];
            }
        }
        if (NArchUsages != 0)
        {
            ArchSuccess = ArchSuccess / NArchUsages;
            NoArchSuccess = NoArchSuccess / (NInds - NArchUsages);
            ArchProbs = ArchSuccess / (ArchSuccess + NoArchSuccess);
            ArchProbs = max(0.1, min(0.9, ArchProbs));
            if (ArchSuccess == 0)
                ArchProbs = 0.5;
        }
        else
            ArchProbs = 0.5;
        int newNInds = round((NIndsMin - NIndsMax) * pow((double(NFEval) / double(MaxFEval)), (1.0 - double(NFEval) / double(MaxFEval))) + NIndsMax);
        if (newNInds < NIndsMin)
            newNInds = NIndsMin;
        if (newNInds > NIndsMax)
            newNInds = NIndsMax;
        int newArchSize = round((NIndsMin - NIndsMax) * pow((double(NFEval) / double(MaxFEval)), (1.0 - double(NFEval) / double(MaxFEval))) + NIndsMax) * ArchiveSizeParam;
        if (newArchSize < NIndsMin)
            newArchSize = NIndsMin;
        ArchiveSize = newArchSize;
        if (CurrentArchiveSize >= ArchiveSize)
            CurrentArchiveSize = ArchiveSize;
        RemoveWorst(NInds, newNInds);
        NInds = newNInds;
        UpdateMemoryCrF();
        SuccessFilled = 0;
        Generation++;
    } while (NFEval < MaxFEval);
}
void Optimizer::Clean()
{
    delete Donor;
    delete Trial;
    delete Rands;
    for (int i = 0; i != NIndsMax; i++)
    {
        delete Popul[i];
        delete PopulTemp[i];
    }
    for (int i = 0; i != NIndsMax * Int_ArchiveSizeParam; i++)
        delete Archive[i];
    delete ArchUsages;
    delete Archive;
    delete Popul;
    delete PopulTemp;
    delete FitMass;
    delete FitMassTemp;
    delete FitMassCopy;
    delete BestInd;
    delete Indexes;
    delete BackIndexes;
    delete tempSuccessCr;
    delete tempSuccessF;
    delete FGenerated;
    delete CrGenerated;
    delete FitDelta;
    delete MemoryCr;
    delete MemoryF;
    delete Weights;
}
int main()
{
    cout << scientific << setprecision(8);
    float f;
    FILE *rs;
    rs = fopen("input_data\\Rand_Seeds.txt", "r+");
    int Seeds[1000];
    // for (int i = 0; i != 1000; i++)
    // {
    //     fscanf(rs, "%f", &f);
    //     Seeds[i] = int(f);
    // }

    for (int GNVarsIter = 0; GNVarsIter != 2; GNVarsIter++)
    {
        if (GNVarsIter == 0)
        {
            GNVars = 10;
            MaxFEval = 200000;
        }
        if (GNVarsIter == 1)
        {
            GNVars = 20;
            MaxFEval = 1000000;
        }
        for (int bench_num = 0; bench_num != 1; bench_num++)
        {
            string benchname;
            switch (bench_num)
            {
            case 0:
            {
                benchname = "000";
                break;
            }
            case 1:
            {
                benchname = "100";
                break;
            }
            case 2:
            {
                benchname = "010";
                break;
            }
            case 3:
            {
                benchname = "001";
                break;
            }
            case 4:
            {
                benchname = "110";
                break;
            }
            case 5:
            {
                benchname = "101";
                break;
            }
            case 6:
            {
                benchname = "011";
                break;
            }
            case 7:
            {
                benchname = "111";
                break;
            }
            }
            for (int RunN = 0; RunN != 30; RunN++)
            {
                for (int func_num = 1; func_num < 13; func_num++)
                {
                    seed1 = Seeds[GNVars * func_num * 3 + RunN - 30];
                    std::mt19937 generator_uni_i(seed1);
                    std::mt19937 generator_uni_r(seed1 + 100);
                    std::mt19937 generator_norm(seed1 + 200);
                    std::mt19937 generator_cachy(seed1 + 300);
                    std::mt19937 generator_uni_i_2(seed1 + 400);
                    //cout<<GNVars*func_num*3+RunN-30<<" Random seeds is: "<<seed1<<endl;
                    switch (func_num)
                    {
                    case 1:
                    {
                        fopt = 300;
                        break;
                    }
                    case 2:
                    {
                        fopt = 400;
                        break;
                    }
                    case 3:
                    {
                        fopt = 600;
                        break;
                    }
                    case 4:
                    {
                        fopt = 800;
                        break;
                    }
                    case 5:
                    {
                        fopt = 900;
                        break;
                    }
                    case 6:
                    {
                        fopt = 1800;
                        break;
                    }
                    case 7:
                    {
                        fopt = 2000;
                        break;
                    }
                    case 8:
                    {
                        fopt = 2200;
                        break;
                    }
                    case 9:
                    {
                        fopt = 2300;
                        break;
                    }
                    case 10:
                    {
                        fopt = 2400;
                        break;
                    }
                    case 11:
                    {
                        fopt = 2600;
                        break;
                    }
                    case 12:
                    {
                        fopt = 2700;
                        break;
                    }
                    }
                    globalbestinit = false;
                    initfinished = false;
                    LastFEcount = 0;
                    NFEval = 0;
                    Optimizer OptZ;
                    ////////////////NInds     NVars  func     Run  bench     memory    arch size
                    OptZ.Initialize(GNVars * 30, GNVars, func_num, RunN, bench_num, 20 * GNVars, 2.1);
                    OptZ.MainCycle();
                    OptZ.Clean();
                    cout << GNVars << "D bench " << benchname << " Run " << RunN + 1 << " f " << func_num << " global best " << globalbest - fopt << endl;
                }
            }

            for (int func_num = 1; func_num < 13; func_num++)
            {
                string filename = "NL-SHADE-RSP_";
                sprintf(buffer, "%d_%d.txt", func_num, GNVars);
                filename = filename + buffer;
                cout << filename << endl;
                ofstream fout(filename.c_str(), ios::trunc);
                fout << scientific << setprecision(6);
                for (int step = 0; step != 16; step++)
                {
                    for (int RunN = 0; RunN != 30; RunN++)
                        fout << ResultsArray[func_num - 1][RunN][step] << "\t";
                    fout << endl;
                }
                for (int RunN = 0; RunN != 30; RunN++)
                    fout << 0 << "\t";
                fout << endl;
            }
        }
    }
    return 0;
}
