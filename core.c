#include <TF1.h>
#include <TMath.h>
#include <Math/GaussIntegrator.h>
#include "Math/WrappedFunction.h"
#include "Math/WrappedTF1.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TRandom2.h>
#include <TMultiGraph.h>
#include <TApplication.h>
#include <iostream>
#include <bitset>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#define NOMINMAX
#include <windows.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include ".\di\back_func.c"

#define DATA_ERROR 0.006 //only for signal, not bkg
#define ERR_ESTIMATE 1
#define MODE_ 1
#define NUMBER_OF_SIMULATIONS 50
#define WRITE_STAT 0
#define MANUAL_LIMITS 1

//for testing integral validity
//class  TF1_EvalWrapper : public ROOT::Math::IGenFunction 
//{
//	 public:
//		    TF1_EvalWrapper(TF1 * f, const Double_t * par, bool useAbsVal, Double_t n = 1, Double_t x0 = 0) :
//			       fFunc(f),
//			       fPar(((par) ? par : f->GetParameters())),
//			       fAbsVal(useAbsVal),
//			       fN(n),
//			       fX0(x0)
//			    {
//			       fFunc->InitArgs(fX, fPar);
//			       if (par) fFunc->SetParameters(par);
//			    }
//		
//			ROOT::Math::IGenFunction * Clone()  const 
//			{
//				   TF1_EvalWrapper * f = new TF1_EvalWrapper(*this);
//			       f->fFunc->InitArgs(f->fX, f->fPar);
//			       return f;
//			}
//			Double_t DoEval(Double_t x) const 
//			{
//			       fX[0] = x;
//			       Double_t fval = fFunc->EvalPar(fX, 0);
//			       if (fAbsVal && fval < 0)  return -fval;
//			       return fval;
//			}
//		    Double_t EvalFirstMom(Double_t x) 
//			{
//				fX[0] = x;
//				return fX[0] * TMath::Abs(fFunc->EvalPar(fX, 0));
//			}
//			Double_t EvalNMom(Double_t x) const  
//			{
//			fX[0] = x;
//			return TMath::Power(fX[0] - fX0, fN) * TMath::Abs(fFunc->EvalPar(fX, 0));
//			}
//		TF1 * fFunc;
//		mutable Double_t fX[1];
//		const double * fPar;
//		Bool_t fAbsVal;
//		Double_t fN;
//		Double_t fX0;
//};

//pars [0] - x, pars[1] - A, pars[2] - offset, pars[3] - sigma, pars[4] - alpha, pars[5] - m
double mult_func(double *t, double *pars) //then integrated by t
{
	if (t[0] <= 0) return 0;
	return pars[1] * TMath::Power(pars[5], pars[4] * t[0])*TMath::Exp(-pars[5])*TMath::Gaus((t[0] - pars[0]+pars[2]), 0, pars[3], 0) / (TMath::Gamma(pars[4] * t[0]));
}

//for testing	 integral validity
//double custom_TF1_integrate(TF1* func, double left, double right, double &status)
//{
//	double result;
//	ROOT::Math::GaussIntegrator iod(1.e-12,1.e-12);
//	TF1_EvalWrapper wtf(func, 0, 1);
//	iod.SetFunction(wtf);
//	if (left != -TMath::Infinity() && right != TMath::Infinity())
//	       result = iod.Integral(left, right);
//	       else if (left == -TMath::Infinity() && right != TMath::Infinity())
//		         result = iod.IntegralLow(right);
//	       else if (left != -TMath::Infinity() && right == TMath::Infinity())
//		          result = iod.IntegralUp(left);
//	       else if (left == -TMath::Infinity() && right == TMath::Infinity())
//		          result = iod.Integral();
//	status = iod.Status();
//	if (status)
//	{
//		std::cout << "Gauss Integral issue" << std::endl;
//	}
//	return result;
//}

//for testing integral validity
// pars[0] - A, pars[1] - offset, pars[2] - sigma, pars[3] - alpha, pars[4] - m, pars[5] - x0(discriminator), pars[6] - a(discriminator)
//void write_pars(double x, double *pars, int place)
//{
//	std::ofstream rr;
//	if (place) rr.open(".\\di\\func_val_1.txt", std::ofstream::app);
//	else rr.open(".\\di\\func_val_2.txt", std::ofstream::app);
//	rr <<x<<"\t"<<pars[0] << "\t" << pars[1] << "\t" << pars[2] << "\t" << pars[3] << "\t" << pars[4] << "\t" << pars[5] << std::endl;
//	rr.close();
//}
//double test(double *t, double *pars) //then integrated by t
//{
//	return 1800*0.5*TMath::Power(0.5*t[0],3.7)*TMath::Exp(-0.5*t[0])/TMath::Gamma(4.7);
//}

//pars[0] - A, pars[1] - offset, pars[2] - sigma, pars[3] - alpha, pars[4] - m, pars[5] - x0 (discriminator), pars[6] - a(discriminator)
double signal_func(double* x, double* pars)
{
	double L1, R1, R2, L2;
	double one = 1 / pars[3];
	double second = x[0] - pars[1]+ pars[2] * pars[2] * pars[3] * TMath::Log(pars[4]);
	double interval1 = 5 / pars[3];
	double interval2 = 5* TMath::Abs(pars[2]);
	L1 = one - interval1;
	L1 = (L1 > 0) ? L1:0;
	R1 = one + interval1;
	R1 = (R1 > 0) ? R1 : 0;
	
	L2 = second - interval2;
	L2 = (L2 > 0) ? L2 : 0;
	R2 = second + interval2;
	R2 = (R2 > 0) ? R2 : 0;
	if ((L2 >= L1) && (L2 < R1))
	{
		R1 = R2;
		L2 = 0;
		R2 = 0;
	}
	if ((L1 >= L2) && (L1 < R2))
	{
		L1 = L2;
		R2 = 0;
		L2 = 0;
	}
	//extra: Gamma(alpha*t)~alpha*t^(alpha*t) < 10^12 - or integral too small
	R2 = (R2>(11 / pars[3])) ? 11 / pars[3] : R2;
	L2 = (L2>(11 / pars[3])) ? 11 / pars[3] : L2;
	double ret1 = 0;
	double ret2 = 0;
	if (R1 > L1)
	{
		TF1 *func = new TF1("signal", mult_func, L1, R1, 6);
		func->SetParameters(x[0], pars[0], pars[1], pars[2], pars[3], pars[4]);
		ret1 = func->Integral(L1, R1);
		delete func;
	}
	if (R2>L2)
	{
		TF1 *func = new TF1("signal", mult_func, L2, R2, 6);
		func->SetParameters(x[0], pars[0], pars[1], pars[2], pars[3], pars[4]);
		ret2 = func->Integral(L2, R2);
		delete func;
	}
	return (ret1+ret2)/(1+TMath::Exp((pars[5]-x[0])/pars[6]));
	//return (ret1 + ret2);
}

//pars[0] - A, pars[1] - offset, pars[2] - sigma, pars[3] - alpha, pars[4] - m
double signal_func_old(double* x, double* pars)
{
	double L1, R1, R2, L2;
	double one = 1 / pars[3];
	double second = x[0] - pars[1] + pars[2] * pars[2] * pars[3] * TMath::Log(pars[4]);
	double interval1 = 5 / pars[3];
	double interval2 = 5 * TMath::Abs(pars[2]);
	L1 = one - interval1;
	L1 = (L1 > 0) ? L1 : 0;
	R1 = one + interval1;
	R1 = (R1 > 0) ? R1 : 0;

	L2 = second - interval2;
	L2 = (L2 > 0) ? L2 : 0;
	R2 = second + interval2;
	R2 = (R2 > 0) ? R2 : 0;
	if ((L2 >= L1) && (L2 < R1))
	{
		R1 = R2;
		L2 = 0;
		R2 = 0;
	}
	if ((L1 >= L2) && (L1 < R2))
	{
		L1 = L2;
		R2 = 0;
		L2 = 0;
	}
	//extra: Gamma(alpha*t)~alpha*t^(alpha*t) < 10^12 - or integral too small
	R2 = (R2>(11 / pars[3])) ? 11 / pars[3] : R2;
	L2 = (L2>(11 / pars[3])) ? 11 / pars[3] : L2;
	double ret1 = 0;
	double ret2 = 0;
	if (R1 > L1)
	{
		TF1 *func = new TF1("signal", mult_func, L1, R1, 6);
		func->SetParameters(x[0], pars[0], pars[1], pars[2], pars[3], pars[4]);
		ret1 = func->Integral(L1, R1);
		delete func;
	}
	if (R2>L2)
	{
		TF1 *func = new TF1("signal", mult_func, L2, R2, 6);
		func->SetParameters(x[0], pars[0], pars[1], pars[2], pars[3], pars[4]);
		ret2 = func->Integral(L2, R2);
		delete func;
	}
	return (ret1 + ret2);
}

double average_signal(double A,double offset, double sigma, double alpha,double m ,double x_min,double x_max)
{
	TF1 *func = new TF1("signal", signal_func_old, x_min, x_max, 5);
	func->SetParameters(A, offset, sigma, alpha, m);
	double ret = func->Mean(x_min,x_max);
	delete func;
	return ret;
}

double average_signal_dis(double A, double offset, double sigma, double alpha, double m, double x0, double a, double x_min, double x_max)
{
	TF1 *func = new TF1("signal", signal_func, x_min, x_max, 7);
	func->SetParameters(A, offset, sigma, alpha, m, x0, a);
	double ret = func->Mean(x_min, x_max);
	delete func;
	return ret;
}

//pars[0,1] - info,pars[2] - offset of background, pars[3] - scale of background pars[4] - x_scale
//pars[5] - A, pars[6] - offset, pars[7] - sigma, pars[8] - alpha, pars[9] - m, pars[10] - x0, par[11]-a
double total_function(double*x, double*pars)
{
	long long temp1, temp2;
	memcpy(&temp1, pars, sizeof(long long)); //size of long long = double =long double= 8 bytes (ON WINDOWS ONLY!!!!!!)
	memcpy(&temp2, pars + 1, sizeof(long long)); //size of long long = double =long double= 8 bytes (ON WINDOWS ONLY!!!!!!)
	temp1 = (temp1 >> 8) & 0xffffffff;
	temp2 = (temp2 >> 8) & 0xffffffff;
	temp1 = temp1 | (temp2<<32);
	extra_info *O_o = reinterpret_cast <extra_info*> (temp1);
	double _y_, _x_,x_scale;
	_x_ = pars[2];
	_y_ = pars[3];
	x_scale = pars[4];
	return approximate(*x, O_o, _x_, _y_,x_scale) +signal_func(x, pars + 5);
	//return signal_func(x, pars + 5);
}

int Read_data(int& N_N, double** x, double** y, char filename_bkg[], char filename_data[], int manual_limits, double x_left, double x_right)
{
	std::ifstream str;
	str.open(filename_bkg, 0);
	if (!str.is_open())
	{
		std::cout << "wrong bkg filename" << std::endl;
		str.close();
		return 1;
	}
	str.close();
	str.open(filename_data, 0);
	if (!str.is_open())
	{
		std::cout << "wrong signal filename" << std::endl;
		return 1;
	}
	int ok = 1;
	double col1, col2;
	do
	{
		str >> col1;
		if (!str.good()) ok = 0;
		str.clear();
		str.ignore(1, 0);

		str >> col2;
		if (!str.good()) ok = 0;
		str.clear();
		if ((ok) && (!manual_limits))
		{
			if (N_N == 0)
			{
				x_right = -col1;
				x_left = -col1;
			}
			else
			{
				x_right = (x_right > -col1) ? x_right : -col1;
				x_left = (x_left > -col1) ? -col1 : x_left;
			}
		}
		if ((x_left > -col1) || (x_right < -col1)) ok = 0;
		str.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	} while (str.good() ? (ok ? (N_N++, 1) : ok = 1) : (str.clear(), N_N += ok, 0));
	N_N++;
	str.seekg(0, str.beg);
	//N_N must be more than 5
	//if (N_N < 6) return;
	if (*x) delete [] *x;
	if (*y) delete [] *y;
	*x = new double[N_N];
	*y = new double[N_N];
	N_N = 0;
	ok = 1;
	do
	{
		str >> col1;
		if (str.fail()) ok = 0;
		str.clear();
		str.ignore(1, 0);

		str >> col2;
		if (str.fail()) ok = 0;
		str.clear();
		if ((x_left > -col1) || (x_right < -col1)) ok = 0;
		str.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	} while (str.good() ? (ok ? ((*x)[N_N] = -col1, (*y)[N_N] = col2, N_N++, 1) : ok = 1)
		: (str.clear(), ok ? ((*x)[N_N] = -col1, (*y)[N_N] = col2, N_N++, 0) : 0, 0));
	std::cout << N_N << std::endl;
	return 0;
}

void Fit_and_Add_to_Draw(extra_info **O_o, TGraph* gr, TMultiGraph* mgr, double x_left, double x_right, double *out_pars) //number of pars is 10+avr=11
{
	TF1 *signal = new TF1("signal", total_function, x_left, x_right, 12);
	double par1, par2;
	long long p_bits1 = reinterpret_cast<long long> (*O_o);
	long long p_bits2 = reinterpret_cast<long long> (*O_o);
	p_bits1 = (p_bits1 & 0x00000000ffffffff) << 8;
	p_bits2 = (p_bits2 & 0xffffffff00000000) >> 24;
	long long tt = 0x3fffff0000000000;
	p_bits1 = p_bits1 | tt;
	p_bits2 = p_bits2 | tt;
	memcpy(&par1, &p_bits1, sizeof(long long));
	memcpy(&par2, &p_bits2, sizeof(long long));
	double A = 2.82e+5;
	double offset = 0.462;
	double sigma = 0.527;
	double alpha = 30; //(>0)
	double m = 1.414;
	double x0 = 0.235;
	double a = 0.0128;
	signal->SetParNames("(reserved)", "(reserved)", "bkg_x_off", "bkg_y_scale", "bkg_x_scale", "A", "signal_off", "sigma", "alpha", "m");
	signal->SetParName(10, "x0");
	signal->SetParName(11, "a");
	signal->SetParameters(par1, par2, 0.001, 0.1, 0.988, A, offset, sigma, alpha, m);
	signal->SetParameter(10, x0);
	signal->SetParameter(11, a);
	signal->FixParameter(0, par1);
	signal->FixParameter(1, par2);

	TGraph* sig = new TGraph(signal);
	sig->SetMarkerColor(3);
	sig->SetMarkerStyle(8);
	sig->SetLineColor(3);
	sig->SetLineWidth(1);
	//if (!out_pars[10])
	//{
		signal->SetParLimits(11, 0.008, 0.09);//a
		signal->SetParLimits(10, 0.15, 0.33); //x0
		signal->SetParLimits(9, 0.1, 12);//m
		signal->SetParLimits(8, 5, 50);//alpha
		signal->SetParLimits(7, 0.06, 3);//sigma
		signal->SetParLimits(6, -0.1, 0.8);//offset
		signal->SetParLimits(5, 50000, 5000000);//A
		signal->SetParLimits(4, 0.95, 2); //background x scale
		signal->SetParLimits(3, 0.01, 5);//background y scale
		signal->SetParLimits(2, -0.1, 0.2);//background x offset

		//signal->SetParLimits(11, 0.005, 0.05);//a
		//signal->SetParLimits(10, 0.17, 0.3); //x0
		//signal->SetParLimits(9, 0.8, 10);//m
		//signal->SetParLimits(8, 1, 80);//alpha
		//signal->SetParLimits(7, 0.02, 5);//sigma
		//signal->SetParLimits(6, 0, 1);//offset
		//signal->SetParLimits(5, 8000, 5000000);//A
		//signal->SetParLimits(4, 0.95, 1.3); //background x scale
		//signal->SetParLimits(3, 0.2, 10);//background y scale
		//signal->SetParLimits(2, -0.1, 0.6);//background x offset
	//}
	//else //min is known, optimization
	//{
	//	signal->SetParameters(par1, par2, out_pars[0], out_pars[1], out_pars[2], out_pars[3], out_pars[4], out_pars[5], out_pars[6], out_pars[7]);
	//	signal->SetParameter(10, out_pars[8]);
	//	signal->SetParameter(11, out_pars[9]);

	//	signal->SetParLimits(11, 0.008, 0.05);//a
	//	signal->SetParLimits(10, 0.10, 0.3); //x0
	//	signal->SetParLimits(9, out_pars[7] * 0.1, out_pars[7] * 10);//m
	//	signal->SetParLimits(8, out_pars[6] * 0.1, out_pars[6] * 10);//alpha
	//	signal->SetParLimits(7, out_pars[5] * 0.1, out_pars[5] * 10);//sigma
	//	signal->SetParLimits(6, out_pars[4] * 0.1, out_pars[4] * 10);//offset
	//	signal->SetParLimits(5, out_pars[3] * 0.1, out_pars[3] * 10);//A
	//	signal->SetParLimits(4, out_pars[2] * 0.1, out_pars[2] * 10); //background x scale
	//	signal->SetParLimits(3, out_pars[1] * 0.1, out_pars[1] * 10);//background y scale
	//	signal->SetParLimits(2, out_pars[0] * 0.1, out_pars[0] * 10);//background x offset
	//}

	std::cout << "starting fitting" << std::endl;
	gr->Fit(signal);

	A = signal->GetParameter(5);
	offset = signal->GetParameter(6);
	sigma = signal->GetParameter(7);
	alpha = signal->GetParameter(8);
	m = signal->GetParameter(9);
	x0 = signal->GetParameter(10);
	a = signal->GetParameter(11);
	double x_avr_left = x_left;
	double x_avr_right = x_right;
	x_avr_left = -3.0;
	x_avr_right = 3.7;
	double avr = average_signal(A, offset, sigma, alpha, m, x_avr_left, x_avr_right);
	out_pars[0] = signal->GetParameter(2);
	out_pars[1] = signal->GetParameter(3);
	out_pars[2] = signal->GetParameter(4);
	out_pars[3] = signal->GetParameter(5);
	out_pars[4] = signal->GetParameter(6);
	out_pars[5] = signal->GetParameter(7);
	out_pars[6] = signal->GetParameter(8);
	out_pars[7] = signal->GetParameter(9);
	out_pars[8] = signal->GetParameter(10);
	out_pars[9] = signal->GetParameter(11);
	out_pars[10] = avr;
	out_pars[11] = signal->GetChisquare();
	out_pars[12] = signal->GetProb();
	//std::cout << "average x (discriminator on): " << average_signal_dis(A, offset, sigma, alpha, m,x0,a, x_left, x_right) << std::endl;
	std::cout << "average x (discriminator off): " << avr<< std::endl;
	std::ofstream strt;
	std::cout << "fitted" << std::endl;
	std::cout << "x_offset: " << signal->GetParameter(2) << std::endl;
	std::cout << "y scale: " << signal->GetParameter(3) << std::endl;
	std::cout << "x scale: " << signal->GetParameter(4) << std::endl;
	std::cout << "A: " << signal->GetParameter(5) << std::endl;
	std::cout << "sig offset: " << signal->GetParameter(6) << std::endl;
	std::cout << "sigma: " << signal->GetParameter(7) << std::endl;
	std::cout << "alpha: " << signal->GetParameter(8) << std::endl;
	std::cout << "m: " << signal->GetParameter(9) << std::endl;
	std::cout << "x0: " << signal->GetParameter(10) << std::endl;
	std::cout << "a: " << signal->GetParameter(11) << std::endl;
	if (mgr)
	{
		mgr->Add(gr);

		TF1 *background = new TF1("background", back_func, x_left, x_right, 5);
		background->SetParameters(par1, par2, signal->GetParameter(2), signal->GetParameter(3), signal->GetParameter(4));
		background->FixParameter(0, par1);
		background->FixParameter(1, par2);
		TGraph* bkg = new TGraph(background);
		bkg->SetMarkerColor(1);
		bkg->SetMarkerStyle(3);
		bkg->SetLineColor(1);
		bkg->SetLineWidth(2);
		mgr->Add(bkg);

		TF1 *pure_signal = new TF1("signal with discriminator", signal_func, x_left, x_right, 7);
		pure_signal->SetParameters(A, offset, sigma, alpha, m, x0, a);
		TGraph* sigga = new TGraph(pure_signal);
		sigga->SetMarkerColor(6);
		sigga->SetMarkerStyle(8);
		sigga->SetLineColor(6);
		sigga->SetLineWidth(2);
		mgr->Add(sigga);

		TF1 *no_discr_sig = new TF1("signal", signal_func_old, x_left, x_right, 5);
		no_discr_sig->SetParameters(A, offset, sigma, alpha, m);
		TGraph* no_dis = new TGraph(no_discr_sig);
		no_dis->SetMarkerColor(3);
		no_dis->SetMarkerStyle(34);
		no_dis->SetLineColor(3);
		no_dis->SetLineWidth(2);
		mgr->Add(no_dis);

		gr->SetName("data");
		bkg->SetName("background");
		sigga->SetName("signal with discriminator");
		no_dis->SetName("signal");
	}
}

void create_bkg_sim(char filename_genuine_bkg[], char filename_bogus_bkg[])
{
	TRandom *rnd_obj = new TRandom2();
	rnd_obj->SetSeed();
	std::ifstream ge;
	std::ofstream bo;
	ge.open(filename_genuine_bkg, 0);
	bo.open(filename_bogus_bkg, std::ofstream::trunc);
	double col1, col2;
	int ok = 1;
	do
	{
		ge >> col1;
		if (!ge.good()) ok = 0;
		ge.clear();
		ge.ignore(1, 0);
		ge >> col2;
		if (!ge.good()) ok = 0;
		ge.clear();
		if (ok)
		{
			bo << col1 << ',' << ((col2>0)? rnd_obj->Poisson(col2):0 )<< std::endl;
		}
		ge.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	} while (ge.good() ? (ok ? 1 : ok = 1) : (ge.clear(), 0));
	ge.close();
	bo.close();
	delete rnd_obj;
}

void test_data_simulation(int step, int N_N, double* signal_sim)
{
	if (!WRITE_STAT) return;
	std::string base(".\\di\\signal_");
	std::string suffix(".txt");
	std::string num = std::to_string(step);
	base.append(num);
	base.append(suffix);
	std::ofstream str;
	str.open(base, std::ofstream::trunc);
	for (int g = 0; g < N_N; g++)
	{
		str <<-signal_sim[g] << "\t" << signal_sim[g + N_N] << std::endl;
	}
	str.close();
}

void test_integral(int N_N, double *x,double x_left,double x_right, extra_info **O_o)
{
	TF1 *signal = new TF1("signal", total_function, x_left, x_right, 12);
	double par1, par2;
	long long p_bits1 = reinterpret_cast<long long> (*O_o);
	long long p_bits2 = reinterpret_cast<long long> (*O_o);
	p_bits1 = (p_bits1 & 0x00000000ffffffff) << 8;
	p_bits2 = (p_bits2 & 0xffffffff00000000) >> 24;
	long long tt = 0x3fffff0000000000;
	p_bits1 = p_bits1 | tt;
	p_bits2 = p_bits2 | tt;
	memcpy(&par1, &p_bits1, sizeof(long long));
	memcpy(&par2, &p_bits2, sizeof(long long));
	double A = 2.8e+5;
	double offset = 0;
	double sigma = 0.02; 
	double alpha = 1; 
	double m = 0.5;
	double x0 = 0.2;
	double a = 0.03;
	for (int JJ = 4; JJ < 116; JJ++) //~300 000 000 calls
	{
		std::cout << "Test completed " << std::setprecision(3) << (100.0*JJ / 116) <<"%"<< std::endl;
		m = 0.5 + JJ*0.1;
		for (int K = 0; K <201 ; K++)
		{
			alpha = 1 + K*0.5;
			for (int KK = 0; KK < 500; KK++)
			{
				sigma = 0.02 + KK*0.02;
				for (int L = 0; L < 21; L++)
				{
					offset = L*0.05;
					signal->SetParameters(par1, par2, 0.00, 1, 1, A, offset, sigma, alpha, m);
					for (int J = 0; J < N_N; J++)
					{
						signal->Eval(x[J]);
					}
				}
			}
		}
	}
}

//IN sizeof (data) = N_N*2 +10, coz data = (x united with y ) + initial fit parameters for the sake of optimization
//OUT size of data = 13 - number of parameters + average + Chi square + fit probability
void core(int do_draw, int n_data, double** data, char filename_bkg[], char filename_data[])
{
	int Number_of_simulations = NUMBER_OF_SIMULATIONS;
	int manual_limits = MANUAL_LIMITS;
	int error_estimation = ERR_ESTIMATE;
	int MODE = MODE_; //0 - only for testing some things (see corresponding code)
	if (!do_draw) MODE = 1;
	if (n_data != 0)
	{
		manual_limits = 1;
		error_estimation = 0;
		MODE = 1;
	}
	//
	double x_left=0, x_right=0;
	if (manual_limits)
	{
		x_left =-0.30;
		x_right = 2.3;
	}//
	double *x = 0, *f = 0; //f - from file or - in the estimation - generated data
	int  N_N = 0;
	if (!n_data)
	{
		if (Read_data(N_N, &x, &f, filename_bkg, filename_data, manual_limits, x_left, x_right)) return;
	}
	else
	{
		std::ifstream str;
		str.open(filename_bkg, 0);
		if (!str.is_open())
		{
			std::cout << "wrong bkg filename" << std::endl;
			str.close();
			return;
		}
		str.close();
		N_N = n_data;
		x = new double[N_N];
		f = new double[N_N];
		for (int kl = 0; kl < N_N; kl++)
		{
			x[kl] = (*data)[kl];
			f[kl] = (*data)[kl+N_N];
		}
	}
	TApplication *theApp = NULL;
	TCanvas *can = NULL;
	TMultiGraph *mgr=NULL;
	if (do_draw)
	{
		theApp = new TApplication("App", 0, 0);//I hope it is deleted automatically// not really in fact
		can = new TCanvas("can", filename_data, 1900 - 1200, 1000 - 900, 1200, 900);
		can->SetGrid();
		can->cd();
		mgr = new TMultiGraph();
	}
	TGraph* gr = new TGraph(N_N, x, f);
	gr->SetMarkerColor(4);
	gr->SetLineColor(4);
	gr->SetLineWidth(1);
	gr->SetMarkerStyle(21);
	gr->SetTitle("data");
	extra_info *O_o = new extra_info;
	O_o->num_of_pars = 5;
	O_o->filename = filename_bkg;
	O_o->last_filename = " ";
	O_o->f_vals = 0;
	O_o->x_vals = 0;
	O_o->Poly = 0;
	O_o->last_x = 0;
	O_o->last_x_ = 0;
	O_o->last_x_scale = 1;
	if (!MODE)
	{
		mgr->Add(gr);
		TF1 *background = new TF1("background", back_func, x_left, x_right, 5);
		double par1, par2;
		long long p_bits1 = reinterpret_cast<long long> (O_o);
		long long p_bits2 = reinterpret_cast<long long> (O_o);
		p_bits1 = (p_bits1 & 0x00000000ffffffff) << 8;
		p_bits2 = (p_bits2 & 0xffffffff00000000) >> 24;
		long long tt = 0x3fffff0000000000;
		p_bits1 = p_bits1 | tt;
		p_bits2 = p_bits2 | tt;
		memcpy(&par1, &p_bits1, sizeof(long long));
		memcpy(&par2, &p_bits2, sizeof(long long));
		background->SetParameters(par1,par2, 0.00, 0.42,1);
		background->FixParameter(0, par1);
		background->FixParameter(1, par2);
		//pars[0] - A, pars[1] - offset, pars[2] - sigma, pars[3] - alpha, pars[4] - m ,pars[5]-x0, pars[6]-a
		TF1 *signal = new TF1("signal", signal_func, x_left, x_right, 7);
		//TF1 *signal = new TF1("signal", signal_func_old, x_left, x_right, 5);
		double A = 2.8e+5;
		double offset = 0;
		double sigma = 0.161;
		double alpha = 20; //(>0)
		double m = 2;
		double x0 = -1.;
		double a = 0.02;

		signal->SetParameters(A, offset, sigma, alpha, m,x0,a);
		//signal->SetParameters(A, offset, sigma, alpha, m);
		std::cout << "average x: " << average_signal(A, offset, sigma, alpha, m, x_left, x_right) << std::endl;
		std::ofstream strt;
		strt.open("test_in.txt", std::ios_base::trunc);
		strt.close();

		TF1 *sum = new TF1("total", total_function, x_left, x_right, 12);
		//TF1 *sum = new TF1("total", total_function, x_left, x_right, 10);
		sum->SetParameters(par1, par2, -0.023, 0.31, 1.025, A, offset, sigma, alpha, m);
		sum->SetParameter(10,x0);
		sum->SetParameter(11,a);
		sum->FixParameter(0, par1);
		sum->FixParameter(1, par2);
		TGraph* tot = new TGraph(sum);
		tot->SetMarkerColor(2);
		tot->SetMarkerStyle(8);
		tot->SetLineColor(2);
		tot->SetLineWidth(2);
		mgr->Add(tot);
		TGraph* sig = new TGraph(signal);
		sig->SetMarkerColor(3);
		sig->SetMarkerStyle(8);
		sig->SetLineColor(3);
		sig->SetLineWidth(2);
		mgr->Add(sig);
		TGraph* bkg = new TGraph(background);
		bkg->SetMarkerColor(1);
		bkg->SetMarkerStyle(3);
		bkg->SetLineColor(1);
		bkg->SetLineWidth(2);
		mgr->Add(bkg);
		gr->SetName("data");
		bkg->SetName("background");
		sig->SetName("signal");
		tot->SetName("signal+background");
	}
	else
	{
		double *result_params=new double [13];
		result_params[10] = n_data;
		//
		result_params[10] = 0;//no "optimization"
		if (n_data)
		{
			for (int jjj = 0; jjj < 10; jjj++)
			{
				result_params[jjj] = (*data)[jjj + 2 * N_N];
			}
		}
		Fit_and_Add_to_Draw(&O_o,gr,mgr,x_left,x_right,result_params);
		if (error_estimation)
		{
			//void test_integral(int N_N, double *x, double x_left, double x_right, extra_info **O_o);
			//test_integral(N_N, x,x_left, x_right,&O_o);
			std::ofstream avr_values;
			avr_values.open(".\\di\\avr_statistic.txt", std::ofstream::trunc);
			avr_values << "avr: " <<"\t"<<"Chi^2"<< "\t" << "x_offset" << "\t" << "y_scale" << "\t" << "x_scale" << "\t" << "A" << "\t" << "signal offset" << "\t" << "sigma" << "\t" 
				<< "alpha" << "\t" << "m" << "\t" << "x0" << "\t" <<"a"<< std::endl;
			avr_values << result_params[10] <<"\t"<<result_params[11]<<"\t" << result_params[0] << "\t" << result_params[1] << "\t" << result_params[2] << "\t" << result_params[3] 
				<< "\t" << result_params[4] << "\t" << result_params[5] << "\t" << result_params[6] << "\t" << result_params[7] << "\t" << result_params[8] << "\t" << result_params[9]
				<< std::endl; 
			
			double *avr_sim = new double[Number_of_simulations];
			double *chi_sim = new double[Number_of_simulations];
			double * orig = new double[N_N * 2];
			for (int jjj = 0; jjj < N_N; jjj++)
			{
				orig[jjj]=x[jjj];
				orig[jjj+N_N] = f[jjj];
			}
			test_data_simulation(0, N_N, orig);
			delete orig;
			for (int Monte = 0; Monte < Number_of_simulations; Monte++)
			{
				double * signal_sim = new double[2 * N_N+10];
				char bkg_sim[40] = ".\\di\\bkg_sim.csv";
				TRandom* rand_obj = new TRandom2();
				rand_obj->SetSeed();
				for (int ooo = 0; ooo < N_N; ooo++)
				{ 
					signal_sim[ooo] = x[ooo];
					//average between 5 points as rough estimate for real value; take distribution with the highest dispersion in that point; and signal positive obv.
					//and gaus works only for signal above disriminator (par[8]), not for background.
					double rrr = ((ooo < 2) || (ooo >(N_N - 3)/*edges - only 1 point*/) || (ooo[x]<(result_params[8] + result_params[9])/*bkg is more precise, no need in smoothing*/)) 
						? f[ooo] : ((abs(((f[ooo - 2] + f[ooo - 1] + f[ooo] + f[ooo + 1] + f[ooo + 2]) / 5) - f[ooo])<1.5*DATA_ERROR*f[ooo]) ?
						((f[ooo - 2] + f[ooo - 1] + f[ooo] + f[ooo + 1] + f[ooo + 2]) / 5):
						((abs((f[ooo - 1] + f[ooo] + f[ooo + 1]) / 3) - f[ooo])<1.5*DATA_ERROR*f[ooo]) ?
						((f[ooo - 1] + f[ooo] + f[ooo + 1]) / 3) : f[ooo]);
					signal_sim[ooo + N_N] = (rrr > 0) ? (((ooo[x] > (result_params[8]+result_params[9])) && (rrr * DATA_ERROR > sqrt(rrr))) ? rand_obj->Gaus(rrr, rrr* DATA_ERROR) 
						: rand_obj->Poisson(rrr)) : 0;
					signal_sim[ooo + N_N] = (signal_sim[ooo + N_N] > 0) ? signal_sim[ooo + N_N] : 0;
				}
				for (int jjj = 0; jjj < 10;jjj++)
				{
					signal_sim[jjj + 2 * N_N] = result_params[jjj];
				}
				delete rand_obj;
				test_data_simulation(Monte+1, N_N, signal_sim); // TEMP
				create_bkg_sim(filename_bkg, bkg_sim);
				core(0, N_N, &signal_sim, bkg_sim, filename_data); //<-----------------------------------------
				avr_values << signal_sim[10] << "\t" << signal_sim[11] << "\t" << signal_sim[0] << "\t" << signal_sim[1] << "\t" << signal_sim[2] << "\t" << signal_sim[3]
					<< "\t" << signal_sim[4] << "\t" << signal_sim[5] << "\t" << signal_sim[6] << "\t" << signal_sim[7] << "\t" << signal_sim[8] << "\t" << signal_sim[9]
					<< std::endl;
				avr_sim[Monte] = signal_sim[10];
				chi_sim[Monte] = signal_sim[11];
				delete signal_sim;
			}
			avr_values.close();
			//calculating error of average
			int N = 0;
			double avr_avr_sim=0;
			double squares_avr_sim = 0;
			for (int tt = 0; tt < Number_of_simulations; tt++)
			{
				if (chi_sim[tt] <11*result_params[11]) //simulation acseptance criterion
				{
					N++;
					avr_avr_sim += avr_sim[tt];
				}
			}
			avr_avr_sim = avr_avr_sim / N;
			for (int tt = 0; tt < Number_of_simulations; tt++)
			{
				if (chi_sim[tt] < 11 * result_params[11]) //simulation acseptance criterion
					squares_avr_sim += (avr_sim[tt] - avr_avr_sim)*(avr_sim[tt] - avr_avr_sim);
			}
			std::cout << "Avrerage: " << result_params[10] << std::endl;
			std::cout << "Avr error: " << sqrt(squares_avr_sim / N)<<std::endl;
			std::cout << "N accepted = " << N << " of " << Number_of_simulations << std::endl;
			//calculating error of average
			delete[] result_params;
			delete avr_sim;
			delete chi_sim;
		}
		else
		{
			if (data)
			{
				delete[] * data;
				(*data) = result_params;
			}
			else delete[] result_params;
		}
	}
	if (do_draw)
	{
		Beep(1500,300);
		Sleep(150);
		Beep(1500, 300);
		Sleep(150);
		Beep(1500, 300);
		std::cout << "draw" << std::endl;
		//mgr->SetTitle("Fitting with background fixed (left slope)");
		//mgr->SetTitle("result of fitting with background fixed (left slope)\n and signal fitted with right slope");
		//mgr->SetTitle("Fitting with arbitrary background and signal");
		//mgr->SetTitle("Manual fitting");
		//mgr->SetTitle("Fitting with signal fitted with right slope");
		mgr->SetTitle("Data fitting");
		//mgr->SetTitle("Data fitting with discriminator");
		mgr->Draw("AP"); //n_data!=0 -> estimation, no need to draw
		can->BuildLegend();
		can->Update();
		can->Paint();
		can->Draw();
		can->Modified();
		theApp->Run();
		can->Close();
	}
	std::cout << "deleting" << std::endl;
	if (O_o->f_vals) delete O_o->f_vals;
	if (O_o->x_vals) delete O_o->x_vals;
	if (O_o->Poly) delete O_o->Poly;
	delete O_o;
	if (mgr) delete mgr;
	if (gr) delete gr;	
	if (can) delete can;
	return;
}
