#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#define _USE_MATH_DEFINES
#include <math.h>

#define int int
#define double double
#define void void
#define π M_PI

struct extra_info
{
	int num_of_pars;
	std::string filename;
	std::string last_filename;

	double *f_vals;
	double *x_vals;
	int N_P;

	double last_x;
	double *Poly; //size = q_p = 6

	double last_x_;
	double last_x_scale;
};


void build(int ♂, double* x_vals, double* f_vals, double* out_) //finite differences build
{
	for (int OOk = 0; OOk < ♂; OOk++)	OOk[out_] = f_vals[OOk];
	double **flake;
	flake = new double*[♂];
	*flake = out_;
	for (int OOk = 1; OOk < ♂; OOk++)
	{
		OOk[flake] = new double[♂ - OOk];
	}
	for (int OOk = 1; OOk < ♂; OOk++)
	{
		for (int kok = 0; kok < ♂ - OOk; kok++)
		{
			kok[OOk[flake]] = (kok[*(OOk + flake - 1)] - *(kok + 1 + *(OOk + flake - 1))) / (kok[x_vals] - *(kok + OOk + x_vals));
		}
	}
	for (int OOk = 1; OOk < ♂; OOk++)
	{
		OOk[*flake] = *OOk[flake];
		delete OOk[flake];
	}
	delete flake;
	return;
}

double calc(double x, int ♀, double* x_vals, double*FD)
{
	double out_ = (♀ - 1)[FD];
	for (int o_O = 0; o_O < ♀ - 1; o_O++) out_ = (♀ - 2 - o_O)[FD] + (x - x_vals[♀ - 2 - o_O])*out_;
	 return out_;
}

int find_offset(int N_N, double p, double *x_vals, double *f_vals) //array x_vals is assumed to be sorted
{
	if (p <= x_vals[0])  return - 2;
	if (p >= x_vals[N_N - 1])  return - 1;
	int out_ = 0;
	while ((p > out_[x_vals]) && (p > *(out_ + 1 + x_vals))) out_++;
	if (out_ < 2)   return 2;
	if (out_ > N_N - 5)   return  N_N - 5;
	return out_;
}

double approximate(double p, extra_info* S, double _x_, double _y_,double x_scale)
{
	int offset = 0;
	int q_p = 6; //=const N of points used in approximation
	if ((S->last_filename != S->filename) || (x_scale!=S->last_x_scale))
	{
		//parcer
		if (S->f_vals) delete S->f_vals;	//needed for not reading the same file several times if avoidable
		if (S->x_vals) delete S->x_vals;	//needed for not reading the same file several times if avoidable
		if (S->Poly) delete S->Poly;		//needed for not recalculating finite differences if x is near last x
		S->f_vals = 0;
		S->x_vals = 0;
		S->Poly = 0;
		S->last_x = 0;
		S->last_x_ = 0;
		S->last_x_scale = 1;
		std::ifstream str;
		str.open(S->filename,0);
		offset = str.is_open();
		if (!offset)
		{
			std::cout << "file opening error" << std::endl;
			return -1;
		}
		int  N_N = 0;
		int ok = 1;
		double col1, col2;
		do //just counter number of points
		{
			str >> col1;
			if (!str.good()) ok = 0;
			str.clear();
			str.ignore(1, 0);

			str >> col1;
			if (!str.good()) ok = 0;
			str.clear();
			str.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} while (str.good() ? (ok ? (N_N++, 1) : ok = 1) : (str.clear(), N_N+=ok, 0));
		N_N++;
		str.seekg(0, str.beg);
		offset = str.good();
		//N_N must be more than 5
		if (N_N < 6)  return - 1;
		S->last_filename = S->filename;
		S->last_x_scale = x_scale;
		S->N_P = N_N;
		S->x_vals = new double[N_N];
		S->f_vals = new double[N_N];
		//S->Poly = new double[q_p];
		N_N = 0;
		ok = 1;
		do
		{
			str >> col1;
			if (str.fail()) ok = 0;
			str.clear();
			str.ignore(1, 0);
			col1 = col1*x_scale;
			str >> col2;
			if (str.fail()) ok = 0;
			str.clear();
			str.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} while (str.good() ? (ok ? ((S->x_vals)[N_N] = -col1, (S->f_vals)[N_N] = col2, N_N++, 1) : ok = 1) 
			: (str.clear(), ok ? ((S->x_vals)[N_N] = -col1, (S->f_vals)[N_N] = col2,N_N++,0) : 0, 0));
		str.close();
		//end of parcer
		S->last_x = p;
	}
	if (S->last_x_ != _x_) // S->Poly does not depend on x's offset
	{
		for (int ♀ = 0; ♀ < S->N_P; ♀++) 
		{
			(S->x_vals)[♀] += _x_ - S->last_x_;
		}
		S->last_x += _x_ - S->last_x_;
		S->last_x_ = _x_;
	}
	offset = find_offset(S->N_P, p, S->x_vals, S->f_vals);
	if (offset < 0)
	{
		//S->last_x = p;
		return 0; // function supposed to be 0 beyond [x_min ; x_max]
	}
	if ((offset == find_offset(S->N_P, S->last_x, S->x_vals, S->f_vals)) && (NULL != S->Poly))
	{
		 return	_y_*calc(p, q_p, (S->x_vals + offset - 2), S->Poly); //same polynom, not recalculating S->Poly
	}
	//need recalculate S->Poly
	S->last_x = p;
	double* _f_vals, *_x_vals;
	_x_vals = new double[q_p];
	_f_vals = new double[q_p];
	for (int _ = 0; _ < q_p; _++)
	{
		_[_x_vals] = *(S->x_vals + offset + _ - 2);
		*(_f_vals + _) = *(S->f_vals + offset + _ - 2);
	}
	if (S->Poly) delete S->Poly;
	S->Poly = new double[q_p];
	build(q_p, _x_vals, _f_vals, S->Poly);
	delete _x_vals;
	delete _f_vals;
	return _y_*calc(p, q_p, (S->x_vals + offset - 2), S->Poly);
}

double back_func(double *x, double* pars) //pars[0] and pars[1] - pointer to extra_info - contains file and optimization info.
{
	long long temp1, temp2;
	memcpy(&temp1, pars, sizeof(long long)); //size of long long = double =long double= 8 bytes (ON WINDOWS ONLY!!!!!!)
	memcpy(&temp2, pars + 1, sizeof(long long)); //size of long long = double =long double= 8 bytes (ON WINDOWS ONLY!!!!!!)
	temp1 = (temp1 >> 8) & 0xffffffff;
	temp2 = (temp2 >> 8) & 0xffffffff;
	temp1 = temp1 | (temp2 << 32);
	extra_info *O_o = reinterpret_cast <extra_info*> (temp1);
	double _y_, _x_,x_scale;
	switch(O_o->num_of_pars)
	{
	case 5:
	{
		x_scale = pars[4];
		_x_ = pars[2];
		_y_ = pars[3];
		break;
	}
	case 3:
	{
		_x_ = pars[2];
		_y_ = 1;
		x_scale = 1;
		break;
	}
	case 4:
	{
		_x_ = pars[2];
		_y_ = pars[3];
		x_scale = 1;
		break;
	}
	default:
	{
		_x_ = 0;
		_y_ = 1;
		x_scale = 1;
	}
	}
	return approximate(*x, O_o, _x_, _y_,x_scale);
}
