#ifndef __TBoltzmannSumulation_H__
#define __TBoltzmannSumulation_H__

#include <assert.h>
#include <iostream>
#include <string>
#include <cmath>
#include "omp.h"

template<class T> 
inline T clamp(T x, T min, T max)
{
	return std::max(min, std::min(x, max));
}

const int NUM_THREADS = 8;

const int NL = 9;
const int dxs[] = { 0, 0, 1, 1, 1, 0,-1,-1,-1 };
const int dys[] = { 0, 1, 1, 0,-1,-1,-1, 0, 1 };
const double weights[] = { 4.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36, 1.0 / 9, 1.0 / 36 };
const int invert[] = { 0, 5,6,7,8,1,2,3,4 };



class BoltzmannSumulation
{
private:
	const int Ny;
	const int Nx;
	const double TAU = 0.58;
	float *new_table;
	float *cur_table;
	float *rho;
	float* ux;
	float* uy;
	bool* isBondary;
	//812
	//703
	//654
	

	BoltzmannSumulation& operator = (const BoltzmannSumulation& other) { throw std::logic_error("MNE LEN' DELAT' ETOT OPERATOR"); }
	

	void init_sumulation()
	{
		memset(isBondary, false, Ny * Nx * sizeof(bool));
		for (int i = 0; i < NL * Ny * Nx; i++)
		{
			cur_table[i] = 0.5 + (abs(rand()) % 1000) / 10000.0; 
			int coor = i % (Nx * Ny);
			int y = coor / Nx;
			int x = coor % Nx;
			if (i / (Nx * Ny) == 3)
				cur_table[i] = 1.3;

			int cylX = x - Nx* 2 / 3;
			int cylY = y - Ny / 2;
			int R = Ny / 10;
			if (cylX * cylX + cylY * cylY <= R * R)
			{
				isBondary[coor] = 1;
			}
		}



	}
	void update_stream()
	{
		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i = 0; i < NL * Ny * Nx; i++)
		{
			int l = i / (Ny * Nx);
			int coor = i % (Ny * Nx);
			int y = coor / Nx;
			int x = coor % Nx;
			int new_y = (y + dys[l] + Ny) % Ny;
			int new_x = (x + dxs[l] + Nx) % Nx;
			new_table[i] = cur_table[l * (Ny * Nx) + new_y * Nx + new_x];
		}

		//std::swap(tmp_table, table);	
	}
	void update_collide()
	{
		memset(rho, 0, Nx * Ny * sizeof(float));
		memset(ux, 0, Nx * Ny * sizeof(float));
		memset(uy, 0, Nx * Ny * sizeof(float));
	
		for (int l = 0; l < NL; l++)
		{
			#pragma omp parallel for  num_threads(NUM_THREADS)
			for (int coor = 0; coor < Ny * Nx; coor++)
			{
				int i = l * Nx * Ny + coor;
				rho[coor] = rho[coor] + cur_table[i];
				ux[coor] = ux[coor] + dxs[l] * cur_table[i];
				uy[coor] = uy[coor] + dys[l] * cur_table[i];
				if (isBondary[coor])
				{
					cur_table[i] = cur_table[invert[l] * Ny * Nx + coor];
				}
			}
		}
		

		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i = 0; i < Nx * Ny; i++)
		{
			ux[i] /= rho[i];
			uy[i] /= rho[i];
			if (isBondary[i])
			{
				ux[i] = 0.0;
				uy[i] = 0.0;
			}
		}

		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i = 0; i < NL * Ny * Nx; i++)
		{
			int l = i / (Ny * Nx);
			int coor = i % (Ny * Nx);
			float tmp = dxs[l] * ux[coor] + dys[l] * uy[coor];
			float sqU = ux[coor] * ux[coor] + uy[coor] * uy[coor];
			double Feq = rho[coor] * weights[l] * (1.0 + 3.0 * tmp + 4.5 * tmp * tmp - 1.5 * sqU);
			cur_table[i] = cur_table[i] - (cur_table[i] - Feq) / TAU;
		}
	}

public:
	
	int get_W() { return Nx; }
	int get_H() { return Ny; }

	float* get_ux() { return ux; }
	float* get_uy() { return uy; }


	BoltzmannSumulation(int _Nx, int _Ny) : Nx(_Nx), Ny(_Ny)
	{
		new_table = new float[NL * Ny * Nx];
		cur_table = new float[NL * Ny * Nx];
		rho = new float[Ny * Nx];
		ux = new float[Ny * Nx];
		uy = new float[Ny * Nx];
		isBondary = new bool[Ny * Nx];
		init_sumulation();
	}

	BoltzmannSumulation(const BoltzmannSumulation& bs) : Nx(bs.Nx), Ny(bs.Ny), TAU(bs.TAU) 
	{ 
		float* tmp_new_table = new float[NL * Ny * Nx];
		float* tmp_cur_table = new float[NL * Ny * Nx];
		float* tmp_rho = new float[Ny * Nx];
		float* tmp_ux = new float[Ny * Nx];
		float* tmp_uy = new float[Ny * Nx];
		bool* tmp_isBondary = new bool[Ny * Nx];

		memcpy(tmp_new_table, bs.new_table, NL * Ny * Nx * sizeof(float));
		memcpy(tmp_cur_table, bs.cur_table, NL * Ny * Nx * sizeof(float));	
		memcpy(tmp_rho, bs.rho, Ny * Nx * sizeof(float));
		memcpy(tmp_ux, bs.ux, Ny * Nx * sizeof(float));
		memcpy(tmp_uy, bs.uy, Ny * Nx * sizeof(float));
		memcpy(tmp_isBondary, bs.isBondary, Ny * Nx * sizeof(bool));

		new_table = tmp_new_table;
		cur_table = tmp_cur_table;
		rho = tmp_rho;
		ux = tmp_ux;
		uy = tmp_uy;
		isBondary = tmp_isBondary;
	}

	~BoltzmannSumulation()
	{
		if (cur_table) delete[] cur_table;
		if (new_table) delete[] new_table;
		if (rho) delete[] rho;
		if (ux) delete[] ux;
		if (uy) delete[] uy;
		if (isBondary) delete[] isBondary;
	}

	void update()
	{
		update_collide();
		update_stream();
		std::swap(cur_table, new_table);
	}



};
#endif