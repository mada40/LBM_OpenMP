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
	double max_speed;
	//812
	//703
	//654
	

	BoltzmannSumulation& operator = (const BoltzmannSumulation& other);
	BoltzmannSumulation(const BoltzmannSumulation& bs): Nx(0), Ny(0), TAU(-1.0) { throw std::logic_error("MNE LEN' DELAT' COPY CONSTRUCTROR"); }

	void init_sumulation()
	{
		memset(isBondary, false, Ny * Nx * sizeof(bool));
		for (int i = 0; i < NL * Ny * Nx; i++)
		{
			cur_table[i] = 1.0 + (abs(rand()) % 1000) / 10000.0; 
			int coor = i % (Nx * Ny);
			int y = coor / Nx;
			int x = coor % Nx;
			if (i / (Nx * Ny) == 3)
				cur_table[i] = 2.5;

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
		

		max_speed = std::sqrt(ux[0] * ux[0] + uy[0] * uy[0]);
		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int i = 0; i < Nx * Ny; i++)
		{
			ux[i] /= rho[i];
			uy[i] /= rho[i];
			max_speed = std::max(max_speed, 1.0*std::sqrt(ux[i] * ux[i] + uy[i] * uy[i]));
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


	void draw(sf::Uint8* pixels, int gird_size)
	{
		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int coor = 0; coor < Nx * Ny * gird_size * gird_size; coor++)
		{
			int row = coor / (gird_size * Nx) / gird_size;
			int col = coor % (gird_size * Nx) / gird_size;
			int i = row * Nx + col;
			double y = uy[i];
			double x = ux[i];
			double speed = std::sqrt(x * x + y * y) * 1250.0;
			int r = -speed * (speed - 384.0) / 128.0;
			int g = -speed * (speed - 256.0) / 64.0;
			int b = -(speed - 256.0) * (speed + 128.0) / 128.0;
			r = clamp(r, 0, 255);
			g = clamp(g, 0, 255);
			b = clamp(b, 0, 255);
			/*r*/pixels[coor * 4 + 0] = r;
			/*g*/pixels[coor * 4 + 1] = g;
			/*b*/pixels[coor * 4 + 2] = b;
			/*a*/pixels[coor * 4 + 3] = 255;
		}
	}


};
#endif