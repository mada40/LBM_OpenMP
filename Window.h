#ifndef __WINDOW__
#define __WINDOW__

#include <SFML/Graphics.hpp>
#include"BoltzmannSimulation.h"
#include "omp.h"
class Window
{

private:
	sf::RenderWindow window;
	int const NUMBER_OF_COLUMNS;
	int const NUMBER_OF_ROWS;
	int const GIRD_SIZE;

	sf::Uint8* pixels;
	sf::Texture texture;
	sf::Sprite sprite; // needed to draw the texture on screen

	BoltzmannSumulation sumulation;

	sf::Clock clock;
	float lastTime;


	void draw()
	{
		int W = NUMBER_OF_COLUMNS * GIRD_SIZE;
		int H = NUMBER_OF_ROWS * GIRD_SIZE;

		const float* ux = sumulation.get_ux();
		const float* uy = sumulation.get_uy();
		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int coor = 0; coor < W*H; coor++)
		{
			int row = (coor / W) / GIRD_SIZE;
			int col = (coor % W) / GIRD_SIZE;
			int i = row * NUMBER_OF_COLUMNS + col;
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

public:
	Window(int gird_size, BoltzmannSumulation bs) : NUMBER_OF_COLUMNS(bs.get_W()), NUMBER_OF_ROWS(bs.get_H()), GIRD_SIZE(gird_size), sumulation(bs)
	{
		int W = NUMBER_OF_COLUMNS * GIRD_SIZE;
		int H = NUMBER_OF_ROWS * GIRD_SIZE;
		window.create(sf::VideoMode(W, H), "LBM!");
		pixels = new sf::Uint8[W * H * 4];
		texture.create(W, H);
		sprite = sf::Sprite(texture);

		lastTime = 0.0F;
	}

	bool isOpen()
	{
		return window.isOpen();
	}

	void check_event()
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
	}

	void update()
	{
		sumulation.update();
	}

	void update_FPS()
	{
		float currentTime = clock.getElapsedTime().asMilliseconds();
		float deltatime = (currentTime - lastTime);
		float fps = 1000.f / deltatime;
		lastTime = currentTime;
		window.setTitle(std::to_string(fps));
	}

	void show()
	{
		draw();
		texture.update(pixels);
		window.draw(sprite);
		window.display();
	}

	void close()
	{
		window.close();
	}

	~Window()
	{
		if(pixels) delete[] pixels;
	}

};
#endif