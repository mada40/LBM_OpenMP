#include <SFML/Graphics.hpp>
#include"BoltzmannSimulation.h"
int const K = 100;
int const H = 9 * K, W = 16 * K;


int main()
{
    sf::RenderWindow window(sf::VideoMode(W, H), "SFML works!");
    //window.setFramerateLimit(60);
    sf::Clock clock;
    float lastTime = 0;

    sf::Uint8* pixels = new sf::Uint8[W * H * 4];
    sf::Texture texture;
    texture.create(W, H);
    sf::Sprite sprite(texture); // needed to draw the texture on screen
    
    BoltzmannSumulation bs(W, H);
    int cnt = 0;
    while (window.isOpen())
    {
        cnt++;
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        bs.update();
        if (cnt % 1 == 0)
        {

            bs.draw(pixels);
            texture.update(pixels);
            window.draw(sprite);
            window.display();
        }

        float currentTime = clock.getElapsedTime().asMilliseconds();
        float deltatime = (currentTime - lastTime);
        float fps = 1000.f / deltatime;
        lastTime = currentTime;
        window.setTitle(std::to_string(fps));
    }

    delete[] pixels;
    return 0;
}