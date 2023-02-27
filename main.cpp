#include"BoltzmannSimulation.h"
#include"Window.h";
#include<iostream>
#include <sstream>
// string = "W H GIRD_SIZE NUMBER_OF_ITERATIONS NEED_IMAGE = false"
int main(int argc, char** argv)
{    

    int W = 600;
    int H = 300;
    int GIRD_SIZE = 2;
    int NUMBER_OF_ITERATIONS = 5000;
    bool NEED_IMAGE = false;

    if (argc != 1 && argc != 5 && argc != 6)
        return -1;


    if (argc > 1)
    {
        W = atoi(argv[1]);
        H = atoi(argv[2]);
        GIRD_SIZE = atoi(argv[3]);
        NUMBER_OF_ITERATIONS = atoi(argv[4]);

        if (argc == 6)
        {
            std::stringstream ss(argv[5]);
            ss >> std::boolalpha >> NEED_IMAGE;
        }
    }



    

    std::cout << "W: " << W << "\n";
    std::cout << "H: " << H << "\n";
    std::cout << "GIRD SIZE: " << GIRD_SIZE << "\n";
    std::cout << "NUMBER OF ITERATIONS: " << NUMBER_OF_ITERATIONS << "\n";
    std::cout << "MODE: " << (NEED_IMAGE ? "Window" : "Console") << "\n";

    std::cout << "\n\t\t--RUN--\n\n";

    if (NEED_IMAGE)
    {
        BoltzmannSumulation bs(W, H);
        Window window(GIRD_SIZE, bs);
        for (int i = 0; i < NUMBER_OF_ITERATIONS && window.isOpen(); i++)
        {
            window.check_event();
            window.update();
            window.show();
            window.update_FPS();
        }

        while (window.isOpen())
        {
            window.check_event();
            window.show();
        }
    }
    else
    {
        BoltzmannSumulation bs(W, H);
        for (int i = 0; i < NUMBER_OF_ITERATIONS; i++)
        {
            bs.update();
            if (i % 1000 == 0)
                std::cout << i <<  "/" << NUMBER_OF_ITERATIONS << " ";
        }
        std::cout << "\n";
    }
    std::cout << "---END---";

    return 0;
}