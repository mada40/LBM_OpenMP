#include"BoltzmannSimulation.h"
#include"Window.h";
int const H = 200, W = 400;
int const GIRD_SIZE = 1;

int main(int argc, char** argv)
{    
    BoltzmannSumulation bs(W, H);

    Window window(GIRD_SIZE, bs);
    
    while (window.isOpen())
    {
        window.show();
    }

    return 0;
}