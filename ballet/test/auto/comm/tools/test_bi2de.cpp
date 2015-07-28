#include <cstdio>
#include <armadillo>
#include <ballet/comm/bi2de.h>

using namespace arma;

int main(int argc, char *argv[])
{
    
    bool test_passed = true;

    // create binary row vector
    // b = [1 1 1 1 1 0 1 1 1 1 1]
    imat b;
    b << 1 << 1  << 1 << 1 << 1 << 0 << 1 
      << 1  << 1 << 1 << 1 << endr;

    // convert binary vector to decimal numbers
    imat d = ballet::bi2de(b);

    // check result
    if (d(0) != 2015)
    {
        test_passed = false;
    }

    if (test_passed)
    {
        printf("test passed\n");
    }
    else
    {
        printf("test failed\n");
    }

    return 0;

}
