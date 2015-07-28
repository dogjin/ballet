#include <cstdio>

#include <armadillo>
#include <ballet/random.h>

using namespace arma;

int main(int argc, char *argv[])
{

    ivec p = ballet::randperm(10);
    p.print();

    return 0;

}
