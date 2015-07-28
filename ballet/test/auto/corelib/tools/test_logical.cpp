#include <cstdio>

#include <armadillo>
#include <ballet/logical.h>

using namespace arma;

int main(int argc, char *argv[])
{

    umat A = ones<umat>(10,1);
    int B = ballet::all(A);

    printf("B = %d\n",B);

    A = zeros<umat>(10,1);
    B = ballet::all(A);

    printf("B = %d\n",B);

    return 0;

}
