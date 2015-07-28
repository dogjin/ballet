#include <cstdio>
#include <ballet/object.h>

int main(int argc, char *argv[])
{
    
    bool test_passed = true;

    ballet::balletObject obj;
    if (obj.isLocked())
    {
        test_passed = false;
    }

    obj.release();

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
