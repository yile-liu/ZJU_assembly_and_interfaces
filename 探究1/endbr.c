#include <stdio.h>

void foo()
{
    puts("inside foo");
    return;
}

int main()
{
    int (*func_p)() = foo;
    (*func_p)();
    // (*(func_p + 0X04))();
    // (*(func_p + 0x05))();

    return 0;
}