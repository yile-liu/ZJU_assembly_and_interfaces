// canary.c
#include <string.h>

void foo(char *str)
{
    char buffer[10];
    // char buffer[0x10];
    strcpy(buffer, str);
}

int main(int argc, char **argv)
{
    char str[100] = {0};
    // char str[100] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    //                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0};

    // char str[100] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    //              11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    //              21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 0};

    foo(str);
    puts("return successfully");
    return 0;
}