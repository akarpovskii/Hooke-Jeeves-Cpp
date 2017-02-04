#include <iostream>
#include "HookeJeevesWrapper.h"

int main() {
#if Debug
    std::cout << HJWrapper::test_exp() << std::endl;
#endif
    return 0;
}