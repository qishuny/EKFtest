
#include <iostream>
// #include "EKFEstimator.h"
#include "EKFEstimator.cpp"


int main()
{    
    EKFEstimator ekf;
    ekf.Init();
    ekf.updateStep();
}
