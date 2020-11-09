#include <stdio.h>
#include <math.h>

/* 待定的二元交互作用参数 */
double beitavij;
double gamavij;
double beitaTij;
double gamaTij;
double Fij;
/* 模型固定参数 */
double modelFixedParameters[10][3] = {
        {-2.4547627E-02, 1, 2},  // nk,dk,tk
        {-2.4120612E-01, 1, 4},
        {-5.1380195E-03, 1, -2},
        {-2.3982483E-02, 2, 1},
        {2.5977234E-01,  3, 4},
        {-1.7201412E-01, 4, 4},
        {4.2949003E-02,  5, 4},
        {-2.0210859E-04, 6, 0},
        {-3.8298423E-03, 6, 4},
        {2.6992331E-06,  8, -2},
};
/* 组分特有的参数 */
double MethaneParameters[14][3] = { //甲烷
        {-2.4547627E-02, 1, 2},  // nk,dk,tk,ck
        {-2.4120612E-01, 1, 4},
        {-5.1380195E-03, 1, -2},
        {-2.3982483E-02, 2, 1},
        {2.5977234E-01,  3, 4},
        {-1.7201412E-01, 4, 4},
        {4.2949003E-02,  5, 4},
        {-2.0210859E-04, 6, 0},
        {-3.8298423E-03, 6, 4},
        {2.6992331E-06,  8, -2},
};
double carbonDioxideParameters[14][3] = { //二氧化碳
        {-2.4547627E-02, 1, 2},  // nk,dk,tk,ck
        {-2.4120612E-01, 1, 4},
        {-5.1380195E-03, 1, -2},
        {-2.3982483E-02, 2, 1},
        {2.5977234E-01,  3, 4},
        {-1.7201412E-01, 4, 4},
        {4.2949003E-02,  5, 4},
        {-2.0210859E-04, 6, 0},
        {-3.8298423E-03, 6, 4},
        {2.6992331E-06,  8, -2},
};

// 二元混合物的 虚拟临界密度的倒数
double mixtureVirtualCriticalDensityReverse(double x1, double criticalDensity1, double criticalDensity2) {
    return (pow(x1, 2) / criticalDensity1 + pow(1 - x1, 2) / criticalDensity2
            + 2 * x1 * (1 - x1) * beitavij * gamavij / (pow(beitavij, 2) * x1 + 1 - x1) * (1 / 8)
              * pow((1 / (pow(criticalDensity1, 1 / 3) + 1 / pow(criticalDensity2, 1 / 3))), 3)
    );
}

// 二元混合物的 虚拟临界温度
double mixtureVirtualCriticalTemperature(double x1, double criticalTemperature1, double criticalTemperature2) {
    return (
            pow(x1, 2) * criticalTemperature1 + pow(1 - x1, 2) * criticalTemperature2
            + 2 * x1 * (1 - x1) * beitaTij * gamaTij / (pow(beitaTij, 2) * x1 + 1 - x1)
              * pow(criticalTemperature1 * criticalTemperature2, 0.5)
    );
};

double sumOfPolynomialOfMixtureReducedHelmholtzExcess(double x1, double Temperature, double density,
                                                      double criticalDensity1, double criticalDensity2,
                                                      double criticalTemperature1, double criticalTemperature2) {
    double sum = 0;
    for (int i = 0; i < sizeof(modelFixedParameters); i++) {
        sum += modelFixedParameters[i][0]
               * pow(density * mixtureVirtualCriticalDensityReverse(x1, criticalDensity1, criticalDensity2),
                     modelFixedParameters[i][1])
               * pow(mixtureVirtualCriticalTemperature(x1, criticalTemperature1, criticalTemperature2) / Temperature,
                     modelFixedParameters[i][2]);
    }
    return sum;
}

// 二元混合物的 超额部分无量度亥姆霍兹自由能
double mixtureReducedHelmholtzExcess(double x1, double Temperature, double density,
                                     double criticalDensity1, double criticalDensity2,
                                     double criticalTemperature1, double criticalTemperature2) {
    return (x1 * (1 - x1) * Fij *
            sumOfPolynomialOfMixtureReducedHelmholtzExcess(x1, Temperature, density,
                                                           criticalDensity1, criticalDensity2,
                                                           criticalTemperature1, criticalTemperature2));
}


// 组分i的 无量度亥姆霍兹自由能的剩余部分
double iReducedHelmholtzResidualPart(char componentName,
                                     double Temperature, double density,
                                     double criticalTemperature, double criticalDensity) {
    double sum = 0;
    if (componentName == 'Methane') {
        for (int i = 0; i < 6; ++i) {
            sum += MethaneParameters[i][0]
                   * pow(density / criticalDensity, MethaneParameters[i][1])
                   * pow(criticalTemperature / Temperature, MethaneParameters[i][2]);
        }
        for (int i = 6; i < 14; ++i) {
            sum += MethaneParameters[i][0]
                   * pow(density / criticalDensity, MethaneParameters[i][1])
                   * pow(criticalTemperature / Temperature, MethaneParameters[i][2])
                   * exp(-pow(density / criticalDensity, MethaneParameters[i][3]));
        }

    } else if (componentName == 'carbonDioxide') {
        for (int i = 0; i < 6; ++i) {
            sum += carbonDioxideParameters[i][0]
                   * pow(density / criticalDensity, carbonDioxideParameters[i][1])
                   * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2]);
        }
        for (int i = 6; i < 14; ++i) {
            sum += carbonDioxideParameters[i][0]
                   * pow(density / criticalDensity, carbonDioxideParameters[i][1])
                   * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2])
                   * exp(-pow(density / criticalDensity, carbonDioxideParameters[i][3]));
        }
    }
    return sum;
}


// 混合物的 无量度亥姆霍兹自由能的剩余部分
double mixtureReducedHelmholtzResidualPart(double x1,
                                           char componentName1, char componentName2,
                                           double Temperature, double density,
                                           double criticalTemperature1, double criticalDensity1,
                                           double criticalTemperature2, double criticalDensity2) {
    return (x1 *
            iReducedHelmholtzResidualPart(componentName1, Temperature, density, criticalTemperature1, criticalDensity1)
            + (1 - x1) *
              iReducedHelmholtzResidualPart(componentName2, Temperature, density, criticalTemperature2,
                                            criticalDensity2)
            + mixtureReducedHelmholtzExcess(x1, Temperature, density, criticalDensity1, criticalDensity2,
                                            criticalTemperature1, criticalTemperature2)
    );
}

// 二元混合物的 超额部分无量度亥姆霍兹自由能 对对比密度的一阶导数
double mixtureReducedHelmholtzExcessDerivativeToContrastDensity(double x1, double Temperature, double density,
                                                                double criticalDensity1, double criticalDensity2,
                                                                double criticalTemperature1,
                                                                double criticalTemperature2) {
    double sum = 0;
    for (int i = 0; i < sizeof(modelFixedParameters); i++) {
        sum += modelFixedParameters[i][0]
               * pow(mixtureVirtualCriticalTemperature(x1, criticalTemperature1, criticalTemperature2) / Temperature,
                     modelFixedParameters[i][2])
               * modelFixedParameters[i][1]
               * pow(density *
                     mixtureVirtualCriticalDensityReverse(x1, criticalDensity1, criticalDensity2),
                     modelFixedParameters[i][1] - 1);
    }
    return sum;
}


// 组分i的 无量度亥姆霍兹自由能的剩余部分 对对比密度的一阶导数
double iReducedHelmholtzResidualPartPartDerivativeToContrastDensity(char componentName,
                                                                    double Temperature, double density,
                                                                    double criticalTemperature,
                                                                    double criticalDensity) {
    double sum = 0;
    if (componentName == 'Methane') {
        for (int i = 0; i < 6; ++i) {
            sum += MethaneParameters[i][0]
                   * MethaneParameters[i][1]
                   * pow(density / criticalDensity, MethaneParameters[i][1] - 1)
                   * pow(criticalTemperature / Temperature, MethaneParameters[i][2]);

        }
        for (int i = 6; i < 14; ++i) {
            sum += MethaneParameters[i][0]
                   * pow(density / criticalDensity, MethaneParameters[i][1] - 1)
                   * (MethaneParameters[i][1] -
                      MethaneParameters[i][3] * pow(density / criticalDensity, MethaneParameters[i][3]))
                   * pow(criticalTemperature / Temperature, MethaneParameters[i][2])
                   * exp(-pow(density / criticalDensity, MethaneParameters[i][3]));
        }
    } else if (componentName == 'carbonDioxide') {
        for (int i = 0; i < 6; ++i) {
            sum += carbonDioxideParameters[i][0]
                   * carbonDioxideParameters[i][1]
                   * pow(density / criticalDensity, carbonDioxideParameters[i][1] - 1)
                   * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2]);

        }
        for (int i = 6; i < 14; ++i) {
            sum += carbonDioxideParameters[i][0]
                   * pow(density / criticalDensity, carbonDioxideParameters[i][1] - 1)
                   * (carbonDioxideParameters[i][1] -
                      carbonDioxideParameters[i][3] * pow(density / criticalDensity, carbonDioxideParameters[i][3]))
                   * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2])
                   * exp(-pow(density / criticalDensity, carbonDioxideParameters[i][3]));
        }
    }
    return sum;
}


// 混合物的 无量度亥姆霍兹自由能的剩余部分 对对比密度的一阶导数
double mixtureReducedHelmholtzResidualPartDerivativeToContrastDensity(double x1,
                                                                      char componentName1, char componentName2,
                                                                      double Temperature, double density,
                                                                      double criticalTemperature1,
                                                                      double criticalDensity1,
                                                                      double criticalTemperature2,
                                                                      double criticalDensity2) {
    return (x1 * iReducedHelmholtzResidualPartPartDerivativeToContrastDensity(componentName1, Temperature, density,
                                                                              criticalTemperature1, criticalDensity1)
            + (1 - x1) *
              iReducedHelmholtzResidualPartPartDerivativeToContrastDensity(componentName2, Temperature, density,
                                                                           criticalTemperature2,
                                                                           criticalDensity2)
            + mixtureReducedHelmholtzExcessDerivativeToContrastDensity(x1, Temperature, density, criticalDensity1,
                                                                       criticalDensity2, criticalTemperature1,
                                                                       criticalTemperature2));
}

// 二元混合物的 超额部分无量度亥姆霍兹自由能 对对比密度的二阶导数
double mixtureReducedHelmholtzExcessSecondDerivativeToContrastDensity(double x1, double Temperature, double density,
                                                                      double criticalDensity1, double criticalDensity2,
                                                                      double criticalTemperature1,
                                                                      double criticalTemperature2) {
    double sum = 0;
    for (int i = 0; i < sizeof(modelFixedParameters); i++) {
        sum += modelFixedParameters[i][0]
               * pow(mixtureVirtualCriticalTemperature(x1, criticalTemperature1, criticalTemperature2) / Temperature,
                     modelFixedParameters[i][2])
               * modelFixedParameters[i][1]
               * (modelFixedParameters[i][1] - 1)
               * pow(density *
                     mixtureVirtualCriticalDensityReverse(x1, criticalDensity1, criticalDensity2),
                     modelFixedParameters[i][1] - 2);
    }
    return sum;
}

// 组分i的 无量度亥姆霍兹自由能的剩余部分 对对比密度的二阶导数
double iReducedHelmholtzResidualPartPartSecondDerivativeToContrastDensity(char componentName,
                                                                          double Temperature, double density,
                                                                          double criticalTemperature,
                                                                          double criticalDensity) {
    double sum = 0;
    if (componentName == 'Methane') {
        for (int i = 0; i < 6; ++i) {
            sum += MethaneParameters[i][0] * MethaneParameters[i][1] * (MethaneParameters[i][1] - 1)
                   * pow(density / criticalDensity, MethaneParameters[i][1] - 2)
                   * pow(criticalTemperature / Temperature, MethaneParameters[i][2]);
        }
        for (int i = 6; i < 14; ++i) {
            MethaneParameters[i][0]
            * pow(density / criticalDensity, MethaneParameters[i][1] - 2)
            *
            ((MethaneParameters[i][1] -
              MethaneParameters[i][3] * pow(density / criticalDensity, MethaneParameters[i][3])) *
             (MethaneParameters[i][1] - 1 -
              MethaneParameters[i][3] * pow(density / criticalDensity, MethaneParameters[i][3])) -
             pow(MethaneParameters[i][3], 2) * pow(density / criticalDensity, MethaneParameters[i][3]))
            * pow(criticalTemperature / Temperature, MethaneParameters[i][2])
            * exp(-pow(density / criticalDensity, MethaneParameters[i][3]));
        }

    } else if (componentName == 'carbonDioxide') {
        for (int i = 0; i < 6; ++i) {
            sum += carbonDioxideParameters[i][0] * carbonDioxideParameters[i][1] * (carbonDioxideParameters[i][1] - 1)
                   * pow(density / criticalDensity, carbonDioxideParameters[i][1] - 2)
                   * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2]);
        }
        for (int i = 6; i < 14; ++i) {
            carbonDioxideParameters[i][0]
            * pow(density / criticalDensity, carbonDioxideParameters[i][1] - 2)
            *
            ((carbonDioxideParameters[i][1] -
              carbonDioxideParameters[i][3] * pow(density / criticalDensity, carbonDioxideParameters[i][3])) *
             (carbonDioxideParameters[i][1] - 1 -
              carbonDioxideParameters[i][3] * pow(density / criticalDensity, carbonDioxideParameters[i][3])) -
             pow(carbonDioxideParameters[i][3], 2) * pow(density / criticalDensity, carbonDioxideParameters[i][3]))
            * pow(criticalTemperature / Temperature, carbonDioxideParameters[i][2])
            * exp(-pow(density / criticalDensity, carbonDioxideParameters[i][3]));
        }

    }
    return sum;

}

// 混合物的 无量度亥姆霍兹自由能的剩余部分 对对比密度的二阶导数
double mixtureReducedHelmholtzResidualPartSecondDerivativeToContrastDensity(double x1,
                                                                            char componentName1, char componentName2,
                                                                            double Temperature, double density,
                                                                            double criticalTemperature1,
                                                                            double criticalDensity1,
                                                                            double criticalTemperature2,
                                                                            double criticalDensity2) {
    return (x1 *
            iReducedHelmholtzResidualPartPartSecondDerivativeToContrastDensity(componentName1, Temperature, density,
                                                                               criticalTemperature1, criticalDensity1)
            + (1 - x1) *
              iReducedHelmholtzResidualPartPartSecondDerivativeToContrastDensity(componentName2, Temperature, density,
                                                                                 criticalTemperature2,
                                                                                 criticalDensity2)
            + mixtureReducedHelmholtzExcessSecondDerivativeToContrastDensity(x1, Temperature, density, criticalDensity1,
                                                                             criticalDensity2, criticalTemperature1,
                                                                             criticalTemperature2));

}


double primitiveFunction(double contrastDensity) {
    double p, T, componentName1, componentName2, x1, density, R,
            criticalDensity1, criticalDensity2, criticalTemperature1, criticalTemperature2;
    return (p - contrastDensity / mixtureVirtualCriticalDensityReverse(x1, criticalDensity1, criticalDensity2) * R * T
                * (1 + contrastDensity *
                       mixtureReducedHelmholtzResidualPartDerivativeToContrastDensity(x1, componentName1,
                                                                                      componentName2,
                                                                                      T, density,
                                                                                      criticalTemperature1,
                                                                                      criticalDensity1,
                                                                                      criticalTemperature2,
                                                                                      criticalDensity2)));
}

double derivativeFunction(double contrastDensity) {
    double T, componentName1, componentName2, x1, density, R,
            criticalDensity1, criticalDensity2, criticalTemperature1, criticalTemperature2;

    return (
            1 / mixtureVirtualCriticalDensityReverse(x1, criticalDensity1, criticalDensity2) * R * T
            * (pow(contrastDensity, 2) *
               mixtureReducedHelmholtzResidualPartSecondDerivativeToContrastDensity(x1, componentName1, componentName2,
                                                                                    T, density, criticalTemperature1,
                                                                                    criticalDensity1,
                                                                                    criticalTemperature2,
                                                                                    criticalDensity2) - 1)
    );
}

int Newton(double *contrastDensity, double precision, int maxCycle) {
    double nextValue, initialValue;
    int k;
    initialValue = *contrastDensity;
    for (k = 0; k < maxCycle; k++) {
        if (derivativeFunction(initialValue) == 0.0)//若通过初值，函数返回值为0
        {
            printf("The derivative is zero during the iteration!\n");
            return 0;
        }
        nextValue = initialValue - primitiveFunction(initialValue) / derivativeFunction(initialValue);//进行牛顿迭代计算
        if (fabs(nextValue - initialValue) < precision || fabs(primitiveFunction(nextValue)) < precision)//达到结束条件
        {
            *contrastDensity = nextValue; //返回结果
            printf("The roots around this value are:%lf\n", contrastDensity);

            return 1;
        } else //未达到结束条件
        {
            initialValue = nextValue; //准备下一次迭代
        }
    }
    printf("The number of iterations exceeds expectations!\n"); //迭代次数达到，仍没有达到精度
    return 0;
}

int main() {
    double contrastDensity, precision;
    int maxCycle;
    printf("Please enter the iteration initialValue:");
    scanf_s("%lf", &contrastDensity);
    printf("Please enter the maximum number of iterations:");
    scanf_s("%d", &maxCycle);
    printf("Please enter the required precision for iteration:");
    scanf_s("%lf", &precision);
    printf("%s", Newton(&contrastDensity, precision, maxCycle) == 1 ? &contrastDensity : "Sorry,Iteration failed!\n");
    return 0;
}