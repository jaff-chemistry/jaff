#ifndef ROSENBROCK_TABLEAU_H
#define ROSENBROCK_TABLEAU_H

namespace rosenbrock {

struct ros2s_tableau {
  static constexpr int stages = 3;
  static constexpr double gamma = 0.292893218813452;

  inline static constexpr double ct(const int i) {
    return i == 2 ? 0.585786437626905 : 1.0;
  }

  inline static constexpr double a(const int i, const int j) {
    if (i == 2 && j == 1) {
      return 2.0000000000000036;
    }
    if (i == 3 && j == 1) {
      return 6.828427124746214;
    }
    if (i == 3 && j == 2) {
      return 3.4142135623731007;
    }
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    if (i == 2 && j == 1) {
      return -6.828427124746214;
    }
    if (i == 3 && j == 1) {
      return -10.949747468305889;
    }
    if (i == 3 && j == 2) {
      return -7.535533905932761;
    }
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    if (i == 1) {
      return 6.828427124746214;
    }
    if (i == 2) {
      return 3.414213562373101;
    }
    return 1.0;
  }

  inline static constexpr double e(const int i) {
    if (i == 1) {
      return -0.23570226039551292;
    }
    if (i == 2) {
      return -0.23570226039551567;
    }
    return -0.13807118745769906;
  }
};

struct ros2_tableau {
  static constexpr int stages = 2;
  static constexpr double gamma = 1.7071067811865475;

  inline static constexpr double ct(const int i) {
    (void)i;
    return 1.0;
  }

  inline static constexpr double a(const int i, const int j) {
    if (i == 2 && j == 1) {
      return 0.585786437626905;
    }
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    if (i == 2 && j == 1) {
      return -1.1715728752538102;
    }
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    if (i == 1) {
      return 0.8786796564403575;
    }
    return 0.2928932188134525;
  }

  inline static constexpr double e(const int i) {
    (void)i;
    return 0.2928932188134525;
  }
};

struct rodas3p_tableau {
  static constexpr int stages = 5;
  static constexpr double gamma = 1.0 / 3.0;

  inline static constexpr double ct(const int i) {
    if (i == 2) {
      return 4.0 / 9.0;
    }
    if (i >= 4) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double a(const int i, const int j) {
    if (i == 2 && j == 1) {
      return 4.0 / 3.0;
    }
    if ((i == 4 || i == 5) && j == 1) {
      return 2.90625;
    }
    if ((i == 4 || i == 5) && j == 2) {
      return 3.375;
    }
    if ((i == 4 || i == 5) && j == 3) {
      return 0.40625;
    }
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    if (i == 2 && j == 1) {
      return -4.0;
    }
    if (i == 3 && j == 1) {
      return 8.25;
    }
    if (i == 3 && j == 2) {
      return 6.75;
    }
    if (i == 4 && j == 1) {
      return 1.21875;
    }
    if (i == 4 && j == 2) {
      return -5.0625;
    }
    if (i == 4 && j == 3) {
      return -1.96875;
    }
    if (i == 5 && j == 1) {
      return 4.03125;
    }
    if (i == 5 && j == 2) {
      return -15.1875;
    }
    if (i == 5 && j == 3) {
      return -4.03125;
    }
    if (i == 5 && j == 4) {
      return 6.0;
    }
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    if (i == 1) {
      return 2.90625;
    }
    if (i == 2) {
      return 3.375;
    }
    if (i == 3) {
      return 0.40625;
    }
    if (i == 5) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double e(const int i) {
    if (i == 4) {
      return -1.0;
    }
    if (i == 5) {
      return 1.0;
    }
    return 0.0;
  }
};

struct rodas5p_tableau {
  static constexpr int stages = 8;
  static constexpr double gamma = 0.21193756319429014;

  inline static constexpr double ct(const int i) {
    if (i == 2) {
      return 0.6358126895828704;
    }
    if (i == 3) {
      return 0.4095798393397535;
    }
    if (i == 4) {
      return 0.9769306725060716;
    }
    if (i == 5) {
      return 0.4288403609558664;
    }
    if (i >= 6) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double a(const int i, const int j) {
    if (i == 2 && j == 1) {
      return 3.0;
    }
    if (i == 3 && j == 1) {
      return 2.849394379747939;
    }
    if (i == 3 && j == 2) {
      return 0.45842242204463923;
    }
    if (i == 4 && j == 1) {
      return -6.954028509809101;
    }
    if (i == 4 && j == 2) {
      return 2.489845061869568;
    }
    if (i == 4 && j == 3) {
      return -10.358996098473584;
    }
    if (i == 5 && j == 1) {
      return 2.8029986275628964;
    }
    if (i == 5 && j == 2) {
      return 0.5072464736228206;
    }
    if (i == 5 && j == 3) {
      return -0.3988312541770524;
    }
    if (i == 5 && j == 4) {
      return -0.04721187230404641;
    }
    if ((i == 6 || i == 7 || i == 8) && j == 1) {
      return -7.502846399306121;
    }
    if ((i == 6 || i == 7 || i == 8) && j == 2) {
      return 2.561846144803919;
    }
    if ((i == 6 || i == 7 || i == 8) && j == 3) {
      return -11.627539656261098;
    }
    if ((i == 6 || i == 7 || i == 8) && j == 4) {
      return -0.18268767659942256;
    }
    if ((i == 6 || i == 7 || i == 8) && j == 5) {
      return 0.030198172008377946;
    }
    if ((i == 7 || i == 8) && j == 6) {
      return 1.0;
    }
    if (i == 8 && j == 7) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    if (i == 2 && j == 1) {
      return -14.155112264123755;
    }
    if (i == 3 && j == 1) {
      return -17.97296035885952;
    }
    if (i == 3 && j == 2) {
      return -2.859693295451294;
    }
    if (i == 4 && j == 1) {
      return 147.12150275711716;
    }
    if (i == 4 && j == 2) {
      return -1.41221402718213;
    }
    if (i == 4 && j == 3) {
      return 71.68940251302358;
    }
    if (i == 5 && j == 1) {
      return 165.43517024871676;
    }
    if (i == 5 && j == 2) {
      return -0.4592823456491126;
    }
    if (i == 5 && j == 3) {
      return 42.90938336958603;
    }
    if (i == 5 && j == 4) {
      return -5.961986721573306;
    }
    if (i == 6 && j == 1) {
      return 24.854864614690072;
    }
    if (i == 6 && j == 2) {
      return -3.0009227002832186;
    }
    if (i == 6 && j == 3) {
      return 47.4931110020768;
    }
    if (i == 6 && j == 4) {
      return 5.5814197821558125;
    }
    if (i == 6 && j == 5) {
      return -0.6610691825249471;
    }
    if (i == 7 && j == 1) {
      return 30.91273214028599;
    }
    if (i == 7 && j == 2) {
      return -3.1208243349937974;
    }
    if (i == 7 && j == 3) {
      return 77.79954646070892;
    }
    if (i == 7 && j == 4) {
      return 34.28646028294783;
    }
    if (i == 7 && j == 5) {
      return -19.097331116725623;
    }
    if (i == 7 && j == 6) {
      return -28.087943162872662;
    }
    if (i == 8 && j == 1) {
      return 37.80277123390563;
    }
    if (i == 8 && j == 2) {
      return -3.2571969029072276;
    }
    if (i == 8 && j == 3) {
      return 112.26918849496327;
    }
    if (i == 8 && j == 4) {
      return 66.9347231244047;
    }
    if (i == 8 && j == 5) {
      return -40.06618937091002;
    }
    if (i == 8 && j == 6) {
      return -54.66780262877968;
    }
    if (i == 8 && j == 7) {
      return -9.48861652309627;
    }
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    if (i == 1) {
      return -7.502846399306121;
    }
    if (i == 2) {
      return 2.561846144803919;
    }
    if (i == 3) {
      return -11.627539656261098;
    }
    if (i == 4) {
      return -0.18268767659942256;
    }
    if (i == 5) {
      return 0.030198172008377946;
    }
    return 1.0;
  }

  inline static constexpr double e(const int i) { return i == 8 ? 1.0 : 0.0; }
};

struct rodas4p_tableau {
  static constexpr int stages = 6;
  static constexpr double gamma = 0.25;

  inline static constexpr double ct(const int i) {
    if (i == 2) {
      return 0.75;
    }
    if (i == 3) {
      return 0.21;
    }
    if (i == 4) {
      return 0.63;
    }
    if (i >= 5) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double a(const int i, const int j) {
    if (i == 2 && j == 1) {
      return 3.0;
    }
    if (i == 3 && j == 1) {
      return 1.831036793486759;
    }
    if (i == 3 && j == 2) {
      return 0.4955183967433795;
    }
    if (i == 4 && j == 1) {
      return 2.304376582692669;
    }
    if (i == 4 && j == 2) {
      return -0.05249275245743001;
    }
    if (i == 4 && j == 3) {
      return -1.176798761832782;
    }
    if ((i == 5 || i == 6) && j == 1) {
      return -7.170454962423024;
    }
    if ((i == 5 || i == 6) && j == 2) {
      return -4.741636671481785;
    }
    if ((i == 5 || i == 6) && j == 3) {
      return -16.31002631330971;
    }
    if ((i == 5 || i == 6) && j == 4) {
      return -1.062004044111401;
    }
    if (i == 6 && j == 5) {
      return 1.0;
    }
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    if (i == 2 && j == 1) {
      return -12.0;
    }
    if (i == 3 && j == 1) {
      return -8.791795173947035;
    }
    if (i == 3 && j == 2) {
      return -2.207865586973518;
    }
    if (i == 4 && j == 1) {
      return 10.81793056857153;
    }
    if (i == 4 && j == 2) {
      return 6.780270611428266;
    }
    if (i == 4 && j == 3) {
      return 19.5348594464241;
    }
    if (i == 5 && j == 1) {
      return 34.19095006749676;
    }
    if (i == 5 && j == 2) {
      return 15.49671153725963;
    }
    if (i == 5 && j == 3) {
      return 54.7476087596413;
    }
    if (i == 5 && j == 4) {
      return 14.16005392148534;
    }
    if (i == 6 && j == 1) {
      return 34.62605830930532;
    }
    if (i == 6 && j == 2) {
      return 15.30084976114473;
    }
    if (i == 6 && j == 3) {
      return 56.99955578662667;
    }
    if (i == 6 && j == 4) {
      return 18.40807009793095;
    }
    if (i == 6 && j == 5) {
      return -5.714285714285717;
    }
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    if (i == 1) {
      return -7.170454962423024;
    }
    if (i == 2) {
      return -4.741636671481785;
    }
    if (i == 3) {
      return -16.31002631330971;
    }
    if (i == 4) {
      return -1.062004044111401;
    }
    return 1.0;
  }

  inline static constexpr double e(const int i) { return i == 6 ? 1.0 : 0.0; }
};

struct rosenbrock_euler_tableau {
  static constexpr int stages = 1;
  static constexpr double gamma = 1.0;

  inline static constexpr double ct(const int i) {
    (void)i;
    return 0.0;
  }

  inline static constexpr double a(const int i, const int j) {
    (void)i;
    (void)j;
    return 0.0;
  }

  inline static constexpr double c(const int i, const int j) {
    (void)i;
    (void)j;
    return 0.0;
  }

  inline static constexpr double b(const int i) {
    (void)i;
    return 1.0;
  }

  inline static constexpr double e(const int i) {
    (void)i;
    return 1.0;
  }
};

using coefficients = rodas5p_tableau;

} // namespace rosenbrock

#endif
