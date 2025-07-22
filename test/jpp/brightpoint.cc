#include <iostream>
#include "JPhysics/JPDF.hh"
#include "JPhysics/KM3NeT.hh"

double  absorptionLengthFactor = 1.0;
double  scatteringLengthFactor = 1.0;


inline double getAbsorptionLength(const double lambda)
{
  return absorptionLengthFactor * NAMESPACE::getAbsorptionLength(lambda);
}


inline double getScatteringLength(const double lambda)
{
  return scatteringLengthFactor * NAMESPACE::getScatteringLength(lambda);
}



int main(int argc, char **argv)
{
  using namespace std;
  using namespace JPP;

  double epsilon = 1.0e-5;
  int numberOfPoints = 5;
  const double P_atm = 240.0; //NAMESPACE::getAmbientPressure();
  const double wmin  = getMinimalWavelength();
  const double wmax  = getMaximalWavelength();


  const JPDF_C
    pdf_c(NAMESPACE::getPhotocathodeArea(),
      NAMESPACE::getQE,
      NAMESPACE::getAngularAcceptance,
      getAbsorptionLength,
      getScatteringLength,
      NAMESPACE::getScatteringProbability,
      P_atm,
      wmin,
      wmax,
      numberOfPoints,
      epsilon);

  double npe;

  std::cout << "R_m,cos_theta,delta_t,dlbright,slbright" << std::endl;
  for (double R_m = 0.4; R_m < 100; R_m *= 2) {
    for (double cos_theta = -1; cos_theta < 1; cos_theta += 0.2) {
      for (double delta_t = -20; delta_t < 50; delta_t += 2) {
        std::cout << R_m << "," << cos_theta << "," << delta_t << ",";
        npe = pdf_c.getDirectLightFromBrightPoint(R_m, cos_theta, delta_t);
        std::cout << npe << ",";
        npe = pdf_c.getScatteredLightFromBrightPoint(R_m, cos_theta, delta_t);
        std::cout << npe;
        std::cout << std::endl;
      }
    }
  }

  return 0;
}
