#include <iostream>
#include "JPhysics/JPDF.hh"
#include "JPhysics/KM3NeT.hh"
#include "JPhysics/JPDFToolkit.hh"

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
  const double P_atm = 240.0;
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

  // Coarser grid: energy-integrated shower is expensive
  std::cout << "E,D,cd,theta,phi,delta_t,npe" << std::endl;
  for (double E = 1; E < 1000; E *= 10) {
    for (double D = 0.4; D < 100; D *= 2) {
      for (double cd = -0.8; cd <= 0.8; cd += 0.4) {
        for (double theta = 0; theta < 2*PI; theta += 1.57) {
          for (double phi = 0; phi < 2*PI; phi += 1.57) {
            for (double delta_t = -40; delta_t < 200; delta_t += 30) {
              std::cout << E << "," << D << "," << cd << "," << theta << "," << phi << "," << delta_t << ",";
              npe = pdf_c.getDirectLightFromEMshower(E, D, cd, theta, phi, delta_t);
              std::cout << npe;
              std::cout << std::endl;
            }
          }
        }
      }
    }
  }

  return 0;
}
