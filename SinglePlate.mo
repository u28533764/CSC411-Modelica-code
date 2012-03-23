within ;
model SinglePlate "Mass and Component balance over a single plate"
  parameter Real V = 1;
  Real V2;
  parameter Real L2 = 1.2;
  Real L;
  Real x;
  parameter Real x2 = 0.07;
  parameter Real y2 = 0.7;
  Real y;
  parameter Real F = 2.8;
  parameter Real z = 0.3;
  Real Psat;
  Real M(start=0);
  parameter Real P = 101.325;
  parameter Real T = 178/2;
  parameter Real A = 16.8958;
  parameter Real B = 3795.17;
  parameter Real C = 230.918;
  parameter Real k = 0.02;
equation
  der(M) = F + L2 - L + V - V2;
  der(M*x) = F*z + L2*x2 - L*x + V*y - V2*y2;
  V = V2;
  L = k*M;
  Psat = exp(A - (B/(T + C)));
  x = P*y/Psat;
  annotation (uses(Modelica(version="3.2")));
end SinglePlate;
