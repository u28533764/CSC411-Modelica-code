within ;
package DistillationWithConnectors
  "A Distillation column using the Connector class"
  connector Vapour "A vapour connector"
    flow Modelica.SIunits.MassFlowRate fv;
    Modelica.SIunits.Temp_C Tv;
    Modelica.SIunits.MassFraction y;
  end Vapour;

  connector Liquid "A liquid connector"
    flow Modelica.SIunits.MassFlowRate fl;
    Modelica.SIunits.Temp_C Tl;
    Modelica.SIunits.MassFraction x;
  end Liquid;

  model Plate "Model of a Single plate"
    parameter Modelica.SIunits.Pressure P = 101.325 "Pressure in the column";
    parameter Real k = 0.02 "Constant for plate hydraulics";
    parameter Modelica.SIunits.Temp_C Tref = 25 "Reference Tempearture";
    parameter Modelica.SIunits.MassFraction z=0.3 "Feed composition";
    parameter Modelica.SIunits.Temp_C T_F=40 "Temperature of the feed";
    parameter Modelica.SIunits.SpecificEnergy Lam_E = 924170.76786;
    parameter Modelica.SIunits.SpecificEnergy Lam_W = 2437190.7246;
    Modelica.SIunits.MassFlowRate F;
    Modelica.SIunits.SpecificEnthalpy h_0;
    Modelica.SIunits.SpecificEnthalpy h_i;
    Modelica.SIunits.SpecificEnthalpy H_0;
    Modelica.SIunits.SpecificEnthalpy H_i;
    input Vapour Vi "Incoming Vapour flow";
    output Vapour V0 "outgoing Vapour flow";
    input Liquid Li "Incoming Liquid flow";
    output Liquid L0 "Outgoing Liquid flow";
    Modelica.SIunits.Mass M "Mass accumulated in the plate";
    Modelica.SIunits.Energy E "Energy Accumulation within the plate";
    Modelica.SIunits.SpecificHeatCapacity Cp_mix
      "Specific Heat Capacity of Mixture";
    Modelica.SIunits.SpecificEnthalpy h_F
      "Specific Enthalpy of the feed mixture";
  initial equation
    der(M) = 0;
    der(E) = 0;
    der(L0.x)=0;
  equation
    der(M) = F + Li.fl + L0.fl;
    der(M*L0.x) = F*z + Li.fl*Li.x + L0.fl*L0.x + Vi.fv*Vi.y + V0.fv*V0.y;
    L0.fl = k*M;
    L0.x*Psat(L0.Tl)= P*V0.y;
    Vi.fv = V0.fv;
    der(E) = F*h_F + Li.fl*h_i + Vi.fv*H_i + V0.fv*H_0 + L0.fl*h_0;
    E = M*Cp_mix*L0.Tl;
    Cp_mix = L0.x*Cp_E(L0.Tl) + (1-L0.x)*Cp_W(L0.Tl);
    h_F = z*Cp_E(T_F)*(T_F - Tref) + (1-z)*Cp_W(T_F)*(T_F - Tref);
    H_0 = V0.y*(Cp_E(V0.Tv)*(V0.Tv-Tref) + Lam_E) + (1-V0.y)*(Cp_W(V0.Tv)*(V0.Tv-Tref) + Lam_W);
    H_i = Vi.y*(Cp_E(Vi.Tv)*(Vi.Tv-Tref) + Lam_E) + (1-Vi.y)*(Cp_W(Vi.Tv)*(Vi.Tv-Tref) + Lam_W);
    h_i = Li.x*Cp_E(Li.Tl)*(Li.Tl-Tref) + (1-Li.x)*Cp_W(Li.Tl)*(Li.Tl-Tref);
    h_0 = L0.x*Cp_E(L0.Tl)*(L0.Tl-Tref) + (1-L0.x)*Cp_W(L0.Tl)*(L0.Tl-Tref);
    L0.Tl = V0.Tv;
  end Plate;

  function Psat "Vapour Pressure of Ethanol"
  protected
    parameter Real A = 16.8958;
    parameter Real B = 3795.17;
    parameter Real C = 230.918;
  public
    output Modelica.SIunits.Pressure P_E;
    input Modelica.SIunits.Temp_C T;
  algorithm
    P_E :=exp(A - (B/(T + C)));
  end Psat;

  function Cp_E "Specific Heat capacity of Ethanol"
  protected
    parameter Real C1 = 1.0264E05;
    parameter Real C2 = -1.3963E02;
    parameter Real C3 = -3.0341E-02;
    parameter Real C4 = 2.0386E-03;
  public
    input Modelica.SIunits.Temp_C T;
    output Modelica.SIunits.SpecificHeatCapacity C;
  algorithm
    C := (C1 + C2*(T+273.15) + C3*((T+273.15)^2) + C4*((T+273.15)^3))/46;
  end Cp_E;

  function Cp_W "Specific Heat Capacity of Water"
  protected
    parameter Real C1 = 2.7637E05;
    parameter Real C2 = -2.0961E03;
    parameter Real C3 = 8.125;
    parameter Real C4 = -1.4116E-02;
    parameter Real C5 = 9.3701E-06;
  public
    output Modelica.SIunits.SpecificHeatCapacity C;
    input Modelica.SIunits.Temp_C T;
  algorithm
    C :=(C1 + C2*(T + 273.15) + C3*((T + 273.15)^2) + C4*((T + 273.15)^3) + C5*
      ((T + 273.15)^4))/18;
  end Cp_W;

  model fullmodel
    Plate plate[2];
    Condenser condenser(x_D=0.8, T_D=78.2);
    Reboiler reboiler(x_B=0.05);
  equation
    connect(plate[1].L0, plate[2].Li);
    connect(plate[1].Vi, plate[2].V0);
    connect(plate[1].Li, condenser.R);
    connect(plate[1].V0, condenser.V);
    connect(plate[2].Vi, reboiler.V);
    connect(plate[2].L0, reboiler.L);

    //reboiler.x_B = 0.05;
    //reboiler.T_B = 96;
    reboiler.B = (plate[1].F+plate[2].F)*((condenser.x_D-plate[1].z)/(condenser.x_D-reboiler.x_B));

    //condenser.x_D = 0.8;
    //condenser.T_D = 78.2;
    condenser.D = (plate[1].F+plate[2].F)*((plate[1].z-reboiler.x_B)/(condenser.x_D-reboiler.x_B));

    plate[1].F = 2.8;
    plate[2].F = 0;
  end fullmodel;

  model Reboiler
    Modelica.SIunits.MassFlowRate B;
    output Vapour V;
    input Liquid L;
    input Modelica.SIunits.MassFraction x_B;
    Modelica.SIunits.Mass M_B;
    Modelica.SIunits.EnergyFlowRate Q_R;
    Modelica.SIunits.Energy E;
    Modelica.SIunits.SpecificEnthalpy H_10;
    Modelica.SIunits.SpecificEnthalpy h_10;
    Modelica.SIunits.SpecificEnthalpy h_B;
    Modelica.SIunits.SpecificHeatCapacity Cp_mix;
    Modelica.SIunits.Temp_C T_B;
    parameter Real k = 0.2;
    parameter Modelica.SIunits.Pressure P = 101.325;
    parameter Modelica.SIunits.Temp_C Tref = 25;
    parameter Modelica.SIunits.SpecificEnergy Lam_E = 924170.76786;
    parameter Modelica.SIunits.SpecificEnergy Lam_W = 2437190.7246;
  initial equation
    der(M_B) = 0;
    der(x_B) = 0;
    der(E) = 0;
  equation
    der(M_B) = L.fl + B + V.fv;
    der(M_B*x_B) = L.fl*L.x + B*x_B + V.fv*V.y;
    der(E) = L.fl*h_10 + B*h_B + V.fv*H_10 + Q_R;
    x_B*Psat(T_B) = P*V.y;
    E = M_B*Cp_mix*(T_B-Tref);
    Cp_mix = x_B*Cp_E(T_B) + (1-x_B)*Cp_W(T_B);
    h_B = x_B*Cp_E(T_B)*(T_B-Tref) + (1-x_B)*Cp_W(T_B)*(T_B-Tref);
    h_10 = L.x*Cp_E(L.Tl)*(L.Tl-Tref) + (1-L.x)*Cp_W(L.Tl)*(L.Tl-Tref);
    H_10 = V.y*(Cp_E(V.Tv)*(V.Tv-Tref)+Lam_E) + (1-V.y)*(Cp_W(V.Tv)*(V.Tv-Tref)+Lam_W);
    V.Tv = T_B;
    //V.fv = B*2;
    B = k*sqrt(M_B);
  end Reboiler;

  model Condenser
    Modelica.SIunits.MassFlowRate D;
    input Vapour V;
    output Liquid R;
    input Modelica.SIunits.MassFraction x_D;
    Modelica.SIunits.Mass M_D;
    Modelica.SIunits.EnergyFlowRate Q_C;
    Modelica.SIunits.Energy E;
    Modelica.SIunits.SpecificEnthalpy H_1;
    Modelica.SIunits.SpecificEnthalpy h_1;
    Modelica.SIunits.SpecificEnthalpy h_D;
    Modelica.SIunits.SpecificHeatCapacity Cp_mix;
    input Modelica.SIunits.Temp_C T_D;
    parameter Real k = 0.2;
    parameter Modelica.SIunits.Pressure P = 101.325;
    parameter Modelica.SIunits.Temp_C Tref = 25;
    parameter Modelica.SIunits.SpecificEnergy Lam_E = 924170.76786;
    parameter Modelica.SIunits.SpecificEnergy Lam_W = 2437190.7246;
  initial equation
    der(M_D) = 0;
    der(x_D) = 0;
    der(E) = 0;
  equation
    der(M_D) = R.fl + D + V.fv;
    der(M_D*x_D) = R.fl*R.x + D*x_D + V.fv*V.y;
    der(E) = R.fl*h_1 + D*h_D + V.fv*H_1 + Q_C;
    E = M_D*Cp_mix*(T_D-Tref);
    Cp_mix = x_D*Cp_E(T_D) + (1-x_D)*Cp_W(T_D);
    h_D = x_D*Cp_E(T_D)*(T_D-Tref) + (1-x_D)*Cp_W(T_D)*(T_D-Tref);
    h_1 = R.x*Cp_E(R.Tl)*(R.Tl-Tref) + (1-R.x)*Cp_W(R.Tl)*(R.Tl-Tref);
    H_1 = V.y*(Cp_E(V.Tv)*(V.Tv-Tref)+Lam_E) + (1-V.y)*(Cp_W(V.Tv)*(V.Tv-Tref)+Lam_W);
    R.Tl = T_D;
    //R.x = x_D;
    //R.fl = 3*D;
    D = k*sqrt(M_D);
  end Condenser;
  annotation (uses(Modelica(version="3.2")));
end DistillationWithConnectors;
