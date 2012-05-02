within ;
package DistillationTenPlates
  "The Model of a Full Distillation Column Without connectors"
  function Psat
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

  model Plate
    parameter Modelica.SIunits.Pressure P = 101.325
      "Pressure in the distillation column";
    parameter Real k = 0.02 "Constant for plate hydrilics";
    parameter Modelica.SIunits.Temp_C Tref = 25 "Reference Tempearture";
    parameter Modelica.SIunits.MassFraction z=0.3 "Feed composition";
    Modelica.SIunits.MassFlowRate F "Feed flow rate";
    Modelica.SIunits.MassFraction x_i;
    Modelica.SIunits.MassFraction x_0(start=0.46);
    Modelica.SIunits.MassFraction y_0(start=0.226);
    Modelica.SIunits.MassFraction y_i;
    Modelica.SIunits.SpecificEnthalpy h_0(start=126840);
    Modelica.SIunits.SpecificEnthalpy h_i;
    Modelica.SIunits.SpecificEnthalpy H_0(start=137600);
    Modelica.SIunits.SpecificEnthalpy H_i;
    Modelica.SIunits.Temp_C t_0(start=61.35);
    Modelica.SIunits.Temp_C t_i;
    Modelica.SIunits.Temp_C T_0(start=61.35);
    Modelica.SIunits.Temp_C T_i;
    Modelica.SIunits.MassFlowRate V_i "Incoming Vapour flow";
    Modelica.SIunits.MassFlowRate V_0(start=1) "outgoing Vapour flow";
    Modelica.SIunits.MassFlowRate L_i "Incoming Liquid flow";
    Modelica.SIunits.MassFlowRate L_0(start=4) "Outgoing Liquid flow";
    Modelica.SIunits.Mass M(start=200) "Mass accumulated in the plate";
    Modelica.SIunits.Energy E(start=42817000)
      "Energy Accumulation within the plate";
    Modelica.SIunits.SpecificHeatCapacity Cp_mix
      "Specific Heat Capacity of the Mixture";
    Modelica.SIunits.SpecificEnthalpy h_F
      "Specific Enthalpy of the feed mixture";
    parameter Modelica.SIunits.Temp_C T_F=40 "Temperature of the feed";
    parameter Modelica.SIunits.SpecificEnergy Lam_E = 924170.76786;
    parameter Modelica.SIunits.SpecificEnergy Lam_W = 2437190.7246;
  initial equation
    der(M) = 0;
    der(E) = 0;
  equation
    der(M) = F + L_i - L_0;
    der(M*x_0) = F*z + L_i*x_i - L_0*x_0 + V_i*y_i - V_0*y_0;
    L_0 = k*M;
    x_0 = P*y_0/Psat(t_0);
    V_i = V_0;
    der(E) = F*h_F + L_i*h_i + V_i*H_i - V_0*H_0 - L_0*h_0;
    E = M*Cp_mix*(t_0-Tref);
    Cp_mix = x_0*Cp_E(t_0) + (1-x_0)*Cp_W(t_0);
    h_F = z*Cp_E(T_F)*(T_F - Tref) + (1-z)*Cp_W(T_F)*(T_F - Tref);
    H_0 = y_0*(Cp_E(T_0)*(T_0-Tref) + Lam_E) + (1-y_0)*(Cp_W(T_0)*(T_0-Tref) + Lam_W);
    H_i = y_i*(Cp_E(T_i)*(T_i-Tref) + Lam_E) + (1-y_i)*(Cp_W(T_i)*(T_i-Tref) + Lam_W);
    h_i = x_i*Cp_E(t_i)*(t_i-Tref) + (1-x_i)*Cp_W(t_i)*(t_i-Tref);
    h_0 = x_0*Cp_E(t_0)*(t_0-Tref) + (1-x_0)*Cp_W(t_0)*(t_0-Tref);
    t_0 = T_0;
    annotation (Diagram(graphics), Icon(graphics));
  end Plate;

  model CombineModels
    Plate plate[10];
    Reboiler reboiler;
    Condenser condenser;
  equation
    plate[10].F = 0;
    reboiler.x_B = 0.05;
    condenser.x_D = 0.8;
    condenser.T_D = 78.2;
    reboiler.T_B = 96;
    reboiler.B = (plate[5].F+plate[8].F)*((condenser.x_D-plate[1].z)/(condenser.x_D-reboiler.x_B));
    condenser.D = (plate[5].F+plate[8].F)*((plate[1].z-reboiler.x_B)/(condenser.x_D-reboiler.x_B));
    for n in 1:9 loop
      plate[n+1].L_i-plate[n].L_0 = 0;
      plate[n].V_i-plate[n+1].V_0 = 0;
      plate[n].y_i-plate[n+1].y_0 = 0;
      plate[n+1].x_i-plate[n].x_0 = 0;
      plate[n+1].t_i-plate[n].t_0 = 0;
      plate[n].T_i-plate[n+1].T_0 = 0;
      plate[n].F = (if n==5 then 2.8 elseif n==8 then 2.8 else 0);
      if n==1 then
        //plate[n].V_0 - condenser.V = 0;
        //plate[n].y_0 - condenser.y = 0;
        //plate[n].T_0 - condenser.T_V = 0;
        plate[n].L_i - condenser.R = 0;
        plate[n].x_i - condenser.x = 0;
        plate[n].t_i - condenser.T_R = 0;
      elseif n==9 then
        plate[n+1].V_i - reboiler.V = 0;
        plate[n+1].y_i - reboiler.y = 0;
        plate[n+1].T_i - reboiler.T_B = 0;
        //plate[n+1].L_0 - reboiler.L = 0;
        //plate[n+1].x_0 - reboiler.x = 0;
        plate[n+1].t_0 - reboiler.T_L = 0;
      end if;
    end for;
                    annotation (Placement(transformation(extent={{-62,46},{-42,66}})),
                Diagram(graphics),
      experiment(StopTime=100),
      __Dymola_experimentSetupOutput);
  end CombineModels;

  function Cp_E "The Specific Heat Capacity of Ethanol"
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

  function Cp_W "The Specific Heat Capacity of Water"
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

  model Reboiler "Partial Reboiler"
    Modelica.SIunits.MassFlowRate B;
    Modelica.SIunits.MassFlowRate V;
    Modelica.SIunits.MassFlowRate L;
    Modelica.SIunits.Temp_C T_L;
    Modelica.SIunits.Temp_C T_V;
    Modelica.SIunits.MassFraction y;
    Modelica.SIunits.MassFraction x;
    Modelica.SIunits.MassFraction x_B;
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
  equation
    der(M_B) = L - B - V;
    der(M_B*x_B) = L*x - B*x_B - V*y;
    der(E) = L*h_10 - B*h_B - V*H_10 + Q_R;
    x_B = P*y/Psat(T_B);
    E = M_B*Cp_mix*(T_B-Tref);
    Cp_mix = x_B*Cp_E(T_B) + (1-x_B)*Cp_W(T_B);
    h_B = x_B*Cp_E(T_B)*(T_B-Tref) + (1-x_B)*Cp_W(T_B)*(T_B-Tref);
    h_10 = x*Cp_E(T_L)*(T_L-Tref) + (1-x)*Cp_W(T_L)*(T_L-Tref);
    H_10 = y*(Cp_E(T_V)*(T_V-Tref)+Lam_E) + (1-y)*(Cp_W(T_V)*(T_V-Tref)+Lam_W);
    T_V = T_B;
    B = k*sqrt(M_B);
    V = B*2;
  end Reboiler;

  model Condenser "Total Condenser"
    Modelica.SIunits.MassFlowRate D;
    Modelica.SIunits.MassFlowRate V;
    Modelica.SIunits.MassFlowRate R;
    Modelica.SIunits.Temp_C T_V;
    Modelica.SIunits.Temp_C T_R;
    Modelica.SIunits.MassFraction y;
    Modelica.SIunits.MassFraction x;
    Modelica.SIunits.MassFraction x_D;
    Modelica.SIunits.Mass M_D;
    Modelica.SIunits.EnergyFlowRate Q_C;
    Modelica.SIunits.Energy E;
    Modelica.SIunits.SpecificEnthalpy H_1;
    Modelica.SIunits.SpecificEnthalpy h_1;
    Modelica.SIunits.SpecificEnthalpy h_D;
    Modelica.SIunits.SpecificHeatCapacity Cp_mix;
    Modelica.SIunits.Temp_C T_D;
    parameter Real k = 0.2;
    parameter Modelica.SIunits.Pressure P = 101.325;
    parameter Modelica.SIunits.Temp_C Tref = 25;
    parameter Modelica.SIunits.SpecificEnergy Lam_E = 924170.76786;
    parameter Modelica.SIunits.SpecificEnergy Lam_W = 2437190.7246;
  equation
    der(M_D) = R + D + V;
    der(M_D*x_D) = R*x + D*x_D + V*y;
    der(E) = R*h_1 + D*h_D + V*H_1 + Q_C;
    E = M_D*Cp_mix*(T_D-Tref);
    Cp_mix = x_D*Cp_E(T_D) + (1-x_D)*Cp_W(T_D);
    h_D = x_D*Cp_E(T_D)*(T_D-Tref) + (1-x_D)*Cp_W(T_D)*(T_D-Tref);
    h_1 = x*Cp_E(T_R)*(T_R-Tref) + (1-x)*Cp_W(T_R)*(T_R-Tref);
    H_1 = y*(Cp_E(T_V)*(T_V-Tref)+Lam_E) + (1-y)*(Cp_W(T_V)*(T_V-Tref)+Lam_W);
    //Q_C = V*H_1 + D*h_D + R*h_1;
    T_R = T_D;
    x = x_D;
    T_V = T_D;
    D = k*sqrt(M_D);
    R = 3*D;
  end Condenser;
  annotation (Diagram(graphics={
        Polygon(
          points={{-66,74},{-70,68},{-62,68},{-66,74}},
          lineColor={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{44,48},{34,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{40,48},{30,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{34,48},{26,40}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{72,40},{72,32}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-56,60},{-16,-54}},
          lineColor={0,0,255},
          radius=10),
        Line(
          points={{-56,40},{-22,40},{-22,44},{-22,36}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-16,30},{-50,30},{-50,34},{-50,26}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,20},{-22,20},{-22,24},{-22,16}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-16,10},{-50,10},{-50,14},{-50,6}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,0},{-22,0},{-22,4},{-22,-4}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-16,-10},{-50,-10},{-50,-6},{-50,-14}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,-20},{-22,-20},{-22,-16},{-22,-24}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-16,-30},{-50,-30},{-50,-26},{-50,-34}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,-40},{-22,-40},{-22,-36},{-22,-44}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-16,50},{-50,50},{-50,54},{-50,46}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,-48},{-16,-48}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,22},{-96,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,-18},{-96,-18}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-56,-18},{-60,-16},{-60,-20},{-56,-18}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Rectangle(
          extent={{-80,-44},{-66,-68}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          radius=10),
        Line(
          points={{-74,-44},{-74,-38},{-56,-38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-72,-80},{-72,-66}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-72,-80},{-38,-80},{-38,-54}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-38,-80},{20,-80}},
          color={0,0,255},
          smooth=Smooth.None),
        Ellipse(extent={{8,-44},{-2,-54}}, lineColor={0,0,255}),
        Line(
          points={{-16,-48},{-2,-48}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-4,-76},{-4,-84},{2,-80},{-4,-76}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{8,-76},{8,-84},{2,-80},{8,-76}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{24,-80},{20,-78},{20,-82},{24,-80}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-56,22},{-60,24},{-60,20},{-56,22}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-56,-38},{-60,-36},{-60,-40},{-56,-38}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-92,-50},{-80,-50}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-80,-62},{-88,-62}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-80,-50},{-84,-48},{-84,-52},{-80,-50}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-92,-62},{-88,-60},{-88,-64},{-92,-62}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Text(
          extent={{4,-48},{4,-48}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          textStyle={TextStyle.Bold},
          textString="%LC",
          fontSize=14),
        Rectangle(
          extent={{-2,-66},{6,-72}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{2,-54},{2,-66}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{2,-72},{2,-80}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-36,60},{-36,80},{-8,80}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-8,86},{24,72}},
          lineColor={0,0,255},
          fillColor={215,215,215},
          fillPattern=FillPattern.Sphere,
          radius=10),
        Line(
          points={{24,78},{36,78},{36,64}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{24,64},{54,38}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255},
          radius=5),
        Line(
          points={{38,38},{38,20},{88,20}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{38,20},{10,20},{10,54},{-16,54}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-18,54},{-14,56},{-14,52},{-18,54}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-8,80},{-12,82},{-12,78},{-8,80}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{90,20},{86,22},{86,18},{90,20}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{24,48},{54,48}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{0,86},{0,94}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{16,86},{16,96}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-2,0},{2,2},{2,-2},{-2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={0,94},
          rotation=-90),
        Polygon(
          points={{2,0},{-2,2},{-2,-2},{2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={16,88},
          rotation=270),
        Ellipse(extent={{78,50},{68,40}}, lineColor={0,0,255}),
        Line(
          points={{54,46},{68,46}},
          color={0,0,255},
          smooth=Smooth.None),
        Text(
          extent={{74,46},{74,46}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          textStyle={TextStyle.Bold},
          textString="%LC",
          fontSize=14),
        Rectangle(
          extent={{68,36},{76,30}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{72,30},{72,20}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{78,24},{78,16},{72,20},{78,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{66,24},{66,16},{72,20},{66,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Text(
          extent={{36,44},{36,44}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="M_D",
          fontSize=14,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-38,-52},{-38,-52}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="M_B",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{24,-76},{38,-84}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="B
x_B",     fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{88,22},{102,14}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="D
x_D"),
        Text(
          extent={{-74,28},{-74,28}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="F_1"),
        Text(
          extent={{-74,-12},{-74,-12}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="F_2"),
        Text(
          extent={{-74,18},{-74,16}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="z_1",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-74,-22},{-74,-24}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="z_2"),
        Text(
          extent={{-4,58},{-4,58}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="R",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-14,84},{-28,84}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="V_1",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Polygon(
          points={{2,0},{-2,2},{-2,-2},{2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={36,66},
          rotation=270),
        Line(
          points={{54,46},{46,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{52,48},{42,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{48,48},{38,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{30,48},{24,42}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{36,22},{40,18}},
          lineColor={0,0,255},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-40,-78},{-36,-82}},
          lineColor={0,0,255},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{30,24},{30,16},{24,20},{30,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{18,24},{18,16},{24,20},{18,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Rectangle(
          extent={{20,34},{28,28}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{24,20},{24,28}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-94,26},{-94,18},{-88,22},{-94,26}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-82,26},{-82,18},{-88,22},{-82,26}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-88,30},{-88,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-92,36},{-84,30}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-92,-14},{-92,-22},{-86,-18},{-92,-14}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-80,-14},{-80,-22},{-86,-18},{-80,-14}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-86,-10},{-86,-18}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-90,-4},{-82,-10}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-96,-18},{-96,78},{-68,78}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-82,78},{-78,80},{-78,76},{-82,78}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Ellipse(
          extent={{-62,78},{-70,70}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-68,74},{-64,76},{-64,72},{-68,74}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-64,74},{-50,74}},
          color={0,0,255},
          smooth=Smooth.None)}),                                  uses(Modelica(
          version="3.2")),
    Icon(graphics={
        Polygon(
          points={{-64,74},{-68,68},{-60,68},{-64,74}},
          lineColor={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{46,48},{36,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{42,48},{32,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{36,48},{28,40}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{74,40},{74,32}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-54,60},{-14,-54}},
          lineColor={0,0,255},
          radius=10),
        Line(
          points={{-54,40},{-20,40},{-20,44},{-20,36}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-14,30},{-48,30},{-48,34},{-48,26}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,20},{-20,20},{-20,24},{-20,16}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-14,10},{-48,10},{-48,14},{-48,6}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,0},{-20,0},{-20,4},{-20,-4}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-14,-10},{-48,-10},{-48,-6},{-48,-14}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,-20},{-20,-20},{-20,-16},{-20,-24}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-14,-30},{-48,-30},{-48,-26},{-48,-34}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,-40},{-20,-40},{-20,-36},{-20,-44}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-14,50},{-48,50},{-48,54},{-48,46}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,-48},{-14,-48}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,22},{-94,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-54,-18},{-94,-18}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-54,-18},{-58,-16},{-58,-20},{-54,-18}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Rectangle(
          extent={{-78,-44},{-64,-68}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          radius=10),
        Line(
          points={{-72,-44},{-72,-38},{-54,-38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-70,-80},{-70,-66}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-70,-80},{-36,-80},{-36,-54}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-36,-80},{22,-80}},
          color={0,0,255},
          smooth=Smooth.None),
        Ellipse(extent={{10,-44},{0,-54}}, lineColor={0,0,255}),
        Line(
          points={{-14,-48},{0,-48}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-2,-76},{-2,-84},{4,-80},{-2,-76}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{10,-76},{10,-84},{4,-80},{10,-76}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{26,-80},{22,-78},{22,-82},{26,-80}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-54,22},{-58,24},{-58,20},{-54,22}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-54,-38},{-58,-36},{-58,-40},{-54,-38}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-90,-50},{-78,-50}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-78,-62},{-86,-62}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-78,-50},{-82,-48},{-82,-52},{-78,-50}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-90,-62},{-86,-60},{-86,-64},{-90,-62}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Text(
          extent={{6,-48},{6,-48}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          textStyle={TextStyle.Bold},
          textString="%LC",
          fontSize=14),
        Rectangle(
          extent={{0,-66},{8,-72}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{4,-54},{4,-66}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{4,-72},{4,-80}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-34,60},{-34,80},{-6,80}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-6,86},{26,72}},
          lineColor={0,0,255},
          fillColor={215,215,215},
          fillPattern=FillPattern.Sphere,
          radius=10),
        Line(
          points={{26,78},{38,78},{38,64}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{26,64},{56,38}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255},
          radius=5),
        Line(
          points={{40,38},{40,20},{90,20}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{40,20},{12,20},{12,54},{-14,54}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-16,54},{-12,56},{-12,52},{-16,54}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-6,80},{-10,82},{-10,78},{-6,80}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{92,20},{88,22},{88,18},{92,20}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{26,48},{56,48}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{2,86},{2,94}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{18,86},{18,96}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-2,0},{2,2},{2,-2},{-2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={2,94},
          rotation=-90),
        Polygon(
          points={{2,0},{-2,2},{-2,-2},{2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={18,88},
          rotation=270),
        Ellipse(extent={{80,50},{70,40}}, lineColor={0,0,255}),
        Line(
          points={{56,46},{70,46}},
          color={0,0,255},
          smooth=Smooth.None),
        Text(
          extent={{76,46},{76,46}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          textStyle={TextStyle.Bold},
          textString="%LC",
          fontSize=14),
        Rectangle(
          extent={{70,36},{78,30}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{74,30},{74,20}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{80,24},{80,16},{74,20},{80,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{68,24},{68,16},{74,20},{68,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Text(
          extent={{38,44},{38,44}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="M_D",
          fontSize=14,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-36,-52},{-36,-52}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="M_B",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{26,-76},{40,-84}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="B
x_B",     fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{88,24},{102,16}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="D
x_D"),
        Text(
          extent={{-72,28},{-72,28}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="F_1"),
        Text(
          extent={{-72,-12},{-72,-12}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="F_2"),
        Text(
          extent={{-72,18},{-72,16}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="z_1",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-72,-22},{-72,-24}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          fontSize=12,
          textStyle={TextStyle.Bold},
          textString="z_2"),
        Text(
          extent={{-2,58},{-2,58}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="R",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Text(
          extent={{-12,84},{-26,84}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="V_1",
          fontSize=12,
          textStyle={TextStyle.Bold}),
        Polygon(
          points={{2,0},{-2,2},{-2,-2},{2,0}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None,
          origin={38,66},
          rotation=270),
        Line(
          points={{56,46},{48,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{54,48},{44,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{50,48},{40,38}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{32,48},{26,42}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{38,22},{42,18}},
          lineColor={0,0,255},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-38,-78},{-34,-82}},
          lineColor={0,0,255},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{32,24},{32,16},{26,20},{32,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{20,24},{20,16},{26,20},{20,24}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Rectangle(
          extent={{22,34},{30,28}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{26,20},{26,28}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-92,26},{-92,18},{-86,22},{-92,26}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-80,26},{-80,18},{-86,22},{-80,26}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-86,30},{-86,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-90,36},{-82,30}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-90,-14},{-90,-22},{-84,-18},{-90,-14}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Polygon(
          points={{-78,-14},{-78,-22},{-84,-18},{-78,-14}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-84,-10},{-84,-18}},
          color={0,0,255},
          smooth=Smooth.None),
        Rectangle(
          extent={{-88,-4},{-80,-10}},
          lineColor={0,0,255},
          radius=90,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-94,-18},{-94,78},{-66,78}},
          color={0,0,255},
          smooth=Smooth.None),
        Polygon(
          points={{-80,78},{-76,80},{-76,76},{-80,78}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Ellipse(
          extent={{-60,78},{-68,70}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-66,74},{-62,76},{-62,72},{-66,74}},
          lineColor={0,0,255},
          fillPattern=FillPattern.Sphere,
          smooth=Smooth.None),
        Line(
          points={{-62,74},{-48,74}},
          color={0,0,255},
          smooth=Smooth.None)}));
end DistillationTenPlates;
