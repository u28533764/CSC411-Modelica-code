within ;
package DistillationColumn
  "The Model of a Full Distillation Column Without connectors"
  function Psat
  protected
    parameter Real A = 16.8958;
    parameter Real B = 3795.17;
    parameter Real C = 230.918;
  public
    output Modelica.SIunits.Pressure P_E;
    output Modelica.SIunits.SpecificHeatCapacity Cp_E;
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
    der(x_0)=0;
  equation
    der(M) = F + L_i - L_0;
    der(M*x_0) = F*z + L_i*x_i - L_0*x_0 + V_i*y_i - V_0*y_0;
    L_0 = k*M;
    x_0 = P*y_0/Psat(t_0);
    V_i = V_0;
    der(E) = F*h_F + L_i*h_i + V_i*H_i - V_0*H_0 - L_0*h_0;
    E = M*Cp_mix*t_0;
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
    Plate plate[2];
  equation
    plate[1].x_i = 0.5;
    plate[1].t_i = 87;
    plate[1].F = 2.8;
    plate[1].L_i = 1.2;
    plate[2].F = 0;
    plate[2].V_i = 1;
    plate[2].y_i = 0.5;
    plate[2].T_i = 89;
    for n in 1:1 loop
      plate[n+1].L_i-plate[n].L_0 = 0;
      plate[n].V_i-plate[n+1].V_0 = 0;
      plate[n].y_i-plate[n+1].y_0 = 0;
      plate[n+1].x_i-plate[n].x_0 = 0;
      plate[n+1].t_i-plate[n].t_0 = 0;
      plate[n].T_i-plate[n+1].T_0 = 0;
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
  annotation (Diagram(graphics={
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
          points={{-56,22},{-84,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-56,-18},{-84,-18}},
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
          fillPattern=FillPattern.Solid)}),                       uses(Modelica(
          version="3.2")));
end DistillationColumn;
