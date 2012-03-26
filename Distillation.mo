within ;
package Distillation
  function PsatEthanol
  protected
    parameter Real A = 16.8958;
    parameter Real B = 3795.17;
    parameter Real C = 230.918;
  public
    output Modelica.SIunits.Pressure Psat;
    input Modelica.SIunits.Temp_C T;
  algorithm
    Psat :=exp(A - (B/(T + C)));
  end PsatEthanol;

  partial model Plate
    Modelica.SIunits.MassFlowRate Fin "Feed flow rate";
    parameter Modelica.SIunits.Pressure P = 101.325
      "Pressure in the distillation column";
    Modelica.SIunits.MassFraction z "Feed composition";
    parameter Real k = 0.02 "Constant for plate hydrilics";
    ConnectFlow V_in "Incoming Vapour flow";
    ConnectFlow V_out "outgoing Vapour flow";
    ConnectFlow L_in "Incoming Liquid flow";
    ConnectFlow L_out "Outgoing Liquid flow";
    Modelica.SIunits.Mass M(start=0) "Mass accumulated in the plate";
  equation
    der(M) = Fin + L_in.F + L_out.F;
    der(M*L_out.C) = Fin*z + L_in.F*L_in.C + L_out.F*L_out.C + V_in.F*V_in.C + V_out.F*V_out.C;
    L_out.F = k*M;
    L_out.C = P*V_out.C/PsatEthanol(87);
    V_in.F = V_out.F;
  end Plate;

  connector ConnectFlow
    flow Modelica.SIunits.MassFlowRate F;
    Modelica.SIunits.MassFraction  C;
  end ConnectFlow;

  model CombineModels
    Plate plate[2];
  equation
    for n in 1:1 loop
      connect(plate[n].V_in,plate[n+1].V_out);
      connect(plate[n].L_out,plate[n+1].L_in);
    end for;
    plate[1].L_in.C = 0.5;
    plate[2].V_out.C = 0.8;
    plate[1].Fin = 2.8;
    plate[1].L_in.F = 1.2;
  end CombineModels;
  annotation (Diagram(graphics={
        Rectangle(extent={{-24,60},{30,-18}}, lineColor={0,0,255}),
        Ellipse(extent={{-24,52},{30,68}}, lineColor={0,0,255}),
        Line(
          points={{-24,22},{-82,22},{-24,22},{-32,26},{-24,22},{-32,18},{-24,22}},
          color={0,0,255},
          smooth=Smooth.None),
        Ellipse(extent={{-24,-12},{30,-26}}, lineColor={0,0,255}),
        Line(
          points={{2,60},{2,86},{64,86},{58,80},{64,86},{58,92}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{2,-26},{2,-52},{58,-52},{50,-46},{58,-52},{50,-58}},
          color={0,0,255},
          smooth=Smooth.None),
        Text(
          extent={{58,-40},{58,-40}},
          lineColor={0,0,255},
          textString="Bottoms",
          textStyle={TextStyle.Bold}),
        Text(
          extent={{58,76},{58,74}},
          lineColor={0,0,255},
          textStyle={TextStyle.Bold},
          textString="Distilate"),
        Text(
          extent={{-2,26},{0,36}},
          lineColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="P_column
T_column"),
        Text(
          extent={{-52,28},{-52,28}},
          lineColor={0,0,255},
          textString="Feed"),
        Ellipse(extent={{64,82},{76,92}}, lineColor={0,0,255})}), uses(Modelica(
          version="3.2")));
end Distillation;
