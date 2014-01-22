within Annex60.Fluid.FMI.BaseClasses;
model Outlet "Model for exposing a fluid outlet to the FMI interface"

  replaceable package Medium =
      Modelica.Media.Interfaces.PartialMedium "Medium model within the source"
     annotation (choicesAllMatching=true);

  Interfaces.Outlet outlet(redeclare final package Medium=Medium)
    "Fluid outlet" annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));

  Modelica.Fluid.Interfaces.FluidPort_a port(redeclare final package Medium=Medium)
    "Fluid port"
                annotation (Placement(
        transformation(extent={{-110,-10},{-90,10}}),
                                                    iconTransformation(extent={{-110,
            -10},{-90,10}})));
equation
  port.p          = outlet.p;
  port.m_flow     = outlet.m_flow;
  port.h_outflow  = outlet.h_outflow;
  port.Xi_outflow = outlet.Xi_outflow;
  port.C_outflow  = outlet.C_outflow;

  inStream(port.h_outflow)  = outlet.h_inflow;
  inStream(port.Xi_outflow) = outlet.Xi_inflow;
  inStream(port.C_outflow)  = outlet.C_inflow;

    annotation (defaultComponentName="boundary",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{60,60},{-60,-60}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Sphere,
          fillColor={0,127,255}),
        Text(
          extent={{-150,110},{150,150}},
          textString="%name",
          lineColor={0,0,255}),
        Line(
          points={{60,0},{100,0}},
          color={0,0,255}),
        Rectangle(
          extent={{-100,20},{-60,-21}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={0,127,255}),
        Text(
          extent={{86,30},{108,12}},
          lineColor={0,0,255},
          textString="outlet")}),
    Documentation(info="<html>
<p>
Defines output values for boundary conditions:
</p>
<ul>
<li> Boundary pressure.</li>
<li> Boundary enthalpy.</li>
<li> Boundary composition (only for multi-substance or trace-substance flow).</li>
</ul>
<p>
This model can be used to send output signals from a fluid flow component
to the FMI interface.
</p>
</html>", revisions="<html>
<ul>
<li>
January 21, 2014 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}), graphics));
end Outlet;