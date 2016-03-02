within Annex60.Fluid.Movers;
model FlowControlled_m_flow
  "Fan or pump with ideally controlled mass flow rate as input signal"
  extends Annex60.Fluid.Movers.BaseClasses.PartialFlowMachine(
    preSou(m_flow_start=m_flow_start),
    final stageInputs(each final unit="kg/s")=massFlowRates,
    final constInput(final unit="kg/s")=constantMassFlowRate,
    filter(
      final y_start=m_flow_start,
      u_nominal=m_flow_nominal,
      u(final unit="kg/s"),
      y(final unit="kg/s")),
    final preVar=Annex60.Fluid.Types.PrescribedVariable.FlowRate,
    eff(exaPowCom=not default_record));

  parameter Modelica.SIunits.MassFlowRate m_flow_start(min=0)=0
    "Initial value of mass flow rate"
    annotation(Dialog(tab="Dynamics", group="Filtered speed"));

  parameter Modelica.SIunits.MassFlowRate constantMassFlowRate=m_flow_nominal
    "Constant pump mass flow rate, used when inputType=Constant"
    annotation(Dialog(enable=inputType == Annex60.Fluid.Types.InputType.Constant));

  parameter Modelica.SIunits.MassFlowRate[:] massFlowRates = m_flow_nominal*{0}
    "Vector of mass flow rate set points, used when inputType=Stage"
    annotation(Dialog(enable=inputType == Annex60.Fluid.Types.InputType.Stages));

  Modelica.Blocks.Interfaces.RealInput m_flow_in(
    final unit="kg/s",
    nominal=m_flow_nominal) if
       inputType == Annex60.Fluid.Types.InputType.Continuous
    "Prescribed mass flow rate"
    annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={0,120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-2,120})));
  Modelica.Blocks.Interfaces.RealOutput m_flow_actual(
    final unit="kg/s",
    nominal=m_flow_nominal) "Actual mass flow rate"
    annotation (Placement(transformation(extent={{100,10},{120,30}}),
        iconTransformation(extent={{100,10},{120,30}})));

equation
  if filteredSpeed then
    connect(filter.y, m_flow_actual) annotation (Line(
      points={{34.7,88},{44,88},{44,20},{110,20}},
      color={0,0,127},
      smooth=Smooth.None));
  else
    connect(inputSwitch.y, preSou.m_flow_in) annotation (Line(
      points={{1,50},{44,50},{44,8}},
      color={0,0,127},
      smooth=Smooth.None));
  end if;
    connect(m_flow_actual, preSou.m_flow_in) annotation (Line(
      points={{110,20},{44,20},{44,8}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(inputSwitch.u, m_flow_in) annotation (Line(
      points={{-22,50},{-26,50},{-26,80},{0,80},{0,120}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (defaultComponentName="fan",
  Documentation(
   info="<html>
<p>
This model describes a fan or pump with prescribed mass flow rate.
The efficiency of the device is computed based
on the efficiency and pressure curves that are defined
in record <code>per</code>, which is of type
<a href=\"modelica://Annex60.Fluid.Movers.SpeedControlled_Nrpm\">
Annex60.Fluid.Movers.SpeedControlled_Nrpm</a>.
</p>
<p>
See the
<a href=\"modelica://Annex60.Fluid.Movers.UsersGuide\">
User's Guide</a> for more information.
</p>
</html>",
      revisions="<html>
<ul>
<li>
March 2, 2016, by Filip Jorissen:<br/>
Refactored model such that it directly extends <code>PartialFlowMachine</code>.
This is for
<a href=\"https://github.com/iea-annex60/modelica-annex60/issues/417\">#417</a>.
</li>
<li>
April 2, 2015, by Filip Jorissen:<br/>
Added code for supporting stage input and constant input.
</li>
<li>
January 6, 2015, by Michael Wetter:<br/>
Revised model for OpenModelica.
</li>
<li>
February 14, 2012, by Michael Wetter:<br/>
Added filter for start-up and shut-down transient.
</li>
<li>
May 25, 2011, by Michael Wetter:<br/>
Revised implementation of energy balance to avoid having to use conditionally removed models.
</li>
<li>
July 27, 2010, by Michael Wetter:<br/>
Redesigned model to fix bug in medium balance.
</li>
<li>March 24, 2010, by Michael Wetter:<br/>
Revised implementation to allow zero flow rate.
</li>
<li>October 1, 2009
    by Michael Wetter:<br/>
       Model added to the Annex60 library.
</ul>
</html>"),
    Icon(graphics={
        Text(
          visible = inputType == Annex60.Fluid.Types.InputType.Continuous,
          extent={{22,146},{114,102}},
          textString="m_flow_in"),
        Line(
          points={{32,50},{100,50}},
          color={0,0,0},
          smooth=Smooth.None),
        Text(extent={{50,68},{100,54}},
          lineColor={0,0,127},
          textString="m_flow"),
        Text(
          visible=inputType == Annex60.Fluid.Types.InputType.Constant,
          extent={{-80,136},{78,102}},
          lineColor={0,0,255},
          textString="%m_flow_nominal")}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}})));
end FlowControlled_m_flow;
