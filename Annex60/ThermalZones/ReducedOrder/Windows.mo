within Annex60.ThermalZones.ReducedOrder;
package Windows "This Package calculates solar gain through windows"
  extends Modelica.Icons.VariantsPackage;
  package BaseClasses "BaseClasses for Windows"
    extends Modelica.Icons.BasesPackage;

    model SkylineShadowing
      "Calculation of the limit elevation angle for shadowing by a skyline (for direct solar irradiation)"
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_northAzimuth;
      import Modelica.Constants.pi;
      parameter Integer n(min = 1) "Number of corner points"
          annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Angle[n] alpha(displayUnit="deg") "Azimuth of corner points, sorted from north to east to south to west,
     azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
          annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Height[n] deltaH
        "Difference between height of corner point and the window centre"
        annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Distance[n] s
        "horizontal distance between corner point and window centre"
        annotation(dialog(group="skyline"));
      parameter Boolean[n-1] gap
        "corner points i and i+1 are gap between buildings: true, else: false"
        annotation(dialog(group="skyline"));

      Modelica.Blocks.Interfaces.RealInput solAzi(
        quantity="Angle",
        unit="rad",
        displayUnit="deg") "Solar azimuth"
         annotation (Placement(transformation(extent={{-140,-30},{-100,10}}),
            iconTransformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Interfaces.RealOutput altLim(
        quantity="Angle",
        unit="rad",
        displayUnit="deg") "limit elevation angle for shadowing by a skyline"
      annotation (Placement(transformation(extent={{100,-14},{128,14}}),
            iconTransformation(extent={{100,-10},{120,10}})));

    protected
      Modelica.SIunits.Angle[n-1] X "Calculation factor to simplify equations";
      Modelica.SIunits.Angle[n-1] Y "Calculation factor to simplify equations";
      Modelica.SIunits.Angle altLimi[n-1](displayUnit="deg")
        "limit elevation angle for shadowing by a skyline for point i and i+1";
      Modelica.SIunits.Angle gamma[n]( min=0,max=pi/2,displayUnit="deg")
        "elevation angle of the obstruction for point i";
    equation
      //Calculating gamma
      for i in 1:n loop
        Modelica.Math.tan(gamma[i])=deltaH[i]/s[i];
      end for;
      //Calculating altLim
      for i in 1:(n-1) loop
        X[i] = pi-Y[i]-(to_northAzimuth(alpha[i+1])-to_northAzimuth(alpha[i]));
        Y[i] = Modelica.Math.atan((Modelica.Math.tan(gamma[i+1])*Modelica.Math.sin(to_northAzimuth(alpha[i+1])-to_northAzimuth(alpha[i])))/
        (Modelica.Math.tan(gamma[i])-Modelica.Math.tan(gamma[i+1])*Modelica.Math.cos(to_northAzimuth(alpha[i+1])-to_northAzimuth(alpha[i]))));
        if gap[i] then
          altLimi[i]=-pi/2;
        else
        altLimi[i] = Modelica.Math.atan(Modelica.Math.sin(gamma[i])/Modelica.Math.cos(gamma[i])*Modelica.Math.sin(pi-X[i]-(to_northAzimuth(solAzi)-to_northAzimuth(alpha[i])))
        /Modelica.Math.sin(X[i]));
        end if;
      end for;
      altLim=max(altLimi);
        annotation(dialog(group="skyline"),
                  Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Bitmap(extent={{-84,-92},{84,90}}, fileName=
                  "modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/SkylineShadowing.png")}),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(revisions="<html>
<ul>
<li>May 23, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>",     info="<html>
<p>This model considers third party shadowing due to horizon vertical exaggeration and/or obstruction for direct radiation based on VDI 6007 part 3. It calculates a limit elevation angle for the sun.</p>
<p><img alt=\"SkylineShadowing\" src=\"modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/SkylineShadowing.png\" height=\"400\"/></p>
<p><img alt=\"SkylineShadowing\" src=\"modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/SkylineShadowing(2).png\" height=\"400\"/></p>
The figures show how the parameter should be set. For every considered cornerpoint of the skyline there should be an azimut, a heightdifference between the cornerpoint and the middle of the window and a distance between the window centre and the cornerpoint. The Boolean gap indicates if there is a gap in the skyline and should be set for every cornerpointpair.
In the example above it should be set {false,false,true,false,false} for the pairs {1+2,2+3,3+4,4+5,5+6}.
<p><b>References</b> </p>
<p>VDI. German Association of Engineers Guideline VDI 6007-3 June 2015. Calculation of transient thermal response of rooms and buildings - modelling of solar radiation.</p>
</html>"));
    end SkylineShadowing;

    block VentilationHeat "heat input due to ventilation with closed sunblind"
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
      parameter Real x_f(min=0,max=1) "percentage of open windowarea"
        annotation(dialog(group="window"));
      parameter Modelica.SIunits.Distance d "distance sunscreen to window"
        annotation(dialog(group="sunscreen"));
      parameter Boolean screen "if screen: true, if blind: false"
        annotation(dialog(group="sunscreen"));
      parameter Modelica.SIunits.TransmissionCoefficient tau_e
        "Transmission coefficient of sunscreen"
         annotation(dialog(group="sunscreen"));
      parameter Modelica.SIunits.ReflectionCoefficient rho_e
        "Reflection Coefficient of sunscreen"
         annotation(dialog(group="sunscreen"));
      parameter Modelica.SIunits.Angle til(displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));

       Modelica.Blocks.Interfaces.BooleanInput sunscreen
        "true: sunscreen closed, false: sunscreen open"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}}),
            iconTransformation(extent={{-120,10},{-100,30}})));
       Modelica.Blocks.Interfaces.RealInput HDirTil(final quantity=
         "RadiantEnergyFluenceRate",final unit="W/m2")
        "direct radiation on a tilted surface"
       annotation (Placement(transformation(extent={{-140,-50},{-100,-10}}),
            iconTransformation(extent={{-120,-30},{-100,-10}})));
       Modelica.Blocks.Interfaces.RealInput HDifTil(final quantity=
       "RadiantEnergyFluenceRate", final unit="W/m2")
        "diffuse radiation on a tilted surface"
       annotation (Placement(transformation(extent={{-140,20},{-100,60}}),
            iconTransformation(extent={{-120,40},{-100,60}})));
       Modelica.Blocks.Interfaces.RealInput HDifHor(final quantity=
         "RadiantEnergyFluenceRate",final unit="W/m2")
        "diffuse radiation on horizontal surface"
       annotation (Placement(transformation(extent={{-140,50},{-100,90}}),
            iconTransformation(extent={{-120,70},{-100,90}})));
       Modelica.Blocks.Interfaces.RealInput HDirNor(final quantity=
         "RadiantEnergyFluenceRate",final unit="W/m2")
        "direct radiation on horizontal surface"
       annotation (Placement(transformation(extent={{-140,-80},{-100,-40}}),
            iconTransformation(extent={{-120,-60},{-100,-40}})));
      Modelica.Blocks.Interfaces.RealOutput HVen(final quantity=
        "RadiantEnergyFluenceRate", final unit="W/m2")
        "heat input due to ventilation"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
            iconTransformation(extent={{100,-12},{124,12}})));
      Modelica.Blocks.Interfaces.RealInput alt(
        final quantity="Angle",
        final unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-140,-110},{-100,-70}}),
            iconTransformation(extent={{-120,-90},{-100,-70}})));
    protected
       parameter Modelica.SIunits.ReflectionCoefficient rho=0.2
        "ground reflection";
       Real factor_gv "calculation factor";

    equation
      //calculating factor_gv
      if d<=0.15 and screen then
         factor_gv=0.35;
      elseif d<=0.4 and d>0.15 and screen then factor_gv=0.17;
      elseif d<=0.15 and screen==false then factor_gv=0.2;
      elseif d>0.15 and d<=0.4 and screen==false then factor_gv=0.1;
      else factor_gv=0;
      end if;
      //calculating Output HVen
      if sunscreen then
        HVen=(HDirTil+HDifTil+(HDifHor+HDirNor*Modelica.Math.sin(alt))*0.5*rho*(1-Modelica.Math.cos(to_surfaceTiltVDI(til))))*(1-rho_e-tau_e)*factor_gv*x_f;
      else
        HVen=0;
      end if;
      annotation (defaultComponentName="ventilationHeat",Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Bitmap(extent={{-94,-96},{96,96}}, fileName=
                  "modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/VentilationHeat.png")}),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
The model considers additional heat input in the event of window ventilation and simultaneously closed external solar protection based on VDI 6007 part 3. 
<br>The closed external solar protection absorbs solar irradiation which is brought into the room by convection.
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>",     revisions="<html>
<ul>
<li>May 5, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end VentilationHeat;

    block SelfShadowing
      "Self-shadowing due to projections for direct radiation"
      parameter Integer n(min = 1) "number of windows"
        annotation(dialog(group="window"));
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_northAzimuth;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
      import Modelica.Constants.pi;
      parameter Modelica.SIunits.Length b[n] "width of window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Height h[n] "height of window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length bLef[n] "window projection left"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length bRig[n] "window projection right"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length dLef[n]
        "distance between projection (left) and window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length dRig[n]
        "distance between projection (right) and window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length bAbo[n] "window projection above"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length bBel[n] "window projection below"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length dAbo[n]
        "distance between projection (above) and window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Length dBel[n]
        "distance between projection (below) and window"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Angle azi[n](displayUnit="degree")
        "Surface azimuth. azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
        annotation (Dialog(group="Window parameter"));
      parameter Modelica.SIunits.Angle til[n](displayUnit="degree")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="Window parameter"));

       Modelica.Blocks.Interfaces.RealInput incAng[n](
        quantity="Angle",
        unit="rad",
        displayUnit="degree")
        "Incidence angle of the sun beam on a tilted surface"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}}),
            iconTransformation(extent={{-120,-80},{-100,-60}})));
      Modelica.Blocks.Interfaces.RealInput solAzi(
        quantity="Angle",
        unit="rad",
        displayUnit="deg") "Solar Azimuth"
         annotation (Placement(transformation(extent={{-140,40},{-100,80}}),
            iconTransformation(extent={{-120,60},{-100,80}})));
      Modelica.Blocks.Interfaces.RealInput alt(
        quantity="Angle",
        unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-140,-30},{-100,10}}),
            iconTransformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Interfaces.RealOutput x_As[n](min=0,
        final unit="1") "not shaded percentage of window area"
        annotation (Placement(transformation(extent={{100,-20},{140,20}}),
            iconTransformation(extent={{100,-8},{118,10}})));

    protected
      Real e_hn[n] "horizontal calculation factor";
      Real e_vn[n] "vertical calculation factor";
      Modelica.SIunits.Distance x1[n] "auxiliary variable for shadow from left";
      Modelica.SIunits.Distance x2[n]
        "auxiliary variable for shadow from right";
      Modelica.SIunits.Distance x3[n]
        "auxiliary variable for shadow from above";
      Modelica.SIunits.Distance x4[n]
        "auxiliary variable for shadow from below";
      Modelica.SIunits.Distance s_h[n] "horizontal reduction of window";
      Modelica.SIunits.Distance s_v[n] "vertical reduction of window";
      Modelica.SIunits.Area A_S[n] "auxiliary variable for effective Area";
      Modelica.SIunits.Area A_s[n] "effective windowarea";
    equation
      for i in 1:n loop
      //Calculating e_hn and e_vn
        if Modelica.Math.cos(incAng[i])<=0 then
           e_hn[i]=10^20;
           e_vn[i]=10^20;
        else
          e_hn[i]=Modelica.Math.sin(to_northAzimuth(azi[i])-to_northAzimuth(solAzi))*Modelica.Math.cos(alt)/Modelica.Math.cos(incAng[i]);
          e_vn[i]=(Modelica.Math.sin(alt)*Modelica.Math.sin(to_surfaceTiltVDI(til[i]))-Modelica.Math.cos(alt)*
          Modelica.Math.cos(to_surfaceTiltVDI(til[i]))*Modelica.Math.cos(to_northAzimuth(azi[i])-to_northAzimuth(solAzi)))/Modelica.Math.cos(incAng[i]);
        end if;

      //Calculating s_h
      x1[i]=e_hn[i]*bLef[i]-dLef[i];
      x2[i]=-e_hn[i]*bRig[i]-dRig[i];
      if e_hn[i]>=0 then
        if x1[i]<0 then
          s_h[i]=0;
        elseif (e_hn[i]*bRig[i]+dRig[i])<0 then
          s_h[i]=x1[i]-e_hn[i]*bRig[i]-dRig[i];
        else
          s_h[i]=x1[i];
        end if;
      else
        if x2[i]<0 then
          s_h[i]=0;
        elseif (-e_hn[i]*bLef[i]+dLef[i])<0 then
          s_h[i]=x2[i]+e_hn[i]*bLef[i]-dLef[i];
        else
          s_h[i]=x2[i];
        end if;
      end if;
      //Calculating s_v
      x3[i]=e_vn[i]*bAbo[i]-dAbo[i];
      x4[i]=-e_vn[i]*bBel[i]-dBel[i];
      if e_vn[i]>=0 then
        if x3[i]<0 then
          s_v[i]=0;
        elseif (e_vn[i]*bBel[i]+dBel[i])<0 then
          s_v[i]=x3[i]-e_hn[i]*bBel[i]-dBel[i];
        else
          s_v[i]=x3[i];
        end if;
      else
        if x4[i]<0 then
          s_v[i]=0;
        elseif (-e_vn[i]*bAbo[i]+dAbo[i])<0 then
          s_v[i]=x4[i]-e_vn[i]*bAbo[i]-dAbo[i];
        else
          s_v[i]=x4[i];
        end if;
      end if;
      //Calculating A_s
      if (b[i]-s_h[i])<0 then
         A_S[i]=0;
      elseif (h[i]-s_v[i])<0 then
         A_S[i]=0;
      else
         A_S[i]=(b[i]-s_h[i])*(h[i]-s_v[i]);
      end if;
      A_s[i] = max(0,A_S[i]);
      x_As[i]=A_s[i]/(b[i]*h[i]);
      end for;
         annotation (defaultComponentName="selfShadowing",
         Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Bitmap(extent={{-90,-98},{94,92}}, fileName="modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/SelfShadowing.png")}),
                                                                        Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This model considers self-shadowing of windows due to projections for direct radiation based on VDI 6007 part 3. It calculates what part of the windowarea is effective. </p>
<p><img alt=\"SelfShadowing\" src=\"modelica://Annex60/Resources/Images/ThermalZones/ReducedOrder/Windows/BaseClasses/SelfShadowing.png\" height=\"400\"/></p>
<p>The image above shows how the parameters should be set and is based on VDI 6007 part 3. Parameters with Index 2 are alligned on the other side</p>
<p>(i.e.: dRig is the distance between the projection and the window on the right handside, dBel is the distance between the projection below and the window). eh and ev are calculated within the model and are shown for demonstration.</p>
<p>The connectors are all solar geometry dimensions and can be calculated by the SolarGeometry package of Annex60. </p>
<p>An Example on how to use this model is the SelfShadowingExample in the Example package. </p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>",     revisions="<html>
<ul>
<li>May 23, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end SelfShadowing;

    model Illumination
      "Determining the activation and deactivation times of the illumination"
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_northAzimuth;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
      import Modelica.SIunits.Conversions.to_deg;
      import Modelica.Constants.pi;
      parameter Real D "daylight quotient";
      parameter Modelica.SIunits.Illuminance e_ILim1
        "internal illumninance required in reference point in the morning and evening";
      parameter Modelica.SIunits.Illuminance e_ILim2
        "internal illumainance required in reference point during working hours";
      parameter Boolean office "if true: room is office";
      final parameter Modelica.SIunits.LuminousEfficacy k_mDifCov=115
        "radiation equivalent for uniformly overcast skies";

      //Window parameter
      parameter Integer n(min=1) "number of windows"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));
      parameter Real r[n] "frame share"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Area A[n] "windowarea"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
        "degree of light transmission"
        annotation (Dialog(group="window"));
      final parameter Modelica.SIunits.ReflectionCoefficient rho=0.2
        "degree of ground reflection";

      Modelica.Blocks.Interfaces.BooleanOutput Illumination
        "If Illumination=true: activation of Illumination"
        annotation (Placement(transformation(extent={{98,-10},{118,10}}),
            iconTransformation(extent={{100,-10},{120,10}})));

      Modelica.Blocks.Interfaces.RealInput HVis[n]( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2")
        "solar energy entering the room in the visible area"
        annotation (Placement(transformation(extent={{-120,50},{-100,70}}),
            iconTransformation(extent={{-120,50},{-100,70}})));
              Modelica.Blocks.Interfaces.RealInput CorTaue_DifCov[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for diffuse irradiation during covered sky"
        annotation (Placement(transformation(extent={{-120,-70},{-100,-50}}),
            iconTransformation(extent={{-120,-70},{-100,-50}})));

      Modelica.Blocks.Interfaces.RealInput CorTaue_Gro[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for ground reflexion radiation"
        annotation (Placement(transformation(extent={{-120,-10},{-100,10}}),
            iconTransformation(extent={{-120,-10},{-100,10}})));
    //protected
      constant Modelica.SIunits.Time day=86400 "Number of seconds in a day";
      constant Modelica.SIunits.Time week=604800 "Number of seconds in a week";

      Modelica.SIunits.Illuminance e_ILim
        "internal illumance in reference point";
      Real r_DifCov[n] "conversion factor";

      Modelica.SIunits.EnergyFlowRate HLimVisi[n] "thresholds within the room";
      Modelica.SIunits.EnergyFlowRate HLimVis "sum of H_LimInsi";

      Modelica.SIunits.EnergyFlowRate HVisi[n]
        "solar energy entering the room in the visible area";
      Modelica.SIunits.EnergyFlowRate HVisSum "sum of HVisi";

    equation
      //pick value for e_ILim
      if (time-integer(time/day)*day)>64800 or (time-integer(time/day)*day)<25200 or (office and time-integer(time/week)*week>432000) then
        e_ILim=0;
      elseif (time-integer(time/day)*day)>28800 and (time-integer(time/day)*day)<57600 then
        e_ILim=e_ILim2;
      else
        e_ILim=e_ILim1;
      end if;
      //calculating H_LimIns
      for i in 1:n loop
        r_DifCov[i]=0.182*(1.178*(1+Modelica.Math.cos(to_surfaceTiltVDI(til[i])))+(pi-to_surfaceTiltVDI(til[i]))*
        Modelica.Math.cos(to_surfaceTiltVDI(til[i]))+Modelica.Math.sin(to_surfaceTiltVDI(til[i])));
        HLimVisi[i]=e_ILim/(D*k_mDifCov)*(r_DifCov[i]*T_L[i]*
        CorTaue_DifCov[i]+0.5*rho*(1-Modelica.Math.cos(to_surfaceTiltVDI(til[i])))*T_L[i]*CorTaue_Gro[i])*
        (1-r[i])*A[i];
        HVisi[i]=HVis[i]*(1-r[i])*A[i];
      end for;
      HLimVis=sum(HLimVisi);
      HVisSum=sum(HVisi);
      //comparing H_RoomL with H_LimIns
      Illumination = (HVisSum<HLimVis);
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Bitmap(extent={{-92,-102},{100,114}}, fileName=
                  "modelica://Annex60/../../../vdi/Resources/Icons/Illumination.png")}),
                                                                     Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(revisions="<html>
<ul>
<li>May 23, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>", info="<html>
This model calculates the activation and deactivation times of the illumination and gives it back as the Boolean \"Illumination\".
It is based on VDI 6007 part 3. <br>
The total solar energy entering the room, which can be calculated by <a href=\"Windows.Window\">Window</a> or <a href=\"Windows.ShadedWindow\">ShadedWindow</a>, is compared to a limit value based on the parameters.
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>"));
    end Illumination;

    model Sunblind "Calculates if sunblind of window is active"
      extends Modelica.Blocks.Icons.Block;
      parameter Modelica.SIunits.RadiantEnergyFluenceRate lim
        "Limit for the sunscreen to become active";

      Modelica.Blocks.Interfaces.RealInput HDifTil
        "Hemispherical diffuse solar irradiation on a tilted surface from the sky"
        annotation (Placement(transformation(extent={{-116,52},{-100,68}}),
            iconTransformation(extent={{-120,48},{-100,68}})));
      Modelica.Blocks.Interfaces.RealInput HDirTil
        "Direct irradition on tilted surface"
        annotation (Placement(transformation(extent={{-114,-64},{-100,-50}}),
            iconTransformation(extent={{-120,-70},{-100,-50}})));
      Modelica.Blocks.Interfaces.BooleanOutput sunscreen
        "if true: sunscreen is closed, else sunscreen is open;"
          annotation (Placement(transformation(extent={{98,-10},{118,10}}),
            iconTransformation(extent={{100,-10},{120,10}})));

    equation
      if (HDifTil+HDirTil)>lim then
        sunscreen=true;
      else
        sunscreen=false;
      end if;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This model computes whether the sunscreen is active or not. Therefor it compares the irradiation on the window with a limit for the sunscreen to be active set as a parameter.</p>

</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end Sunblind;

    model HeatIllumination "heating energy due to Illumination"
      extends Modelica.Blocks.Icons.Block;
      parameter Modelica.SIunits.EnergyFlowRate HIll1
        "Energyoutput of Illumination in the morning and evening";
      parameter Modelica.SIunits.EnergyFlowRate HIll2
        "Energyoutput of Illumination during daytime";
      Modelica.Blocks.Interfaces.BooleanInput Illumination
        annotation (Placement(transformation(extent={{-140,-30},{-100,10}}),
            iconTransformation(extent={{-120,-10},{-100,10}})));
      Modelica.Blocks.Interfaces.RealOutput HIll(final quantity=
        "EnergyFlowRate", final unit="W") "Energoutput of Illumination"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
            iconTransformation(extent={{100,-12},{124,12}})));
    protected
      constant Modelica.SIunits.Time day=86400;
    equation
      if Illumination==false then
        HIll=0;
      else
        if time-integer(time/day)*day < 25200 or time-integer(time/day)*day > 68400 then
          HIll=0;
        elseif time-integer(time/day)*day > 25200 and time-integer(time/day)*day < 64800 then
          HIll=HIll2;
        else
          HIll=HIll1;
        end if;
      end if;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This model calculates the heat input in the room due to illumination.</p>
</html>"));
    end HeatIllumination;

    model HVisible
      "Calculates the solar energy entering the room in the visible area"
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
      import Modelica.SIunits.Conversions.to_deg;
      parameter Integer n(min=1) "number of windows"
        annotation (Dialog(group="window"));

      parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
        "degree of light transmission"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDir[n]
        "degree of light transmission for direct irradiation, with sunscreen"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDif[n]
        "degree of light transmission for diffuse irradiation, with sunscreen"
        annotation (Dialog(group="window"));

      parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));
      final parameter Modelica.SIunits.ReflectionCoefficient rho=0.2
        "degree of ground reflection";

      Modelica.Blocks.Interfaces.BooleanInput sunscreen[n]
        "true: sunscreen closed, false: sunscreen open"
        annotation (Placement(transformation(extent={{-120,-40},{-80,0}}),
            iconTransformation(extent={{-114,-6},{-100,8}})));
      Modelica.Blocks.Interfaces.RealInput CorTaue_Dir[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for direct irradiation"
        annotation (Placement(transformation(extent={{-128,78},{-108,98}}),
            iconTransformation(extent={{-114,-106},{-100,-92}})));

      Modelica.Blocks.Interfaces.RealInput CorTaue_DifCle[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for diffuse irradiation during clear sky"
        annotation (Placement(transformation(extent={{-120,-92},{-100,-72}}),
            iconTransformation(extent={{-114,-86},{-100,-72}})));

      Modelica.Blocks.Interfaces.RealInput CorTaue_DifCov[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for diffuse irradiation during covered sky"
        annotation (Placement(transformation(extent={{-120,-72},{-100,-52}}),
            iconTransformation(extent={{-114,-66},{-100,-52}})));

      Modelica.Blocks.Interfaces.RealInput CorTaue_Gro[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Correction value for transluance for ground reflexion radiation"
        annotation (Placement(transformation(extent={{-120,-52},{-100,-32}}),
            iconTransformation(extent={{-114,-46},{-100,-32}})));
         Modelica.Blocks.Interfaces.RealInput alt(
        final quantity="Angle",
        final unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-120,-32},{-100,-12}}),
            iconTransformation(extent={{-114,-26},{-100,-12}})));

       Modelica.Blocks.Interfaces.RealInput HDirTil[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct irradition on tilted surface"
        annotation (Placement(transformation(extent={{-114,76},{-100,90}}),
            iconTransformation(extent={{-114,76},{-100,90}})));

       Modelica.Blocks.Interfaces.RealInput HDirNor(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct normal radiation."
        annotation (Placement(transformation(extent={{-122,84},{-100,106}}),
            iconTransformation(extent={{-114,92},{-100,106}})));

     Modelica.Blocks.Interfaces.RealInput HDifHorCov(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Horizontal diffuse solar radiation at covered sky."
        annotation (Placement(transformation(extent={{-116,26},{-100,42}}),
            iconTransformation(extent={{-114,28},{-100,42}})));

      Modelica.Blocks.Interfaces.RealInput HDifHorCle(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Horizontal diffuse solar radiation at clear sky."
        annotation (Placement(transformation(extent={{-116,10},{-100,26}}),
            iconTransformation(extent={{-114,12},{-100,26}})));

      Modelica.Blocks.Interfaces.RealInput HDifTilCov[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface at covered sky"
        annotation (Placement(transformation(extent={{-116,42},{-100,58}}),
            iconTransformation(extent={{-114,44},{-100,58}})));
      Modelica.Blocks.Interfaces.RealInput HDifTilCle[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface at clear sky"
        annotation (Placement(transformation(extent={{-116,58},{-100,74}}),
            iconTransformation(extent={{-114,60},{-100,74}})));

      Modelica.Blocks.Interfaces.RealOutput HVis[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "solar energy entering the room in the visible area"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
            iconTransformation(extent={{100,-10},{120,10}})));

    protected
      parameter Real Cor_KMDifCov=1
        "correction factor for diffuse irradiation at uniformly overcast skies according to DIN 5034-2";
      Real Cor_KMDir
        "correction factor for direct solar irradiation according to DIN 5034-2";
      Real Cor_KMDifCle
        "correction factor for diffuse irradiation at cloudless clear skies according to DIN 5034-2";
      Modelica.SIunits.TransmissionCoefficient T_LDifx[n]
        "calculation variable for the degree of light transmission for diffuse irradiation";
      Modelica.SIunits.TransmissionCoefficient T_LDirx[n]
        "calculation variable for the degree of light transmission for direct irradiation";
      Modelica.SIunits.EnergyFlowRate H_EvaHor[n]
        "evaluated solar irradiation onto the horizontal for determining the ground reflexion radiation";

    equation
      //calculating H_RoomL
      Cor_KMDir=(17.72+4.4585*to_deg(alt)-0.087563*to_deg(alt)^2+7.39487*10^(-4)*to_deg(alt)^3-2.167*10^(-6)*
      to_deg(alt)^4-8.4132*10^(-10)*to_deg(alt)^5)/115;
      Cor_KMDifCle=(15.1+3.1076*to_deg(alt)+0.0048*to_deg(alt)^2-0.0014*to_deg(alt)^3+2.04*10^(-5)*
      to_deg(alt)^4-8.91*10^(-8)*to_deg(alt)^5)/115;

      for i in 1:n loop
        if sunscreen[i] then
          T_LDifx[i]=T_LTotDif[i];
          T_LDirx[i]=T_LTotDir[i];
        else
          T_LDifx[i]=T_L[i];
          T_LDirx[i]=T_L[i];
        end if;
        H_EvaHor[i]=(HDirNor*Cor_KMDir*Modelica.Math.sin(alt)+HDifHorCle*Cor_KMDifCle+HDifHorCov*
        Cor_KMDifCov)*T_LDifx[i];
        HVis[i]=(HDirTil[i]*T_LDirx[i]*CorTaue_Dir[i]*Cor_KMDir+HDifTilCle[i]*
        T_LDifx[i]*CorTaue_DifCle[i]*Cor_KMDifCle+HDifTilCov[i]*T_LDifx[i]*
        CorTaue_DifCov[i]*Cor_KMDifCov+H_EvaHor[i]*0.5*rho*(1-Modelica.Math.cos(
        to_surfaceTiltVDI(til[i])))*CorTaue_Gro[i]);
      end for;
        annotation (defaultComponentName="HVis",Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics),
        Documentation(revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end HVisible;

    model HWindow "Calculates the solar heat input through the window"
      extends Modelica.Blocks.Icons.Block;
      import
        Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
      parameter Integer n(min=1) "number of windows"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g[n]
        "Total energy transmittance of windows"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDir[n]
        "Total energy transmittance of windows with closed sunscreen for direct radiation"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDif[n]
        "Total energy transmittance of windows with closed sunscreen for diffuse radiation"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));
      final parameter Modelica.SIunits.ReflectionCoefficient rho=0.2
        "degree of ground reflection";
      Modelica.Blocks.Interfaces.RealInput alt(
        final quantity="Angle",
        final unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-120,8},{-100,28}}),
            iconTransformation(extent={{-114,14},{-100,28}})));
      Modelica.Blocks.Interfaces.RealInput HDifTilCov[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface at covered sky"
        annotation (Placement(transformation(extent={{-116,-68},{-100,-52}}),
            iconTransformation(extent={{-114,-66},{-100,-52}})));
      Modelica.Blocks.Interfaces.RealInput HDifTilCle[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface at clear sky"
        annotation (Placement(transformation(extent={{-116,-88},{-100,-72}}),
            iconTransformation(extent={{-114,-86},{-100,-72}})));
      Modelica.Blocks.Interfaces.RealInput HDirNor(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct normal radiation."
        annotation (Placement(transformation(extent={{-122,-56},{-100,-34}}),
            iconTransformation(extent={{-114,-48},{-100,-34}})));
      Modelica.Blocks.Interfaces.RealInput HDifHor(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Horizontal diffuse solar radiation."
        annotation (Placement(transformation(extent={{-116,-108},{-100,-92}}),
            iconTransformation(extent={{-114,-106},{-100,-92}})));
      Modelica.Blocks.Interfaces.RealInput HDirTil[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct irradition on tilted surface"
        annotation (Placement(transformation(extent={{-120,-26},{-106,-12}}),
            iconTransformation(extent={{-114,-26},{-100,-12}})));
      Modelica.Blocks.Interfaces.BooleanInput sunscreen[n]
        "true: sunscreen closed, false: sunscreen open"
        annotation (Placement(transformation(extent={{-120,-40},{-80,0}}),
            iconTransformation(extent={{-114,-6},{-100,8}})));
      Modelica.Blocks.Interfaces.RealInput CorG_Dir[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Transmission coefficient correction factor for direct radiation"
        annotation (Placement(transformation(extent={{-120,86},{-100,106}}),
        iconTransformation(extent={{-114,92},{-100,106}})));
      Modelica.Blocks.Interfaces.RealInput CorG_DifCle[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Transmission coefficient correction factor for diffuse radiation while clear sky"
        annotation (Placement(transformation(extent={{-130,68},{-110,88}}),
        iconTransformation(extent={{-114,74},{-100,88}})));
      Modelica.Blocks.Interfaces.RealInput CorG_DifCov[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Transmission coefficient correction factor for diffuse radiation while covered sky"
        annotation (Placement(transformation(extent={{-120,48},{-100,68}}),
        iconTransformation(extent={{-114,54},{-100,68}})));
      Modelica.Blocks.Interfaces.RealInput CorG_Gro[n](
        final quantity="TransmissionCoefficient",
        final unit="1")
        "Transmission coefficient correction factor for ground reflection radiation"
        annotation (Placement(transformation(extent={{-120,28},{-100,48}}),
        iconTransformation(extent={{-114,34},{-100,48}})));
      Modelica.Blocks.Interfaces.RealOutput HWin[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Solar radiation transmitted through aggregated window"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
            iconTransformation(extent={{100,-10},{120,10}})));
    protected
      Modelica.SIunits.TransmissionCoefficient g_Dirx[n]
        "calculation variable to determine the active total energy transmittance of windows for direct radiation";
      Modelica.SIunits.TransmissionCoefficient g_Difx[n]
        "calculation variable to determine the active total energy transmittance of windows for diffuse radiation";
    equation
      for i in 1:n loop
       if sunscreen[i] then
          g_Dirx[i]=g_TotDir[i];
          g_Difx[i]=g_TotDif[i];
        else
          g_Dirx[i]=g[i];
          g_Difx[i]=g[i];
        end if;
        HWin[i]=(HDirTil[i]*g_Dirx[i]*CorG_Dir[i]+HDifTilCov[i]*g_Difx[i]*CorG_DifCov[i]+HDifTilCle[i]*g_Difx[i]*CorG_DifCle[i]+0.5*rho*(1-Modelica.Math.cos(to_surfaceTiltVDI(til[i])))
          *(HDirNor*Modelica.Math.sin(alt)+HDifHor)*g_Difx[i]*CorG_Gro[i]);
      end for
      annotation (defaultComponentName="HWin",Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
      annotation (Documentation(revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end HWindow;

    model Window "Calculation of solar energy transmitted through windows"
      import vdi6007 = Annex60.ThermalZones.ReducedOrder.Windows;
      parameter Integer n(min = 1) "number of windows"
        annotation(dialog(group="window"));
      parameter Modelica.SIunits.CoefficientOfHeatTransfer UWin
        "Thermal transmission coefficient of whole window"
        annotation(dialog(group="window"));
       parameter Modelica.SIunits.TransmissionCoefficient g[n]
        "Total energy transmittance of windows"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDir[n]
        "Total energy transmittance of windows with closed sunscreen for direct radiation"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDif[n]
        "Total energy transmittance of windows with closed sunscreen for diffuse radiation"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
        "degree of light transmission for direct irradiation"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDir[n]
        "degree of light transmission for direct irradiation, with sunscreen"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDif[n]
        "degree of light transmission for diffuse irradiation, with sunscreen"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.RadiantEnergyFluenceRate lim
        "Limit for the sunscreen to become active"
        annotation(dialog(group="sunscreen"));
      parameter Modelica.SIunits.Angle xi(  displayUnit="degree")= 0
        "elevation angle";
      parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));

      Modelica.Blocks.Interfaces.RealInput incAng[n](
        final quantity="Angle",
        final unit="rad",
        displayUnit="degree")
        "Incidence angles of the sun beam on a tilted surface"
        annotation (Placement(transformation(extent={{-124,68},{-100,92}}),
        iconTransformation(extent={{-120,72},{-100,92}})));
       Modelica.Blocks.Interfaces.RealInput alt(
        final quantity="Angle",
        final unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-124,12},{-100,36}}),
            iconTransformation(extent={{-120,16},{-100,36}})));

       Modelica.Blocks.Interfaces.RealInput nTot( min=0,max=1,
        final unit="1") "Total sky Cover"
        annotation (Placement(transformation(extent={{-124,42},{-100,66}}),
            iconTransformation(extent={{-120,46},{-100,66}})));
       Modelica.Blocks.Interfaces.RealInput HDirTil[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct irradition on tilted surface"
        annotation (Placement(transformation(extent={{-124,-92},{-100,-68}}),
            iconTransformation(extent={{-120,-88},{-100,-68}})));

        Modelica.Blocks.Interfaces.RealInput HDirNor(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Direct normal radiation."
        annotation (Placement(transformation(extent={{-124,-64},{-100,-40}}),
            iconTransformation(extent={{-120,-60},{-100,-40}})));

       Modelica.Blocks.Interfaces.RealInput HDifHor(
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "Horizontal diffuse solar radiation."
        annotation (Placement(transformation(extent={{-124,-36},{-100,-12}}),
            iconTransformation(extent={{-120,-32},{-100,-12}})));

        Modelica.Blocks.Interfaces.RealInput HDifTil[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface from the sky"
        annotation (Placement(transformation(extent={{-124,-12},{-100,12}}),
            iconTransformation(extent={{-120,-8},{-100,12}})));

        Modelica.Blocks.Interfaces.RealOutput HVis[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2") "solar energy entering the room in the visible area"
        annotation (Placement(transformation(extent={{100,30},{120,50}}),
            iconTransformation(extent={{100,30},{120,50}})));

        Modelica.Blocks.Interfaces.RealOutput HWin[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Solar radiation transmitted through aggregated window"
        annotation (Placement(transformation(extent={{100,-50},{120,-30}}),
            iconTransformation(extent={{100,-50},{120,-30}})));

      vdi6007.BaseClasses.Sunblind sunscreen[n](lim=lim)
        annotation (Placement(transformation(extent={{-68,-42},{-56,-30}})));
      vdi6007.SolarGain.CorrectionGTaueDoublePane CorGTaue(
        final n=n,
        final UWin=UWin,
        final xi=xi,
        final til=til)
        annotation (Placement(transformation(extent={{-30,-22},{-6,8}})));
      vdi6007.BaseClasses.HVisible HVisible(
        final n=n,
        final T_LTotDir=T_LTotDir,
        final T_LTotDif=T_LTotDif,
        final til=til,
        final T_L=T_L)
        annotation (Placement(transformation(extent={{56,18},{100,62}})));
      vdi6007.BaseClasses.Conversions.HDif_toClearCovered hDif_toClearCovered[n]
        annotation (Placement(transformation(extent={{-74,-18},{-54,2}})));
      vdi6007.BaseClasses.HWindow HWindow(
        final g=g,
        final til=til,
        final n=n,
        final g_TotDir=g_TotDir,
        final g_TotDif=g_TotDif)
        annotation (Placement(transformation(extent={{54,-62},{100,-18}})));

    equation
      connect(HDifTil, sunscreen.HDifTil) annotation (Line(points={{-112,0},{-92,0},
              {-92,-32.52},{-68.6,-32.52}},
                                          color={0,0,127}));
      connect(HDirTil, sunscreen.HDirTil) annotation (Line(points={{-112,-80},{-81.5,
              -80},{-81.5,-39.6},{-68.6,-39.6}},
                                               color={0,0,127}));
      connect(sunscreen.sunscreen, CorGTaue.sunscreen) annotation (Line(points={{-55.4,
              -36},{-38,-36},{-38,-10},{-28.8,-10}},
                                                color={255,0,255}));
      connect(HVisible.HVis, HVis) annotation (Line(points={{102.2,40},{102.2,
              40},{110,40}}, color={0,0,127}));
      connect(sunscreen.sunscreen, HVisible.sunscreen) annotation (Line(points=
              {{-55.4,-36},{-38,-36},{-38,40.22},{54.46,40.22}}, color={255,0,
              255}));
      connect(CorGTaue.CorTaue_Dir, HVisible.CorTaue_Dir) annotation (Line(
            points={{-7.2,-4},{22,-4},{22,18.22},{54.46,18.22}}, color={0,0,127}));
      connect(CorGTaue.CorTaue_DifCle, HVisible.CorTaue_DifCle) annotation (
          Line(points={{-7.2,-1},{18,-1},{18,22.62},{54.46,22.62}}, color={0,0,
              127}));
      connect(CorGTaue.CorTaue_DifCov, HVisible.CorTaue_DifCov) annotation (
          Line(points={{-7.2,2},{14,2},{14,27.02},{54.46,27.02}}, color={0,0,
              127}));
      connect(CorGTaue.CorTaue_Gro, HVisible.CorTaue_Gro) annotation (Line(
            points={{-7.2,5},{10,5},{10,31.42},{54.46,31.42}}, color={0,0,127}));
      connect(alt, HVisible.alt) annotation (Line(points={{-112,24},{-30,24},{-30,
              35.82},{54.46,35.82}}, color={0,0,127}));
      connect(hDif_toClearCovered.HDifTilCle, HVisible.HDifTilCle) annotation (
          Line(points={{-53,-2},{-48,-2},{-48,54.74},{54.46,54.74}}, color={0,0,
              127}));
      connect(hDif_toClearCovered.HDifTilCov, HVisible.HDifTilCov) annotation (
          Line(points={{-53,-6},{-46,-6},{-46,51.22},{54.46,51.22}}, color={0,0,
              127}));
      connect(hDif_toClearCovered[1].HDifHorCov, HVisible.HDifHorCov)
        annotation (Line(points={{-53,-10},{-44,-10},{-44,47.7},{54.46,47.7}},
            color={0,0,127}));
      connect(hDif_toClearCovered[1].HDifHorCle, HVisible.HDifHorCle)
        annotation (Line(points={{-53,-14},{-42,-14},{-42,44.18},{54.46,44.18}},
            color={0,0,127}));
      connect(HDifTil, hDif_toClearCovered.HDifTil) annotation (Line(points={{-112,0},
              {-92,0},{-92,-2},{-75,-2}},  color={0,0,127}));

      connect(HWindow.HWin, HWin) annotation (Line(points={{102.3,-40},{102.3,-40},
              {110,-40}}, color={0,0,127}));
      connect(sunscreen.sunscreen, HWindow.sunscreen) annotation (Line(points={
              {-55.4,-36},{-2,-36},{-2,-39.78},{52.39,-39.78}}, color={255,0,
              255}));
      connect(CorGTaue.CorG_Dir, HWindow.CorG_Dir) annotation (Line(points={{-7.2,
              -10},{22,-10},{22,-18.22},{52.39,-18.22}}, color={0,0,127}));
      connect(CorGTaue.CorG_DifCle, HWindow.CorG_DifCle) annotation (Line(
            points={{-7.2,-13},{18.4,-13},{18.4,-22.18},{52.39,-22.18}}, color=
              {0,0,127}));
      connect(CorGTaue.CorG_DifCov, HWindow.CorG_DifCov) annotation (Line(
            points={{-7.2,-16},{16,-16},{16,-26.58},{52.39,-26.58}}, color={0,0,
              127}));
      connect(CorGTaue.CorG_Gro, HWindow.CorG_Gro) annotation (Line(points={{-7.2,
              -19},{12.4,-19},{12.4,-30.98},{52.39,-30.98}}, color={0,0,127}));
      connect(HDifHor, HWindow.HDifHor) annotation (Line(points={{-112,-24},{-94,
              -24},{-94,-72},{44,-72},{44,-61.78},{48,-61.78},{52.39,-61.78}},
            color={0,0,127}));
      connect(hDif_toClearCovered.HDifTilCle, HWindow.HDifTilCle) annotation (
          Line(points={{-53,-2},{-48,-2},{-48,-57.38},{52.39,-57.38}}, color={0,
              0,127}));
      connect(hDif_toClearCovered.HDifTilCov, HWindow.HDifTilCov) annotation (
          Line(points={{-53,-6},{-50,-6},{-50,-52.98},{52.39,-52.98}}, color={0,
              0,127}));
      connect(HDirNor, HWindow.HDirNor) annotation (Line(points={{-112,-52},{-85.5,
              -52},{-85.5,-49.02},{52.39,-49.02}}, color={0,0,127}));
      connect(HDirTil, HWindow.HDirTil) annotation (Line(points={{-112,-80},{
              38.5,-80},{38.5,-44.18},{52.39,-44.18}}, color={0,0,127}));
      connect(incAng, CorGTaue.incAng) annotation (Line(points={{-112,80},{-32,80},{
              -32,-4},{-28.8,-4}},  color={0,0,127}));
      connect(HDirTil, HVisible.HDirTil) annotation (Line(points={{-112,-80},{
              29.5,-80},{29.5,58.26},{54.46,58.26}}, color={0,0,127}));
      connect(HDirNor, HVisible.HDirNor) annotation (Line(points={{-112,-52},{
              27.5,-52},{27.5,61.78},{54.46,61.78}}, color={0,0,127}));
      connect(alt, HWindow.alt) annotation (Line(points={{-112,24},{4,24},{4,-35.38},
              {52.39,-35.38}}, color={0,0,127}));
      for i in 1:n loop
          connect(nTot, hDif_toClearCovered[i].nTot) annotation (Line(points={{-112,54},
              {-82,54},{-82,-8},{-75,-8}}, color={0,0,127}));
      connect(HDifHor, hDif_toClearCovered[i].HDifHor) annotation (Line(points={{-112,
              -24},{-94,-24},{-94,-14},{-75,-14}}, color={0,0,127}));
      end for;
       annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-96,96},{96,-96}},
              pattern=LinePattern.None,
              lineColor={0,0,0},
              fillColor={156,232,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-2,96},{2,-96}},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}),
            Rectangle(
              extent={{-96,2},{96,-2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-34,116},{40,82}},
              lineColor={28,108,200},
              fillColor={170,255,255},
              fillPattern=FillPattern.Solid,
              textString="%name
")}),                                                                Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
        Documentation(revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>",     info="<html>
<p>This model calculates the input of heat and visible light into the room due to solar irradiation. Therefore it uses the calculations of VDI 6007 part 3.  It considers  the correction values for non-vertical and non-parallel radiation incidence.</p>

  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>"));
    end Window;

    model ShadedWindow
      "Calculation of solar energy transmitted through windows considering shadowing."
      parameter Integer n(min = 1) "number of windows"
        annotation(dialog(group="window"));
      parameter Modelica.SIunits.CoefficientOfHeatTransfer UWin
        "Thermal transmission coefficient of whole window"
        annotation(dialog(group="window"));
       parameter Modelica.SIunits.TransmissionCoefficient g[n]
        "Total energy transmittance of windows"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
        "degree of light transmission for direct irradiation"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDir[n]
        "degree of light transmission for direct irradiation, with sunscreen"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient T_LTotDif[n]
        "degree of light transmission for diffuse irradiation, with sunscreen"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.RadiantEnergyFluenceRate lim
        "Limit for the sunscreen to become active"
        annotation(dialog(group="sunscreen"));
      parameter Modelica.SIunits.Angle xi(  displayUnit="degree")= 0
        "elevation angle";
      parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
        "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length b[n] "width of window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Height h[n] "height of window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length bLef[n] "window projection left"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length bRig[n] "window projection right"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length dLef[n]
        "distance between projection (left) and window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length dRig[n]
        "distance between projection (right) and window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length bAbo[n] "window projection above"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length bBel[n] "window projection below"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length dAbo[n]
        "distance between projection (above) and window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Length dBel[n]
        "distance between projection (below) and window"
        annotation (Dialog(group="window"));
      parameter Modelica.SIunits.Angle azi[n](displayUnit="degree")
        "Surface azimuth. azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
        annotation (Dialog(group="window"));
      parameter Integer nCorPoi(min = 1) "Number of corner points"
          annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Angle[nCorPoi] alpha(displayUnit="deg") "Azimuth of corner points, sorted from north to east to south to west,
     azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
          annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Height[nCorPoi] deltaH
        "Difference between height of corner point and the window centre"
        annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.Distance[nCorPoi] s
        "horizontal distance between corner point and window centre"
        annotation(dialog(group="skyline"));
      parameter Boolean[nCorPoi-1] gap
        "corner points i and i+1 are gap between buildings: true, else: false"
        annotation(dialog(group="skyline"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDir[n]
        "Total energy transmittance of windows with closed sunscreen for direct radiation"
        annotation(Dialog(group="window"));
      parameter Modelica.SIunits.TransmissionCoefficient g_TotDif[n]
        "Total energy transmittance of windows with closed sunscreen for diffuse radiation"
        annotation(Dialog(group="window"));
        Modelica.Blocks.Interfaces.RealInput incAng[n](
        final quantity="Angle",
        final unit="rad",
        displayUnit="degree")
        "Incidence angles of the sun beam on a tilted surface"
        annotation (Placement(transformation(extent={{-124,84},{-100,108}}),
        iconTransformation(extent={{-120,88},{-100,108}})));
       Modelica.Blocks.Interfaces.RealInput solAzi(quantity="Angle", unit="rad")
        "Solar azimuth angle"
        annotation (Placement(transformation(extent={{-124,-58},{-100,-34}}),
            iconTransformation(extent={{-120,-54},{-100,-34}})));
       Modelica.Blocks.Interfaces.RealInput alt(
        quantity="Angle",
        unit="rad",
        displayUnit="deg") "Solar altitude angle"
        annotation (Placement(transformation(extent={{-124,-30},{-100,-6}}),
            iconTransformation(extent={{-120,-26},{-100,-6}})));

       Modelica.Blocks.Interfaces.RealInput nTot( min=0,max=1,
        final unit="1") "Total sky Cover"
        annotation (Placement(transformation(extent={{-124,60},{-100,84}}),
            iconTransformation(extent={{-120,64},{-100,84}})));
       Modelica.Blocks.Interfaces.RealInput HDirTil[n]( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2")
        "Direct irradition on tilted surface"
        annotation (Placement(transformation(extent={{-124,30},{-100,54}}),
            iconTransformation(extent={{-120,34},{-100,54}})));

        Modelica.Blocks.Interfaces.RealInput HDirNor( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2") "Direct normal radiation."
        annotation (Placement(transformation(extent={{-124,0},{-100,24}}),
            iconTransformation(extent={{-120,4},{-100,24}})));

       Modelica.Blocks.Interfaces.RealInput HDifHor( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2")
        "Horizontal diffuse solar radiation."
        annotation (Placement(transformation(extent={{-124,-110},{-100,-86}}),
            iconTransformation(extent={{-120,-106},{-100,-86}})));

        Modelica.Blocks.Interfaces.RealInput HDifTil[n]( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2")
        "Hemispherical diffuse solar irradiation on a tilted surface from the sky"
        annotation (Placement(transformation(extent={{-124,-88},{-100,-64}}),
            iconTransformation(extent={{-120,-84},{-100,-64}})));
       Modelica.Blocks.Interfaces.RealOutput HVis[n]( quantity=
        "RadiantEnergyFluenceRate", unit="W/m2")
        "solar energy entering the room in the visible area"
        annotation (Placement(transformation(extent={{100,30},{120,50}}),
            iconTransformation(extent={{100,30},{120,50}})));

       Modelica.Blocks.Interfaces.RealOutput HWin[n](
        final quantity="RadiantEnergyFluenceRate",
        final unit="W/m2")
        "Solar radiation transmitted through aggregated window"
        annotation (Placement(transformation(extent={{100,-50},{120,-30}}),
            iconTransformation(extent={{100,-50},{120,-30}})));

      BaseClasses.Window window(
        final n=n,
        final UWin=UWin,
        final g=g,
        final T_L=T_L,
        final T_LTotDir=T_LTotDir,
        final T_LTotDif=T_LTotDif,
        final lim=lim,
        final til=til,
        final g_TotDir=g_TotDir,
        final g_TotDif=g_TotDif)
        annotation (Placement(transformation(extent={{40,-2},{88,58}})));
      BaseClasses.SelfShadowing selfShadowing(
        final b=b,
        final h=h,
        final bLef=bLef,
        final bRig=bRig,
        final dLef=dLef,
        final dRig=dRig,
        final bAbo=bAbo,
        final bBel=bBel,
        final dAbo=dAbo,
        final dBel=dBel,
        final azi=azi,
        final til=til,
        final n=n) annotation (Placement(transformation(extent={{-34,-62},{-14,-42}})));
      Modelica.Blocks.Math.Product productSelfShadowing[n]
        annotation (Placement(transformation(extent={{4,-6},{18,8}})));
      SkylineShadowing skylineShadow(
        final alpha=alpha,
        final deltaH=deltaH,
        final s=s,
        final gap=gap,
        final n=nCorPoi)
        annotation (Placement(transformation(extent={{-88,-40},{-68,-20}})));
      Modelica.Blocks.Math.Product productSkylineShadowing[n]
        annotation (Placement(transformation(extent={{-30,-6},{-16,8}})));
      Modelica.Blocks.Math.BooleanToReal boolToReal
        annotation (Placement(transformation(extent={{-48,-18},{-38,-8}})));
      Modelica.Blocks.Logical.Greater greater
        annotation (Placement(transformation(extent={{-60,-18},{-50,-8}})));
    equation
      for i in 1:n loop

        connect(boolToReal.y, productSkylineShadowing[i].u2) annotation (Line(
              points={{-37.5,-13},{-37.5,-12},{-34,-12},{-34,-4},{-31.4,-4},{
                -31.4,-3.2}}, color={0,0,127}));
      end for;
      connect(solAzi, selfShadowing.solAzi) annotation (Line(points={{-112,-46},
              {-92,-46},{-92,-45},{-35,-45}},color={0,0,127}));
      connect(alt, selfShadowing.alt) annotation (Line(points={{-112,-18},{-92,
              -18},{-92,-52},{-35,-52}},
                                    color={0,0,127}));
      connect(window.HVis, HVis)
        annotation (Line(points={{90.4,40},{90.4,40},{110,40}}, color={0,0,127}));
      connect(window.HWin, HWin) annotation (Line(points={{90.4,16},{96,16},{96,
              -40},{110,-40}},
                          color={0,0,127}));
      connect(productSkylineShadowing.y, productSelfShadowing.u1) annotation (
          Line(points={{-15.3,1},{-9.75,1},{-9.75,5.2},{2.6,5.2}}, color={0,0,
              127}));
      connect(productSelfShadowing.y, window.HDirTil) annotation (Line(points={
              {18.7,1},{24,1},{24,4},{36,4},{36,4.6},{37.6,4.6}}, color={0,0,
              127}));
      connect(HDirTil, productSkylineShadowing.u1) annotation (Line(points={{
              -112,42},{-68,42},{-68,5.2},{-31.4,5.2}}, color={0,0,127}));
      connect(alt, greater.u1) annotation (Line(points={{-112,-18},{-92,-18},{
              -92,-14},{-76,-14},{-76,-13},{-61,-13}},
                          color={0,0,127}));
      connect(skylineShadow.altLim, greater.u2) annotation (Line(points={{-67,-30},
              {-66,-30},{-66,-18},{-64,-18},{-64,-17},{-61,-17}},
                                              color={0,0,127}));
      connect(greater.y, boolToReal.u)
        annotation (Line(points={{-49.5,-13},{-49,-13}}, color={255,0,255}));
      connect(incAng, window.incAng) annotation (Line(points={{-112,96},{-44,96},
              {-44,52.6},{37.6,52.6}},
                                  color={0,0,127}));
      connect(incAng, selfShadowing.incAng) annotation (Line(points={{-112,96},
              {-94,96},{-94,-59},{-35,-59}},
                                        color={0,0,127}));
      connect(HDirNor, window.HDirNor) annotation (Line(points={{-112,12},{-46,
              12},{-46,13},{37.6,13}},
                                  color={0,0,127}));
      connect(HDifHor, window.HDifHor) annotation (Line(points={{-112,-98},{32,
              -98},{32,21.4},{37.6,21.4}},
                                      color={0,0,127}));
      connect(HDifTil, window.HDifTil) annotation (Line(points={{-112,-76},{-76,
              -76},{-76,-72},{28,-72},{28,28.6},{37.6,28.6}},
                                                         color={0,0,127}));
      connect(nTot, window.nTot) annotation (Line(points={{-112,72},{-48,72},{
              -48,44.8},{37.6,44.8}},
                            color={0,0,127}));
      connect(alt, window.alt) annotation (Line(points={{-112,-18},{-92,-18},{
              -92,35.8},{37.6,35.8}},
                            color={0,0,127}));
      connect(solAzi, skylineShadow.solAzi) annotation (Line(points={{-112,-46},
              {-96,-46},{-96,-30},{-89,-30}},
                                           color={0,0,127}));

      connect(selfShadowing.x_As, productSelfShadowing.u2) annotation (Line(
            points={{-13.1,-51.9},{-9.1,-51.9},{-9.1,-12},{-6,-12},{-6,-8},{2.6,
              -8},{2.6,-3.2}}, color={0,0,127}));
      annotation (defaultComponentName="shadedWindow",Icon(coordinateSystem(preserveAspectRatio=false),
            graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-96,96},{96,-96}},
              fillColor={175,175,175},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Rectangle(
              extent={{-2,96},{2,-96}},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}),
            Rectangle(
              extent={{-96,2},{96,-2}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-34,116},{40,82}},
              lineColor={28,108,200},
              fillColor={170,255,255},
              fillPattern=FillPattern.Solid,
              textString="%name
")}),                                                                                                    Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
<p>This model calculates the input of heat and visible light into the room due to solar irradiation. This model calculates the input of heat and visible light into the room due to solar irradiation. Therefore it uses the calculations of VDI 6007 part 3.  It considers  the correction values for non-vertical and non-parallel radiation incidence.</p>
<p> Additionaly to the <a href=\"vdi6007.Window\">Window</a> model it includes the formation of shades because of the window itself and because of the surrounding skyline.  </p>

  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>",   revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end ShadedWindow;

    package Conversions
      "This package provides conversions to translate dimensions from Annex60 definition to the definition of VDI 6007."

      function to_northAzimuth
        "Conversion from azimuth based on Annex60 to north based azimuth"
        extends Modelica.SIunits.Icons.Conversion;
        import Modelica.Constants.pi;
        input Modelica.SIunits.Angle azi;
        output Modelica.SIunits.Angle alpha;
      algorithm
        alpha:=pi+azi;

        annotation (Documentation(info="<html>
<p>This Function converts the azimuth based on <a href=\"Annex60\">Annex60</a> to the north based definition.</p>
</html>",     revisions="<html>
<ul>
<li>June 07, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
      end to_northAzimuth;

      function to_surfaceTiltVDI
        "Conversion of Annex60 surface tilt to surface tilt according to VDI 6007"
        extends Modelica.SIunits.Icons.Conversion;
        import Modelica.Constants.pi;
        input Modelica.SIunits.Angle til;
        output Modelica.SIunits.Angle gamma_F;
      algorithm
        gamma_F:=pi-til;

        annotation (Documentation(info="<html>
This function converts the inclination of a surface from the <a href=\"Annex60\">Annex60</a> definition to the definition of the VDI 6007 part 3.
</html>",     revisions="<html>
<ul>
<li>June 07, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
      end to_surfaceTiltVDI;

      model to_HDirNor
        "Converts the direct irradiation onto a horizontal surface to direct irradiation on a normal surface"
        extends Modelica.Blocks.Icons.Block;

        Modelica.Blocks.Interfaces.RealInput alt(
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg") "Solar altitude angle"
          annotation (Placement(transformation(extent={{-140,-60},{-100,-20}}),
              iconTransformation(extent={{-140,-60},{-100,-20}})));
        Modelica.Blocks.Interfaces.RealInput HDirHor(quantity=
              "RadiantEnergyFluenceRate", unit="W/m2")
          "Direct normal radiation"
          annotation (Placement(transformation(extent={{-140,20},{-100,60}}),
              iconTransformation(extent={{-140,20},{-100,60}})));
        Modelica.Blocks.Interfaces.RealOutput HDirNor(quantity=
              "RadiantEnergyFluenceRate", unit="W/m2")
          "Direct normal radiation"
          annotation (Placement(transformation(extent={{100,-20},{140,20}})));

      equation
        HDirNor=HDirHor/Modelica.Math.sin(alt);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Line(points={{-78,0},{42,0}}, color={191,0,0}),
              Polygon(
                points={{94,0},{34,20},{34,-20},{94,0}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<p>This model converts the direct irradiation on a horizontal surface to the direct irradiation on a normal surface. Therefore it needs the solar altitude angle.</p>
</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
      end to_HDirNor;

      model HDif_toClearCovered
        "Splits the total diffuse irradiation in diffuse irradiation at clear and covered sky"
        extends Modelica.Blocks.Icons.Block;

        Modelica.Blocks.Interfaces.RealInput HDifHor( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Horizontal diffuse solar radiation."
           annotation (Placement(transformation(extent={{-116,-66},{-100,-50}}),
              iconTransformation(extent={{-120,-70},{-100,-50}})));

         Modelica.Blocks.Interfaces.RealInput HDifTil( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Hemispherical diffuse solar irradiation on a tilted surface from the sky"
          annotation (Placement(transformation(extent={{-116,54},{-100,70}}),
              iconTransformation(extent={{-120,50},{-100,70}})));
         Modelica.Blocks.Interfaces.RealOutput HDifHorCov( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Horizontal diffuse solar radiation at covered sky."
          annotation (Placement(transformation(extent={{104,-26},{120,-10}}),
              iconTransformation(extent={{100,-30},{120,-10}})));

        Modelica.Blocks.Interfaces.RealOutput HDifHorCle( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Horizontal diffuse solar radiation at clear sky."
          annotation (Placement(transformation(extent={{104,-66},{120,-50}}),
              iconTransformation(extent={{100,-70},{120,-50}})));

        Modelica.Blocks.Interfaces.RealOutput HDifTilCov( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Hemispherical diffuse solar irradiation on a tilted surface at covered sky"
          annotation (Placement(transformation(extent={{104,14},{120,30}}),
              iconTransformation(extent={{100,10},{120,30}})));
        Modelica.Blocks.Interfaces.RealOutput HDifTilCle( quantity=
          "RadiantEnergyFluenceRate", unit="W/m2")
          "Hemispherical diffuse solar irradiation on a tilted surface at clear sky"
          annotation (Placement(transformation(extent={{104,54},{120,70}}),
              iconTransformation(extent={{100,50},{120,70}})));
        Modelica.Blocks.Interfaces.RealInput nTot( min=0,max=1,
          final unit="1") "Total sky Cover"
          annotation (Placement(transformation(extent={{-120,-10},{-100,10}}),
              iconTransformation(extent={{-120,-10},{-100,10}})));

      equation
        HDifTilCle=HDifTil*(1-nTot);
        HDifTilCov=HDifTil*nTot;
        HDifHorCle=HDifHor*(1-nTot);
        HDifHorCov=HDifHor*nTot;
        annotation (defaultComponentName="HDif_toCleCov",
          Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Line(points={{-78,0},{42,0}}, color={191,0,0}),
              Polygon(
                points={{94,0},{34,20},{34,-20},{94,0}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid)}),                      Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<p>This model calculates the diffuse irradiation at clear and covered sky out of the total diffuse irradiation. Therefore it uses the total sky cover.</p>
</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
      end HDif_toClearCovered;
      annotation (Documentation(info="<html>
<p>This package contains the conversions needed to use the equations of VDI6007 with inputs and parameters with <code>Annex60</code> definition.</p>
</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
    end Conversions;
    annotation (Documentation(info="<html>
<p>The BaseClasses package provides the basic models for the window calculations.<\\p>
</html>", revisions="<html>
<ul>
<li>July 1 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
  end BaseClasses;

  package Examples "Examples for the windows submodels"
    extends Modelica.Icons.ExamplesPackage;

    model Illumination "Testmodel for Illumination"
      import vdi6007 = Annex60.ThermalZones.ReducedOrder.Windows;
      extends Modelica.Icons.Example;

      Windows.BaseClasses.HeatIllumination heatIllumination(HIll1=120, HIll2=240)
        "Heat input into the room due to the illumination"
        annotation (Placement(transformation(extent={{76,-10},{96,10}})));
      Windows.BaseClasses.Illumination illumination(
        D=0.27,
        e_ILim1=250,
        e_ILim2=500,
        office=true,
        n=1,
        r={0.21},
        A={5.13},
        T_L={0.72},
        til={1.5707963267949})
        "Determining the switch times for the illumination in the room"
        annotation (Placement(transformation(extent={{52,-10},{72,10}})));
      Windows.SolarGain.CorrectionGTaueDoublePane CorGTaue(
        n=1,
        UWin=1.4,
        til={90}) "Correction values for non-vertical non-parallel irradiation"
                  annotation (Placement(transformation(extent={{18,-42},{38,-22}})));
      Windows.BaseClasses.Sunblind sunblind(lim=200)
        "Calculates if the sunblind of the window is active"
        annotation (Placement(transformation(extent={{-16,-64},{-6,-54}})));
      BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil(
        azi=0,
        til=1.5707963267949,
        lat=0.86393797973719) "Diffuse irradiation on the window"
        annotation (Placement(transformation(extent={{-60,-18},{-40,2}})));
      Windows.Window window(
        n=1,
        UWin=1.4,
        g={0.64},
        g_TotDir={0.08},
        g_TotDif={0.27},
        T_L={0.72},
        T_LTotDir={0.08},
        T_LTotDif={0.32},
        lim=200,
        azi={0},
        lat=0.86393797973719,
        til={1.5707963267949}) "Window facing the south in a wall"
        annotation (Placement(transformation(extent={{-10,24},{10,44}})));
      BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
        filNam="modelica://Annex60/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
        annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil(
        azi=0,
        til=1.5707963267949,
        lat=0.86393797973719) "Direct irradiation on the surface"
        annotation (Placement(transformation(extent={{-66,-64},{-46,-44}})));
    equation
      connect(illumination.Illumination, heatIllumination.Illumination)
        annotation (Line(points={{73,0},{75,0}}, color={255,0,255}));
      connect(CorGTaue.CorTaue_Gro, illumination.CorTaue_Gro) annotation (Line(
            points={{37,-24},{40,-24},{40,0},{51,0}}, color={0,0,127}));
      connect(CorGTaue.CorTaue_DifCov, illumination.CorTaue_DifCov) annotation (
         Line(points={{37,-26},{44,-26},{44,-6},{51,-6}}, color={0,0,127}));
      connect(HDifTil.H, sunblind.HDifTil) annotation (Line(points={{-39,-8},{-32,
              -8},{-32,-42},{-32,-44},{-32,-56},{-24,-56},{-20,-56},{-20,-56.1},{-16.5,
              -56.1}},                         color={0,0,127}));
      connect(sunblind.sunscreen, CorGTaue.sunscreen[1]) annotation (Line(
            points={{-5.5,-59},{-5.5,-58},{-6,-58},{-4,-58},{16,-58},{16,-34},{
              19,-34}}, color={255,0,255}));
      connect(weaDat.weaBus, HDifTil.weaBus) annotation (Line(
          points={{-80,0},{-70,0},{-70,-8},{-60,-8}},
          color={255,204,51},
          thickness=0.5));
      connect(weaDat.weaBus, window.weaBus) annotation (Line(
          points={{-80,0},{-64,0},{-64,34},{-9.8,34}},
          color={255,204,51},
          thickness=0.5));
      connect(window.HVis, illumination.HVis)
        annotation (Line(points={{11,38},{32,38},{32,6},{51,6}}, color={0,0,127}));
      connect(weaDat.weaBus, HDirTil.weaBus) annotation (Line(
          points={{-80,0},{-74,0},{-74,-54},{-66,-54}},
          color={255,204,51},
          thickness=0.5));
      connect(HDirTil.H, sunblind.HDirTil) annotation (Line(points={{-45,-54},{
              -36,-54},{-36,-62},{-16.5,-62}}, color={0,0,127}));
      connect(HDirTil.inc, CorGTaue.incAng[1]) annotation (Line(points={{-45,-58},{-40,
              -58},{-40,-30},{19,-30}}, color={0,0,127}));
      annotation (experiment(StartTime=0,StopTime=31536000),Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}})),
        Documentation(info="<html>
  <p>This example shows the application of
  <a href=\"Windows.BaseClasses.Illumination\">Illumination</a>.
   For solar radiation, the example relies on the standard
  weather file in Annex60.</p>
  <p>The idea of the example is to show a typical application of all
  sub-models and to use the example in unit tests. The results are
  reasonable, but not related to any real use case or measurement
  data.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>",     revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
    end Illumination;

    model Window "Testmodel for Window"
      import vdi6007 = Annex60.ThermalZones.ReducedOrder.Windows;
        extends Modelica.Icons.Example;

      BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
        filNam="modelica://Annex60/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
        annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));
      Windows.Window window(
        n=2,
        UWin=1.4,
        g={0.64,0.64},
        g_TotDir={0.08,0.08},
        g_TotDif={0.27,0.27},
        T_L={0.72,0.72},
        T_LTotDir={0.08,0.08},
        T_LTotDif={0.32,0.32},
        lim=200,
        lat=0.86393797973719,
        til={1.5707963267949,1.5707963267949},
        azi={0,1.5707963267949})
        "Two windows: One facing south, one facing west"
        annotation (Placement(transformation(extent={{12,-10},{32,10}})));
    equation
      connect(weaDat.weaBus, window.weaBus) annotation (Line(
          points={{-22,0},{12.2,0}},
          color={255,204,51},
          thickness=0.5));
      annotation (experiment(StartTime=0,StopTime=31536000),Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
  <p>This example shows the application of
  <a href=\"vdi6007.Window\">Window</a>.
   For solar radiation, the example relies on the standard
  weather file in Annex60.</p>
  <p>The idea of the example is to show a typical application of all
  sub-models and to use the example in unit tests. The results are
  reasonable, but not related to any real use case or measurement
  data.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>", revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
    end Window;

    model ShadedWindow "Testmodel for ShadedWindow"
      import Annex60.ThermalZones.ReducedOrder.Windows;
      extends Modelica.Icons.Example;

      BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
        filNam="modelica://Annex60/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
        annotation (Placement(transformation(extent={{-42,-10},{-22,10}})));
      Windows.ShadedWindow shadedWindow(
        n=2,
        UWin=1.4,
        g={0.64,0.64},
        g_TotDir={0.08,0.08},
        g_TotDif={0.27,0.27},
        T_L={0.72,0.72},
        T_LTotDir={0.08,0.08},
        T_LTotDif={0.32,0.32},
        lim=200,
        nCorPoi=4,
        b={2,2},
        h={2,2},
        bLef={0.5,0.5},
        bRig={0.5,0.5},
        dLef={0.1,0.1},
        dRig={0.1,0.1},
        bAbo={1,1},
        bBel={0,0},
        dAbo={0.5,0.5},
        dBel={0,0},
        azi(displayUnit="deg") = {0,1.5707963267949},
        deltaH={5,5,10,10},
        s={100,100,20,20},
        gap={false,true,false},
        lat=0.86393797973719,
        til={1.5707963267949,1.5707963267949},
        alpha={-0.34906585039887,0.34906585039887,1.7453292519943,
            2.0943951023932})
        "Two shaded windows: One facing south, one facing west."
        annotation (Placement(transformation(extent={{20,-12},{48,14}})));
    equation
      connect(weaDat.weaBus, shadedWindow.weaBus) annotation (Line(
          points={{-22,0},{-2,0},{-2,1},{20,1}},
          color={255,204,51},
          thickness=0.5));
      annotation (experiment(StartTime=0,StopTime=31536000),Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
  <p>This example shows the application of
  <a href=\"Windows.ShadedWindow\">ShadedWindow</a>.
   For solar radiation, the example relies on the standard
  weather file in Annex60.</p>
  <p>The idea of the example is to show a typical application of all
  sub-models and to use the example in unit tests. The results are
  reasonable, but not related to any real use case or measurement
  data.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>", revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
    end ShadedWindow;

    model VentilationHeat "Testmodel for VentilationHeat"
      extends Modelica.Icons.Example;

      BaseClasses.VentilationHeat ventilationHeat(
        d=0.1,
        screen=true,
        x_f=0.8,
        tau_e=0.1,
        rho_e=0.7125,
        til=1.5707963267949)
        "Heat input due to ventilation with closed sunblind"
        annotation (Placement(transformation(extent={{56,-10},{76,10}})));
      BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
        filNam="modelica://Annex60/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")
        annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil(
        azi=0,
        til=1.5707963267949,
        lat=0.86393797973719) "Diffuse irradiation on the window"
               annotation (Placement(transformation(extent={{-42,20},{-26,36}})));
      BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil(
        azi=0,
        til=1.5707963267949,
        lat=0.86393797973719) "Direct solar irradiation on the window"
               annotation (Placement(transformation(extent={{-42,-8},{-26,8}})));
      BaseClasses.Sunblind sunblind(lim=200) "Sunblind of the window"
        annotation (Placement(transformation(extent={{-16,8},{-2,22}})));
      BoundaryConditions.SolarGeometry.BaseClasses.AltitudeAngle altAng
        "Solar altitude angle"
        annotation (Placement(transformation(extent={{10,-20},{20,-10}})));
      BoundaryConditions.SolarGeometry.ZenithAngle zen(lat=0.86393797973719)
        "Solar zenith angle"
        annotation (Placement(transformation(extent={{-2,-20},{8,-10}})));
      BoundaryConditions.WeatherData.Bus weaBus "Weather bus"
                                                annotation (Placement(
            transformation(extent={{-70,38},{-38,68}}), iconTransformation(extent={{
                -196,50},{-176,70}})));
      ThermalZones.ReducedOrder.Windows.Window window(
        n=1,
        UWin=1.4,
        g={0.64},
        g_TotDir={0.08},
        g_TotDif={0.27},
        T_L={0.72},
        T_LTotDir={0.08},
        T_LTotDif={0.32},
        lim=200,
        azi={0},
        lat=0.86393797973719,
        til={1.5707963267949}) "Window facing south"
        annotation (Placement(transformation(extent={{56,-54},{76,-34}})));
      Modelica.Blocks.Math.Add add
        "Total solar energy entering the room through the window"
        annotation (Placement(transformation(extent={{92,-4},{100,4}})));
    equation
      connect(altAng.zen, zen.y)
        annotation (Line(points={{9,-15},{9,-15},{8.5,-15}}, color={0,0,127}));
      connect(altAng.alt, ventilationHeat.alt) annotation (Line(points={{20.5,-15},{
              37.25,-15},{37.25,-8},{55,-8}}, color={0,0,127}));
      connect(HDifTil.H, ventilationHeat.HDifTil) annotation (Line(points={{-25.2,28},
              {14,28},{14,5},{55,5}}, color={0,0,127}));
      connect(weaDat.weaBus, HDifTil.weaBus) annotation (Line(
          points={{-80,0},{-62,0},{-62,28},{-42,28}},
          color={255,204,51},
          thickness=0.5));
      connect(HDirTil.H, ventilationHeat.HDirTil) annotation (Line(points={{-25.2,0},
              {16,0},{16,-2},{55,-2}}, color={0,0,127}));
      connect(HDifTil.H, sunblind.HDifTil) annotation (Line(points={{-25.2,28},
              {-22,28},{-22,19.06},{-16.7,19.06}}, color={0,0,127}));
      connect(HDirTil.H, sunblind.HDirTil) annotation (Line(points={{-25.2,0},{
              -22,0},{-22,10.8},{-16.7,10.8}}, color={0,0,127}));
      connect(sunblind.sunscreen, ventilationHeat.sunscreen) annotation (Line(
            points={{-1.3,15},{4,15},{4,16},{10,16},{10,2},{55,2}}, color={255,
              0,255}));
      connect(weaDat.weaBus, weaBus) annotation (Line(
          points={{-80,0},{-68,0},{-68,53},{-54,53}},
          color={255,204,51},
          thickness=0.5), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      connect(weaBus.HDifHor, ventilationHeat.HDifHor) annotation (Line(
          points={{-54,53},{32,53},{32,8},{55,8}},
          color={255,204,51},
          thickness=0.5), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(weaBus.HDirNor, ventilationHeat.HDirNor) annotation (Line(
          points={{-54,53},{32,53},{32,-5},{55,-5}},
          color={255,204,51},
          thickness=0.5), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(weaDat.weaBus, HDirTil.weaBus) annotation (Line(
          points={{-80,0},{-62,0},{-42,0}},
          color={255,204,51},
          thickness=0.5));
      connect(weaDat.weaBus, zen.weaBus) annotation (Line(
          points={{-80,0},{-56,0},{-56,-15},{-2,-15}},
          color={255,204,51},
          thickness=0.5));
      connect(weaDat.weaBus, window.weaBus) annotation (Line(
          points={{-80,0},{-68,0},{-68,-44},{56.2,-44}},
          color={255,204,51},
          thickness=0.5));
      connect(ventilationHeat.HVen, add.u1) annotation (Line(points={{77.2,0},{86,0},
              {86,2.4},{91.2,2.4}}, color={0,0,127}));
      connect(window.HWin[1], add.u2) annotation (Line(points={{77,-48},{88,-48},{88,
              -2.4},{91.2,-2.4}}, color={0,0,127}));
      annotation (experiment(StartTime=0,StopTime=31536000),Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        Documentation(info="<html>
  <p>This example shows the application of
  <a href=\"Windows.BaseClasses.VentilationHeat\">VentilationHeat</a>.
   For solar radiation, the example relies on the standard
  weather file in Annex60.</p>
  <p>The idea of the example is to show a typical application of all
  sub-models and to use the example in unit tests. The results are
  reasonable, but not related to any real use case or measurement
  data.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>", revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
    end VentilationHeat;
  end Examples;

  package Validation "Collecion of validation models"
    extends Modelica.Icons.ExamplesPackage;

    model SkylineShadowingTest
      import Annex60;
      extends Modelica.Icons.Example;

      Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.SkylineShadowing
        skylineShadow(
        n=4,
        gap={false,true,false},
        deltaH={5,5,100,100},
        s={20,20,20,20},
        alpha={1.3962634015955,1.7453292519943,-1.3962634015955,-1.7453292519943})
        "Shadow due to buildings on the west and east side"
        annotation (Placement(transformation(extent={{28,-10},{48,10}})));

      Modelica.Blocks.Sources.Sine solAziSine(freqHz=1, amplitude=Modelica.Constants.pi)
        "Solar azimuth input generated as sine"
        annotation (Placement(transformation(extent={{-48,-10},{-28,10}})));
    equation
      connect(solAziSine.y, skylineShadow.solAzi)
        annotation (Line(points={{-27,0},{-27,0},{27,0}},color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics),
        Documentation(revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>", info="<html>
This is an example for the <a href=\"Windows.BaseClasses.SkylineShadowing\">SkylineShadowing</a> model. It simulates two buildings which shade the window. One smaller building is on the east side and one bigger building on the west side. Between the buildings there is a gap.<\\p>
</html>"));
    end SkylineShadowingTest;

    package BaseClasses
        extends Modelica.Icons.BasesPackage;

      block IncidenceAngleVDI6007
        "The solar incidence angle on a tilted surface"
        extends Modelica.Blocks.Icons.Block;
        import
          Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
        import
          Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_northAzimuth;
        parameter Modelica.SIunits.Angle azi(displayUnit="degree")
          "Surface azimuth. azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south";
        parameter Modelica.SIunits.Angle til(displayUnit="degree")
          "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof";
        Modelica.Blocks.Interfaces.RealInput solAzi(quantity="Angle", unit="rad")
          "Solar azimuth angle"
          annotation (Placement(transformation(extent={{-140,-68},{-100,-28}})));
        Modelica.Blocks.Interfaces.RealInput alt(quantity="Angle", unit="rad")
          "solar altitude angle"
          annotation (Placement(transformation(extent={{-142,34},{-102,74}})));
        Modelica.Blocks.Interfaces.RealOutput incAng(
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg") "Incidence angle on a tilted surface"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));

      equation
        incAng = max(0, Modelica.Math.acos(Modelica.Math.cos(to_surfaceTiltVDI(til))*Modelica.Math.sin(alt) +
        Modelica.Math.sin(to_surfaceTiltVDI(til))*Modelica.Math.cos(alt)*Modelica.Math.cos((abs(to_northAzimuth(azi)-to_northAzimuth(solAzi))))));

        annotation (
          defaultComponentName="incAng",
          Documentation(info="<html>
<p>
This component computes the solar incidence angle on a tilted surface using the solar hour angle and the declination angle as input.
</p>
</html>",       revisions="<html>
<ul>
<li>
Dec 7, 2010, by Michael Wetter:<br/>
Rewrote equation in explicit form to avoid nonlinear equations in room model.
</li>
<li>
May 19, 2010, by Wangda Zuo:<br/>
First implementation.
</li>
</ul>
</html>"),Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={
              Text(
                extent={{-150,110},{150,150}},
                textString="%name",
                lineColor={0,0,255}),
              Bitmap(extent={{-90,92},{90,-94}}, fileName=
                    "modelica://Annex60/Resources/Images/BoundaryConditions/SolarGeometry/BaseClasses/IncidenceAngle.png")}));
      end IncidenceAngleVDI6007;

      model SolarDeclinationAngleVDI6007 "Calculates the solar azimuth angle based on the equations of VDI 6007 part 3.
  It is modelled to test other models based on VDI2078. "
        extends Modelica.Blocks.Icons.Block;
        import Modelica.SIunits.Conversions.from_deg;

        Modelica.Blocks.Interfaces.RealOutput decAng(
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg") "Solar declination angle"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));
      protected
          constant Modelica.SIunits.Time day=86400 "Number of seconds in a day";
          Modelica.SIunits.Angle J_;
      equation
        J_=from_deg(360*105/365);
        decAng=from_deg(0.3948-23.2559*Modelica.Math.cos(J_+from_deg(9.1))-0.3915*Modelica.Math.cos(2*J_+from_deg(5.4))-0.1764*Modelica.Math.cos(3*J_+from_deg(26)));

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
This model computes the solar declination angle for test case 1 and 3 of the VDI2078 in April.
</html>"));
      end SolarDeclinationAngleVDI6007;

      model SolarHourAngleVDI6007 "Calculates the solar hour angle every hour based on the equations of VDI 6007 part 3.
  It is modelled to test other Models based on VDI2078. It doesn't consider summer time"

        extends Modelica.Blocks.Icons.Block;
        import Modelica.SIunits.Conversions.from_deg;
        parameter Modelica.SIunits.Angle lon "Longitude";
        import Modelica.SIunits.Conversions.to_deg;
        Modelica.SIunits.Angle J_;
        Modelica.SIunits.Time zgl "time equation";
        Modelica.SIunits.Time woz "true time";
        Modelica.Blocks.Interfaces.RealOutput solHouAng(
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg") "Solar hour angle"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));
      protected
        constant Modelica.SIunits.Time day=86400 "Number of seconds in a day";
        Real clock;
      equation
        clock=(integer(time/3600)-0.5-integer(time/day)*24);
        J_=from_deg(360*105/365);
        zgl=0.0066+7.3525*Modelica.Math.cos(J_+from_deg(85.9))+9.9359*Modelica.Math.cos(2*J_+from_deg(108.9))+0.3387*Modelica.Math.cos(3*J_+from_deg(105.2));
        woz=(integer(time/3600)-0.5-integer(time/day)*24)-4*(15-to_deg(lon))/60+zgl/60;
        solHouAng=(12-woz)*2*Modelica.Constants.pi/24;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
This model computes the solar hour angle for test case 1 and 3 of VDI2078 in April.
</html>"));
      end SolarHourAngleVDI6007;
      annotation (Documentation(info="<html>
<p> This package includes BaseClasses that are only used for validation causes. </p>
</html>"));
    end BaseClasses;

    package VDI2078
        extends Modelica.Icons.ExamplesPackage;
      model TestCase1_Illumination
          extends Modelica.Icons.Example;

        Modelica.Blocks.Sources.CombiTimeTable HDirHor(
          columns={2},
          tableName="HDirHor",
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              6.97804373; 28800,45.26143833; 32400,110.4153032; 36000,
              186.2726158; 39600,255.5241448; 43200,302.9030765; 46800,
              318.0561637; 50400,297.6769067; 54000,246.2124597; 57600,
              174.9154955; 61200,99.54837583; 64800,37.47098217; 68400,
              4.241288132; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,6.97804373;
              115200,45.26143833; 118800,110.4153032; 122400,186.2726158;
              126000,255.5241448; 129600,302.9030765; 133200,318.0561637;
              136800,297.6769067; 140400,246.2124597; 144000,174.9154955;
              147600,99.54837583; 151200,37.47098217; 154800,4.241288132;
              158400,0; 162000,0; 165600,0; 169200,0; 172800,0; 176400,0;
              180000,0; 183600,0; 187200,0; 190800,0; 194400,0; 198000,
              6.97804373; 201600,45.26143833; 205200,110.4153032; 208800,
              186.2726158; 212400,255.5241448; 216000,302.9030765; 219600,
              318.0561637; 223200,297.6769067; 226800,246.2124597; 230400,
              174.9154955; 234000,99.54837583; 237600,37.47098217; 241200,
              4.241288132; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,6.97804373; 288000,45.26143833; 291600,110.4153032; 295200,
              186.2726158; 298800,255.5241448; 302400,302.9030765; 306000,
              318.0561637; 309600,297.6769067; 313200,246.2124597; 316800,
              174.9154955; 320400,99.54837583; 324000,37.47098217; 327600,
              4.241288132; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,6.97804373; 374400,45.26143833; 378000,110.4153032; 381600,
              186.2726158; 385200,255.5241448; 388800,302.9030765; 392400,
              318.0561637; 396000,297.6769067; 399600,246.2124597; 403200,
              174.9154955; 406800,99.54837583; 410400,37.47098217; 414000,
              4.241288132; 417600,0; 421200,0; 424800,0; 428400,0; 432000,0;
              435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,6.97804373; 460800,45.26143833; 464400,110.4153032; 468000,
              186.2726158; 471600,255.5241448; 475200,302.9030765; 478800,
              318.0561637; 482400,297.6769067; 486000,246.2124597; 489600,
              174.9154955; 493200,99.54837583; 496800,37.47098217; 500400,
              4.241288132; 504000,0; 507600,0; 511200,0; 514800,0; 518400,0;
              522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,6.97804373; 547200,45.26143833; 550800,110.4153032; 554400,
              186.2726158; 558000,255.5241448; 561600,302.9030765; 565200,
              318.0561637; 568800,297.6769067; 572400,246.2124597; 576000,
              174.9154955; 579600,99.54837583; 583200,37.47098217; 586800,
              4.241288132; 590400,0; 594000,0; 597600,0; 601200,0; 604800,0;
              608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,6.97804373; 633600,45.26143833; 637200,110.4153032; 640800,
              186.2726158; 644400,255.5241448; 648000,302.9030765; 651600,
              318.0561637; 655200,297.6769067; 658800,246.2124597; 662400,
              174.9154955; 666000,99.54837583; 669600,37.47098217; 673200,
              4.241288132; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,6.97804373; 720000,45.26143833; 723600,110.4153032; 727200,
              186.2726158; 730800,255.5241448; 734400,302.9030765; 738000,
              318.0561637; 741600,297.6769067; 745200,246.2124597; 748800,
              174.9154955; 752400,99.54837583; 756000,37.47098217; 759600,
              4.241288132; 763200,0; 766800,0; 770400,0; 774000,0; 777600,0;
              781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,6.97804373; 806400,45.26143833; 810000,110.4153032; 813600,
              186.2726158; 817200,255.5241448; 820800,302.9030765; 824400,
              318.0561637; 828000,297.6769067; 831600,246.2124597; 835200,
              174.9154955; 838800,99.54837583; 842400,37.47098217; 846000,
              4.241288132; 849600,0; 853200,0; 856800,0; 860400,0; 864000,0;
              867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,6.97804373; 892800,45.26143833; 896400,110.4153032; 900000,
              186.2726158; 903600,255.5241448; 907200,302.9030765; 910800,
              318.0561637; 914400,297.6769067; 918000,246.2124597; 921600,
              174.9154955; 925200,99.54837583; 928800,37.47098217; 932400,
              4.241288132; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,6.97804373; 979200,45.26143833; 982800,110.4153032; 986400,
              186.2726158; 990000,255.5241448; 993600,302.9030765; 997200,
              318.0561637; 1000800,297.6769067; 1004400,246.2124597; 1008000,
              174.9154955; 1011600,99.54837583; 1015200,37.47098217; 1018800,
              4.241288132; 1022400,0; 1026000,0; 1029600,0; 1033200,0; 1036800,
              0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0; 1058400,
              0; 1062000,6.97804373; 1065600,45.26143833; 1069200,110.4153032;
              1072800,186.2726158; 1076400,255.5241448; 1080000,302.9030765;
              1083600,318.0561637; 1087200,297.6769067; 1090800,246.2124597;
              1094400,174.9154955; 1098000,99.54837583; 1101600,37.47098217;
              1105200,4.241288132; 1108800,0; 1112400,0; 1116000,0; 1119600,0;
              1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0; 1141200,0;
              1144800,0; 1148400,6.97804373; 1152000,45.26143833; 1155600,
              110.4153032; 1159200,186.2726158; 1162800,255.5241448; 1166400,
              302.9030765; 1170000,318.0561637; 1173600,297.6769067; 1177200,
              246.2124597; 1180800,174.9154955; 1184400,99.54837583; 1188000,
              37.47098217; 1191600,4.241288132; 1195200,0; 1198800,0; 1202400,0;
              1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0; 1224000,0;
              1227600,0; 1231200,0; 1234800,39.07398846; 1238400,164.3396228;
              1242000,319.1618766; 1245600,466.1278073; 1249200,583.2867887;
              1252800,657.0058834; 1256400,679.6830794; 1260000,649.0893148;
              1263600,568.2428437; 1267200,445.5861318; 1270800,295.7855335;
              1274400,142.4655924; 1278000,26.47655754; 1281600,0; 1285200,0;
              1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0; 1306800,0;
              1310400,0; 1314000,0; 1317600,0; 1321200,39.07398846; 1324800,
              164.3396228; 1328400,319.1618766; 1332000,466.1278073; 1335600,
              583.2867887; 1339200,657.0058834; 1342800,679.6830794; 1346400,
              649.0893148; 1350000,568.2428437; 1353600,445.5861318; 1357200,
              295.7855335; 1360800,142.4655924; 1364400,26.47655754; 1368000,0;
              1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0; 1389600,0;
              1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,39.07398846;
              1411200,164.3396228; 1414800,319.1618766; 1418400,466.1278073;
              1422000,583.2867887; 1425600,657.0058834; 1429200,679.6830794;
              1432800,649.0893148; 1436400,568.2428437; 1440000,445.5861318;
              1443600,295.7855335; 1447200,142.4655924; 1450800,26.47655754;
              1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,0; 1472400,0;
              1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,0; 1494000,
              39.07398846; 1497600,164.3396228; 1501200,319.1618766; 1504800,
              466.1278073; 1508400,583.2867887; 1512000,657.0058834; 1515600,
              679.6830794; 1519200,649.0893148; 1522800,568.2428437; 1526400,
              445.5861318; 1530000,295.7855335; 1533600,142.4655924; 1537200,
              26.47655754; 1540800,0; 1544400,0; 1548000,0; 1551600,0; 1555200,
              0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0; 1576800,
              0; 1580400,39.07398846; 1584000,164.3396228; 1587600,319.1618766;
              1591200,466.1278073; 1594800,583.2867887; 1598400,657.0058834;
              1602000,679.6830794; 1605600,649.0893148; 1609200,568.2428437;
              1612800,445.5861318; 1616400,295.7855335; 1620000,142.4655924;
              1623600,26.47655754; 1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Direct horizontal irradiation"
          annotation (Placement(transformation(extent={{-102,94},{-82,114}})));

        Modelica.Blocks.Sources.CombiTimeTable alt(
          columns={2},
          tableName="alt",
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableOnFile=false,
          table=[0,-0.524971501; 3600,-0.540504674; 7200,-0.511513997; 10800,-0.436196642;
              14400,-0.323752929; 18000,-0.184550738; 21600,-0.027807643; 25200,
              0.138810415; 28800,0.308550932; 32400,0.474473156; 36000,
              0.627965928; 39600,0.756893419; 43200,0.844457705; 46800,
              0.872876752; 50400,0.834722167; 54000,0.739755433; 57600,
              0.606187916; 61200,0.450168115; 64800,0.28318588; 68400,
              0.113490388; 72000,-0.052089885; 75600,-0.206727187; 79200,-0.342547809;
              82800,-0.450098665; 86400,-0.524971501; 90000,-0.540504674; 93600,
              -0.511513997; 97200,-0.436196642; 100800,-0.323752929; 104400,-0.184550738;
              108000,-0.027807643; 111600,0.138810415; 115200,0.308550932;
              118800,0.474473156; 122400,0.627965928; 126000,0.756893419;
              129600,0.844457705; 133200,0.872876752; 136800,0.834722167;
              140400,0.739755433; 144000,0.606187916; 147600,0.450168115;
              151200,0.28318588; 154800,0.113490388; 158400,-0.052089885;
              162000,-0.206727187; 165600,-0.342547809; 169200,-0.450098665;
              172800,-0.524971501; 176400,-0.540504674; 180000,-0.511513997;
              183600,-0.436196642; 187200,-0.323752929; 190800,-0.184550738;
              194400,-0.027807643; 198000,0.138810415; 201600,0.308550932;
              205200,0.474473156; 208800,0.627965928; 212400,0.756893419;
              216000,0.844457705; 219600,0.872876752; 223200,0.834722167;
              226800,0.739755433; 230400,0.606187916; 234000,0.450168115;
              237600,0.28318588; 241200,0.113490388; 244800,-0.052089885;
              248400,-0.206727187; 252000,-0.342547809; 255600,-0.450098665;
              259200,-0.524971501; 262800,-0.540504674; 266400,-0.511513997;
              270000,-0.436196642; 273600,-0.323752929; 277200,-0.184550738;
              280800,-0.027807643; 284400,0.138810415; 288000,0.308550932;
              291600,0.474473156; 295200,0.627965928; 298800,0.756893419;
              302400,0.844457705; 306000,0.872876752; 309600,0.834722167;
              313200,0.739755433; 316800,0.606187916; 320400,0.450168115;
              324000,0.28318588; 327600,0.113490388; 331200,-0.052089885;
              334800,-0.206727187; 338400,-0.342547809; 342000,-0.450098665;
              345600,-0.524971501; 349200,-0.540504674; 352800,-0.511513997;
              356400,-0.436196642; 360000,-0.323752929; 363600,-0.184550738;
              367200,-0.027807643; 370800,0.138810415; 374400,0.308550932;
              378000,0.474473156; 381600,0.627965928; 385200,0.756893419;
              388800,0.844457705; 392400,0.872876752; 396000,0.834722167;
              399600,0.739755433; 403200,0.606187916; 406800,0.450168115;
              410400,0.28318588; 414000,0.113490388; 417600,-0.052089885;
              421200,-0.206727187; 424800,-0.342547809; 428400,-0.450098665;
              432000,-0.524971501; 435600,-0.540504674; 439200,-0.511513997;
              442800,-0.436196642; 446400,-0.323752929; 450000,-0.184550738;
              453600,-0.027807643; 457200,0.138810415; 460800,0.308550932;
              464400,0.474473156; 468000,0.627965928; 471600,0.756893419;
              475200,0.844457705; 478800,0.872876752; 482400,0.834722167;
              486000,0.739755433; 489600,0.606187916; 493200,0.450168115;
              496800,0.28318588; 500400,0.113490388; 504000,-0.052089885;
              507600,-0.206727187; 511200,-0.342547809; 514800,-0.450098665;
              518400,-0.524971501; 522000,-0.540504674; 525600,-0.511513997;
              529200,-0.436196642; 532800,-0.323752929; 536400,-0.184550738;
              540000,-0.027807643; 543600,0.138810415; 547200,0.308550932;
              550800,0.474473156; 554400,0.627965928; 558000,0.756893419;
              561600,0.844457705; 565200,0.872876752; 568800,0.834722167;
              572400,0.739755433; 576000,0.606187916; 579600,0.450168115;
              583200,0.28318588; 586800,0.113490388; 590400,-0.052089885;
              594000,-0.206727187; 597600,-0.342547809; 601200,-0.450098665;
              604800,-0.524971501; 608400,-0.540504674; 612000,-0.511513997;
              615600,-0.436196642; 619200,-0.323752929; 622800,-0.184550738;
              626400,-0.027807643; 630000,0.138810415; 633600,0.308550932;
              637200,0.474473156; 640800,0.627965928; 644400,0.756893419;
              648000,0.844457705; 651600,0.872876752; 655200,0.834722167;
              658800,0.739755433; 662400,0.606187916; 666000,0.450168115;
              669600,0.28318588; 673200,0.113490388; 676800,-0.052089885;
              680400,-0.206727187; 684000,-0.342547809; 687600,-0.450098665;
              691200,-0.524971501; 694800,-0.540504674; 698400,-0.511513997;
              702000,-0.436196642; 705600,-0.323752929; 709200,-0.184550738;
              712800,-0.027807643; 716400,0.138810415; 720000,0.308550932;
              723600,0.474473156; 727200,0.627965928; 730800,0.756893419;
              734400,0.844457705; 738000,0.872876752; 741600,0.834722167;
              745200,0.739755433; 748800,0.606187916; 752400,0.450168115;
              756000,0.28318588; 759600,0.113490388; 763200,-0.052089885;
              766800,-0.206727187; 770400,-0.342547809; 774000,-0.450098665;
              777600,-0.524971501; 781200,-0.540504674; 784800,-0.511513997;
              788400,-0.436196642; 792000,-0.323752929; 795600,-0.184550738;
              799200,-0.027807643; 802800,0.138810415; 806400,0.308550932;
              810000,0.474473156; 813600,0.627965928; 817200,0.756893419;
              820800,0.844457705; 824400,0.872876752; 828000,0.834722167;
              831600,0.739755433; 835200,0.606187916; 838800,0.450168115;
              842400,0.28318588; 846000,0.113490388; 849600,-0.052089885;
              853200,-0.206727187; 856800,-0.342547809; 860400,-0.450098665;
              864000,-0.524971501; 867600,-0.540504674; 871200,-0.511513997;
              874800,-0.436196642; 878400,-0.323752929; 882000,-0.184550738;
              885600,-0.027807643; 889200,0.138810415; 892800,0.308550932;
              896400,0.474473156; 900000,0.627965928; 903600,0.756893419;
              907200,0.844457705; 910800,0.872876752; 914400,0.834722167;
              918000,0.739755433; 921600,0.606187916; 925200,0.450168115;
              928800,0.28318588; 932400,0.113490388; 936000,-0.052089885;
              939600,-0.206727187; 943200,-0.342547809; 946800,-0.450098665;
              950400,-0.524971501; 954000,-0.540504674; 957600,-0.511513997;
              961200,-0.436196642; 964800,-0.323752929; 968400,-0.184550738;
              972000,-0.027807643; 975600,0.138810415; 979200,0.308550932;
              982800,0.474473156; 986400,0.627965928; 990000,0.756893419;
              993600,0.844457705; 997200,0.872876752; 1000800,0.834722167;
              1004400,0.739755433; 1008000,0.606187916; 1011600,0.450168115;
              1015200,0.28318588; 1018800,0.113490388; 1022400,-0.052089885;
              1026000,-0.206727187; 1029600,-0.342547809; 1033200,-0.450098665;
              1036800,-0.524971501; 1040400,-0.540504674; 1044000,-0.511513997;
              1047600,-0.436196642; 1051200,-0.323752929; 1054800,-0.184550738;
              1058400,-0.027807643; 1062000,0.138810415; 1065600,0.308550932;
              1069200,0.474473156; 1072800,0.627965928; 1076400,0.756893419;
              1080000,0.844457705; 1083600,0.872876752; 1087200,0.834722167;
              1090800,0.739755433; 1094400,0.606187916; 1098000,0.450168115;
              1101600,0.28318588; 1105200,0.113490388; 1108800,-0.052089885;
              1112400,-0.206727187; 1116000,-0.342547809; 1119600,-0.450098665;
              1123200,-0.524971501; 1126800,-0.540504674; 1130400,-0.511513997;
              1134000,-0.436196642; 1137600,-0.323752929; 1141200,-0.184550738;
              1144800,-0.027807643; 1148400,0.138810415; 1152000,0.308550932;
              1155600,0.474473156; 1159200,0.627965928; 1162800,0.756893419;
              1166400,0.844457705; 1170000,0.872876752; 1173600,0.834722167;
              1177200,0.739755433; 1180800,0.606187916; 1184400,0.450168115;
              1188000,0.28318588; 1191600,0.113490388; 1195200,-0.052089885;
              1198800,-0.206727187; 1202400,-0.342547809; 1206000,-0.450098665;
              1209600,-0.524971501; 1213200,-0.540504674; 1216800,-0.511513997;
              1220400,-0.436196642; 1224000,-0.323752929; 1227600,-0.184550738;
              1231200,-0.027807643; 1234800,0.138810415; 1238400,0.308550932;
              1242000,0.474473156; 1245600,0.627965928; 1249200,0.756893419;
              1252800,0.844457705; 1256400,0.872876752; 1260000,0.834722167;
              1263600,0.739755433; 1267200,0.606187916; 1270800,0.450168115;
              1274400,0.28318588; 1278000,0.113490388; 1281600,-0.052089885;
              1285200,-0.206727187; 1288800,-0.342547809; 1292400,-0.450098665;
              1296000,-0.524971501; 1299600,-0.540504674; 1303200,-0.511513997;
              1306800,-0.436196642; 1310400,-0.323752929; 1314000,-0.184550738;
              1317600,-0.027807643; 1321200,0.138810415; 1324800,0.308550932;
              1328400,0.474473156; 1332000,0.627965928; 1335600,0.756893419;
              1339200,0.844457705; 1342800,0.872876752; 1346400,0.834722167;
              1350000,0.739755433; 1353600,0.606187916; 1357200,0.450168115;
              1360800,0.28318588; 1364400,0.113490388; 1368000,-0.052089885;
              1371600,-0.206727187; 1375200,-0.342547809; 1378800,-0.450098665;
              1382400,-0.524971501; 1386000,-0.540504674; 1389600,-0.511513997;
              1393200,-0.436196642; 1396800,-0.323752929; 1400400,-0.184550738;
              1404000,-0.027807643; 1407600,0.138810415; 1411200,0.308550932;
              1414800,0.474473156; 1418400,0.627965928; 1422000,0.756893419;
              1425600,0.844457705; 1429200,0.872876752; 1432800,0.834722167;
              1436400,0.739755433; 1440000,0.606187916; 1443600,0.450168115;
              1447200,0.28318588; 1450800,0.113490388; 1454400,-0.052089885;
              1458000,-0.206727187; 1461600,-0.342547809; 1465200,-0.450098665;
              1468800,-0.524971501; 1472400,-0.540504674; 1476000,-0.511513997;
              1479600,-0.436196642; 1483200,-0.323752929; 1486800,-0.184550738;
              1490400,-0.027807643; 1494000,0.138810415; 1497600,0.308550932;
              1501200,0.474473156; 1504800,0.627965928; 1508400,0.756893419;
              1512000,0.844457705; 1515600,0.872876752; 1519200,0.834722167;
              1522800,0.739755433; 1526400,0.606187916; 1530000,0.450168115;
              1533600,0.28318588; 1537200,0.113490388; 1540800,-0.052089885;
              1544400,-0.206727187; 1548000,-0.342547809; 1551600,-0.450098665;
              1555200,-0.524971501; 1558800,-0.540504674; 1562400,-0.511513997;
              1566000,-0.436196642; 1569600,-0.323752929; 1573200,-0.184550738;
              1576800,-0.027807643; 1580400,0.138810415; 1584000,0.308550932;
              1587600,0.474473156; 1591200,0.627965928; 1594800,0.756893419;
              1598400,0.844457705; 1602000,0.872876752; 1605600,0.834722167;
              1609200,0.739755433; 1612800,0.606187916; 1616400,0.450168115;
              1620000,0.28318588; 1623600,0.113490388; 1627200,-0.052089885;
              1630800,-0.206727187; 1634400,-0.342547809; 1638000,-0.450098665])
          "solar altitude angle"
          annotation (Placement(transformation(extent={{-102,-70},{-82,-50}})));

        SolarGain.CorrectionGTaueDoublePane CorGTaue(
          n=1,
          UWin=1.4,
          xi=0,
          til(displayUnit="deg") = {1.5707963267949})
          "Correction values for non-parallel and non-vertical irradiation for VDI2078 test case 1"
          annotation (Placement(transformation(extent={{-52,46},{-32,66}})));
        Windows.BaseClasses.Sunblind sunblind(lim=200)
          "Calculates if the sunblind is active"
          annotation (Placement(transformation(extent={{-20,-46},{-6,-36}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifHorCle(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifHorCle",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              21.04231965; 28800,38.76402419; 32400,49.79122826; 36000,
              56.90519903; 39600,61.37224652; 43200,63.79161511; 46800,
              64.48257797; 50400,63.54466077; 54000,60.84410045; 57600,
              56.02824899; 61200,48.43452029; 64800,36.62912553; 68400,
              17.58508841; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              21.04231965; 115200,38.76402419; 118800,49.79122826; 122400,
              56.90519903; 126000,61.37224652; 129600,63.79161511; 133200,
              64.48257797; 136800,63.54466077; 140400,60.84410045; 144000,
              56.02824899; 147600,48.43452029; 151200,36.62912553; 154800,
              17.58508841; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,21.04231965; 201600,38.76402419; 205200,49.79122826;
              208800,56.90519903; 212400,61.37224652; 216000,63.79161511;
              219600,64.48257797; 223200,63.54466077; 226800,60.84410045;
              230400,56.02824899; 234000,48.43452029; 237600,36.62912553;
              241200,17.58508841; 244800,0; 248400,0; 252000,0; 255600,0;
              259200,0; 262800,0; 266400,0; 270000,0; 273600,0; 277200,0;
              280800,0; 284400,21.04231965; 288000,38.76402419; 291600,
              49.79122826; 295200,56.90519903; 298800,61.37224652; 302400,
              63.79161511; 306000,64.48257797; 309600,63.54466077; 313200,
              60.84410045; 316800,56.02824899; 320400,48.43452029; 324000,
              36.62912553; 327600,17.58508841; 331200,0; 334800,0; 338400,0;
              342000,0; 345600,0; 349200,0; 352800,0; 356400,0; 360000,0;
              363600,0; 367200,0; 370800,21.04231965; 374400,38.76402419;
              378000,49.79122826; 381600,56.90519903; 385200,61.37224652;
              388800,63.79161511; 392400,64.48257797; 396000,63.54466077;
              399600,60.84410045; 403200,56.02824899; 406800,48.43452029;
              410400,36.62912553; 414000,17.58508841; 417600,0; 421200,0;
              424800,0; 428400,0; 432000,0; 435600,0; 439200,0; 442800,0;
              446400,0; 450000,0; 453600,0; 457200,21.04231965; 460800,
              38.76402419; 464400,49.79122826; 468000,56.90519903; 471600,
              61.37224652; 475200,63.79161511; 478800,64.48257797; 482400,
              63.54466077; 486000,60.84410045; 489600,56.02824899; 493200,
              48.43452029; 496800,36.62912553; 500400,17.58508841; 504000,0;
              507600,0; 511200,0; 514800,0; 518400,0; 522000,0; 525600,0;
              529200,0; 532800,0; 536400,0; 540000,0; 543600,21.04231965;
              547200,38.76402419; 550800,49.79122826; 554400,56.90519903;
              558000,61.37224652; 561600,63.79161511; 565200,64.48257797;
              568800,63.54466077; 572400,60.84410045; 576000,56.02824899;
              579600,48.43452029; 583200,36.62912553; 586800,17.58508841;
              590400,0; 594000,0; 597600,0; 601200,0; 604800,0; 608400,0;
              612000,0; 615600,0; 619200,0; 622800,0; 626400,0; 630000,
              21.04231965; 633600,38.76402419; 637200,49.79122826; 640800,
              56.90519903; 644400,61.37224652; 648000,63.79161511; 651600,
              64.48257797; 655200,63.54466077; 658800,60.84410045; 662400,
              56.02824899; 666000,48.43452029; 669600,36.62912553; 673200,
              17.58508841; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,21.04231965; 720000,38.76402419; 723600,49.79122826;
              727200,56.90519903; 730800,61.37224652; 734400,63.79161511;
              738000,64.48257797; 741600,63.54466077; 745200,60.84410045;
              748800,56.02824899; 752400,48.43452029; 756000,36.62912553;
              759600,17.58508841; 763200,0; 766800,0; 770400,0; 774000,0;
              777600,0; 781200,0; 784800,0; 788400,0; 792000,0; 795600,0;
              799200,0; 802800,21.04231965; 806400,38.76402419; 810000,
              49.79122826; 813600,56.90519903; 817200,61.37224652; 820800,
              63.79161511; 824400,64.48257797; 828000,63.54466077; 831600,
              60.84410045; 835200,56.02824899; 838800,48.43452029; 842400,
              36.62912553; 846000,17.58508841; 849600,0; 853200,0; 856800,0;
              860400,0; 864000,0; 867600,0; 871200,0; 874800,0; 878400,0;
              882000,0; 885600,0; 889200,21.04231965; 892800,38.76402419;
              896400,49.79122826; 900000,56.90519903; 903600,61.37224652;
              907200,63.79161511; 910800,64.48257797; 914400,63.54466077;
              918000,60.84410045; 921600,56.02824899; 925200,48.43452029;
              928800,36.62912553; 932400,17.58508841; 936000,0; 939600,0;
              943200,0; 946800,0; 950400,0; 954000,0; 957600,0; 961200,0;
              964800,0; 968400,0; 972000,0; 975600,21.04231965; 979200,
              38.76402419; 982800,49.79122826; 986400,56.90519903; 990000,
              61.37224652; 993600,63.79161511; 997200,64.48257797; 1000800,
              63.54466077; 1004400,60.84410045; 1008000,56.02824899; 1011600,
              48.43452029; 1015200,36.62912553; 1018800,17.58508841; 1022400,0;
              1026000,0; 1029600,0; 1033200,0; 1036800,0; 1040400,0; 1044000,0;
              1047600,0; 1051200,0; 1054800,0; 1058400,0; 1062000,21.04231965;
              1065600,38.76402419; 1069200,49.79122826; 1072800,56.90519903;
              1076400,61.37224652; 1080000,63.79161511; 1083600,64.48257797;
              1087200,63.54466077; 1090800,60.84410045; 1094400,56.02824899;
              1098000,48.43452029; 1101600,36.62912553; 1105200,17.58508841;
              1108800,0; 1112400,0; 1116000,0; 1119600,0; 1123200,0; 1126800,0;
              1130400,0; 1134000,0; 1137600,0; 1141200,0; 1144800,0; 1148400,
              21.04231965; 1152000,38.76402419; 1155600,49.79122826; 1159200,
              56.90519903; 1162800,61.37224652; 1166400,63.79161511; 1170000,
              64.48257797; 1173600,63.54466077; 1177200,60.84410045; 1180800,
              56.02824899; 1184400,48.43452029; 1188000,36.62912553; 1191600,
              17.58508841; 1195200,0; 1198800,0; 1202400,0; 1206000,0; 1209600,
              0; 1213200,0; 1216800,0; 1220400,0; 1224000,0; 1227600,0; 1231200,
              0; 1234800,40.47459946; 1238400,67.00103172; 1242000,82.68760638;
              1245600,93.44518302; 1249200,100.7065715; 1252800,104.8421528;
              1256400,106.0514354; 1260000,104.4130318; 1263600,99.82357411;
              1267200,92.07110863; 1270800,80.71803866; 1274400,63.97542723;
              1278000,34.69852159; 1281600,0; 1285200,0; 1288800,0; 1292400,0;
              1296000,0; 1299600,0; 1303200,0; 1306800,0; 1310400,0; 1314000,0;
              1317600,0; 1321200,40.47459946; 1324800,67.00103172; 1328400,
              82.68760638; 1332000,93.44518302; 1335600,100.7065715; 1339200,
              104.8421528; 1342800,106.0514354; 1346400,104.4130318; 1350000,
              99.82357411; 1353600,92.07110863; 1357200,80.71803866; 1360800,
              63.97542723; 1364400,34.69852159; 1368000,0; 1371600,0; 1375200,0;
              1378800,0; 1382400,0; 1386000,0; 1389600,0; 1393200,0; 1396800,0;
              1400400,0; 1404000,0; 1407600,40.47459946; 1411200,67.00103172;
              1414800,82.68760638; 1418400,93.44518302; 1422000,100.7065715;
              1425600,104.8421528; 1429200,106.0514354; 1432800,104.4130318;
              1436400,99.82357411; 1440000,92.07110863; 1443600,80.71803866;
              1447200,63.97542723; 1450800,34.69852159; 1454400,0; 1458000,0;
              1461600,0; 1465200,0; 1468800,0; 1472400,0; 1476000,0; 1479600,0;
              1483200,0; 1486800,0; 1490400,0; 1494000,40.47459946; 1497600,
              67.00103172; 1501200,82.68760638; 1504800,93.44518302; 1508400,
              100.7065715; 1512000,104.8421528; 1515600,106.0514354; 1519200,
              104.4130318; 1522800,99.82357411; 1526400,92.07110863; 1530000,
              80.71803866; 1533600,63.97542723; 1537200,34.69852159; 1540800,0;
              1544400,0; 1548000,0; 1551600,0; 1555200,0; 1558800,0; 1562400,0;
              1566000,0; 1569600,0; 1573200,0; 1576800,0; 1580400,40.47459946;
              1584000,67.00103172; 1587600,82.68760638; 1591200,93.44518302;
              1594800,100.7065715; 1598400,104.8421528; 1602000,106.0514354;
              1605600,104.4130318; 1609200,99.82357411; 1612800,92.07110863;
              1616400,80.71803866; 1620000,63.97542723; 1623600,34.69852159;
              1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Diffuse irradiation at clear sky on horizontal surface"
          annotation (Placement(transformation(extent={{-102,38},{-82,58}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifHorCov(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifHorCov",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              24.22614885; 28800,69.2322226; 32400,121.6772143; 36000,
              171.3048325; 39600,210.9338328; 43200,235.8864705; 46800,
              243.5618298; 50400,233.20683; 54000,205.8425596; 57600,
              164.3628129; 61200,113.7867437; 64800,61.72452797; 68400,
              18.86655303; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              24.22614885; 115200,69.2322226; 118800,121.6772143; 122400,
              171.3048325; 126000,210.9338328; 129600,235.8864705; 133200,
              243.5618298; 136800,233.20683; 140400,205.8425596; 144000,
              164.3628129; 147600,113.7867437; 151200,61.72452797; 154800,
              18.86655303; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,24.22614885; 201600,69.2322226; 205200,121.6772143; 208800,
              171.3048325; 212400,210.9338328; 216000,235.8864705; 219600,
              243.5618298; 223200,233.20683; 226800,205.8425596; 230400,
              164.3628129; 234000,113.7867437; 237600,61.72452797; 241200,
              18.86655303; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,24.22614885; 288000,69.2322226; 291600,121.6772143; 295200,
              171.3048325; 298800,210.9338328; 302400,235.8864705; 306000,
              243.5618298; 309600,233.20683; 313200,205.8425596; 316800,
              164.3628129; 320400,113.7867437; 324000,61.72452797; 327600,
              18.86655303; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,24.22614885; 374400,69.2322226; 378000,121.6772143; 381600,
              171.3048325; 385200,210.9338328; 388800,235.8864705; 392400,
              243.5618298; 396000,233.20683; 399600,205.8425596; 403200,
              164.3628129; 406800,113.7867437; 410400,61.72452797; 414000,
              18.86655303; 417600,0; 421200,0; 424800,0; 428400,0; 432000,0;
              435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,24.22614885; 460800,69.2322226; 464400,121.6772143; 468000,
              171.3048325; 471600,210.9338328; 475200,235.8864705; 478800,
              243.5618298; 482400,233.20683; 486000,205.8425596; 489600,
              164.3628129; 493200,113.7867437; 496800,61.72452797; 500400,
              18.86655303; 504000,0; 507600,0; 511200,0; 514800,0; 518400,0;
              522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,24.22614885; 547200,69.2322226; 550800,121.6772143; 554400,
              171.3048325; 558000,210.9338328; 561600,235.8864705; 565200,
              243.5618298; 568800,233.20683; 572400,205.8425596; 576000,
              164.3628129; 579600,113.7867437; 583200,61.72452797; 586800,
              18.86655303; 590400,0; 594000,0; 597600,0; 601200,0; 604800,0;
              608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,24.22614885; 633600,69.2322226; 637200,121.6772143; 640800,
              171.3048325; 644400,210.9338328; 648000,235.8864705; 651600,
              243.5618298; 655200,233.20683; 658800,205.8425596; 662400,
              164.3628129; 666000,113.7867437; 669600,61.72452797; 673200,
              18.86655303; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,24.22614885; 720000,69.2322226; 723600,121.6772143; 727200,
              171.3048325; 730800,210.9338328; 734400,235.8864705; 738000,
              243.5618298; 741600,233.20683; 745200,205.8425596; 748800,
              164.3628129; 752400,113.7867437; 756000,61.72452797; 759600,
              18.86655303; 763200,0; 766800,0; 770400,0; 774000,0; 777600,0;
              781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,24.22614885; 806400,69.2322226; 810000,121.6772143; 813600,
              171.3048325; 817200,210.9338328; 820800,235.8864705; 824400,
              243.5618298; 828000,233.20683; 831600,205.8425596; 835200,
              164.3628129; 838800,113.7867437; 842400,61.72452797; 846000,
              18.86655303; 849600,0; 853200,0; 856800,0; 860400,0; 864000,0;
              867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,24.22614885; 892800,69.2322226; 896400,121.6772143; 900000,
              171.3048325; 903600,210.9338328; 907200,235.8864705; 910800,
              243.5618298; 914400,233.20683; 918000,205.8425596; 921600,
              164.3628129; 925200,113.7867437; 928800,61.72452797; 932400,
              18.86655303; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,24.22614885; 979200,69.2322226; 982800,121.6772143; 986400,
              171.3048325; 990000,210.9338328; 993600,235.8864705; 997200,
              243.5618298; 1000800,233.20683; 1004400,205.8425596; 1008000,
              164.3628129; 1011600,113.7867437; 1015200,61.72452797; 1018800,
              18.86655303; 1022400,0; 1026000,0; 1029600,0; 1033200,0; 1036800,
              0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0; 1058400,
              0; 1062000,24.22614885; 1065600,69.2322226; 1069200,121.6772143;
              1072800,171.3048325; 1076400,210.9338328; 1080000,235.8864705;
              1083600,243.5618298; 1087200,233.20683; 1090800,205.8425596;
              1094400,164.3628129; 1098000,113.7867437; 1101600,61.72452797;
              1105200,18.86655303; 1108800,0; 1112400,0; 1116000,0; 1119600,0;
              1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0; 1141200,0;
              1144800,0; 1148400,24.22614885; 1152000,69.2322226; 1155600,
              121.6772143; 1159200,171.3048325; 1162800,210.9338328; 1166400,
              235.8864705; 1170000,243.5618298; 1173600,233.20683; 1177200,
              205.8425596; 1180800,164.3628129; 1184400,113.7867437; 1188000,
              61.72452797; 1191600,18.86655303; 1195200,0; 1198800,0; 1202400,0;
              1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0; 1224000,0;
              1227600,0; 1231200,0; 1234800,5.860931292; 1238400,16.65168637;
              1242000,28.66366397; 1245600,39.74586572; 1249200,48.47858361;
              1252800,53.93982313; 1256400,55.6149166; 1260000,53.35451152;
              1263600,47.36107656; 1267200,38.20692416; 1270800,26.88090607;
              1274400,14.891867; 1278000,4.534899719; 1281600,0; 1285200,0;
              1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0; 1306800,0;
              1310400,0; 1314000,0; 1317600,0; 1321200,5.860931292; 1324800,
              16.65168637; 1328400,28.66366397; 1332000,39.74586572; 1335600,
              48.47858361; 1339200,53.93982313; 1342800,55.6149166; 1346400,
              53.35451152; 1350000,47.36107656; 1353600,38.20692416; 1357200,
              26.88090607; 1360800,14.891867; 1364400,4.534899719; 1368000,0;
              1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0; 1389600,0;
              1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,5.860931292;
              1411200,16.65168637; 1414800,28.66366397; 1418400,39.74586572;
              1422000,48.47858361; 1425600,53.93982313; 1429200,55.6149166;
              1432800,53.35451152; 1436400,47.36107656; 1440000,38.20692416;
              1443600,26.88090607; 1447200,14.891867; 1450800,4.534899719;
              1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,0; 1472400,0;
              1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,0; 1494000,
              5.860931292; 1497600,16.65168637; 1501200,28.66366397; 1504800,
              39.74586572; 1508400,48.47858361; 1512000,53.93982313; 1515600,
              55.6149166; 1519200,53.35451152; 1522800,47.36107656; 1526400,
              38.20692416; 1530000,26.88090607; 1533600,14.891867; 1537200,
              4.534899719; 1540800,0; 1544400,0; 1548000,0; 1551600,0; 1555200,
              0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0; 1576800,
              0; 1580400,5.860931292; 1584000,16.65168637; 1587600,28.66366397;
              1591200,39.74586572; 1594800,48.47858361; 1598400,53.93982313;
              1602000,55.6149166; 1605600,53.35451152; 1609200,47.36107656;
              1612800,38.20692416; 1616400,26.88090607; 1620000,14.891867;
              1623600,4.534899719; 1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Diffuse irradiation at covered sky on horizontal surface"
          annotation (Placement(transformation(extent={{-102,10},{-82,30}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifTilCle(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifTilCle",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              12.09644867; 28800,27.25095027; 32400,38.27385402; 36000,45.82451;
              39600,51.57958573; 43200,54.91750562; 46800,55.89716878; 50400,
              54.57011009; 54000,50.87097992; 57600,44.75882584; 61200,
              37.07302936; 64800,25.04319434; 68400,9.799046246; 72000,0; 75600,
              0; 79200,0; 82800,0; 86400,0; 90000,0; 93600,0; 97200,0; 100800,0;
              104400,0; 108000,0; 111600,12.09644867; 115200,27.25095027;
              118800,38.27385402; 122400,45.82451; 126000,51.57958573; 129600,
              54.91750562; 133200,55.89716878; 136800,54.57011009; 140400,
              50.87097992; 144000,44.75882584; 147600,37.07302936; 151200,
              25.04319434; 154800,9.799046246; 158400,0; 162000,0; 165600,0;
              169200,0; 172800,0; 176400,0; 180000,0; 183600,0; 187200,0;
              190800,0; 194400,0; 198000,12.09644867; 201600,27.25095027;
              205200,38.27385402; 208800,45.82451; 212400,51.57958573; 216000,
              54.91750562; 219600,55.89716878; 223200,54.57011009; 226800,
              50.87097992; 230400,44.75882584; 234000,37.07302936; 237600,
              25.04319434; 241200,9.799046246; 244800,0; 248400,0; 252000,0;
              255600,0; 259200,0; 262800,0; 266400,0; 270000,0; 273600,0;
              277200,0; 280800,0; 284400,12.09644867; 288000,27.25095027;
              291600,38.27385402; 295200,45.82451; 298800,51.57958573; 302400,
              54.91750562; 306000,55.89716878; 309600,54.57011009; 313200,
              50.87097992; 316800,44.75882584; 320400,37.07302936; 324000,
              25.04319434; 327600,9.799046246; 331200,0; 334800,0; 338400,0;
              342000,0; 345600,0; 349200,0; 352800,0; 356400,0; 360000,0;
              363600,0; 367200,0; 370800,12.09644867; 374400,27.25095027;
              378000,38.27385402; 381600,45.82451; 385200,51.57958573; 388800,
              54.91750562; 392400,55.89716878; 396000,54.57011009; 399600,
              50.87097992; 403200,44.75882584; 406800,37.07302936; 410400,
              25.04319434; 414000,9.799046246; 417600,0; 421200,0; 424800,0;
              428400,0; 432000,0; 435600,0; 439200,0; 442800,0; 446400,0;
              450000,0; 453600,0; 457200,12.09644867; 460800,27.25095027;
              464400,38.27385402; 468000,45.82451; 471600,51.57958573; 475200,
              54.91750562; 478800,55.89716878; 482400,54.57011009; 486000,
              50.87097992; 489600,44.75882584; 493200,37.07302936; 496800,
              25.04319434; 500400,9.799046246; 504000,0; 507600,0; 511200,0;
              514800,0; 518400,0; 522000,0; 525600,0; 529200,0; 532800,0;
              536400,0; 540000,0; 543600,12.09644867; 547200,27.25095027;
              550800,38.27385402; 554400,45.82451; 558000,51.57958573; 561600,
              54.91750562; 565200,55.89716878; 568800,54.57011009; 572400,
              50.87097992; 576000,44.75882584; 579600,37.07302936; 583200,
              25.04319434; 586800,9.799046246; 590400,0; 594000,0; 597600,0;
              601200,0; 604800,0; 608400,0; 612000,0; 615600,0; 619200,0;
              622800,0; 626400,0; 630000,12.09644867; 633600,27.25095027;
              637200,38.27385402; 640800,45.82451; 644400,51.57958573; 648000,
              54.91750562; 651600,55.89716878; 655200,54.57011009; 658800,
              50.87097992; 662400,44.75882584; 666000,37.07302936; 669600,
              25.04319434; 673200,9.799046246; 676800,0; 680400,0; 684000,0;
              687600,0; 691200,0; 694800,0; 698400,0; 702000,0; 705600,0;
              709200,0; 712800,0; 716400,12.09644867; 720000,27.25095027;
              723600,38.27385402; 727200,45.82451; 730800,51.57958573; 734400,
              54.91750562; 738000,55.89716878; 741600,54.57011009; 745200,
              50.87097992; 748800,44.75882584; 752400,37.07302936; 756000,
              25.04319434; 759600,9.799046246; 763200,0; 766800,0; 770400,0;
              774000,0; 777600,0; 781200,0; 784800,0; 788400,0; 792000,0;
              795600,0; 799200,0; 802800,12.09644867; 806400,27.25095027;
              810000,38.27385402; 813600,45.82451; 817200,51.57958573; 820800,
              54.91750562; 824400,55.89716878; 828000,54.57011009; 831600,
              50.87097992; 835200,44.75882584; 838800,37.07302936; 842400,
              25.04319434; 846000,9.799046246; 849600,0; 853200,0; 856800,0;
              860400,0; 864000,0; 867600,0; 871200,0; 874800,0; 878400,0;
              882000,0; 885600,0; 889200,12.09644867; 892800,27.25095027;
              896400,38.27385402; 900000,45.82451; 903600,51.57958573; 907200,
              54.91750562; 910800,55.89716878; 914400,54.57011009; 918000,
              50.87097992; 921600,44.75882584; 925200,37.07302936; 928800,
              25.04319434; 932400,9.799046246; 936000,0; 939600,0; 943200,0;
              946800,0; 950400,0; 954000,0; 957600,0; 961200,0; 964800,0;
              968400,0; 972000,0; 975600,12.09644867; 979200,27.25095027;
              982800,38.27385402; 986400,45.82451; 990000,51.57958573; 993600,
              54.91750562; 997200,55.89716878; 1000800,54.57011009; 1004400,
              50.87097992; 1008000,44.75882584; 1011600,37.07302936; 1015200,
              25.04319434; 1018800,9.799046246; 1022400,0; 1026000,0; 1029600,0;
              1033200,0; 1036800,0; 1040400,0; 1044000,0; 1047600,0; 1051200,0;
              1054800,0; 1058400,0; 1062000,12.09644867; 1065600,27.25095027;
              1069200,38.27385402; 1072800,45.82451; 1076400,51.57958573;
              1080000,54.91750562; 1083600,55.89716878; 1087200,54.57011009;
              1090800,50.87097992; 1094400,44.75882584; 1098000,37.07302936;
              1101600,25.04319434; 1105200,9.799046246; 1108800,0; 1112400,0;
              1116000,0; 1119600,0; 1123200,0; 1126800,0; 1130400,0; 1134000,0;
              1137600,0; 1141200,0; 1144800,0; 1148400,12.09644867; 1152000,
              27.25095027; 1155600,38.27385402; 1159200,45.82451; 1162800,
              51.57958573; 1166400,54.91750562; 1170000,55.89716878; 1173600,
              54.57011009; 1177200,50.87097992; 1180800,44.75882584; 1184400,
              37.07302936; 1188000,25.04319434; 1191600,9.799046246; 1195200,0;
              1198800,0; 1202400,0; 1206000,0; 1209600,0; 1213200,0; 1216800,0;
              1220400,0; 1224000,0; 1227600,0; 1231200,0; 1234800,23.26734519;
              1238400,47.10145092; 1242000,63.56086174; 1245600,75.24935853;
              1249200,84.63765841; 1252800,90.25746571; 1256400,91.93142037;
              1260000,89.6665522; 1263600,83.46122296; 1267200,73.55208829;
              1270800,61.78366585; 1274400,43.73975719; 1278000,19.33526916;
              1281600,0; 1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0;
              1303200,0; 1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,
              23.26734519; 1324800,47.10145092; 1328400,63.56086174; 1332000,
              75.24935853; 1335600,84.63765841; 1339200,90.25746571; 1342800,
              91.93142037; 1346400,89.6665522; 1350000,83.46122296; 1353600,
              73.55208829; 1357200,61.78366585; 1360800,43.73975719; 1364400,
              19.33526916; 1368000,0; 1371600,0; 1375200,0; 1378800,0; 1382400,
              0; 1386000,0; 1389600,0; 1393200,0; 1396800,0; 1400400,0; 1404000,
              0; 1407600,23.26734519; 1411200,47.10145092; 1414800,63.56086174;
              1418400,75.24935853; 1422000,84.63765841; 1425600,90.25746571;
              1429200,91.93142037; 1432800,89.6665522; 1436400,83.46122296;
              1440000,73.55208829; 1443600,61.78366585; 1447200,43.73975719;
              1450800,19.33526916; 1454400,0; 1458000,0; 1461600,0; 1465200,0;
              1468800,0; 1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0;
              1490400,0; 1494000,23.26734519; 1497600,47.10145092; 1501200,
              63.56086174; 1504800,75.24935853; 1508400,84.63765841; 1512000,
              90.25746571; 1515600,91.93142037; 1519200,89.6665522; 1522800,
              83.46122296; 1526400,73.55208829; 1530000,61.78366585; 1533600,
              43.73975719; 1537200,19.33526916; 1540800,0; 1544400,0; 1548000,0;
              1551600,0; 1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0;
              1573200,0; 1576800,0; 1580400,23.26734519; 1584000,47.10145092;
              1587600,63.56086174; 1591200,75.24935853; 1594800,84.63765841;
              1598400,90.25746571; 1602000,91.93142037; 1605600,89.6665522;
              1609200,83.46122296; 1612800,73.55208829; 1616400,61.78366585;
              1620000,43.73975719; 1623600,19.33526916; 1627200,0; 1630800,0;
              1634400,0; 1638000,0])
          "Diffuse irradiation at clear sky on the tilted window"
          annotation (Placement(transformation(extent={{-102,-18},{-82,2}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifTilCov(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifTilCov",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              9.603148499; 28800,27.44337611; 32400,48.23236103; 36000,
              67.90455038; 39600,83.61332758; 43200,93.50445335; 46800,
              96.54693507; 50400,92.4422546; 54000,81.59516725; 57600,
              65.15276159; 61200,45.10461005; 64800,24.46735599; 68400,
              7.478626154; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              9.603148499; 115200,27.44337611; 118800,48.23236103; 122400,
              67.90455038; 126000,83.61332758; 129600,93.50445335; 133200,
              96.54693507; 136800,92.4422546; 140400,81.59516725; 144000,
              65.15276159; 147600,45.10461005; 151200,24.46735599; 154800,
              7.478626154; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,9.603148499; 201600,27.44337611; 205200,48.23236103;
              208800,67.90455038; 212400,83.61332758; 216000,93.50445335;
              219600,96.54693507; 223200,92.4422546; 226800,81.59516725; 230400,
              65.15276159; 234000,45.10461005; 237600,24.46735599; 241200,
              7.478626154; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,9.603148499; 288000,27.44337611; 291600,48.23236103;
              295200,67.90455038; 298800,83.61332758; 302400,93.50445335;
              306000,96.54693507; 309600,92.4422546; 313200,81.59516725; 316800,
              65.15276159; 320400,45.10461005; 324000,24.46735599; 327600,
              7.478626154; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,9.603148499; 374400,27.44337611; 378000,48.23236103;
              381600,67.90455038; 385200,83.61332758; 388800,93.50445335;
              392400,96.54693507; 396000,92.4422546; 399600,81.59516725; 403200,
              65.15276159; 406800,45.10461005; 410400,24.46735599; 414000,
              7.478626154; 417600,0; 421200,0; 424800,0; 428400,0; 432000,0;
              435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,9.603148499; 460800,27.44337611; 464400,48.23236103;
              468000,67.90455038; 471600,83.61332758; 475200,93.50445335;
              478800,96.54693507; 482400,92.4422546; 486000,81.59516725; 489600,
              65.15276159; 493200,45.10461005; 496800,24.46735599; 500400,
              7.478626154; 504000,0; 507600,0; 511200,0; 514800,0; 518400,0;
              522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,9.603148499; 547200,27.44337611; 550800,48.23236103;
              554400,67.90455038; 558000,83.61332758; 561600,93.50445335;
              565200,96.54693507; 568800,92.4422546; 572400,81.59516725; 576000,
              65.15276159; 579600,45.10461005; 583200,24.46735599; 586800,
              7.478626154; 590400,0; 594000,0; 597600,0; 601200,0; 604800,0;
              608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,9.603148499; 633600,27.44337611; 637200,48.23236103;
              640800,67.90455038; 644400,83.61332758; 648000,93.50445335;
              651600,96.54693507; 655200,92.4422546; 658800,81.59516725; 662400,
              65.15276159; 666000,45.10461005; 669600,24.46735599; 673200,
              7.478626154; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,9.603148499; 720000,27.44337611; 723600,48.23236103;
              727200,67.90455038; 730800,83.61332758; 734400,93.50445335;
              738000,96.54693507; 741600,92.4422546; 745200,81.59516725; 748800,
              65.15276159; 752400,45.10461005; 756000,24.46735599; 759600,
              7.478626154; 763200,0; 766800,0; 770400,0; 774000,0; 777600,0;
              781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,9.603148499; 806400,27.44337611; 810000,48.23236103;
              813600,67.90455038; 817200,83.61332758; 820800,93.50445335;
              824400,96.54693507; 828000,92.4422546; 831600,81.59516725; 835200,
              65.15276159; 838800,45.10461005; 842400,24.46735599; 846000,
              7.478626154; 849600,0; 853200,0; 856800,0; 860400,0; 864000,0;
              867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,9.603148499; 892800,27.44337611; 896400,48.23236103;
              900000,67.90455038; 903600,83.61332758; 907200,93.50445335;
              910800,96.54693507; 914400,92.4422546; 918000,81.59516725; 921600,
              65.15276159; 925200,45.10461005; 928800,24.46735599; 932400,
              7.478626154; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,9.603148499; 979200,27.44337611; 982800,48.23236103;
              986400,67.90455038; 990000,83.61332758; 993600,93.50445335;
              997200,96.54693507; 1000800,92.4422546; 1004400,81.59516725;
              1008000,65.15276159; 1011600,45.10461005; 1015200,24.46735599;
              1018800,7.478626154; 1022400,0; 1026000,0; 1029600,0; 1033200,0;
              1036800,0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0;
              1058400,0; 1062000,9.603148499; 1065600,27.44337611; 1069200,
              48.23236103; 1072800,67.90455038; 1076400,83.61332758; 1080000,
              93.50445335; 1083600,96.54693507; 1087200,92.4422546; 1090800,
              81.59516725; 1094400,65.15276159; 1098000,45.10461005; 1101600,
              24.46735599; 1105200,7.478626154; 1108800,0; 1112400,0; 1116000,0;
              1119600,0; 1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0;
              1141200,0; 1144800,0; 1148400,9.603148499; 1152000,27.44337611;
              1155600,48.23236103; 1159200,67.90455038; 1162800,83.61332758;
              1166400,93.50445335; 1170000,96.54693507; 1173600,92.4422546;
              1177200,81.59516725; 1180800,65.15276159; 1184400,45.10461005;
              1188000,24.46735599; 1191600,7.478626154; 1195200,0; 1198800,0;
              1202400,0; 1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0;
              1224000,0; 1227600,0; 1231200,0; 1234800,2.32324972; 1238400,
              6.600661871; 1242000,11.36216174; 1245600,15.75510219; 1249200,
              19.21671663; 1252800,21.38153013; 1256400,22.04553048; 1260000,
              21.14951495; 1263600,18.7737413; 1267200,15.14507191; 1270800,
              10.65548364; 1274400,5.90307651; 1278000,1.797616109; 1281600,0;
              1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0;
              1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,2.32324972;
              1324800,6.600661871; 1328400,11.36216174; 1332000,15.75510219;
              1335600,19.21671663; 1339200,21.38153013; 1342800,22.04553048;
              1346400,21.14951495; 1350000,18.7737413; 1353600,15.14507191;
              1357200,10.65548364; 1360800,5.90307651; 1364400,1.797616109;
              1368000,0; 1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0;
              1389600,0; 1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,
              2.32324972; 1411200,6.600661871; 1414800,11.36216174; 1418400,
              15.75510219; 1422000,19.21671663; 1425600,21.38153013; 1429200,
              22.04553048; 1432800,21.14951495; 1436400,18.7737413; 1440000,
              15.14507191; 1443600,10.65548364; 1447200,5.90307651; 1450800,
              1.797616109; 1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,
              0; 1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,
              0; 1494000,2.32324972; 1497600,6.600661871; 1501200,11.36216174;
              1504800,15.75510219; 1508400,19.21671663; 1512000,21.38153013;
              1515600,22.04553048; 1519200,21.14951495; 1522800,18.7737413;
              1526400,15.14507191; 1530000,10.65548364; 1533600,5.90307651;
              1537200,1.797616109; 1540800,0; 1544400,0; 1548000,0; 1551600,0;
              1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0;
              1576800,0; 1580400,2.32324972; 1584000,6.600661871; 1587600,
              11.36216174; 1591200,15.75510219; 1594800,19.21671663; 1598400,
              21.38153013; 1602000,22.04553048; 1605600,21.14951495; 1609200,
              18.7737413; 1612800,15.14507191; 1616400,10.65548364; 1620000,
              5.90307651; 1623600,1.797616109; 1627200,0; 1630800,0; 1634400,0;
              1638000,0])
          "Diffuse irradiation at covered sky on tilted surface"
          annotation (Placement(transformation(extent={{-102,-44},{-82,-24}})));
        BoundaryConditions.SolarIrradiation.BaseClasses.DirectTiltedSurface
          HDirTil "Direct irradiation on the tilted window"
          annotation (Placement(transformation(extent={{-26,-96},{-6,-76}})));
        Windows.BaseClasses.Conversions.to_HDirNor to_HDirNor
          "Convertion of the horizontal direct irradiation to the normal direct irradiation"
          annotation (Placement(transformation(extent={{-58,-74},{-48,-64}})));
        Windows.BaseClasses.Illumination illumination(
          e_ILim1=250,
          e_ILim2=500,
          office=true,
          n=1,
          r={0.21},
          A={5.13},
          T_L={0.72},
          D=0.027,
          til={1.5707963267949})
          "determining the switch moments for VDI2078 test case 1"
          annotation (Placement(transformation(extent={{62,-8},{82,12}})));
        Windows.Validation.BaseClasses.SolarHourAngleVDI6007 solarHourAngleVDI(lon=
              0.15009831567151)
          "Solar hour angle based on the calculations of VDI6007"
          annotation (Placement(transformation(extent={{-76,-98},{-68,-90}})));
        BoundaryConditions.SolarGeometry.BaseClasses.IncidenceAngle
          incAng(
          azi(displayUnit="deg") = 0,
          til(displayUnit="deg") = 1.5707963267949,
          lat=0.86393797973719) "Solar incidence angle on the tilted window"
          annotation (Placement(transformation(extent={{-54,-94},{-44,-84}})));
        Windows.Validation.BaseClasses.SolarDeclinationAngleVDI6007 solarDeclinationAngleVDI
          "Solar declination angle based on the calculations of VDI6007"
          annotation (Placement(transformation(extent={{-76,-84},{-68,-76}})));
        Windows.BaseClasses.HVisible HVis(
          n=1,
          T_LTotDir={0.08},
          T_LTotDif={0.32},
          T_L={0.72},
          til={1.5707963267949}) "visible light entering the room"
          annotation (Placement(transformation(extent={{40,24},{90,84}})));
        Modelica.Blocks.Sources.CombiTimeTable HVisSum_VDI2078(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HVis",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              41.97104168; 28800,138.4544724; 32400,318.6762614; 36000,
              188.5098051; 39600,243.2590501; 43200,278.7124998; 46800,
              289.7422498; 50400,274.8750443; 54000,236.1080045; 57600,
              179.1680283; 61200,284.9573681; 64800,119.9320428; 68400,
              31.8869081; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              41.97104168; 115200,138.4544724; 118800,318.6762614; 122400,
              188.5098051; 126000,243.2590501; 129600,278.7124998; 133200,
              289.7422498; 136800,274.8750443; 140400,236.1080045; 144000,
              179.1680283; 147600,284.9573681; 151200,119.9320428; 154800,
              31.8869081; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,41.97104168; 201600,138.4544724; 205200,318.6762614;
              208800,188.5098051; 212400,243.2590501; 216000,278.7124998;
              219600,289.7422498; 223200,274.8750443; 226800,236.1080045;
              230400,179.1680283; 234000,284.9573681; 237600,119.9320428;
              241200,31.8869081; 244800,0; 248400,0; 252000,0; 255600,0; 259200,
              0; 262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,41.97104168; 288000,138.4544724; 291600,318.6762614;
              295200,188.5098051; 298800,243.2590501; 302400,278.7124998;
              306000,289.7422498; 309600,274.8750443; 313200,236.1080045;
              316800,179.1680283; 320400,284.9573681; 324000,119.9320428;
              327600,31.8869081; 331200,0; 334800,0; 338400,0; 342000,0; 345600,
              0; 349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,41.97104168; 374400,138.4544724; 378000,318.6762614;
              381600,188.5098051; 385200,243.2590501; 388800,278.7124998;
              392400,289.7422498; 396000,274.8750443; 399600,236.1080045;
              403200,179.1680283; 406800,284.9573681; 410400,119.9320428;
              414000,31.8869081; 417600,0; 421200,0; 424800,0; 428400,0; 432000,
              0; 435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,41.97104168; 460800,138.4544724; 464400,318.6762614;
              468000,188.5098051; 471600,243.2590501; 475200,278.7124998;
              478800,289.7422498; 482400,274.8750443; 486000,236.1080045;
              489600,179.1680283; 493200,284.9573681; 496800,119.9320428;
              500400,31.8869081; 504000,0; 507600,0; 511200,0; 514800,0; 518400,
              0; 522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,41.97104168; 547200,138.4544724; 550800,318.6762614;
              554400,188.5098051; 558000,243.2590501; 561600,278.7124998;
              565200,289.7422498; 568800,274.8750443; 572400,236.1080045;
              576000,179.1680283; 579600,284.9573681; 583200,119.9320428;
              586800,31.8869081; 590400,0; 594000,0; 597600,0; 601200,0; 604800,
              0; 608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,41.97104168; 633600,138.4544724; 637200,318.6762614;
              640800,188.5098051; 644400,243.2590501; 648000,278.7124998;
              651600,289.7422498; 655200,274.8750443; 658800,236.1080045;
              662400,179.1680283; 666000,284.9573681; 669600,119.9320428;
              673200,31.8869081; 676800,0; 680400,0; 684000,0; 687600,0; 691200,
              0; 694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,41.97104168; 720000,138.4544724; 723600,318.6762614;
              727200,188.5098051; 730800,243.2590501; 734400,278.7124998;
              738000,289.7422498; 741600,274.8750443; 745200,236.1080045;
              748800,179.1680283; 752400,284.9573681; 756000,119.9320428;
              759600,31.8869081; 763200,0; 766800,0; 770400,0; 774000,0; 777600,
              0; 781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,41.97104168; 806400,138.4544724; 810000,318.6762614;
              813600,188.5098051; 817200,243.2590501; 820800,278.7124998;
              824400,289.7422498; 828000,274.8750443; 831600,236.1080045;
              835200,179.1680283; 838800,284.9573681; 842400,119.9320428;
              846000,31.8869081; 849600,0; 853200,0; 856800,0; 860400,0; 864000,
              0; 867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,41.97104168; 892800,138.4544724; 896400,318.6762614;
              900000,188.5098051; 903600,243.2590501; 907200,278.7124998;
              910800,289.7422498; 914400,274.8750443; 918000,236.1080045;
              921600,179.1680283; 925200,284.9573681; 928800,119.9320428;
              932400,31.8869081; 936000,0; 939600,0; 943200,0; 946800,0; 950400,
              0; 954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,41.97104168; 979200,138.4544724; 982800,318.6762614;
              986400,188.5098051; 990000,243.2590501; 993600,278.7124998;
              997200,289.7422498; 1000800,274.8750443; 1004400,236.1080045;
              1008000,179.1680283; 1011600,284.9573681; 1015200,119.9320428;
              1018800,31.8869081; 1022400,0; 1026000,0; 1029600,0; 1033200,0;
              1036800,0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0;
              1058400,0; 1062000,41.97104168; 1065600,138.4544724; 1069200,
              318.6762614; 1072800,188.5098051; 1076400,243.2590501; 1080000,
              278.7124998; 1083600,289.7422498; 1087200,274.8750443; 1090800,
              236.1080045; 1094400,179.1680283; 1098000,284.9573681; 1101600,
              119.9320428; 1105200,31.8869081; 1108800,0; 1112400,0; 1116000,0;
              1119600,0; 1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0;
              1141200,0; 1144800,0; 1148400,41.97104168; 1152000,138.4544724;
              1155600,318.6762614; 1159200,188.5098051; 1162800,243.2590501;
              1166400,278.7124998; 1170000,289.7422498; 1173600,274.8750443;
              1177200,236.1080045; 1180800,179.1680283; 1184400,284.9573681;
              1188000,119.9320428; 1191600,31.8869081; 1195200,0; 1198800,0;
              1202400,0; 1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0;
              1224000,0; 1227600,0; 1231200,0; 1234800,33.60423886; 1238400,
              134.7689486; 1242000,141.0670302; 1245600,219.8126796; 1249200,
              283.8509415; 1252800,323.4889321; 1256400,335.5236891; 1260000,
              319.269018; 1263600,275.6788729; 1267200,208.5610811; 1270800,
              129.0805857; 1274400,111.0152963; 1278000,24.65905213; 1281600,0;
              1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0;
              1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,33.60423886;
              1324800,134.7689486; 1328400,141.0670302; 1332000,219.8126796;
              1335600,283.8509415; 1339200,323.4889321; 1342800,335.5236891;
              1346400,319.269018; 1350000,275.6788729; 1353600,208.5610811;
              1357200,129.0805857; 1360800,111.0152963; 1364400,24.65905213;
              1368000,0; 1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0;
              1389600,0; 1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,
              33.60423886; 1411200,134.7689486; 1414800,141.0670302; 1418400,
              219.8126796; 1422000,283.8509415; 1425600,323.4889321; 1429200,
              335.5236891; 1432800,319.269018; 1436400,275.6788729; 1440000,
              208.5610811; 1443600,129.0805857; 1447200,111.0152963; 1450800,
              24.65905213; 1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,
              0; 1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,
              0; 1494000,33.60423886; 1497600,134.7689486; 1501200,141.0670302;
              1504800,219.8126796; 1508400,283.8509415; 1512000,323.4889321;
              1515600,335.5236891; 1519200,319.269018; 1522800,275.6788729;
              1526400,208.5610811; 1530000,129.0805857; 1533600,111.0152963;
              1537200,24.65905213; 1540800,0; 1544400,0; 1548000,0; 1551600,0;
              1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0;
              1576800,0; 1580400,33.60423886; 1584000,134.7689486; 1587600,
              141.0670302; 1591200,219.8126796; 1594800,283.8509415; 1598400,
              323.4889321; 1602000,335.5236891; 1605600,319.269018; 1609200,
              275.6788729; 1612800,208.5610811; 1616400,129.0805857; 1620000,
              111.0152963; 1623600,24.65905213; 1627200,0; 1630800,0; 1634400,0;
              1638000,0]) "Comparison for  the entering visible light"
          annotation (Placement(transformation(extent={{36,-70},{56,-50}})));
        Modelica.Blocks.Sources.CombiTimeTable HLimVis_VDI2078(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HLimVis",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              98.34992957; 28800,196.6998591; 32400,196.6998591; 36000,
              196.6998591; 39600,196.6998591; 43200,196.6998591; 46800,
              196.6998591; 50400,196.6998591; 54000,196.6998591; 57600,
              98.34992957; 61200,98.34992957; 64800,0; 68400,0; 72000,0; 75600,
              0; 79200,0; 82800,0; 86400,0; 90000,0; 93600,0; 97200,0; 100800,0;
              104400,0; 108000,0; 111600,98.34992957; 115200,196.6998591;
              118800,196.6998591; 122400,196.6998591; 126000,196.6998591;
              129600,196.6998591; 133200,196.6998591; 136800,196.6998591;
              140400,196.6998591; 144000,98.34992957; 147600,98.34992957;
              151200,0; 154800,0; 158400,0; 162000,0; 165600,0; 169200,0;
              172800,0; 176400,0; 180000,0; 183600,0; 187200,0; 190800,0;
              194400,0; 198000,98.34992957; 201600,196.6998591; 205200,
              196.6998591; 208800,196.6998591; 212400,196.6998591; 216000,
              196.6998591; 219600,196.6998591; 223200,196.6998591; 226800,
              196.6998591; 230400,98.34992957; 234000,98.34992957; 237600,0;
              241200,0; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,98.34992957; 288000,196.6998591; 291600,196.6998591;
              295200,196.6998591; 298800,196.6998591; 302400,196.6998591;
              306000,196.6998591; 309600,196.6998591; 313200,196.6998591;
              316800,98.34992957; 320400,98.34992957; 324000,0; 327600,0;
              331200,0; 334800,0; 338400,0; 342000,0; 345600,0; 349200,0;
              352800,0; 356400,0; 360000,0; 363600,0; 367200,0; 370800,
              98.34992957; 374400,196.6998591; 378000,196.6998591; 381600,
              196.6998591; 385200,196.6998591; 388800,196.6998591; 392400,
              196.6998591; 396000,196.6998591; 399600,196.6998591; 403200,
              98.34992957; 406800,98.34992957; 410400,0; 414000,0; 417600,0;
              421200,0; 424800,0; 428400,0; 432000,0; 435600,0; 439200,0;
              442800,0; 446400,0; 450000,0; 453600,0; 457200,0; 460800,0;
              464400,0; 468000,0; 471600,0; 475200,0; 478800,0; 482400,0;
              486000,0; 489600,0; 493200,0; 496800,0; 500400,0; 504000,0;
              507600,0; 511200,0; 514800,0; 518400,0; 522000,0; 525600,0;
              529200,0; 532800,0; 536400,0; 540000,0; 543600,0; 547200,0;
              550800,0; 554400,0; 558000,0; 561600,0; 565200,0; 568800,0;
              572400,0; 576000,0; 579600,0; 583200,0; 586800,0; 590400,0;
              594000,0; 597600,0; 601200,0; 604800,0; 608400,0; 612000,0;
              615600,0; 619200,0; 622800,0; 626400,0; 630000,98.34992957;
              633600,196.6998591; 637200,196.6998591; 640800,196.6998591;
              644400,196.6998591; 648000,196.6998591; 651600,196.6998591;
              655200,196.6998591; 658800,196.6998591; 662400,98.34992957;
              666000,98.34992957; 669600,0; 673200,0; 676800,0; 680400,0;
              684000,0; 687600,0; 691200,0; 694800,0; 698400,0; 702000,0;
              705600,0; 709200,0; 712800,0; 716400,98.34992957; 720000,
              196.6998591; 723600,196.6998591; 727200,196.6998591; 730800,
              196.6998591; 734400,196.6998591; 738000,196.6998591; 741600,
              196.6998591; 745200,196.6998591; 748800,98.34992957; 752400,
              98.34992957; 756000,0; 759600,0; 763200,0; 766800,0; 770400,0;
              774000,0; 777600,0; 781200,0; 784800,0; 788400,0; 792000,0;
              795600,0; 799200,0; 802800,98.34992957; 806400,196.6998591;
              810000,196.6998591; 813600,196.6998591; 817200,196.6998591;
              820800,196.6998591; 824400,196.6998591; 828000,196.6998591;
              831600,196.6998591; 835200,98.34992957; 838800,98.34992957;
              842400,0; 846000,0; 849600,0; 853200,0; 856800,0; 860400,0;
              864000,0; 867600,0; 871200,0; 874800,0; 878400,0; 882000,0;
              885600,0; 889200,98.34992957; 892800,196.6998591; 896400,
              196.6998591; 900000,196.6998591; 903600,196.6998591; 907200,
              196.6998591; 910800,196.6998591; 914400,196.6998591; 918000,
              196.6998591; 921600,98.34992957; 925200,98.34992957; 928800,0;
              932400,0; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,98.34992957; 979200,196.6998591; 982800,196.6998591;
              986400,196.6998591; 990000,196.6998591; 993600,196.6998591;
              997200,196.6998591; 1000800,196.6998591; 1004400,196.6998591;
              1008000,98.34992957; 1011600,98.34992957; 1015200,0; 1018800,0;
              1022400,0; 1026000,0; 1029600,0; 1033200,0; 1036800,0; 1040400,0;
              1044000,0; 1047600,0; 1051200,0; 1054800,0; 1058400,0; 1062000,0;
              1065600,0; 1069200,0; 1072800,0; 1076400,0; 1080000,0; 1083600,0;
              1087200,0; 1090800,0; 1094400,0; 1098000,0; 1101600,0; 1105200,0;
              1108800,0; 1112400,0; 1116000,0; 1119600,0; 1123200,0; 1126800,0;
              1130400,0; 1134000,0; 1137600,0; 1141200,0; 1144800,0; 1148400,0;
              1152000,0; 1155600,0; 1159200,0; 1162800,0; 1166400,0; 1170000,0;
              1173600,0; 1177200,0; 1180800,0; 1184400,0; 1188000,0; 1191600,0;
              1195200,0; 1198800,0; 1202400,0; 1206000,0; 1209600,0; 1213200,0;
              1216800,0; 1220400,0; 1224000,0; 1227600,0; 1231200,0; 1234800,
              98.34992957; 1238400,196.6998591; 1242000,196.6998591; 1245600,
              196.6998591; 1249200,196.6998591; 1252800,196.6998591; 1256400,
              196.6998591; 1260000,196.6998591; 1263600,196.6998591; 1267200,
              98.34992957; 1270800,98.34992957; 1274400,0; 1278000,0; 1281600,0;
              1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0;
              1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,98.34992957;
              1324800,196.6998591; 1328400,196.6998591; 1332000,196.6998591;
              1335600,196.6998591; 1339200,196.6998591; 1342800,196.6998591;
              1346400,196.6998591; 1350000,196.6998591; 1353600,98.34992957;
              1357200,98.34992957; 1360800,0; 1364400,0; 1368000,0; 1371600,0;
              1375200,0; 1378800,0; 1382400,0; 1386000,0; 1389600,0; 1393200,0;
              1396800,0; 1400400,0; 1404000,0; 1407600,98.34992957; 1411200,
              196.6998591; 1414800,196.6998591; 1418400,196.6998591; 1422000,
              196.6998591; 1425600,196.6998591; 1429200,196.6998591; 1432800,
              196.6998591; 1436400,196.6998591; 1440000,98.34992957; 1443600,
              98.34992957; 1447200,0; 1450800,0; 1454400,0; 1458000,0; 1461600,
              0; 1465200,0; 1468800,0; 1472400,0; 1476000,0; 1479600,0; 1483200,
              0; 1486800,0; 1490400,0; 1494000,98.34992957; 1497600,196.6998591;
              1501200,196.6998591; 1504800,196.6998591; 1508400,196.6998591;
              1512000,196.6998591; 1515600,196.6998591; 1519200,196.6998591;
              1522800,196.6998591; 1526400,98.34992957; 1530000,98.34992957;
              1533600,0; 1537200,0; 1540800,0; 1544400,0; 1548000,0; 1551600,0;
              1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0;
              1576800,0; 1580400,98.34992957; 1584000,196.6998591; 1587600,
              196.6998591; 1591200,196.6998591; 1594800,196.6998591; 1598400,
              196.6998591; 1602000,196.6998591; 1605600,196.6998591; 1609200,
              196.6998591; 1612800,98.34992957; 1616400,98.34992957; 1620000,0;
              1623600,0; 1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Comparison for the limit to activate the illumination"
          annotation (Placement(transformation(extent={{36,-42},{56,-22}})));
        Modelica.Blocks.Logical.Greater illumination_VDI2078
          "Comparison for the illumination boolean"
          annotation (Placement(transformation(extent={{82,-42},{102,-22}})));
        Modelica.Blocks.Math.Add HDifTil
          annotation (Placement(transformation(extent={{-54,-2},{-38,14}})));
      equation
        connect(sunblind.sunscreen, CorGTaue.sunscreen[1]) annotation (Line(
              points={{-5.3,-41},{-5.3,20},{-52,20},{-56,20},{-56,42},{-56,54},
                {-51,54}}, color={255,0,255}));
        connect(to_HDirNor.HDirNor, HDirTil.HDirNor) annotation (Line(points={{-47,-69},
                {-34,-69},{-34,-80},{-28,-80}},  color={0,0,127}));
        connect(HDirHor.y[1], to_HDirNor.HDirHor) annotation (Line(points={{-81,104},
                {-70,104},{-70,-67},{-59,-67}},      color={0,0,127}));
        connect(alt.y[1], to_HDirNor.alt) annotation (Line(points={{-81,-60},{
                -59,-60},{-59,-71}}, color={0,0,127}));
        connect(HDirTil.HDirTil, sunblind.HDirTil) annotation (Line(
            points={{-5,-86},{0,-86},{0,-56},{-34,-56},{-34,-44},{-20.7,-44}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solarHourAngleVDI.solHouAng, incAng.solHouAng) annotation (Line(
              points={{-67.6,-94},{-62,-94},{-62,-91.4},{-55,-91.4}}, color={0,0,
                127}));
        connect(solarDeclinationAngleVDI.decAng, incAng.decAng) annotation (Line(
              points={{-67.6,-80},{-63.8,-80},{-63.8,-86.3},{-55.1,-86.3}}, color=
               {0,0,127}));
        connect(incAng.incAng, HDirTil.incAng) annotation (Line(points={{-43.5,
                -89},{-35.75,-89},{-35.75,-92},{-28,-92}}, color={0,0,127}));
        connect(incAng.incAng, CorGTaue.incAng[1]) annotation (Line(points={{-43.5,
                -89},{-43.5,-18},{-60,-18},{-60,58},{-51,58}},       color={0,0,
                127}));
        connect(CorGTaue.CorTaue_Gro, illumination.CorTaue_Gro) annotation (Line(
              points={{-33,64},{0,64},{0,54},{14,54},{14,46},{12,46},{12,2},{61,
                2}},
              color={0,0,127}));
        connect(CorGTaue.CorTaue_DifCov, illumination.CorTaue_DifCov) annotation (
           Line(points={{-33,62},{24,62},{24,-4},{61,-4}}, color={0,0,127}));
        connect(HVis.HVis, illumination.HVis) annotation (Line(points={{92.5,54},
                {98,54},{98,52},{98,20},{56,20},{56,8},{61,8}},         color={0,
                0,127}));
        connect(to_HDirNor.HDirNor, HVis.HDirNor) annotation (Line(points={{-47,-69},
                {-2,-69},{-2,83.7},{38.25,83.7}},        color={0,0,127}));
        connect(HDirTil.HDirTil, HVis.HDirTil[1]) annotation (Line(points={{-5,-86},
                {2,-86},{2,78.9},{38.25,78.9}},        color={0,0,127}));
        connect(HDifHorCov.y[1], HVis.HDifHorCov) annotation (Line(points={{-81,
                20},{-22,20},{-22,64.5},{38.25,64.5}}, color={0,0,127}));
        connect(HDifHorCle.y[1], HVis.HDifHorCle) annotation (Line(points={{-81,
                48},{-54,48},{-54,42},{-28,42},{-28,59.7},{38.25,59.7}}, color=
                {0,0,127}));
        connect(sunblind.sunscreen, HVis.sunscreen[1]) annotation (Line(points=
                {{-5.3,-41},{-5.3,6},{18,6},{18,64},{38.25,64},{38.25,54.3}},
              color={255,0,255}));
        connect(alt.y[1], HVis.alt) annotation (Line(points={{-81,-60},{-30,-60},
                {-30,38},{4,38},{4,48.3},{38.25,48.3}}, color={0,0,127}));
        connect(CorGTaue.CorTaue_Gro, HVis.CorTaue_Gro) annotation (Line(points={
                {-33,64},{2,64},{2,42.3},{38.25,42.3}}, color={0,0,127}));
        connect(CorGTaue.CorTaue_DifCov, HVis.CorTaue_DifCov) annotation (Line(
              points={{-33,62},{2,62},{2,36.3},{38.25,36.3}}, color={0,0,127}));
        connect(CorGTaue.CorTaue_DifCle, HVis.CorTaue_DifCle) annotation (Line(
              points={{-33,60},{2,60},{2,30.3},{38.25,30.3}}, color={0,0,127}));
        connect(CorGTaue.CorTaue_Dir, HVis.CorTaue_Dir) annotation (Line(points={
                {-33,58},{2,58},{2,24.3},{38.25,24.3}}, color={0,0,127}));
        connect(HDifTilCov.y, HVis.HDifTilCov) annotation (Line(points={{-81,
                -34},{-22,-34},{-22,69.3},{38.25,69.3}}, color={0,0,127}));
        connect(HDifTilCle.y, HVis.HDifTilCle) annotation (Line(points={{-81,-8},
                {-22,-8},{-22,74.1},{38.25,74.1}}, color={0,0,127}));
        connect(HLimVis_VDI2078.y[1], illumination_VDI2078.u1) annotation (Line(
              points={{57,-32},{66,-32},{80,-32}}, color={0,0,127}));
        connect(HVisSum_VDI2078.y[1], illumination_VDI2078.u2) annotation (Line(
              points={{57,-60},{68,-60},{68,-40},{80,-40}}, color={0,0,127}));
        connect(HDifTilCle.y[1], HDifTil.u1) annotation (Line(
            points={{-81,-8},{-68,-8},{-68,10.8},{-55.6,10.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(HDifTilCov.y[1], HDifTil.u2) annotation (Line(
            points={{-81,-34},{-64,-34},{-64,0},{-60,0},{-60,1.2},{-55.6,1.2}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(HDifTil.y, sunblind.HDifTil) annotation (Line(
            points={{-37.2,6},{-28,6},{-28,-38.1},{-20.7,-38.1}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (experiment(StartTime=0,StopTime=1638000),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
                  -100},{120,120}}), graphics),
          Icon(coordinateSystem(extent={{-120,-100},{120,120}})),
          Documentation(info="<html>
<p>This model simulates parts of VDI2078 test case 1. The solar irradiation is treated as an input. To calculate the boundary conditions <a href=\"Annex60.BoundaryConditions\">Annex60</a> models are used. </p>
</html>"));
      end TestCase1_Illumination;

      model TestCase3_VentilationHeat
        extends Modelica.Icons.Example;

        Windows.BaseClasses.VentilationHeat ventilationHeat(
          x_f=0.8,
          d=0.1,
          screen=false,
          tau_e=0,
          rho_e=0.8125,
          til=1.5707963267949)
          "Calculates the heat input due to ventialtion for test case 3 of VDI2078"
          annotation (Placement(transformation(extent={{80,-10},{100,10}})));
        Windows.BaseClasses.Sunblind sunblind(lim=200)
          "Calculates if the sunblind of the window is active"
          annotation (Placement(transformation(extent={{48,-2},{56,6}})));
        Modelica.Blocks.Sources.CombiTimeTable HDirHor(
          columns={2},
          tableName="HDirHor",
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableOnFile=true,
          table=[0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0; 0.0,0.0;
              0.0,0.0],
          fileName=
              "C:/Users/Stanley/Documents/BA/VDI6007-3/Resources/WeatherData/Trimble/Trimble_Wetterdaten_April_HDirHor.txt")
          "Direct irradiation on horizontal surface"
          annotation (Placement(transformation(extent={{-96,-26},{-82,-12}})));
        Modelica.Blocks.Sources.CombiTimeTable alt(
          columns={2},
          tableName="alt",
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableOnFile=false,
          table=[0,-0.524971501; 3600,-0.540504674; 7200,-0.511513997; 10800,-0.436196642;
              14400,-0.323752929; 18000,-0.184550738; 21600,-0.027807643; 25200,
              0.138810415; 28800,0.308550932; 32400,0.474473156; 36000,
              0.627965928; 39600,0.756893419; 43200,0.844457705; 46800,
              0.872876752; 50400,0.834722167; 54000,0.739755433; 57600,
              0.606187916; 61200,0.450168115; 64800,0.28318588; 68400,
              0.113490388; 72000,-0.052089885; 75600,-0.206727187; 79200,-0.342547809;
              82800,-0.450098665; 86400,-0.524971501; 90000,-0.540504674; 93600,
              -0.511513997; 97200,-0.436196642; 100800,-0.323752929; 104400,-0.184550738;
              108000,-0.027807643; 111600,0.138810415; 115200,0.308550932;
              118800,0.474473156; 122400,0.627965928; 126000,0.756893419;
              129600,0.844457705; 133200,0.872876752; 136800,0.834722167;
              140400,0.739755433; 144000,0.606187916; 147600,0.450168115;
              151200,0.28318588; 154800,0.113490388; 158400,-0.052089885;
              162000,-0.206727187; 165600,-0.342547809; 169200,-0.450098665;
              172800,-0.524971501; 176400,-0.540504674; 180000,-0.511513997;
              183600,-0.436196642; 187200,-0.323752929; 190800,-0.184550738;
              194400,-0.027807643; 198000,0.138810415; 201600,0.308550932;
              205200,0.474473156; 208800,0.627965928; 212400,0.756893419;
              216000,0.844457705; 219600,0.872876752; 223200,0.834722167;
              226800,0.739755433; 230400,0.606187916; 234000,0.450168115;
              237600,0.28318588; 241200,0.113490388; 244800,-0.052089885;
              248400,-0.206727187; 252000,-0.342547809; 255600,-0.450098665;
              259200,-0.524971501; 262800,-0.540504674; 266400,-0.511513997;
              270000,-0.436196642; 273600,-0.323752929; 277200,-0.184550738;
              280800,-0.027807643; 284400,0.138810415; 288000,0.308550932;
              291600,0.474473156; 295200,0.627965928; 298800,0.756893419;
              302400,0.844457705; 306000,0.872876752; 309600,0.834722167;
              313200,0.739755433; 316800,0.606187916; 320400,0.450168115;
              324000,0.28318588; 327600,0.113490388; 331200,-0.052089885;
              334800,-0.206727187; 338400,-0.342547809; 342000,-0.450098665;
              345600,-0.524971501; 349200,-0.540504674; 352800,-0.511513997;
              356400,-0.436196642; 360000,-0.323752929; 363600,-0.184550738;
              367200,-0.027807643; 370800,0.138810415; 374400,0.308550932;
              378000,0.474473156; 381600,0.627965928; 385200,0.756893419;
              388800,0.844457705; 392400,0.872876752; 396000,0.834722167;
              399600,0.739755433; 403200,0.606187916; 406800,0.450168115;
              410400,0.28318588; 414000,0.113490388; 417600,-0.052089885;
              421200,-0.206727187; 424800,-0.342547809; 428400,-0.450098665;
              432000,-0.524971501; 435600,-0.540504674; 439200,-0.511513997;
              442800,-0.436196642; 446400,-0.323752929; 450000,-0.184550738;
              453600,-0.027807643; 457200,0.138810415; 460800,0.308550932;
              464400,0.474473156; 468000,0.627965928; 471600,0.756893419;
              475200,0.844457705; 478800,0.872876752; 482400,0.834722167;
              486000,0.739755433; 489600,0.606187916; 493200,0.450168115;
              496800,0.28318588; 500400,0.113490388; 504000,-0.052089885;
              507600,-0.206727187; 511200,-0.342547809; 514800,-0.450098665;
              518400,-0.524971501; 522000,-0.540504674; 525600,-0.511513997;
              529200,-0.436196642; 532800,-0.323752929; 536400,-0.184550738;
              540000,-0.027807643; 543600,0.138810415; 547200,0.308550932;
              550800,0.474473156; 554400,0.627965928; 558000,0.756893419;
              561600,0.844457705; 565200,0.872876752; 568800,0.834722167;
              572400,0.739755433; 576000,0.606187916; 579600,0.450168115;
              583200,0.28318588; 586800,0.113490388; 590400,-0.052089885;
              594000,-0.206727187; 597600,-0.342547809; 601200,-0.450098665;
              604800,-0.524971501; 608400,-0.540504674; 612000,-0.511513997;
              615600,-0.436196642; 619200,-0.323752929; 622800,-0.184550738;
              626400,-0.027807643; 630000,0.138810415; 633600,0.308550932;
              637200,0.474473156; 640800,0.627965928; 644400,0.756893419;
              648000,0.844457705; 651600,0.872876752; 655200,0.834722167;
              658800,0.739755433; 662400,0.606187916; 666000,0.450168115;
              669600,0.28318588; 673200,0.113490388; 676800,-0.052089885;
              680400,-0.206727187; 684000,-0.342547809; 687600,-0.450098665;
              691200,-0.524971501; 694800,-0.540504674; 698400,-0.511513997;
              702000,-0.436196642; 705600,-0.323752929; 709200,-0.184550738;
              712800,-0.027807643; 716400,0.138810415; 720000,0.308550932;
              723600,0.474473156; 727200,0.627965928; 730800,0.756893419;
              734400,0.844457705; 738000,0.872876752; 741600,0.834722167;
              745200,0.739755433; 748800,0.606187916; 752400,0.450168115;
              756000,0.28318588; 759600,0.113490388; 763200,-0.052089885;
              766800,-0.206727187; 770400,-0.342547809; 774000,-0.450098665;
              777600,-0.524971501; 781200,-0.540504674; 784800,-0.511513997;
              788400,-0.436196642; 792000,-0.323752929; 795600,-0.184550738;
              799200,-0.027807643; 802800,0.138810415; 806400,0.308550932;
              810000,0.474473156; 813600,0.627965928; 817200,0.756893419;
              820800,0.844457705; 824400,0.872876752; 828000,0.834722167;
              831600,0.739755433; 835200,0.606187916; 838800,0.450168115;
              842400,0.28318588; 846000,0.113490388; 849600,-0.052089885;
              853200,-0.206727187; 856800,-0.342547809; 860400,-0.450098665;
              864000,-0.524971501; 867600,-0.540504674; 871200,-0.511513997;
              874800,-0.436196642; 878400,-0.323752929; 882000,-0.184550738;
              885600,-0.027807643; 889200,0.138810415; 892800,0.308550932;
              896400,0.474473156; 900000,0.627965928; 903600,0.756893419;
              907200,0.844457705; 910800,0.872876752; 914400,0.834722167;
              918000,0.739755433; 921600,0.606187916; 925200,0.450168115;
              928800,0.28318588; 932400,0.113490388; 936000,-0.052089885;
              939600,-0.206727187; 943200,-0.342547809; 946800,-0.450098665;
              950400,-0.524971501; 954000,-0.540504674; 957600,-0.511513997;
              961200,-0.436196642; 964800,-0.323752929; 968400,-0.184550738;
              972000,-0.027807643; 975600,0.138810415; 979200,0.308550932;
              982800,0.474473156; 986400,0.627965928; 990000,0.756893419;
              993600,0.844457705; 997200,0.872876752; 1000800,0.834722167;
              1004400,0.739755433; 1008000,0.606187916; 1011600,0.450168115;
              1015200,0.28318588; 1018800,0.113490388; 1022400,-0.052089885;
              1026000,-0.206727187; 1029600,-0.342547809; 1033200,-0.450098665;
              1036800,-0.524971501; 1040400,-0.540504674; 1044000,-0.511513997;
              1047600,-0.436196642; 1051200,-0.323752929; 1054800,-0.184550738;
              1058400,-0.027807643; 1062000,0.138810415; 1065600,0.308550932;
              1069200,0.474473156; 1072800,0.627965928; 1076400,0.756893419;
              1080000,0.844457705; 1083600,0.872876752; 1087200,0.834722167;
              1090800,0.739755433; 1094400,0.606187916; 1098000,0.450168115;
              1101600,0.28318588; 1105200,0.113490388; 1108800,-0.052089885;
              1112400,-0.206727187; 1116000,-0.342547809; 1119600,-0.450098665;
              1123200,-0.524971501; 1126800,-0.540504674; 1130400,-0.511513997;
              1134000,-0.436196642; 1137600,-0.323752929; 1141200,-0.184550738;
              1144800,-0.027807643; 1148400,0.138810415; 1152000,0.308550932;
              1155600,0.474473156; 1159200,0.627965928; 1162800,0.756893419;
              1166400,0.844457705; 1170000,0.872876752; 1173600,0.834722167;
              1177200,0.739755433; 1180800,0.606187916; 1184400,0.450168115;
              1188000,0.28318588; 1191600,0.113490388; 1195200,-0.052089885;
              1198800,-0.206727187; 1202400,-0.342547809; 1206000,-0.450098665;
              1209600,-0.524971501; 1213200,-0.540504674; 1216800,-0.511513997;
              1220400,-0.436196642; 1224000,-0.323752929; 1227600,-0.184550738;
              1231200,-0.027807643; 1234800,0.138810415; 1238400,0.308550932;
              1242000,0.474473156; 1245600,0.627965928; 1249200,0.756893419;
              1252800,0.844457705; 1256400,0.872876752; 1260000,0.834722167;
              1263600,0.739755433; 1267200,0.606187916; 1270800,0.450168115;
              1274400,0.28318588; 1278000,0.113490388; 1281600,-0.052089885;
              1285200,-0.206727187; 1288800,-0.342547809; 1292400,-0.450098665;
              1296000,-0.524971501; 1299600,-0.540504674; 1303200,-0.511513997;
              1306800,-0.436196642; 1310400,-0.323752929; 1314000,-0.184550738;
              1317600,-0.027807643; 1321200,0.138810415; 1324800,0.308550932;
              1328400,0.474473156; 1332000,0.627965928; 1335600,0.756893419;
              1339200,0.844457705; 1342800,0.872876752; 1346400,0.834722167;
              1350000,0.739755433; 1353600,0.606187916; 1357200,0.450168115;
              1360800,0.28318588; 1364400,0.113490388; 1368000,-0.052089885;
              1371600,-0.206727187; 1375200,-0.342547809; 1378800,-0.450098665;
              1382400,-0.524971501; 1386000,-0.540504674; 1389600,-0.511513997;
              1393200,-0.436196642; 1396800,-0.323752929; 1400400,-0.184550738;
              1404000,-0.027807643; 1407600,0.138810415; 1411200,0.308550932;
              1414800,0.474473156; 1418400,0.627965928; 1422000,0.756893419;
              1425600,0.844457705; 1429200,0.872876752; 1432800,0.834722167;
              1436400,0.739755433; 1440000,0.606187916; 1443600,0.450168115;
              1447200,0.28318588; 1450800,0.113490388; 1454400,-0.052089885;
              1458000,-0.206727187; 1461600,-0.342547809; 1465200,-0.450098665;
              1468800,-0.524971501; 1472400,-0.540504674; 1476000,-0.511513997;
              1479600,-0.436196642; 1483200,-0.323752929; 1486800,-0.184550738;
              1490400,-0.027807643; 1494000,0.138810415; 1497600,0.308550932;
              1501200,0.474473156; 1504800,0.627965928; 1508400,0.756893419;
              1512000,0.844457705; 1515600,0.872876752; 1519200,0.834722167;
              1522800,0.739755433; 1526400,0.606187916; 1530000,0.450168115;
              1533600,0.28318588; 1537200,0.113490388; 1540800,-0.052089885;
              1544400,-0.206727187; 1548000,-0.342547809; 1551600,-0.450098665;
              1555200,-0.524971501; 1558800,-0.540504674; 1562400,-0.511513997;
              1566000,-0.436196642; 1569600,-0.323752929; 1573200,-0.184550738;
              1576800,-0.027807643; 1580400,0.138810415; 1584000,0.308550932;
              1587600,0.474473156; 1591200,0.627965928; 1594800,0.756893419;
              1598400,0.844457705; 1602000,0.872876752; 1605600,0.834722167;
              1609200,0.739755433; 1612800,0.606187916; 1616400,0.450168115;
              1620000,0.28318588; 1623600,0.113490388; 1627200,-0.052089885;
              1630800,-0.206727187; 1634400,-0.342547809; 1638000,-0.450098665])
          "Solar altitude angle"
          annotation (Placement(transformation(extent={{-96,-52},{-82,-38}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifTilCle(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifTilCle",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              12.09644867; 28800,27.25095027; 32400,38.27385402; 36000,45.82451;
              39600,51.57958573; 43200,54.91750562; 46800,55.89716878; 50400,
              54.57011009; 54000,50.87097992; 57600,44.75882584; 61200,
              37.07302936; 64800,25.04319434; 68400,9.799046246; 72000,0; 75600,
              0; 79200,0; 82800,0; 86400,0; 90000,0; 93600,0; 97200,0; 100800,0;
              104400,0; 108000,0; 111600,12.09644867; 115200,27.25095027;
              118800,38.27385402; 122400,45.82451; 126000,51.57958573; 129600,
              54.91750562; 133200,55.89716878; 136800,54.57011009; 140400,
              50.87097992; 144000,44.75882584; 147600,37.07302936; 151200,
              25.04319434; 154800,9.799046246; 158400,0; 162000,0; 165600,0;
              169200,0; 172800,0; 176400,0; 180000,0; 183600,0; 187200,0;
              190800,0; 194400,0; 198000,12.09644867; 201600,27.25095027;
              205200,38.27385402; 208800,45.82451; 212400,51.57958573; 216000,
              54.91750562; 219600,55.89716878; 223200,54.57011009; 226800,
              50.87097992; 230400,44.75882584; 234000,37.07302936; 237600,
              25.04319434; 241200,9.799046246; 244800,0; 248400,0; 252000,0;
              255600,0; 259200,0; 262800,0; 266400,0; 270000,0; 273600,0;
              277200,0; 280800,0; 284400,12.09644867; 288000,27.25095027;
              291600,38.27385402; 295200,45.82451; 298800,51.57958573; 302400,
              54.91750562; 306000,55.89716878; 309600,54.57011009; 313200,
              50.87097992; 316800,44.75882584; 320400,37.07302936; 324000,
              25.04319434; 327600,9.799046246; 331200,0; 334800,0; 338400,0;
              342000,0; 345600,0; 349200,0; 352800,0; 356400,0; 360000,0;
              363600,0; 367200,0; 370800,12.09644867; 374400,27.25095027;
              378000,38.27385402; 381600,45.82451; 385200,51.57958573; 388800,
              54.91750562; 392400,55.89716878; 396000,54.57011009; 399600,
              50.87097992; 403200,44.75882584; 406800,37.07302936; 410400,
              25.04319434; 414000,9.799046246; 417600,0; 421200,0; 424800,0;
              428400,0; 432000,0; 435600,0; 439200,0; 442800,0; 446400,0;
              450000,0; 453600,0; 457200,12.09644867; 460800,27.25095027;
              464400,38.27385402; 468000,45.82451; 471600,51.57958573; 475200,
              54.91750562; 478800,55.89716878; 482400,54.57011009; 486000,
              50.87097992; 489600,44.75882584; 493200,37.07302936; 496800,
              25.04319434; 500400,9.799046246; 504000,0; 507600,0; 511200,0;
              514800,0; 518400,0; 522000,0; 525600,0; 529200,0; 532800,0;
              536400,0; 540000,0; 543600,12.09644867; 547200,27.25095027;
              550800,38.27385402; 554400,45.82451; 558000,51.57958573; 561600,
              54.91750562; 565200,55.89716878; 568800,54.57011009; 572400,
              50.87097992; 576000,44.75882584; 579600,37.07302936; 583200,
              25.04319434; 586800,9.799046246; 590400,0; 594000,0; 597600,0;
              601200,0; 604800,0; 608400,0; 612000,0; 615600,0; 619200,0;
              622800,0; 626400,0; 630000,12.09644867; 633600,27.25095027;
              637200,38.27385402; 640800,45.82451; 644400,51.57958573; 648000,
              54.91750562; 651600,55.89716878; 655200,54.57011009; 658800,
              50.87097992; 662400,44.75882584; 666000,37.07302936; 669600,
              25.04319434; 673200,9.799046246; 676800,0; 680400,0; 684000,0;
              687600,0; 691200,0; 694800,0; 698400,0; 702000,0; 705600,0;
              709200,0; 712800,0; 716400,12.09644867; 720000,27.25095027;
              723600,38.27385402; 727200,45.82451; 730800,51.57958573; 734400,
              54.91750562; 738000,55.89716878; 741600,54.57011009; 745200,
              50.87097992; 748800,44.75882584; 752400,37.07302936; 756000,
              25.04319434; 759600,9.799046246; 763200,0; 766800,0; 770400,0;
              774000,0; 777600,0; 781200,0; 784800,0; 788400,0; 792000,0;
              795600,0; 799200,0; 802800,12.09644867; 806400,27.25095027;
              810000,38.27385402; 813600,45.82451; 817200,51.57958573; 820800,
              54.91750562; 824400,55.89716878; 828000,54.57011009; 831600,
              50.87097992; 835200,44.75882584; 838800,37.07302936; 842400,
              25.04319434; 846000,9.799046246; 849600,0; 853200,0; 856800,0;
              860400,0; 864000,0; 867600,0; 871200,0; 874800,0; 878400,0;
              882000,0; 885600,0; 889200,12.09644867; 892800,27.25095027;
              896400,38.27385402; 900000,45.82451; 903600,51.57958573; 907200,
              54.91750562; 910800,55.89716878; 914400,54.57011009; 918000,
              50.87097992; 921600,44.75882584; 925200,37.07302936; 928800,
              25.04319434; 932400,9.799046246; 936000,0; 939600,0; 943200,0;
              946800,0; 950400,0; 954000,0; 957600,0; 961200,0; 964800,0;
              968400,0; 972000,0; 975600,12.09644867; 979200,27.25095027;
              982800,38.27385402; 986400,45.82451; 990000,51.57958573; 993600,
              54.91750562; 997200,55.89716878; 1000800,54.57011009; 1004400,
              50.87097992; 1008000,44.75882584; 1011600,37.07302936; 1015200,
              25.04319434; 1018800,9.799046246; 1022400,0; 1026000,0; 1029600,0;
              1033200,0; 1036800,0; 1040400,0; 1044000,0; 1047600,0; 1051200,0;
              1054800,0; 1058400,0; 1062000,12.09644867; 1065600,27.25095027;
              1069200,38.27385402; 1072800,45.82451; 1076400,51.57958573;
              1080000,54.91750562; 1083600,55.89716878; 1087200,54.57011009;
              1090800,50.87097992; 1094400,44.75882584; 1098000,37.07302936;
              1101600,25.04319434; 1105200,9.799046246; 1108800,0; 1112400,0;
              1116000,0; 1119600,0; 1123200,0; 1126800,0; 1130400,0; 1134000,0;
              1137600,0; 1141200,0; 1144800,0; 1148400,12.09644867; 1152000,
              27.25095027; 1155600,38.27385402; 1159200,45.82451; 1162800,
              51.57958573; 1166400,54.91750562; 1170000,55.89716878; 1173600,
              54.57011009; 1177200,50.87097992; 1180800,44.75882584; 1184400,
              37.07302936; 1188000,25.04319434; 1191600,9.799046246; 1195200,0;
              1198800,0; 1202400,0; 1206000,0; 1209600,0; 1213200,0; 1216800,0;
              1220400,0; 1224000,0; 1227600,0; 1231200,0; 1234800,23.26734519;
              1238400,47.10145092; 1242000,63.56086174; 1245600,75.24935853;
              1249200,84.63765841; 1252800,90.25746571; 1256400,91.93142037;
              1260000,89.6665522; 1263600,83.46122296; 1267200,73.55208829;
              1270800,61.78366585; 1274400,43.73975719; 1278000,19.33526916;
              1281600,0; 1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0;
              1303200,0; 1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,
              23.26734519; 1324800,47.10145092; 1328400,63.56086174; 1332000,
              75.24935853; 1335600,84.63765841; 1339200,90.25746571; 1342800,
              91.93142037; 1346400,89.6665522; 1350000,83.46122296; 1353600,
              73.55208829; 1357200,61.78366585; 1360800,43.73975719; 1364400,
              19.33526916; 1368000,0; 1371600,0; 1375200,0; 1378800,0; 1382400,
              0; 1386000,0; 1389600,0; 1393200,0; 1396800,0; 1400400,0; 1404000,
              0; 1407600,23.26734519; 1411200,47.10145092; 1414800,63.56086174;
              1418400,75.24935853; 1422000,84.63765841; 1425600,90.25746571;
              1429200,91.93142037; 1432800,89.6665522; 1436400,83.46122296;
              1440000,73.55208829; 1443600,61.78366585; 1447200,43.73975719;
              1450800,19.33526916; 1454400,0; 1458000,0; 1461600,0; 1465200,0;
              1468800,0; 1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0;
              1490400,0; 1494000,23.26734519; 1497600,47.10145092; 1501200,
              63.56086174; 1504800,75.24935853; 1508400,84.63765841; 1512000,
              90.25746571; 1515600,91.93142037; 1519200,89.6665522; 1522800,
              83.46122296; 1526400,73.55208829; 1530000,61.78366585; 1533600,
              43.73975719; 1537200,19.33526916; 1540800,0; 1544400,0; 1548000,0;
              1551600,0; 1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0;
              1573200,0; 1576800,0; 1580400,23.26734519; 1584000,47.10145092;
              1587600,63.56086174; 1591200,75.24935853; 1594800,84.63765841;
              1598400,90.25746571; 1602000,91.93142037; 1605600,89.6665522;
              1609200,83.46122296; 1612800,73.55208829; 1616400,61.78366585;
              1620000,43.73975719; 1623600,19.33526916; 1627200,0; 1630800,0;
              1634400,0; 1638000,0])
          "Diffuse irradiation on tilted surface at clear sky"
          annotation (Placement(transformation(extent={{-96,16},{-82,30}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifTilCov(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifTilCov",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              9.603148499; 28800,27.44337611; 32400,48.23236103; 36000,
              67.90455038; 39600,83.61332758; 43200,93.50445335; 46800,
              96.54693507; 50400,92.4422546; 54000,81.59516725; 57600,
              65.15276159; 61200,45.10461005; 64800,24.46735599; 68400,
              7.478626154; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              9.603148499; 115200,27.44337611; 118800,48.23236103; 122400,
              67.90455038; 126000,83.61332758; 129600,93.50445335; 133200,
              96.54693507; 136800,92.4422546; 140400,81.59516725; 144000,
              65.15276159; 147600,45.10461005; 151200,24.46735599; 154800,
              7.478626154; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,9.603148499; 201600,27.44337611; 205200,48.23236103;
              208800,67.90455038; 212400,83.61332758; 216000,93.50445335;
              219600,96.54693507; 223200,92.4422546; 226800,81.59516725; 230400,
              65.15276159; 234000,45.10461005; 237600,24.46735599; 241200,
              7.478626154; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,9.603148499; 288000,27.44337611; 291600,48.23236103;
              295200,67.90455038; 298800,83.61332758; 302400,93.50445335;
              306000,96.54693507; 309600,92.4422546; 313200,81.59516725; 316800,
              65.15276159; 320400,45.10461005; 324000,24.46735599; 327600,
              7.478626154; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,9.603148499; 374400,27.44337611; 378000,48.23236103;
              381600,67.90455038; 385200,83.61332758; 388800,93.50445335;
              392400,96.54693507; 396000,92.4422546; 399600,81.59516725; 403200,
              65.15276159; 406800,45.10461005; 410400,24.46735599; 414000,
              7.478626154; 417600,0; 421200,0; 424800,0; 428400,0; 432000,0;
              435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,9.603148499; 460800,27.44337611; 464400,48.23236103;
              468000,67.90455038; 471600,83.61332758; 475200,93.50445335;
              478800,96.54693507; 482400,92.4422546; 486000,81.59516725; 489600,
              65.15276159; 493200,45.10461005; 496800,24.46735599; 500400,
              7.478626154; 504000,0; 507600,0; 511200,0; 514800,0; 518400,0;
              522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,9.603148499; 547200,27.44337611; 550800,48.23236103;
              554400,67.90455038; 558000,83.61332758; 561600,93.50445335;
              565200,96.54693507; 568800,92.4422546; 572400,81.59516725; 576000,
              65.15276159; 579600,45.10461005; 583200,24.46735599; 586800,
              7.478626154; 590400,0; 594000,0; 597600,0; 601200,0; 604800,0;
              608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,9.603148499; 633600,27.44337611; 637200,48.23236103;
              640800,67.90455038; 644400,83.61332758; 648000,93.50445335;
              651600,96.54693507; 655200,92.4422546; 658800,81.59516725; 662400,
              65.15276159; 666000,45.10461005; 669600,24.46735599; 673200,
              7.478626154; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,9.603148499; 720000,27.44337611; 723600,48.23236103;
              727200,67.90455038; 730800,83.61332758; 734400,93.50445335;
              738000,96.54693507; 741600,92.4422546; 745200,81.59516725; 748800,
              65.15276159; 752400,45.10461005; 756000,24.46735599; 759600,
              7.478626154; 763200,0; 766800,0; 770400,0; 774000,0; 777600,0;
              781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,9.603148499; 806400,27.44337611; 810000,48.23236103;
              813600,67.90455038; 817200,83.61332758; 820800,93.50445335;
              824400,96.54693507; 828000,92.4422546; 831600,81.59516725; 835200,
              65.15276159; 838800,45.10461005; 842400,24.46735599; 846000,
              7.478626154; 849600,0; 853200,0; 856800,0; 860400,0; 864000,0;
              867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,9.603148499; 892800,27.44337611; 896400,48.23236103;
              900000,67.90455038; 903600,83.61332758; 907200,93.50445335;
              910800,96.54693507; 914400,92.4422546; 918000,81.59516725; 921600,
              65.15276159; 925200,45.10461005; 928800,24.46735599; 932400,
              7.478626154; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,9.603148499; 979200,27.44337611; 982800,48.23236103;
              986400,67.90455038; 990000,83.61332758; 993600,93.50445335;
              997200,96.54693507; 1000800,92.4422546; 1004400,81.59516725;
              1008000,65.15276159; 1011600,45.10461005; 1015200,24.46735599;
              1018800,7.478626154; 1022400,0; 1026000,0; 1029600,0; 1033200,0;
              1036800,0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0;
              1058400,0; 1062000,9.603148499; 1065600,27.44337611; 1069200,
              48.23236103; 1072800,67.90455038; 1076400,83.61332758; 1080000,
              93.50445335; 1083600,96.54693507; 1087200,92.4422546; 1090800,
              81.59516725; 1094400,65.15276159; 1098000,45.10461005; 1101600,
              24.46735599; 1105200,7.478626154; 1108800,0; 1112400,0; 1116000,0;
              1119600,0; 1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0;
              1141200,0; 1144800,0; 1148400,9.603148499; 1152000,27.44337611;
              1155600,48.23236103; 1159200,67.90455038; 1162800,83.61332758;
              1166400,93.50445335; 1170000,96.54693507; 1173600,92.4422546;
              1177200,81.59516725; 1180800,65.15276159; 1184400,45.10461005;
              1188000,24.46735599; 1191600,7.478626154; 1195200,0; 1198800,0;
              1202400,0; 1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0;
              1224000,0; 1227600,0; 1231200,0; 1234800,2.32324972; 1238400,
              6.600661871; 1242000,11.36216174; 1245600,15.75510219; 1249200,
              19.21671663; 1252800,21.38153013; 1256400,22.04553048; 1260000,
              21.14951495; 1263600,18.7737413; 1267200,15.14507191; 1270800,
              10.65548364; 1274400,5.90307651; 1278000,1.797616109; 1281600,0;
              1285200,0; 1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0;
              1306800,0; 1310400,0; 1314000,0; 1317600,0; 1321200,2.32324972;
              1324800,6.600661871; 1328400,11.36216174; 1332000,15.75510219;
              1335600,19.21671663; 1339200,21.38153013; 1342800,22.04553048;
              1346400,21.14951495; 1350000,18.7737413; 1353600,15.14507191;
              1357200,10.65548364; 1360800,5.90307651; 1364400,1.797616109;
              1368000,0; 1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0;
              1389600,0; 1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,
              2.32324972; 1411200,6.600661871; 1414800,11.36216174; 1418400,
              15.75510219; 1422000,19.21671663; 1425600,21.38153013; 1429200,
              22.04553048; 1432800,21.14951495; 1436400,18.7737413; 1440000,
              15.14507191; 1443600,10.65548364; 1447200,5.90307651; 1450800,
              1.797616109; 1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,
              0; 1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,
              0; 1494000,2.32324972; 1497600,6.600661871; 1501200,11.36216174;
              1504800,15.75510219; 1508400,19.21671663; 1512000,21.38153013;
              1515600,22.04553048; 1519200,21.14951495; 1522800,18.7737413;
              1526400,15.14507191; 1530000,10.65548364; 1533600,5.90307651;
              1537200,1.797616109; 1540800,0; 1544400,0; 1548000,0; 1551600,0;
              1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0;
              1576800,0; 1580400,2.32324972; 1584000,6.600661871; 1587600,
              11.36216174; 1591200,15.75510219; 1594800,19.21671663; 1598400,
              21.38153013; 1602000,22.04553048; 1605600,21.14951495; 1609200,
              18.7737413; 1612800,15.14507191; 1616400,10.65548364; 1620000,
              5.90307651; 1623600,1.797616109; 1627200,0; 1630800,0; 1634400,0;
              1638000,0]) "Diffuse irradiation on tilted surface"
          annotation (Placement(transformation(extent={{-96,-6},{-82,8}})));
        Modelica.Blocks.Math.Add HDifTil
          "Diffuse irradiation on a tilted surface"
          annotation (Placement(transformation(extent={{-32,8},{-22,18}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifHorCov(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifHorCov",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              24.22614885; 28800,69.2322226; 32400,121.6772143; 36000,
              171.3048325; 39600,210.9338328; 43200,235.8864705; 46800,
              243.5618298; 50400,233.20683; 54000,205.8425596; 57600,
              164.3628129; 61200,113.7867437; 64800,61.72452797; 68400,
              18.86655303; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              24.22614885; 115200,69.2322226; 118800,121.6772143; 122400,
              171.3048325; 126000,210.9338328; 129600,235.8864705; 133200,
              243.5618298; 136800,233.20683; 140400,205.8425596; 144000,
              164.3628129; 147600,113.7867437; 151200,61.72452797; 154800,
              18.86655303; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,24.22614885; 201600,69.2322226; 205200,121.6772143; 208800,
              171.3048325; 212400,210.9338328; 216000,235.8864705; 219600,
              243.5618298; 223200,233.20683; 226800,205.8425596; 230400,
              164.3628129; 234000,113.7867437; 237600,61.72452797; 241200,
              18.86655303; 244800,0; 248400,0; 252000,0; 255600,0; 259200,0;
              262800,0; 266400,0; 270000,0; 273600,0; 277200,0; 280800,0;
              284400,24.22614885; 288000,69.2322226; 291600,121.6772143; 295200,
              171.3048325; 298800,210.9338328; 302400,235.8864705; 306000,
              243.5618298; 309600,233.20683; 313200,205.8425596; 316800,
              164.3628129; 320400,113.7867437; 324000,61.72452797; 327600,
              18.86655303; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,24.22614885; 374400,69.2322226; 378000,121.6772143; 381600,
              171.3048325; 385200,210.9338328; 388800,235.8864705; 392400,
              243.5618298; 396000,233.20683; 399600,205.8425596; 403200,
              164.3628129; 406800,113.7867437; 410400,61.72452797; 414000,
              18.86655303; 417600,0; 421200,0; 424800,0; 428400,0; 432000,0;
              435600,0; 439200,0; 442800,0; 446400,0; 450000,0; 453600,0;
              457200,24.22614885; 460800,69.2322226; 464400,121.6772143; 468000,
              171.3048325; 471600,210.9338328; 475200,235.8864705; 478800,
              243.5618298; 482400,233.20683; 486000,205.8425596; 489600,
              164.3628129; 493200,113.7867437; 496800,61.72452797; 500400,
              18.86655303; 504000,0; 507600,0; 511200,0; 514800,0; 518400,0;
              522000,0; 525600,0; 529200,0; 532800,0; 536400,0; 540000,0;
              543600,24.22614885; 547200,69.2322226; 550800,121.6772143; 554400,
              171.3048325; 558000,210.9338328; 561600,235.8864705; 565200,
              243.5618298; 568800,233.20683; 572400,205.8425596; 576000,
              164.3628129; 579600,113.7867437; 583200,61.72452797; 586800,
              18.86655303; 590400,0; 594000,0; 597600,0; 601200,0; 604800,0;
              608400,0; 612000,0; 615600,0; 619200,0; 622800,0; 626400,0;
              630000,24.22614885; 633600,69.2322226; 637200,121.6772143; 640800,
              171.3048325; 644400,210.9338328; 648000,235.8864705; 651600,
              243.5618298; 655200,233.20683; 658800,205.8425596; 662400,
              164.3628129; 666000,113.7867437; 669600,61.72452797; 673200,
              18.86655303; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,24.22614885; 720000,69.2322226; 723600,121.6772143; 727200,
              171.3048325; 730800,210.9338328; 734400,235.8864705; 738000,
              243.5618298; 741600,233.20683; 745200,205.8425596; 748800,
              164.3628129; 752400,113.7867437; 756000,61.72452797; 759600,
              18.86655303; 763200,0; 766800,0; 770400,0; 774000,0; 777600,0;
              781200,0; 784800,0; 788400,0; 792000,0; 795600,0; 799200,0;
              802800,24.22614885; 806400,69.2322226; 810000,121.6772143; 813600,
              171.3048325; 817200,210.9338328; 820800,235.8864705; 824400,
              243.5618298; 828000,233.20683; 831600,205.8425596; 835200,
              164.3628129; 838800,113.7867437; 842400,61.72452797; 846000,
              18.86655303; 849600,0; 853200,0; 856800,0; 860400,0; 864000,0;
              867600,0; 871200,0; 874800,0; 878400,0; 882000,0; 885600,0;
              889200,24.22614885; 892800,69.2322226; 896400,121.6772143; 900000,
              171.3048325; 903600,210.9338328; 907200,235.8864705; 910800,
              243.5618298; 914400,233.20683; 918000,205.8425596; 921600,
              164.3628129; 925200,113.7867437; 928800,61.72452797; 932400,
              18.86655303; 936000,0; 939600,0; 943200,0; 946800,0; 950400,0;
              954000,0; 957600,0; 961200,0; 964800,0; 968400,0; 972000,0;
              975600,24.22614885; 979200,69.2322226; 982800,121.6772143; 986400,
              171.3048325; 990000,210.9338328; 993600,235.8864705; 997200,
              243.5618298; 1000800,233.20683; 1004400,205.8425596; 1008000,
              164.3628129; 1011600,113.7867437; 1015200,61.72452797; 1018800,
              18.86655303; 1022400,0; 1026000,0; 1029600,0; 1033200,0; 1036800,
              0; 1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0; 1058400,
              0; 1062000,24.22614885; 1065600,69.2322226; 1069200,121.6772143;
              1072800,171.3048325; 1076400,210.9338328; 1080000,235.8864705;
              1083600,243.5618298; 1087200,233.20683; 1090800,205.8425596;
              1094400,164.3628129; 1098000,113.7867437; 1101600,61.72452797;
              1105200,18.86655303; 1108800,0; 1112400,0; 1116000,0; 1119600,0;
              1123200,0; 1126800,0; 1130400,0; 1134000,0; 1137600,0; 1141200,0;
              1144800,0; 1148400,24.22614885; 1152000,69.2322226; 1155600,
              121.6772143; 1159200,171.3048325; 1162800,210.9338328; 1166400,
              235.8864705; 1170000,243.5618298; 1173600,233.20683; 1177200,
              205.8425596; 1180800,164.3628129; 1184400,113.7867437; 1188000,
              61.72452797; 1191600,18.86655303; 1195200,0; 1198800,0; 1202400,0;
              1206000,0; 1209600,0; 1213200,0; 1216800,0; 1220400,0; 1224000,0;
              1227600,0; 1231200,0; 1234800,5.860931292; 1238400,16.65168637;
              1242000,28.66366397; 1245600,39.74586572; 1249200,48.47858361;
              1252800,53.93982313; 1256400,55.6149166; 1260000,53.35451152;
              1263600,47.36107656; 1267200,38.20692416; 1270800,26.88090607;
              1274400,14.891867; 1278000,4.534899719; 1281600,0; 1285200,0;
              1288800,0; 1292400,0; 1296000,0; 1299600,0; 1303200,0; 1306800,0;
              1310400,0; 1314000,0; 1317600,0; 1321200,5.860931292; 1324800,
              16.65168637; 1328400,28.66366397; 1332000,39.74586572; 1335600,
              48.47858361; 1339200,53.93982313; 1342800,55.6149166; 1346400,
              53.35451152; 1350000,47.36107656; 1353600,38.20692416; 1357200,
              26.88090607; 1360800,14.891867; 1364400,4.534899719; 1368000,0;
              1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0; 1389600,0;
              1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,5.860931292;
              1411200,16.65168637; 1414800,28.66366397; 1418400,39.74586572;
              1422000,48.47858361; 1425600,53.93982313; 1429200,55.6149166;
              1432800,53.35451152; 1436400,47.36107656; 1440000,38.20692416;
              1443600,26.88090607; 1447200,14.891867; 1450800,4.534899719;
              1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,0; 1472400,0;
              1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,0; 1494000,
              5.860931292; 1497600,16.65168637; 1501200,28.66366397; 1504800,
              39.74586572; 1508400,48.47858361; 1512000,53.93982313; 1515600,
              55.6149166; 1519200,53.35451152; 1522800,47.36107656; 1526400,
              38.20692416; 1530000,26.88090607; 1533600,14.891867; 1537200,
              4.534899719; 1540800,0; 1544400,0; 1548000,0; 1551600,0; 1555200,
              0; 1558800,0; 1562400,0; 1566000,0; 1569600,0; 1573200,0; 1576800,
              0; 1580400,5.860931292; 1584000,16.65168637; 1587600,28.66366397;
              1591200,39.74586572; 1594800,48.47858361; 1598400,53.93982313;
              1602000,55.6149166; 1605600,53.35451152; 1609200,47.36107656;
              1612800,38.20692416; 1616400,26.88090607; 1620000,14.891867;
              1623600,4.534899719; 1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Diffuse irradiation on a horizontal surface at covered sky"
          annotation (Placement(transformation(extent={{-96,56},{-82,70}})));
        Modelica.Blocks.Math.Add HDifHor
          "Diffuse irradiation on a horizontal surface"
          annotation (Placement(transformation(extent={{-32,32},{-22,42}})));
        Modelica.Blocks.Sources.CombiTimeTable HDifHorCle(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HDifHorCle",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              21.04231965; 28800,38.76402419; 32400,49.79122826; 36000,
              56.90519903; 39600,61.37224652; 43200,63.79161511; 46800,
              64.48257797; 50400,63.54466077; 54000,60.84410045; 57600,
              56.02824899; 61200,48.43452029; 64800,36.62912553; 68400,
              17.58508841; 72000,0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0;
              93600,0; 97200,0; 100800,0; 104400,0; 108000,0; 111600,
              21.04231965; 115200,38.76402419; 118800,49.79122826; 122400,
              56.90519903; 126000,61.37224652; 129600,63.79161511; 133200,
              64.48257797; 136800,63.54466077; 140400,60.84410045; 144000,
              56.02824899; 147600,48.43452029; 151200,36.62912553; 154800,
              17.58508841; 158400,0; 162000,0; 165600,0; 169200,0; 172800,0;
              176400,0; 180000,0; 183600,0; 187200,0; 190800,0; 194400,0;
              198000,21.04231965; 201600,38.76402419; 205200,49.79122826;
              208800,56.90519903; 212400,61.37224652; 216000,63.79161511;
              219600,64.48257797; 223200,63.54466077; 226800,60.84410045;
              230400,56.02824899; 234000,48.43452029; 237600,36.62912553;
              241200,17.58508841; 244800,0; 248400,0; 252000,0; 255600,0;
              259200,0; 262800,0; 266400,0; 270000,0; 273600,0; 277200,0;
              280800,0; 284400,21.04231965; 288000,38.76402419; 291600,
              49.79122826; 295200,56.90519903; 298800,61.37224652; 302400,
              63.79161511; 306000,64.48257797; 309600,63.54466077; 313200,
              60.84410045; 316800,56.02824899; 320400,48.43452029; 324000,
              36.62912553; 327600,17.58508841; 331200,0; 334800,0; 338400,0;
              342000,0; 345600,0; 349200,0; 352800,0; 356400,0; 360000,0;
              363600,0; 367200,0; 370800,21.04231965; 374400,38.76402419;
              378000,49.79122826; 381600,56.90519903; 385200,61.37224652;
              388800,63.79161511; 392400,64.48257797; 396000,63.54466077;
              399600,60.84410045; 403200,56.02824899; 406800,48.43452029;
              410400,36.62912553; 414000,17.58508841; 417600,0; 421200,0;
              424800,0; 428400,0; 432000,0; 435600,0; 439200,0; 442800,0;
              446400,0; 450000,0; 453600,0; 457200,21.04231965; 460800,
              38.76402419; 464400,49.79122826; 468000,56.90519903; 471600,
              61.37224652; 475200,63.79161511; 478800,64.48257797; 482400,
              63.54466077; 486000,60.84410045; 489600,56.02824899; 493200,
              48.43452029; 496800,36.62912553; 500400,17.58508841; 504000,0;
              507600,0; 511200,0; 514800,0; 518400,0; 522000,0; 525600,0;
              529200,0; 532800,0; 536400,0; 540000,0; 543600,21.04231965;
              547200,38.76402419; 550800,49.79122826; 554400,56.90519903;
              558000,61.37224652; 561600,63.79161511; 565200,64.48257797;
              568800,63.54466077; 572400,60.84410045; 576000,56.02824899;
              579600,48.43452029; 583200,36.62912553; 586800,17.58508841;
              590400,0; 594000,0; 597600,0; 601200,0; 604800,0; 608400,0;
              612000,0; 615600,0; 619200,0; 622800,0; 626400,0; 630000,
              21.04231965; 633600,38.76402419; 637200,49.79122826; 640800,
              56.90519903; 644400,61.37224652; 648000,63.79161511; 651600,
              64.48257797; 655200,63.54466077; 658800,60.84410045; 662400,
              56.02824899; 666000,48.43452029; 669600,36.62912553; 673200,
              17.58508841; 676800,0; 680400,0; 684000,0; 687600,0; 691200,0;
              694800,0; 698400,0; 702000,0; 705600,0; 709200,0; 712800,0;
              716400,21.04231965; 720000,38.76402419; 723600,49.79122826;
              727200,56.90519903; 730800,61.37224652; 734400,63.79161511;
              738000,64.48257797; 741600,63.54466077; 745200,60.84410045;
              748800,56.02824899; 752400,48.43452029; 756000,36.62912553;
              759600,17.58508841; 763200,0; 766800,0; 770400,0; 774000,0;
              777600,0; 781200,0; 784800,0; 788400,0; 792000,0; 795600,0;
              799200,0; 802800,21.04231965; 806400,38.76402419; 810000,
              49.79122826; 813600,56.90519903; 817200,61.37224652; 820800,
              63.79161511; 824400,64.48257797; 828000,63.54466077; 831600,
              60.84410045; 835200,56.02824899; 838800,48.43452029; 842400,
              36.62912553; 846000,17.58508841; 849600,0; 853200,0; 856800,0;
              860400,0; 864000,0; 867600,0; 871200,0; 874800,0; 878400,0;
              882000,0; 885600,0; 889200,21.04231965; 892800,38.76402419;
              896400,49.79122826; 900000,56.90519903; 903600,61.37224652;
              907200,63.79161511; 910800,64.48257797; 914400,63.54466077;
              918000,60.84410045; 921600,56.02824899; 925200,48.43452029;
              928800,36.62912553; 932400,17.58508841; 936000,0; 939600,0;
              943200,0; 946800,0; 950400,0; 954000,0; 957600,0; 961200,0;
              964800,0; 968400,0; 972000,0; 975600,21.04231965; 979200,
              38.76402419; 982800,49.79122826; 986400,56.90519903; 990000,
              61.37224652; 993600,63.79161511; 997200,64.48257797; 1000800,
              63.54466077; 1004400,60.84410045; 1008000,56.02824899; 1011600,
              48.43452029; 1015200,36.62912553; 1018800,17.58508841; 1022400,0;
              1026000,0; 1029600,0; 1033200,0; 1036800,0; 1040400,0; 1044000,0;
              1047600,0; 1051200,0; 1054800,0; 1058400,0; 1062000,21.04231965;
              1065600,38.76402419; 1069200,49.79122826; 1072800,56.90519903;
              1076400,61.37224652; 1080000,63.79161511; 1083600,64.48257797;
              1087200,63.54466077; 1090800,60.84410045; 1094400,56.02824899;
              1098000,48.43452029; 1101600,36.62912553; 1105200,17.58508841;
              1108800,0; 1112400,0; 1116000,0; 1119600,0; 1123200,0; 1126800,0;
              1130400,0; 1134000,0; 1137600,0; 1141200,0; 1144800,0; 1148400,
              21.04231965; 1152000,38.76402419; 1155600,49.79122826; 1159200,
              56.90519903; 1162800,61.37224652; 1166400,63.79161511; 1170000,
              64.48257797; 1173600,63.54466077; 1177200,60.84410045; 1180800,
              56.02824899; 1184400,48.43452029; 1188000,36.62912553; 1191600,
              17.58508841; 1195200,0; 1198800,0; 1202400,0; 1206000,0; 1209600,
              0; 1213200,0; 1216800,0; 1220400,0; 1224000,0; 1227600,0; 1231200,
              0; 1234800,40.47459946; 1238400,67.00103172; 1242000,82.68760638;
              1245600,93.44518302; 1249200,100.7065715; 1252800,104.8421528;
              1256400,106.0514354; 1260000,104.4130318; 1263600,99.82357411;
              1267200,92.07110863; 1270800,80.71803866; 1274400,63.97542723;
              1278000,34.69852159; 1281600,0; 1285200,0; 1288800,0; 1292400,0;
              1296000,0; 1299600,0; 1303200,0; 1306800,0; 1310400,0; 1314000,0;
              1317600,0; 1321200,40.47459946; 1324800,67.00103172; 1328400,
              82.68760638; 1332000,93.44518302; 1335600,100.7065715; 1339200,
              104.8421528; 1342800,106.0514354; 1346400,104.4130318; 1350000,
              99.82357411; 1353600,92.07110863; 1357200,80.71803866; 1360800,
              63.97542723; 1364400,34.69852159; 1368000,0; 1371600,0; 1375200,0;
              1378800,0; 1382400,0; 1386000,0; 1389600,0; 1393200,0; 1396800,0;
              1400400,0; 1404000,0; 1407600,40.47459946; 1411200,67.00103172;
              1414800,82.68760638; 1418400,93.44518302; 1422000,100.7065715;
              1425600,104.8421528; 1429200,106.0514354; 1432800,104.4130318;
              1436400,99.82357411; 1440000,92.07110863; 1443600,80.71803866;
              1447200,63.97542723; 1450800,34.69852159; 1454400,0; 1458000,0;
              1461600,0; 1465200,0; 1468800,0; 1472400,0; 1476000,0; 1479600,0;
              1483200,0; 1486800,0; 1490400,0; 1494000,40.47459946; 1497600,
              67.00103172; 1501200,82.68760638; 1504800,93.44518302; 1508400,
              100.7065715; 1512000,104.8421528; 1515600,106.0514354; 1519200,
              104.4130318; 1522800,99.82357411; 1526400,92.07110863; 1530000,
              80.71803866; 1533600,63.97542723; 1537200,34.69852159; 1540800,0;
              1544400,0; 1548000,0; 1551600,0; 1555200,0; 1558800,0; 1562400,0;
              1566000,0; 1569600,0; 1573200,0; 1576800,0; 1580400,40.47459946;
              1584000,67.00103172; 1587600,82.68760638; 1591200,93.44518302;
              1594800,100.7065715; 1598400,104.8421528; 1602000,106.0514354;
              1605600,104.4130318; 1609200,99.82357411; 1612800,92.07110863;
              1616400,80.71803866; 1620000,63.97542723; 1623600,34.69852159;
              1627200,0; 1630800,0; 1634400,0; 1638000,0])
          "Diffuse irradiation on a horizontal surface at clear sky"
          annotation (Placement(transformation(extent={{-96,36},{-82,50}})));
        BoundaryConditions.SolarGeometry.BaseClasses.IncidenceAngle
          incAng(
          azi(displayUnit="deg") = 0,
          til(displayUnit="deg") = 1.5707963267949,
          lat=0.86393797973719) "Incidence angle for the window"
          annotation (Placement(transformation(extent={{-56,-40},{-44,-28}})));
        BoundaryConditions.SolarIrradiation.BaseClasses.DirectTiltedSurface
          HDirTil "Direct solar irradiation on the window"
          annotation (Placement(transformation(extent={{-30,-26},{-14,-10}})));
        Windows.BaseClasses.Conversions.to_HDirNor to_HDirNor
          "Converts the direct irradiation onto a horizontal surface to direct irradiation on a normal surface"
          annotation (Placement(transformation(extent={{-58,-8},{-46,4}})));
        Windows.Validation.BaseClasses.SolarHourAngleVDI6007 SolHouAng(lon(displayUnit="deg") = 0.15009831567151)
          "Solar hour angle based on calculation of VDI 6007"
          annotation (Placement(transformation(extent={{-64,-40},{-58,-34}})));
        Windows.Validation.BaseClasses.SolarDeclinationAngleVDI6007 decAng
          "Declination angle  based on the calculations of VDI6007"
          annotation (Placement(transformation(extent={{-64,-34},{-58,-28}})));
        Modelica.Blocks.Sources.CombiTimeTable HVen_VDI2078(
          columns={2},
          smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          tableName="HVen",
          tableOnFile=false,
          table=[0,0; 3600,0; 7200,0; 10800,0; 14400,0; 18000,0; 21600,0; 25200,
              0; 28800,0; 32400,0; 36000,8.775122167; 39600,11.77075859; 43200,
              13.80359492; 46800,14.45097624; 50400,13.58000818; 54000,
              11.36966839; 57600,8.280950332; 61200,0; 64800,0; 68400,0; 72000,
              0; 75600,0; 79200,0; 82800,0; 86400,0; 90000,0; 93600,0; 97200,0;
              100800,0; 104400,0; 108000,0; 111600,0; 115200,0; 118800,0;
              122400,8.775122167; 126000,11.77075859; 129600,13.80359492;
              133200,14.45097624; 136800,13.58000818; 140400,11.36966839;
              144000,8.280950332; 147600,0; 151200,0; 154800,0; 158400,0;
              162000,0; 165600,0; 169200,0; 172800,0; 176400,0; 180000,0;
              183600,0; 187200,0; 190800,0; 194400,0; 198000,0; 201600,0;
              205200,0; 208800,8.775122167; 212400,11.77075859; 216000,
              13.80359492; 219600,14.45097624; 223200,13.58000818; 226800,
              11.36966839; 230400,8.280950332; 234000,0; 237600,0; 241200,0;
              244800,0; 248400,0; 252000,0; 255600,0; 259200,0; 262800,0;
              266400,0; 270000,0; 273600,0; 277200,0; 280800,0; 284400,0;
              288000,0; 291600,0; 295200,8.775122167; 298800,11.77075859;
              302400,13.80359492; 306000,14.45097624; 309600,13.58000818;
              313200,11.36966839; 316800,8.280950332; 320400,0; 324000,0;
              327600,0; 331200,0; 334800,0; 338400,0; 342000,0; 345600,0;
              349200,0; 352800,0; 356400,0; 360000,0; 363600,0; 367200,0;
              370800,0; 374400,0; 378000,0; 381600,8.775122167; 385200,
              11.77075859; 388800,13.80359492; 392400,14.45097624; 396000,
              13.58000818; 399600,11.36966839; 403200,8.280950332; 406800,0;
              410400,0; 414000,0; 417600,0; 421200,0; 424800,0; 428400,0;
              432000,0; 435600,0; 439200,0; 442800,0; 446400,0; 450000,0;
              453600,0; 457200,0; 460800,0; 464400,0; 468000,0; 471600,0;
              475200,0; 478800,0; 482400,0; 486000,0; 489600,0; 493200,0;
              496800,0; 500400,0; 504000,0; 507600,0; 511200,0; 514800,0;
              518400,0; 522000,0; 525600,0; 529200,0; 532800,0; 536400,0;
              540000,0; 543600,0; 547200,0; 550800,0; 554400,0; 558000,0;
              561600,0; 565200,0; 568800,0; 572400,0; 576000,0; 579600,0;
              583200,0; 586800,0; 590400,0; 594000,0; 597600,0; 601200,0;
              604800,0; 608400,0; 612000,0; 615600,0; 619200,0; 622800,0;
              626400,0; 630000,0; 633600,0; 637200,0; 640800,8.775122167;
              644400,11.77075859; 648000,13.80359492; 651600,14.45097624;
              655200,13.58000818; 658800,11.36966839; 662400,8.280950332;
              666000,0; 669600,0; 673200,0; 676800,0; 680400,0; 684000,0;
              687600,0; 691200,0; 694800,0; 698400,0; 702000,0; 705600,0;
              709200,0; 712800,0; 716400,0; 720000,0; 723600,0; 727200,
              8.775122167; 730800,11.77075859; 734400,13.80359492; 738000,
              14.45097624; 741600,13.58000818; 745200,11.36966839; 748800,
              8.280950332; 752400,0; 756000,0; 759600,0; 763200,0; 766800,0;
              770400,0; 774000,0; 777600,0; 781200,0; 784800,0; 788400,0;
              792000,0; 795600,0; 799200,0; 802800,0; 806400,0; 810000,0;
              813600,8.775122167; 817200,11.77075859; 820800,13.80359492;
              824400,14.45097624; 828000,13.58000818; 831600,11.36966839;
              835200,8.280950332; 838800,0; 842400,0; 846000,0; 849600,0;
              853200,0; 856800,0; 860400,0; 864000,0; 867600,0; 871200,0;
              874800,0; 878400,0; 882000,0; 885600,0; 889200,0; 892800,0;
              896400,0; 900000,8.775122167; 903600,11.77075859; 907200,
              13.80359492; 910800,14.45097624; 914400,13.58000818; 918000,
              11.36966839; 921600,8.280950332; 925200,0; 928800,0; 932400,0;
              936000,0; 939600,0; 943200,0; 946800,0; 950400,0; 954000,0;
              957600,0; 961200,0; 964800,0; 968400,0; 972000,0; 975600,0;
              979200,0; 982800,0; 986400,8.775122167; 990000,11.77075859;
              993600,13.80359492; 997200,14.45097624; 1000800,13.58000818;
              1004400,11.36966839; 1008000,8.280950332; 1011600,0; 1015200,0;
              1018800,0; 1022400,0; 1026000,0; 1029600,0; 1033200,0; 1036800,0;
              1040400,0; 1044000,0; 1047600,0; 1051200,0; 1054800,0; 1058400,0;
              1062000,0; 1065600,0; 1069200,0; 1072800,0; 1076400,0; 1080000,0;
              1083600,0; 1087200,0; 1090800,0; 1094400,0; 1098000,0; 1101600,0;
              1105200,0; 1108800,0; 1112400,0; 1116000,0; 1119600,0; 1123200,0;
              1126800,0; 1130400,0; 1134000,0; 1137600,0; 1141200,0; 1144800,0;
              1148400,0; 1152000,0; 1155600,0; 1159200,0; 1162800,0; 1166400,0;
              1170000,0; 1173600,0; 1177200,0; 1180800,0; 1184400,0; 1188000,0;
              1191600,0; 1195200,0; 1198800,0; 1202400,0; 1206000,0; 1209600,0;
              1213200,0; 1216800,0; 1220400,0; 1224000,0; 1227600,0; 1231200,0;
              1234800,0; 1238400,0; 1242000,9.411054257; 1245600,14.83746656;
              1249200,19.30942141; 1252800,22.15794523; 1256400,23.03789677;
              1260000,21.85111281; 1263600,18.73088748; 1267200,14.06286399;
              1270800,0; 1274400,0; 1278000,0; 1281600,0; 1285200,0; 1288800,0;
              1292400,0; 1296000,0; 1299600,0; 1303200,0; 1306800,0; 1310400,0;
              1314000,0; 1317600,0; 1321200,0; 1324800,0; 1328400,9.411054257;
              1332000,14.83746656; 1335600,19.30942141; 1339200,22.15794523;
              1342800,23.03789677; 1346400,21.85111281; 1350000,18.73088748;
              1353600,14.06286399; 1357200,0; 1360800,0; 1364400,0; 1368000,0;
              1371600,0; 1375200,0; 1378800,0; 1382400,0; 1386000,0; 1389600,0;
              1393200,0; 1396800,0; 1400400,0; 1404000,0; 1407600,0; 1411200,0;
              1414800,9.411054257; 1418400,14.83746656; 1422000,19.30942141;
              1425600,22.15794523; 1429200,23.03789677; 1432800,21.85111281;
              1436400,18.73088748; 1440000,14.06286399; 1443600,0; 1447200,0;
              1450800,0; 1454400,0; 1458000,0; 1461600,0; 1465200,0; 1468800,0;
              1472400,0; 1476000,0; 1479600,0; 1483200,0; 1486800,0; 1490400,0;
              1494000,0; 1497600,0; 1501200,9.411054257; 1504800,14.83746656;
              1508400,19.30942141; 1512000,22.15794523; 1515600,23.03789677;
              1519200,21.85111281; 1522800,18.73088748; 1526400,14.06286399;
              1530000,0; 1533600,0; 1537200,0; 1540800,0; 1544400,0; 1548000,0;
              1551600,0; 1555200,0; 1558800,0; 1562400,0; 1566000,0; 1569600,0;
              1573200,0; 1576800,0; 1580400,0; 1584000,0; 1587600,9.411054257;
              1591200,14.83746656; 1594800,19.30942141; 1598400,22.15794523;
              1602000,23.03789677; 1605600,21.85111281; 1609200,18.73088748;
              1612800,14.06286399; 1616400,0; 1620000,0; 1623600,0; 1627200,0;
              1630800,0; 1634400,0; 1638000,0])
          "Comparison for the heat input due to ventilation with active sunblind"
          annotation (Placement(transformation(extent={{80,-38},{100,-18}})));
      equation
        connect(HDifTilCle.y[1], HDifTil.u1) annotation (Line(points={{-81.3,23},
                {-62,23},{-62,16},{-42,16},{-34,16},{-33,16}}, color={0,0,127}));
        connect(HDifTilCov.y[1], HDifTil.u2) annotation (Line(points={{-81.3,1},
                {-64,1},{-64,10},{-33,10}}, color={0,0,127}));
        connect(HDifTil.y, sunblind.HDifTil) annotation (Line(points={{-21.5,13},
                {16,13},{16,4},{32,4},{32,4.32},{47.6,4.32}}, color={0,0,127}));
        connect(HDifTil.y, ventilationHeat.HDifTil) annotation (Line(points={{-21.5,
                13},{60,13},{60,5},{79,5}},             color={0,0,127}));
        connect(HDifHorCov.y[1], HDifHor.u1) annotation (Line(points={{-81.3,63},
                {-64.75,63},{-64.75,40},{-33,40}}, color={0,0,127}));
        connect(HDifHorCle.y[1], HDifHor.u2) annotation (Line(points={{-81.3,43},
                {-57.75,43},{-57.75,34},{-33,34}}, color={0,0,127}));
        connect(HDifHor.y, ventilationHeat.HDifHor) annotation (Line(points={{-21.5,
                37},{-12,37},{-12,18},{20,18},{66,18},{66,8},{79,8}},
                           color={0,0,127}));
        connect(sunblind.sunscreen, ventilationHeat.sunscreen) annotation (Line(
              points={{56.4,2},{56.4,2},{79,2}}, color={255,0,255}));
        connect(incAng.incAng, HDirTil.incAng) annotation (Line(points={{-43.4,
                -34},{-42,-34},{-42,-22.8},{-31.6,-22.8}}, color={0,0,127}));
        connect(to_HDirNor.HDirNor, HDirTil.HDirNor) annotation (Line(points={{-44.8,
                -2},{-38,-2},{-38,-13.2},{-31.6,-13.2}},         color={0,0,127}));
        connect(HDirHor.y[1], to_HDirNor.HDirHor) annotation (Line(points={{-81.3,
                -19},{-72,-19},{-72,-4},{-66,-4},{-66,0.4},{-59.2,0.4}},
                                                                   color={0,0,
                127}));
        connect(alt.y[1], to_HDirNor.alt) annotation (Line(points={{-81.3,-45},
                {-81.3,-28},{-80,-28},{-76,-28},{-68,-28},{-68,-16},{-68,-6},{
                -62,-6},{-62,-4.4},{-59.2,-4.4}},
                               color={0,0,127}));
        connect(HDirTil.HDirTil, ventilationHeat.HDirTil) annotation (Line(points={{-13.2,
                -18},{64,-18},{64,-2},{79,-2}},        color={0,0,127}));
        connect(HDirTil.HDirTil, sunblind.HDirTil) annotation (Line(points={{-13.2,
                -18},{16,-18},{16,-0.4},{47.6,-0.4}}, color={0,0,127}));
        connect(to_HDirNor.HDirNor, ventilationHeat.HDirNor) annotation (Line(
              points={{-44.8,-2},{4,-2},{4,-32},{68,-32},{68,-4},{72,-4},{76,-4},
                {76,-5},{79,-5}}, color={0,0,127}));
        connect(alt.y[1], ventilationHeat.alt) annotation (Line(points={{-81.3,
                -45},{-80,-45},{-80,-32},{-76,-32},{-76,-56},{6,-56},{72,-56},{
                72,-8},{79,-8}}, color={0,0,127}));
        connect(incAng.decAng, decAng.decAng) annotation (Line(points={{-57.32,
                -30.76},{-57.51,-30.76},{-57.51,-31},{-57.7,-31}}, color={0,0,127}));
        connect(incAng.solHouAng, SolHouAng.solHouAng) annotation (Line(points={{
                -57.2,-36.88},{-57.45,-36.88},{-57.45,-37},{-57.7,-37}}, color={0,
                0,127}));
        annotation (experiment(StartTime=0,StopTime=1638000),Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
<p>This model simulates parts of VDI2078 test case 3. The solar irradiation is treated as an input. To calculate the boundary conditions <a href=\"Annex60.BoundaryConditions\">Annex60</a> models are used.</p>
</html>"));
      end TestCase3_VentilationHeat;
      annotation (Documentation(info="<html>

<p>This package contains parts of VDI2078 test cases 1 and 3. In test case 1 the <a href=\"Windows.BaseClasses.Illumination\">Illumination</a> model is tested. In test case 3 the <a href=\"Windows.BaseClasses.VentilationHeat\">VentilationHeat</a> model is tested.</p>
</html>"));
    end VDI2078;

    package SelfShadowing
      extends Modelica.Icons.ExamplesPackage;
      model SelfShadowingTestAbove
        extends Modelica.Icons.Example;

        Windows.BaseClasses.SelfShadowing selfShadowingAbove(
          final bRig={0},
          final n=1,
          final b={1},
          final h={1},
          final bLef={0},
          final dLef={0},
          final dRig={0},
          final bAbo={1},
          final bBel={0},
          final dAbo={0.01},
          final dBel={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949})
          "Shadowing due to a projection above the window"
                annotation (Placement(transformation(extent={{56,46},{88,74}})));
        BaseClasses.IncidenceAngleVDI6007 incAng1(azi=0, til=90)
          "Incidence angle for the window"
          annotation (Placement(transformation(extent={{-26,40},{-6,60}})));
        Modelica.Blocks.Sources.Constant solAzi(k=0)
          "Constant soar azimuth angle (North)"
          annotation (Placement(transformation(extent={{-88,24},{-68,44}})));
        Modelica.Blocks.Sources.Constant alt(k=Modelica.Constants.pi/6)
          "Constant altitude angle"
          annotation (Placement(transformation(extent={{-88,-20},{-68,0}})));
        Modelica.Blocks.Sources.Sine altSine(freqHz=1, amplitude=Modelica.Constants.pi
              /3) "Altitude angle generated as a sine"
          annotation (Placement(transformation(extent={{-88,56},{-68,76}})));
        Windows.BaseClasses.SelfShadowing selfShadowingAboveSin(
          final n=1,
          final b={1},
          final h={1},
          final bLef={0},
          final bRig={0},
          final dLef={0},
          final dRig={0},
          final bAbo={1},
          final bBel={0},
          final dAbo={0.01},
          final dBel={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949})
          "Shadowing due to a projection above the window"
                annotation (Placement(transformation(extent={{56,-32},{88,-4}})));
        BaseClasses.IncidenceAngleVDI6007 incAng2(azi=0, til=90)
          "Incidence angle for the window"
          annotation (Placement(transformation(extent={{-26,-38},{-6,-18}})));
        Modelica.Blocks.Sources.Sine solAziSine(freqHz=0.25, amplitude=2*Modelica.Constants.pi)
          "Solar azimuth generated as a sine"
          annotation (Placement(transformation(extent={{-88,-52},{-68,-32}})));
      equation
        connect(altSine.y, incAng1.alt) annotation (Line(
            points={{-67,66},{-44,66},{-44,55.4},{-28.2,55.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(altSine.y, selfShadowingAbove.alt) annotation (Line(
            points={{-67,66},{2,66},{2,60},{54.4,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi.y, incAng1.solAzi) annotation (Line(
            points={{-67,34},{-44,34},{-44,45.2},{-28,45.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi.y, selfShadowingAbove.solAzi) annotation (Line(
            points={{-67,34},{38,34},{38,69.8},{54.4,69.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAziSine.y, incAng2.solAzi) annotation (Line(points={{-67,-42},
                {-48,-42},{-48,-32.8},{-28,-32.8}}, color={0,0,127}));
        connect(solAziSine.y, selfShadowingAboveSin.solAzi) annotation (Line(
              points={{-67,-42},{24,-42},{24,-8.2},{54.4,-8.2}}, color={0,0,127}));
        connect(alt.y, incAng2.alt) annotation (Line(points={{-67,-10},{-48,-10},
                {-48,-22.6},{-28.2,-22.6}}, color={0,0,127}));
        connect(alt.y, selfShadowingAboveSin.alt) annotation (Line(points={{-67,
                -10},{-26,-10},{-26,-8},{14,-8},{14,-18},{54.4,-18}}, color={0,0,
                127}));
        connect(incAng1.incAng, selfShadowingAbove.incAng[1]) annotation (Line(
              points={{-5,50},{54.4,50},{54.4,50.2}}, color={0,0,127}));
        connect(incAng2.incAng, selfShadowingAboveSin.incAng[1]) annotation (Line(
              points={{-5,-28},{54.4,-28},{54.4,-27.8}}, color={0,0,127}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}})),
          Documentation(info="<html>
<p>This model simulates a projection above the window.</p>
</html>"));
      end SelfShadowingTestAbove;

      model SelfShadowingTestBelow
        extends Modelica.Icons.Example;

        Windows.BaseClasses.SelfShadowing selfShadowingBelow(
          final bRig={0},
          final n=1,
          final b={1},
          final h={1},
          final bLef={0},
          final dRig={0},
          final dLef={0},
          final bAbo={0},
          final dAbo={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949},
          final dBel={0.01},
          final bBel={1}) "Shadowing due to a projection below"
          annotation (Placement(transformation(extent={{56,46},{88,74}})));
        BaseClasses.IncidenceAngleVDI6007 incAng1(azi=0, til=90)
          "Incidence Angle for the window"
          annotation (Placement(transformation(extent={{-26,40},{-6,60}})));
        Modelica.Blocks.Sources.Constant solAzi(k=0)
          "Constant solar azimuth angle (north)"
          annotation (Placement(transformation(extent={{-88,24},{-68,44}})));
        Modelica.Blocks.Sources.Sine altSine(freqHz=1, amplitude=Modelica.Constants.pi
              /3) "Solar altitude angle generated as a sine"
          annotation (Placement(transformation(extent={{-88,56},{-68,76}})));
        Windows.BaseClasses.SelfShadowing selfShadowingBalkony(
          final bRig={0},
          final n=1,
          final b={1},
          final h={1},
          final bLef={0},
          final dRig={0},
          final dLef={0},
          final bAbo={0},
          final dAbo={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949},
          final bBel={1},
          final dBel={-0.2}) "Shadowing due to a balkony"
          annotation (Placement(transformation(extent={{56,-40},{88,-12}})));
        BaseClasses.IncidenceAngleVDI6007 incAng2(azi=0, til=90)
          "Incidence angle for the window"
          annotation (Placement(transformation(extent={{-26,-46},{-6,-26}})));
        Modelica.Blocks.Sources.Constant solAzi1(k=0)
          "Constant solar azimuth angle (north)"
          annotation (Placement(transformation(extent={{-88,-62},{-68,-42}})));
        Modelica.Blocks.Sources.Sine altSine1(
                                             freqHz=1, amplitude=Modelica.Constants.pi
              /3) "Solar altitude angle generated as a sine"
          annotation (Placement(transformation(extent={{-88,-30},{-68,-10}})));
      equation
        connect(incAng1.incAng, selfShadowingBelow.incAng[1]) annotation (Line(
              points={{-5,50},{28,50},{28,50.2},{54.4,50.2}}, color={0,0,127}));
        connect(altSine.y, incAng1.alt) annotation (Line(
            points={{-67,66},{-44,66},{-44,55.4},{-28.2,55.4}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(altSine.y, selfShadowingBelow.alt) annotation (Line(
            points={{-67,66},{2,66},{2,60},{54.4,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi.y, incAng1.solAzi) annotation (Line(
            points={{-67,34},{-44,34},{-44,45.2},{-28,45.2}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi.y, selfShadowingBelow.solAzi) annotation (Line(
            points={{-67,34},{38,34},{38,69.8},{54.4,69.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(incAng2.incAng, selfShadowingBalkony.incAng[1]) annotation (Line(points={
                {-5,-36},{28,-36},{28,-35.8},{54.4,-35.8}}, color={0,0,127}));
        connect(altSine1.y, incAng2.alt) annotation (Line(
            points={{-67,-20},{-44,-20},{-44,-30.6},{-28.2,-30.6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(altSine1.y, selfShadowingBalkony.alt) annotation (Line(
            points={{-67,-20},{2,-20},{2,-26},{54.4,-26}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi1.y, incAng2.solAzi) annotation (Line(
            points={{-67,-52},{-44,-52},{-44,-40.8},{-28,-40.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAzi1.y, selfShadowingBalkony.solAzi) annotation (Line(
            points={{-67,-52},{38,-52},{38,-16.2},{54.4,-16.2}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}})),
          Documentation(info="<html>
<p>This model simulates a projection below the window and a balcony.</p>
</html>"));
      end SelfShadowingTestBelow;

      model SelfShadowingTestLeft
        extends Modelica.Icons.Example;

        Windows.BaseClasses.SelfShadowing selfShadowingLeft(
          final bRig={0},
          final n=1,
          final b={1},
          final h={1},
          final dRig={0},
          final bAbo={0},
          final bBel={0},
          final dAbo={0},
          final dBel={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949},
          final bLef={1},
          final dLef={0.01})
          "Shadowing due to a projection on the left-hand side of the window"
          annotation (Placement(transformation(extent={{62,-4},{94,24}})));
        BaseClasses.IncidenceAngleVDI6007 incAng1(azi=0, til=90)
          "Incidence angle for the window"
          annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
        Modelica.Blocks.Sources.Constant alt(k=0.3490658504)
          "Constant altitude angle"
          annotation (Placement(transformation(extent={{-76,8},{-56,28}})));
        Modelica.Blocks.Sources.Sine solAziSine(amplitude=Modelica.Constants.pi,
            freqHz=1) "Solar azimuth angle generated as a sine"
          annotation (Placement(transformation(extent={{-76,-26},{-56,-6}})));
      equation
        connect(incAng1.incAng, selfShadowingLeft.incAng[1]) annotation (Line(
              points={{11,-2},{36,-2},{36,0.2},{60.4,0.2}}, color={0,0,127}));
        connect(alt.y, selfShadowingLeft.alt) annotation (Line(points={{-55,18},{-18,
                18},{-18,16},{16,16},{30,16},{36,16},{36,10},{60.4,10}}, color={0,0,
                127}));
        connect(alt.y, incAng1.alt) annotation (Line(points={{-55,18},{-34,18},{-34,
                3.4},{-12.2,3.4}},     color={0,0,127}));
        connect(solAziSine.y, incAng1.solAzi) annotation (Line(
            points={{-55,-16},{-32,-16},{-32,-6.8},{-12,-6.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAziSine.y, selfShadowingLeft.solAzi) annotation (Line(
            points={{-55,-16},{44,-16},{44,19.8},{60.4,19.8}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}})),
          Documentation(info="<html>
<p>This model simulates a projection on the left-hand side of the window.</p>
</html>"));
      end SelfShadowingTestLeft;

      model SelfShadowingTestRight
        extends Modelica.Icons.Example;

        Windows.BaseClasses.SelfShadowing selfShadowingRight(
          n=1,
          final b={1},
          final h={1},
          final bLef={0},
          final dLef={0},
          final bBel={0},
          final dBel={0},
          final azi(displayUnit="deg") = {0},
          final til(displayUnit="deg") = {1.5707963267949},
          final bRig={1},
          final dRig={0.01},
          final bAbo={0},
          final dAbo={0})
          "Shadowing due to a projection on the right-hand side"
          annotation (Placement(transformation(extent={{60,-4},{92,24}})));
        BaseClasses.IncidenceAngleVDI6007 incAng1(azi=0, til=90)
          "Incidence Angle for the window"
          annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
        Modelica.Blocks.Sources.Constant alt(k=0.3490658504)
          "Constant altitude angle"
          annotation (Placement(transformation(extent={{-74,8},{-54,28}})));
        Modelica.Blocks.Sources.Sine solAziSine(amplitude=Modelica.Constants.pi,
            freqHz=1) "Solar azimuth angle generated as a sine"
          annotation (Placement(transformation(extent={{-76,-26},{-56,-6}})));
      equation
        connect(incAng1.incAng, selfShadowingRight.incAng[1]) annotation (Line(
              points={{11,-2},{36,-2},{36,0.2},{58.4,0.2}}, color={0,0,127}));
        connect(alt.y, selfShadowingRight.alt) annotation (Line(points={{-53,18},{-18,
                18},{-18,16},{16,16},{30,16},{36,16},{36,10},{58.4,10}}, color={0,0,
                127}));
        connect(alt.y, incAng1.alt) annotation (Line(points={{-53,18},{-34,18},{
                -34,3.4},{-12.2,3.4}}, color={0,0,127}));
        connect(solAziSine.y, incAng1.solAzi) annotation (Line(
            points={{-55,-16},{-32,-16},{-32,-6.8},{-12,-6.8}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(solAziSine.y, selfShadowingRight.solAzi) annotation (Line(
            points={{-55,-16},{44,-16},{44,19.8},{58.4,19.8}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics),
          Documentation(info="<html>
<p>This model simulates a projection on the right-hand side of the window.</p>
</html>"));
      end SelfShadowingTestRight;

      annotation (Documentation(info="<html>
This package contains four tests for the <a href=\"Windows.BaseClasses.SelfShadowing\">SelfShadowing</a> model. 
</html>"));
    end SelfShadowing;
    annotation (Documentation(info="<html>
<p>This package contains models for validation of Windows models. 
<\\p>

</html>", revisions="<html>
<ul>
<li>July 13, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
<ul>
</html>"));
  end Validation;

  model Window "Calculation of solar energy transmitted through windows"
    parameter Modelica.SIunits.Angle lat "Latitude";
    parameter Integer n(min = 1) "number of windows"
      annotation(dialog(group="window"));
    parameter Modelica.SIunits.CoefficientOfHeatTransfer UWin
      "Thermal transmission coefficient of whole window"
      annotation(dialog(group="window"));
     parameter Modelica.SIunits.TransmissionCoefficient g[n]
      "Total energy transmittance of windows"
      annotation(Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient g_TotDir[n]
      "Total energy transmittance of windows with closed sunscreen for direct radiation"
      annotation(Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient g_TotDif[n]
      "Total energy transmittance of windows with closed sunscreen for diffuse radiation"
      annotation(Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
      "degree of light transmission for direct irradiation"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_LTotDir[n]
      "degree of light transmission for direct irradiation, with sunscreen"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_LTotDif[n]
      "degree of light transmission for diffuse irradiation, with sunscreen"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.RadiantEnergyFluenceRate lim
      "Limit for the sunscreen to become active"
      annotation(dialog(group="sunscreen"));
    parameter Modelica.SIunits.Angle xi(displayUnit="degree")= 0
      "elevation angle";
    parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
      "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Angle azi[n] "Surface azimuth"
      annotation (Dialog(group="window"));
    extends Modelica.Blocks.Icons.Block;
    Modelica.Blocks.Interfaces.RealOutput HVis[n](
      final quantity="RadiantEnergyFluenceRate",
      final unit="W/m2") "solar energy entering the room in the visible area"
      annotation (Placement(transformation(extent={{100,30},{120,50}}),
          iconTransformation(extent={{100,30},{120,50}})));

    Modelica.Blocks.Interfaces.RealOutput HWin[n](
      final quantity="RadiantEnergyFluenceRate",
      final unit="W/m2")
      "Solar radiation transmitted through aggregated window"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}}),
          iconTransformation(extent={{100,-50},{120,-30}})));

    BaseClasses.Window window(
      final n=n,
      final UWin=UWin,
      final g=g,
      final g_TotDir=g_TotDir,
      final g_TotDif=g_TotDif,
      final T_L=T_L,
      final T_LTotDir=T_LTotDir,
      final T_LTotDif=T_LTotDif,
      final lim=lim,
      final til=til) annotation (Placement(transformation(extent={{20,-32},{86,30}})));
    Annex60.BoundaryConditions.SolarGeometry.IncidenceAngle incAng[n](
      each lat=lat,
      final azi=azi,
      final til=til) annotation (Placement(transformation(extent={{-40,74},{-30,84}})));
    Annex60.BoundaryConditions.WeatherData.Bus weaBus "Weather data"
      annotation (Placement(transformation(extent={{-108,-10},{-88,10}}),
          iconTransformation(extent={{-108,-10},{-88,10}})));
    Annex60.BoundaryConditions.SolarGeometry.BaseClasses.AltitudeAngle altAng
      annotation (Placement(transformation(extent={{-18,56},{-10,64}})));
    Annex60.BoundaryConditions.SolarGeometry.ZenithAngle zen(each lat=lat)
      annotation (Placement(transformation(extent={{-28,56},{-20,64}})));
    Annex60.BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil[n](
      final til=til,
      each lat=lat,
      final azi=azi) annotation (Placement(transformation(extent={{-32,-6},{-20,6}})));
    Annex60.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil[n](
      final til=til,
      each lat=lat,
      final azi=azi)
      annotation (Placement(transformation(extent={{-60,-54},{-40,-34}})));
  equation
    connect(window.HVis, HVis) annotation (Line(points={{89.3,11.4},{94,11.4},{94,
            40},{110,40}}, color={0,0,127}));
    connect(window.HWin, HWin) annotation (Line(points={{89.3,-13.4},{96,-13.4},{96,
            -40},{110,-40}}, color={0,0,127}));
    connect(incAng.y, window.incAng) annotation (Line(points={{-29.5,79},{-10,79},
            {-10,80},{8,80},{8,24.42},{16.7,24.42}}, color={0,0,127}));
    connect(weaBus.nTot, window.nTot) annotation (Line(
        points={{-98,0},{-90,0},{-90,72},{2,72},{2,18},{10,18},{10,16.36},{16.7,
            16.36}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(altAng.alt, window.alt) annotation (Line(points={{-9.6,60},{-8,60},{-8,
            7.06},{16.7,7.06}}, color={0,0,127}));
    connect(altAng.zen, zen.y)   annotation (Line(points={{-18.8,60},{-19.6,60}}, color={0,0,127}));
    connect(weaBus, zen.weaBus) annotation (Line(
        points={{-98,0},{-64,0},{-64,60},{-28,60}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(HDifTil.H, window.HDifTil) annotation (Line(points={{-19.4,0},{-16,0},
            {-16,-0.38},{16.7,-0.38}}, color={0,0,127}));
    connect(weaBus.HDifHor, window.HDifHor) annotation (Line(
        points={{-98,0},{-42,0},{-42,-7.82},{16.7,-7.82}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(weaBus.HDirNor, window.HDirNor) annotation (Line(
        points={{-98,0},{-42,0},{-42,-16.5},{16.7,-16.5}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(HDirTil.H, window.HDirTil) annotation (Line(points={{-39,-44},{-12,-44},
            {-12,-25.18},{16.7,-25.18}}, color={0,0,127}));
        for i in 1:n loop
            connect(weaBus, HDifTil[i].weaBus) annotation (Line(
        points={{-98,0},{-66,0},{-66,0},{-32,0}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
            connect(incAng[i].weaBus, weaBus) annotation (Line(
        points={{-40,79},{-90,79},{-90,0},{-98,0}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}}));
            connect(weaBus, HDirTil[i].weaBus) annotation (Line(
        points={{-98,0},{-80,0},{-80,-44},{-60,-44}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
        end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-96,96},{96,-96}},
            pattern=LinePattern.None,
            lineColor={0,0,0},
            fillColor={156,232,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-2,96},{2,-96}},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            lineColor={0,0,0}),
          Rectangle(
            extent={{-96,2},{96,-2}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),                      Diagram(
          coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
      Documentation(info="<html>
<p>This model calculates the input of heat and visible light into the room due to solar irradiation. Therefore it uses the calculations of VDI 6007 part 3.  It considers  the correction values for non-vertical and non-parallel radiation incidence.</p>
<p>To calculate the solar irradiation and the solar geometry it uses the models of the  <a href=\"Annex60.BoundaryConditions\">BoundaryConditions</a> package.</p>
<p>An example on how this model should be used is  <a href=\"vdi6007.Examples.Window\">Window</a>. To consider the additional heat input in case of ventilation with the solar protection the  <a href=\"vdi6007.BaseClasses.VentilationHeat\">VentilationHeat</a> model can be used.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
  end Window;

  model ShadedWindow
    "Calculation of solar energy transmitted through windows considering shadowing."
    parameter Modelica.SIunits.Angle lat "Latitude";
    parameter Integer n(min = 1) "number of windows"
      annotation(dialog(group="window"));
    parameter Modelica.SIunits.CoefficientOfHeatTransfer UWin
      "Thermal transmission coefficient of whole window"
      annotation(dialog(group="window"));
     parameter Modelica.SIunits.TransmissionCoefficient g[n]
      "Total energy transmittance of windows"
      annotation(Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_L[n]
      "degree of light transmission for direct irradiation"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_LTotDir[n]
      "degree of light transmission for direct irradiation, with sunscreen"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient T_LTotDif[n]
      "degree of light transmission for diffuse irradiation, with sunscreen"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.RadiantEnergyFluenceRate lim
      "Limit for the sunscreen to become active"
      annotation(dialog(group="sunscreen"));
    parameter Modelica.SIunits.Angle xi(  displayUnit="degree")= 0
      "elevation angle";
    parameter Modelica.SIunits.Angle til[n](displayUnit="deg")
      "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length b[n] "width of window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Height h[n] "height of window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length bLef[n] "window projection left"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length bRig[n] "window projection right"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length dLef[n]
      "distance between projection (left) and window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length dRig[n]
      "distance between projection (right) and window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length bAbo[n] "window projection above"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length bBel[n] "window projection below"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length dAbo[n]
      "distance between projection (above) and window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Length dBel[n]
      "distance between projection (below) and window"
      annotation (Dialog(group="window"));
    parameter Modelica.SIunits.Angle azi[n](displayUnit="degree")
      "Surface azimuth. azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
      annotation (Dialog(group="window"));
    parameter Integer nCorPoi(min = 1) "Number of corner points"
        annotation(dialog(group="skyline"));
    parameter Modelica.SIunits.Angle[nCorPoi] alpha(displayUnit="deg") "Azimuth of corner points, sorted from north to east to south to west,
     azi=-90 degree if surface outward unit normal points toward east; azi=0 if it points toward south"
        annotation(dialog(group="skyline"));
    parameter Modelica.SIunits.Height[nCorPoi] deltaH
      "Difference between height of corner point and the window centre"
      annotation(dialog(group="skyline"));
    parameter Modelica.SIunits.Distance[nCorPoi] s
      "horizontal distance between corner point and window centre"
      annotation(dialog(group="skyline"));
    parameter Boolean[nCorPoi-1] gap
      "corner points i and i+1 are gap between buildings: true, else: false"
      annotation(dialog(group="skyline"));
    parameter Modelica.SIunits.TransmissionCoefficient g_TotDir[n]
      "Total energy transmittance of windows with closed sunscreen for direct radiation"
      annotation(Dialog(group="window"));
    parameter Modelica.SIunits.TransmissionCoefficient g_TotDif[n]
      "Total energy transmittance of windows with closed sunscreen for diffuse radiation"
      annotation(Dialog(group="window"));

      Modelica.Blocks.Interfaces.RealOutput HVis[n](
      final quantity="RadiantEnergyFluenceRate",
      final unit="W/m2") "solar energy entering the room in the visible area"
      annotation (Placement(transformation(extent={{100,30},{120,50}}),
          iconTransformation(extent={{100,30},{120,50}})));

      Modelica.Blocks.Interfaces.RealOutput HWin[n](
      final quantity="RadiantEnergyFluenceRate",
      final unit="W/m2")
      "Solar radiation transmitted through aggregated window"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}}),
          iconTransformation(extent={{100,-50},{120,-30}})));
    BaseClasses.ShadedWindow shadedWindow(
      final n=n,
      final UWin=UWin,
      final g=g,
      final T_L=T_L,
      final T_LTotDir=T_LTotDir,
      final T_LTotDif=T_LTotDif,
      final lim=lim,
      final til=til,
      final b=b,
      final h=h,
      final bLef=bLef,
      final bRig=bRig,
      final dLef=dLef,
      final dRig=dRig,
      final bAbo=bAbo,
      final bBel=bBel,
      final dAbo=dAbo,
      final dBel=dBel,
      final azi=azi,
      final nCorPoi=nCorPoi,
      final alpha=alpha,
      final deltaH=deltaH,
      final s=s,
      final gap=gap,
      final g_TotDir=g_TotDir,
      final g_TotDif=g_TotDif)
      annotation (Placement(transformation(extent={{12,-34},{86,36}})));
    Annex60.BoundaryConditions.SolarGeometry.IncidenceAngle incAng[n](
      each lat=lat,
      final azi=azi,
      final til=til) annotation (Placement(transformation(extent={{-40,70},{-30,80}})));
    Annex60.BoundaryConditions.SolarGeometry.BaseClasses.AltitudeAngle altAng
      annotation (Placement(transformation(extent={{-30,-8},{-22,0}})));
    Annex60.BoundaryConditions.SolarGeometry.ZenithAngle zen(each lat=lat)
      annotation (Placement(transformation(extent={{-62,-12},{-54,-4}})));
    Annex60.BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil[n](
      final til=til,
      each lat=lat,
      final azi=azi) annotation (Placement(transformation(extent={{-38,-78},{-26,-66}})));
    Annex60.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil[n](
      final til=til,
      each lat=lat,
      final azi=azi)
      annotation (Placement(transformation(extent={{-66,40},{-52,54}})));
    Annex60.BoundaryConditions.WeatherData.Bus weaBus "Weather data"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    Annex60.BoundaryConditions.SolarGeometry.BaseClasses.SolarAzimuth solAzi(lat=lat)
      annotation (Placement(transformation(extent={{-40,-34},{-32,-26}})));
    Annex60.BoundaryConditions.SolarGeometry.BaseClasses.Declination decAng
      annotation (Placement(transformation(extent={{-66,-36},{-60,-30}})));
  equation
    connect(shadedWindow.HVis, HVis) annotation (Line(points={{89.7,15},{95.4,15},
            {95.4,40},{110,40}}, color={0,0,127}));
    connect(shadedWindow.HWin, HWin) annotation (Line(points={{89.7,-13},{95.4,-13},
            {95.4,-40},{110,-40}}, color={0,0,127}));
    connect(altAng.zen,zen. y)   annotation (Line(points={{-30.8,-4},{-30,-4},{-30,
            -4},{-36,-4},{-44,-4},{-44,-8},{-53.6,-8}},                           color={0,0,127}));
    connect(weaBus,zen. weaBus) annotation (Line(
        points={{-100,0},{-84,0},{-84,-6},{-68,-6},{-68,-8},{-62,-8},{-62,-8}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));

    connect(weaBus.HDifHor, shadedWindow.HDifHor) annotation (Line(
        points={{-100,0},{-84,0},{-84,-86},{0,-86},{0,-32.6},{8.3,-32.6}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(HDifTil.H, shadedWindow.HDifTil) annotation (Line(points={{-25.4,-72},
            {-20,-72},{-20,-24.9},{8.3,-24.9}}, color={0,0,127}));
    connect(altAng.alt, shadedWindow.alt) annotation (Line(points={{-21.6,-4},{-18,
            -4},{-18,-4.6},{8.3,-4.6}}, color={0,0,127}));
    connect(zen.y, solAzi.zen) annotation (Line(points={{-53.6,-8},{-50,-8},{-50,-27.6},
            {-40.8,-27.6}}, color={0,0,127}));
    connect(solAzi.solAzi, shadedWindow.solAzi) annotation (Line(points={{-31.6,-30},
            {-26,-30},{-26,-14.4},{8.3,-14.4}}, color={0,0,127}));
    connect(decAng.decAng, solAzi.decAng) annotation (Line(points={{-59.7,-33},{-45.85,
            -33},{-45.85,-30},{-40.8,-30}}, color={0,0,127}));
    connect(weaBus.cloTim, decAng.nDay) annotation (Line(
        points={{-100,0},{-84,0},{-84,-33},{-66.6,-33}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(weaBus.solTim, solAzi.solTim) annotation (Line(
        points={{-100,0},{-84,0},{-84,-44},{-84,-72},{-56,-72},{-56,-44},{-48,-44},
            {-42,-44},{-42,-32.4},{-40.8,-32.4}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(incAng.y, shadedWindow.incAng) annotation (Line(points={{-29.5,75},
            {-11.75,75},{-11.75,35.3},{8.3,35.3}},
                                       color={0,0,127}));
    connect(weaBus.nTot, shadedWindow.nTot) annotation (Line(
        points={{-100,0},{-84,0},{-84,64},{-18,64},{-18,26.9},{8.3,26.9}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    connect(HDirTil.H, shadedWindow.HDirTil) annotation (Line(points={{-51.3,47},{
            -20.65,47},{-20.65,16.4},{8.3,16.4}}, color={0,0,127}));
    for i in 1:n loop
                      connect(weaBus,HDifTil [i].weaBus) annotation (Line(
        points={{-100,0},{-84,0},{-84,-72},{-38,-72}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
            connect(incAng[i].weaBus,weaBus)  annotation (Line(
        points={{-40,75},{-84,75},{-84,0},{-100,0}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}}));
            connect(weaBus,HDirTil [i].weaBus) annotation (Line(
        points={{-100,0},{-84,0},{-84,47},{-66,47}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    end for;
    connect(weaBus.HDirNor, shadedWindow.HDirNor) annotation (Line(
        points={{-100,0},{-48,0},{-48,5.9},{8.3,5.9}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-6,3},{-6,3}}));
    annotation (defaultComponentName="shadedWindow",Icon(coordinateSystem(preserveAspectRatio=false),
          graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-96,96},{96,-96}},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Rectangle(
            extent={{-2,96},{2,-96}},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            lineColor={0,0,0}),
          Rectangle(
            extent={{-96,2},{96,-2}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-34,122},{40,88}},
            lineColor={28,108,200},
            fillColor={170,255,255},
            fillPattern=FillPattern.Solid,
            textString="%name
")}),                                                                                                  Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
<p>This model calculates the input of heat and visible light into the room due to solar irradiation. This model calculates the input of heat and visible light into the room due to solar irradiation. Therefore it uses the calculations of VDI 6007 part 3.  It considers  the correction values for non-vertical and non-parallel radiation incidence.</p>
<p> Additionaly to the <a href=\"Windows.Window\">Window</a> model it includes the formation of shades because of the window itself and because of the surrounding skyline.  </p>
<p>An example on how this model should be used is  <a href=\"Windows.Examples.ShadedWindow\">ShadedWindow</a>. To consider the additional heat input in case of ventilation with the solar protection the  <a href=\"Windows.BaseClasses.VentilationHeat\">VentilationHeat</a> model can be used.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
</html>", revisions="<html>
<ul>
<li>June 30, 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>"));
  end ShadedWindow;

  package SolarGain
    "Package with models for solar gain corrections according to VDI 6007 Part 3"
    extends Modelica.Icons.VariantsPackage;

    package BaseClasses "Package with base classes for SolarGain"
      extends Modelica.Icons.BasesPackage;

      partial model PartialCorrectionGTaue
        "Partial model for correction of the solar gain factor and for the transluence"
        import Modelica.Constants.pi;
        import Modelica.SIunits.Conversions.to_deg;
        import
          Annex60.ThermalZones.ReducedOrder.Windows.BaseClasses.Conversions.to_surfaceTiltVDI;
        parameter Integer n(min = 1) "number of windows"
          annotation(dialog(group="window"));
        parameter Modelica.SIunits.CoefficientOfHeatTransfer UWin
          "Thermal transmission coefficient of whole window"
          annotation(dialog(group="window"));
        parameter Modelica.SIunits.Angle xi( displayUnit="degree")=0
          "elevation angle";
        parameter Modelica.SIunits.Angle[n] til(displayUnit="degree")
          "Surface tilt. til=90 degree for walls; til=0 for ceilings; til=180 for roof"
          annotation(dialog(group="window"));

        Modelica.Blocks.Interfaces.RealOutput[n] CorG_Dir(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Transmission coefficient correction factor for direct radiation"
          annotation (Placement(transformation(extent={{80,-30},{100,-10}}),
          iconTransformation(extent={{80,-30},{100,-10}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorG_DifCle(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Transmission coefficient correction factor for diffuse radiation while clear sky"
          annotation (Placement(transformation(extent={{80,-50},{100,-30}}),
          iconTransformation(extent={{80,-50},{100,-30}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorG_DifCov(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Transmission coefficient correction factor for diffuse radiation while covered sky"
          annotation (Placement(transformation(extent={{80,-70},{100,-50}}),
          iconTransformation(extent={{80,-70},{100,-50}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorG_Gro(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Transmission coefficient correction factor for ground reflection radiation"
          annotation (Placement(transformation(extent={{80,-90},{100,-70}}),
          iconTransformation(extent={{80,-90},{100,-70}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorTaue_Dir(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Correction value for transluance for direct irradiation"
          annotation (Placement(transformation(extent={{80,10},{100,30}}),
          iconTransformation(extent={{80,10},{100,30}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorTaue_DifCle(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Correction value for transluance for diffuse irradiation during clear sky"
          annotation (Placement(transformation(extent={{80,30},{100,50}}),
          iconTransformation(extent={{80,30},{100,50}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorTaue_DifCov(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Correction value for transluance for diffuse irradiation during covered sky"
          annotation (Placement(transformation(extent={{80,50},{100,70}}),
          iconTransformation(extent={{80,50},{100,70}})));
        Modelica.Blocks.Interfaces.RealOutput[n] CorTaue_Gro(
          final quantity="TransmissionCoefficient",
          final unit="1")
          "Correction value for transluance for ground reflexion radiation"
          annotation (Placement(transformation(extent={{80,70},{100,90}}),
          iconTransformation(extent={{80,70},{100,90}})));

        Modelica.Blocks.Interfaces.RealInput incAng[n](
          final quantity="Angle",
          final unit="rad",
          displayUnit="degree")
          "Incidence angles of the sun beam on a tilted surface"
          annotation (Placement(transformation(extent={{-120,-10},{-80,30}}),
          iconTransformation(extent={{-100,10},{-80,30}})));
        Modelica.Blocks.Interfaces.BooleanInput sunscreen[n]
          "true: sunscreen closed, false: sunscreen open"
          annotation (Placement(transformation(extent={{-120,-40},{-80,0}}),
              iconTransformation(extent={{-100,-30},{-80,-10}})));

      protected
        parameter Real A0=0.918
          "Constant 0 to calculate reference transmission";
        parameter Real A1=2.21*10^(-4)
          "Constant 1 to calculate reference transmission";
        parameter Real A2=-2.75*10^(-5)
          "Constant 2 to calculate reference transmission";
        parameter Real A3=-3.82*10^(-7)
          "Constant 3 to calculate reference transmission";
        parameter Real A4=5.83*10^(-8)
          "Constant 4 to calculate reference transmission";
        parameter Real A5=-1.15*10^(-9)
          "Constant 5 to calculate reference transmission";
        parameter Real A6=4.74*10^(-12)
          "Constant 6 to calculate reference transmission";

        parameter Modelica.SIunits.TransmissionCoefficient tau_1DifCov= tau_DifCov*tau_iDif
          "Degreee of transmission for single pane window";
        parameter Modelica.SIunits.ReflectionCoefficient rho_T1DifCov=1-(tau_DifCov)
          "Part of degree of transmission for single pane window related to tau_1DifCov";
        parameter Modelica.SIunits.ReflectionCoefficient rho_11DifCov=rho_T1DifCov/
          (2-(rho_T1DifCov)) "Part of degree of transmission for single pane window related to
    rho_T1_diff";
        parameter Modelica.SIunits.ReflectionCoefficient rho_1DifCov= rho_11DifCov+
          (((1-rho_11DifCov)*tau_iDif)^2*rho_11DifCov)/(1-(rho_11DifCov*tau_iDif)^2)
          "Degree of reflection for single pane window";
        parameter Modelica.SIunits.TransmissionCoefficient tau_DifCov=0.84
          "Energetic degree of transmission for diffuse radiation for uniformly overcast sky";

        parameter Modelica.SIunits.TransmissionCoefficient tau_iDif=0.903
          "Pure degree of transmission for diffuse radiation";
        Modelica.SIunits.Angle[n] gamma_x
          "calculation factor for ground reflexion radiation";
        Modelica.SIunits.TransmissionCoefficient[n] tau_Dir
          "Energetic degree of transmission for direct radiation";
        Modelica.SIunits.TransmissionCoefficient[n] taui_Dir
          "Pure degree of transmission for direct radiation";
        Modelica.SIunits.TransmissionCoefficient[n] tau_1Dir
          "Pure degree of transmission for single pane window";
        Modelica.SIunits.ReflectionCoefficient[n] rho_T1Dir
          "Part of degree of transmission for single pane window related to tau_1Dir";
        Modelica.SIunits.ReflectionCoefficient[n] rho_11Dir
          "Part of degree of transmission for single pane window related to T1_Dir";
        Modelica.SIunits.ReflectionCoefficient[n] rho_1Dir
          "Degree of reflection for single pane window";
        Modelica.SIunits.Emissivity[n] a_1Dir
          "Degree of absorption for single pane window";

        Modelica.SIunits.TransmissionCoefficient[n] tau_DifCle
          "Energetic degree of transmission for diffuse radiation for clear sky";
        Modelica.SIunits.TransmissionCoefficient[n] tau_1DifCle
          "Degreee of transmission for single pane window";
        Modelica.SIunits.ReflectionCoefficient[n] rho_T1DifCle
          "Part of degree of transmission for single pane window related to tau_1DifCle";
        Modelica.SIunits.ReflectionCoefficient[n] rho_11DifCle
          "Part of degree of transmission for single pane window related to T1_DifCle";
        Modelica.SIunits.ReflectionCoefficient[n] rho_1DifCle
          "Degree of reflection for single pane window";
        Modelica.SIunits.Emissivity[n] a_1DifCle
          "Degree of absorption for single pane window";
        Modelica.SIunits.TransmissionCoefficient[n] tau_Gro
          "Energetic degree of transmission for ground reflexion radiation";
        Modelica.SIunits.TransmissionCoefficient[n] tau_1Gro
          "Degreee of transmission for single pane window";
        Modelica.SIunits.ReflectionCoefficient[n] rho_T1Gro
          "Part of degree of transmission for single pane window related to tau_1Gro";
        Modelica.SIunits.ReflectionCoefficient[n] rho_11Gro
          "Part of degree of transmission for single pane window related to T1_gr";
        Modelica.SIunits.ReflectionCoefficient[n] rho_1Gro
          "Degree of reflection for single pane window";
        Modelica.SIunits.Emissivity[n] a_1Gro
          "Degree of absorption for single pane window";
      equation

        for i in 1:n loop
        //Calculating variables for direct irradiation
        taui_Dir[i]= 0.907^(1/sqrt(1-(Modelica.Math.sin(incAng[i])/1.515)^2));
        if (((((A6*to_deg(incAng[i])+A5)*to_deg(incAng[i])+A4)*to_deg(incAng[i])+A3)*
        to_deg(incAng[i])+A2)*to_deg(incAng[i])+A1)*to_deg(incAng[i])+A0 <0 then
        tau_Dir[i]=0;
        else
        tau_Dir[i]= (((((A6*to_deg(incAng[i])+A5)*to_deg(incAng[i])+A4)*to_deg(incAng[i])+A3)*
        to_deg(incAng[i])+A2)*to_deg(incAng[i])+A1)*to_deg(incAng[i])+A0;
        end if;
        tau_1Dir[i]= tau_Dir[i]*taui_Dir[i];
        rho_T1Dir[i]= 1-tau_Dir[i];
        rho_11Dir[i]= rho_T1Dir[i]/(2-rho_T1Dir[i]);
        rho_1Dir[i]=rho_11Dir[i]+(((1-rho_11Dir[i])*taui_Dir[i])^2*rho_11Dir[i])/
        (1-(rho_11Dir[i]*taui_Dir[i])^2);
        a_1Dir[i]= 1-tau_1Dir[i]-rho_1Dir[i];
        //Calculating variables for diffuse, clear irradiation
        if 0.83-0.075*(to_deg(to_surfaceTiltVDI(til[i]))/70-1)^2+(0.052+0.033*(to_deg(to_surfaceTiltVDI(til[i]))/90-1)^2)
        *(Modelica.Math.cos(incAng[i])+0.15)^2 < 0 then
        tau_DifCle[i] = 0;
        else
        tau_DifCle[i]=0.83-0.075*(to_deg(to_surfaceTiltVDI(til[i]))/70-1)^2+(0.052+0.033*(to_deg(to_surfaceTiltVDI(til[i]))/90-1)^2)
        *(Modelica.Math.cos(incAng[i])+0.15)^2;
        end if;
        tau_1DifCle[i]= tau_DifCle[i]*tau_iDif;
        rho_T1DifCle[i]= 1-tau_DifCle[i];
        rho_11DifCle[i]= rho_T1DifCle[i]/(2-rho_T1DifCle[i]);
        rho_1DifCle[i]=rho_11DifCle[i]+(((1-rho_11DifCle[i])*tau_iDif)^2
        *rho_11DifCle[i])/(1-(rho_11DifCle[i]*tau_iDif)^2);
        a_1DifCle[i]=1-tau_1DifCle[i]-rho_1DifCle[i];
        //Calculating variables for ground reflexion radiation
        if (xi+to_surfaceTiltVDI(til[i]))<0 then
          gamma_x[i]=0;
        elseif (xi+to_surfaceTiltVDI(til[i]))>pi/2 then
          gamma_x[i]=pi/2;
        else
          gamma_x[i]=xi+to_surfaceTiltVDI(til[i]);
        end if;
        tau_Gro[i] = 0.84*Modelica.Math.sin(gamma_x[i])^(0.88*(1-0.5*abs(Modelica.Math.sin(2*gamma_x[i]))));
        tau_1Gro[i]=tau_Gro[i]*tau_iDif;
        rho_T1Gro[i]= 1-tau_Gro[i];
        rho_11Gro[i]= rho_T1Gro[i]/(2-rho_T1Gro[i]);
        rho_1Gro[i]=rho_11Gro[i]+(((1-rho_11Gro[i])*tau_iDif)^2*
        rho_11Gro[i])/(1-(rho_11Gro[i]*tau_iDif)^2);
        a_1Gro[i]=1-tau_1Gro[i]-rho_1Gro[i];
        end for;
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
        -100,-100},{100,100}})),Icon(coordinateSystem(
        preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
        extent={{-80,80},{80,-80}},
        lineColor={0,0,0},
        fillColor={215,215,215},
        fillPattern=FillPattern.Solid), Text(
        extent={{-52,24},{62,-16}},
        lineColor={0,0,0},
        textString="%name")}),
        Documentation(info="<html>
  <p>Partial model for correction factors for transmitted solar radiation
  through a transparent element.</p> 
<p><a href=\"vdi6007.BaseClasses.CorrrectionGTaueDoublePane\">CorrectionGTaueDoublePane</a> uses this model to calculate the correction values for double pane windows. This model can be used as a partial model to calculate the correction values for single pane windows and triple pane windows according to the VDI Guideline <\\p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3
  June 2015. Calculation of transient thermal response of rooms
  and buildings - modelling of solar radiation.</p>
  </html>",       revisions="<html>
  <p><i>February 24, 2014</i> by Reza Tavakoli:</p>
  <p>Implemented. </p>
<p><i>May 25, 2016 </i>by Stanley Risch:</p>
<p>Added the correction of the transluence factor according to VDI6007 Part 3</p>
  </html>"));
      end PartialCorrectionGTaue;
    annotation (Documentation(info="<html>
<p>
This package contains base classes to calculate solar gain through windows.
</p>
</html>"));
    end BaseClasses;

    model CorrectionGTaueDoublePane "Correction of the solar gain factor and the transluance factor according to 
  VDI6007 Part 3"
      extends BaseClasses.PartialCorrectionGTaue;
      import Modelica.SIunits.Conversions.to_deg;
      // Parameters for the transmission correction factor based on VDI 6007 Part 3
      // A0 to A6 are experimental constants VDI 6007 Part 3 page 20

      //Calculating the correction factor for direct solar radiation
      Modelica.SIunits.ReflectionCoefficient[n] XN2_Dir
        "Calculation factor to simplify equations";
      Modelica.SIunits.TransmissionCoefficient[n] tau_2Dir
        "Energetic dregree of transmission for second pane";
      Real[n] Q21_Dir
        "Coefficient of heat transfer for exterior pane of double pane window";
      Real[n] Q22_Dir
        "Coefficient of heat transfer for interior pane of double pane window";
      Real[n] Qsek2_Dir
        "Overall coefficient of heat transfer for double pane window";

      //diffuse clear
      Modelica.SIunits.ReflectionCoefficient[n] XN2_DifCle
        "Calculation factor to simplify equations";
      Modelica.SIunits.TransmissionCoefficient[n] tau_2DifCle
        "Energetic dregree of transmission for second pane";
      Real[n] Q21_DifCle
        "Coefficient of heat transfer for exterior pane of double pane window";
      Real[n] Q22_DifCle
        "Coefficient of heat transfer for interior pane of double pane window";
      Real[n] Qsek2_DifCle
        "Overall coefficient of heat transfer for double pane window";

      //ground
      Modelica.SIunits.ReflectionCoefficient[n] XN2_Gro
        "Calculation factor to simplify equations";
      Modelica.SIunits.TransmissionCoefficient[n] tau_2Gro
        "Energetic dregree of transmission for second pane";
      Real[n] Q21_Gro
        "Coefficient of heat transfer for exterior pane of double pane window";
      Real[n] Q22_Gro
        "Coefficient of heat transfer for interior pane of double pane window";
      Real[n] Qsek2_Gro
        "Overall coefficient of heat transfer for double pane window";
    protected
      parameter Modelica.SIunits.TransmissionCoefficient g_Dir0=taue_Dir0+Q210+Q220 "Reference vertical parallel transmission coefficient for direct radiation
    for double pane window";
      parameter Modelica.SIunits.TransmissionCoefficient Q210=(1-rho_1Dir0-0.907*A0)*(1+(0.907*A0*rho_1Dir0/(1-rho_1Dir0^2)))*UWin/25;
      parameter Modelica.SIunits.TransmissionCoefficient Q220=(1-rho_1Dir0-0.907*A0)*(0.907*A0/(1-rho_1Dir0^2))*(1-UWin/7.7);
      parameter Modelica.SIunits.TransmissionCoefficient taue_Dir0=(A0*0.907)^2/(1-rho_1Dir0^2)
        "Reference vertical parallel transmission coefficient for direct radiation";
      parameter Modelica.SIunits.ReflectionCoefficient rho_1Dir0=rho_11Dir0+(((1-rho_11Dir0)*0.907)^2*rho_11Dir0)/
      (1-(rho_11Dir0*0.907)^2);
      parameter Modelica.SIunits.ReflectionCoefficient rho_11Dir0=(1-A0)/(2-(1-A0));
      parameter Modelica.SIunits.ReflectionCoefficient XN2_DifCov=1-rho_1DifCov^2
        "Calculation factor to simplify equations";
      parameter Modelica.SIunits.TransmissionCoefficient tau_2DifCov=(tau_1DifCov^2)/
        XN2_DifCov "Energetic dregree of transmission for second pane";
      parameter Modelica.SIunits.Emissivity a_1DifCov=1-tau_1DifCov-rho_1DifCov
        "Degree of absorption for single pane window";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer Q21_DifCov=
        a_1DifCov*(1+(tau_1DifCov*rho_1DifCov/XN2_DifCov))*UWin/25
        "Coefficient of heat transfer for exterior pane of double pane window";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer Q22_DifCov=
        a_1DifCov*(tau_1DifCov/XN2_DifCov)*(1-(UWin/7.7))
        "Coefficient of heat transfer for interior pane of double pane window";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer Qsek2_DifCov=
        Q21_DifCov+Q22_DifCov
        "Overall coefficient of heat transfer for double pane window";

    equation
      for i in 1:n loop
        //Calculating variables for the overall degree of energy passage for direct irradiation
        if (1-rho_1Dir[i]^2)==0 then
          XN2_Dir[i]=10^(-20);
        else
          XN2_Dir[i]= 1-rho_1Dir[i]^2;
        end if;
        Q21_Dir[i]=a_1Dir[i]*(1+(tau_1Dir[i]*rho_1Dir[i]/XN2_Dir[i]))*UWin/25;
        Q22_Dir[i]= a_1Dir[i]*(tau_1Dir[i]/XN2_Dir[i])*(1-(UWin/7.7));
        Qsek2_Dir[i]=Q21_Dir[i]+Q22_Dir[i];
        tau_2Dir[i]= tau_1Dir[i]^2/XN2_Dir[i];

        //Calculating variables for diffuse irradiation at clear sky
        if (1-rho_1DifCle[i]^2)==0 then
          XN2_DifCle[i]=10^(-20);
        else
          XN2_DifCle[i]= 1-rho_1DifCle[i]^2;
        end if;
        Q21_DifCle[i]=a_1DifCle[i]*(1+(tau_1DifCle[i]*rho_1DifCle[i]/XN2_DifCle[i]))*UWin/25;
        Q22_DifCle[i]= a_1DifCle[i]*(tau_1DifCle[i]/XN2_DifCle[i])*(1-(UWin/7.7));
        Qsek2_DifCle[i]=Q21_DifCle[i]+Q22_DifCle[i];
        tau_2DifCle[i]= tau_1DifCle[i]^2/XN2_DifCle[i];

        //Calculating variables for the overall degree of energy passage for ground reflexion radiation
        if (1-rho_1Gro[i]^2)==0 then
          XN2_Gro[i]=10^(-20);
        else
          XN2_Gro[i]= 1-rho_1Gro[i]^2;
        end if;
        Q21_Gro[i]=a_1Gro[i]*(1+(tau_1Gro[i]*rho_1Gro[i]/XN2_Gro[i]))*UWin/25;
        Q22_Gro[i]= a_1Gro[i]*(tau_1Gro[i]/XN2_Gro[i])*(1-(UWin/7.7));
        Qsek2_Gro[i]=Q21_Gro[i]+Q22_Gro[i];
        tau_2Gro[i]= tau_1Gro[i]^2/XN2_Gro[i];

        //Calculating correction values
        CorG_DifCov[i]=(tau_2DifCov+Qsek2_DifCov)/g_Dir0;
        CorTaue_DifCov[i]=tau_2DifCov/taue_Dir0;
         if sunscreen[i] then
          CorTaue_DifCle[i]=CorTaue_DifCov[i];
          CorTaue_Gro[i]=CorTaue_DifCov[i];
          CorTaue_Dir[i]=CorTaue_DifCov[i];
          CorG_Dir[i]=CorG_DifCov[i];
          CorG_DifCle[i]=CorG_DifCle[i];
          CorG_Gro[i]=CorG_DifCle[i];
         else
          CorTaue_DifCle[i]=tau_2DifCle[i]/taue_Dir0;
          CorTaue_Gro[i]=tau_2Gro[i]/taue_Dir0;
          CorTaue_Dir[i]=tau_2Dir[i]/taue_Dir0;
          CorG_Dir[i]= (tau_2Dir[i]+Qsek2_Dir[i])/g_Dir0;
          CorG_DifCle[i]= (tau_2DifCle[i]+Qsek2_DifCle[i])/g_Dir0;
          CorG_Gro[i]= (tau_2Gro[i]+Qsek2_Gro[i])/g_Dir0;
         end if;
      end for;

      annotation (defaultComponentName="CorGTaue",
      Diagram(coordinateSystem(
      preserveAspectRatio=false,
      extent={{-100,-100},{100,100}},
      grid={2,2})),
      Icon(coordinateSystem(
      preserveAspectRatio=false,
      extent={{-100,-100},{100,100}},
      grid={2,2})),
      Documentation(info="<html>
  <p><a href=\"vdi6007.BaseClasses.CorrrectionGTaueDoublePane\">CorrectionGTaueDoublePane</a> computes 
  transmission correction factors for the g-factor and the transluence. Transmission properties of transparent
  elements are in general dependent on the solar incidence angle. To take this
  dependency into account, correction factors can multiplied with the solar
  radiation. These factors should not be mistaken as calculation of solar
  radiation on tilted surfaces, calculation of g-value or consideration of
  sunblinds, it is an additional step. The implemented calculations are
  defined in the German Guideline VDI 6007 Part 3 (VDI, 2015). The given model
  is only valid for double pane windows. The guideline describes also
  calculations for single pane and triple pane windows.</p>
  <h4>References</h4>
  <p>VDI. German Association of Engineers Guideline VDI 6007-3 June 2015.
  Calculation of transient thermal response of rooms and buildings -
  modelling of solar radiation.</p>
  </html>",
      revisions="<html>
<p><i>February 24, 2014</i> by Reza Tavakoli: </p>
<p>Implemented. </p>
<p><i>September 12, 2015 </i>by Moritz Lauster: </p>
<p>Adapted to Annex 60 requirements. </p>
<p><i>May 25, 2016 </i>by Stanley Risch:</p>
<p>Added the correction of the transluence factor according to VDI6007 Part 3</p>
</html>"));
    end CorrectionGTaueDoublePane;
  annotation (Documentation(info="<html>
<p>
This package contains models to compute solar heat gains.
</p>
</html>"));
  end SolarGain;

  package UsersGuide "User's Guide"
      extends Modelica.Icons.Information;
    annotation (Documentation(info="<html>
<p>The Windows package contains the models to simulate transparent objects in building simulations. It is based on the modelling of solar radiation of VDI 6007. The upper models <a href=\"Windows.Window\">Window</a> and <a href=\"Windows.ShadedWindow\">ShadedWindow</a> calculate the entering solar energy into the room. They use the <a href=\"Windows.SolarGain.CorrectionGTaueDoublePane\">CorrectionGTaueDoublePane</a> to calculate correction values for non-vertical and non-parrallel irradiation. To respect the heat input  at closed sunscreen and open window the <a href=\"Windows.BaseClasses.VentilationHeat\">VentilationHeat</a><\\p> model can be used. The <a href=\"Windows.ShadedWindow\">ShadedWindow</a> additionally considers shading because of window projections through the <a href=\"Windows.BaseClasses.SelfShadowing\">SelfShadowing</a> model and shading because of surrounding buildings through the <a href=\"Windows.BaseClasses.SkylineShadowing\">SkylineShadowing</a>  model. The entering visible light is also calculated by the upper classes. It can be used to determine the switch moment of the lighting with the <a href=\"Windows.BaseClasses.Illumination\">Illumination</a> model. The information sections of the individual models give extra information on the calculations.<\\p> 
<p>The <a href=\"Windows.Examples\">Examples</a> package contains examples on how the models should be integrated.<\\p>
<p>The <a href=\"Windows.Validation\">Validation</a> is splitted in <a href=\"Windows.Validation.VDI2078\">VDI2078</a> and <a href=\"Windows.Validation.Shadowing\">Shadowing</a>. The <a href=\"Windows.BaseClasses.Illumination\">Illumination</a> and the <a href=\"Windows.BaseClasses.VentilationHeat\">VentilationHeat</a> model were tested with parts of test case 1 and parts of test case 3 of VDI2078.  <a href=\"Windows.SolarGain.CorrectionGTaueDoublePane\">CorrectionGTaueDoublePane</a> is tested within the test cases. The shadowing models are not included in the validation of VDI2078. Therfore the models were tested on plausibility with simple examples.<\\p>
</html>"));
  end UsersGuide;
  annotation (Documentation(revisions="<html>
<ul>
<li>July 17 2016,&nbsp; by Stanley Risch:<br>Implemented. </li>
</html>", info="<html>
<p>This package provides two models which calculate the solar heat transmitted through windows into the room. <a href=\"Windows.Window\">Window</a> considers correction values for non-parallel and non-vertical radiation. <a href=\"Windows.ShadedWindow\">ShadedWindow</a> additionally includes shadowing due to the window itself and surrounding buildings. The <a href=\"Windows.BaseClasses\">BaseClasses</a>-package contains an <a href=\"Windows.Illumination\">Illumination</a> model which calculates the activation and deactivation of the illumination <\\p>
</html>"));
end Windows;
