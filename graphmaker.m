#!/usr/local/bin/MathematicaScript -script


list1 = ToExpression[ $ScriptCommandLine[[2]] ]

plot = ListPlot[list1, Joined -> True, PlotMarkers ->{Automatic, Large}, PlotRange -> {{15, 25}, {15, 25}}, PlotStyle -> Thick] 
plot2 = ListPlot3D[list10, PlotMarkers -> {Automatic, Large}, PlotRange -> {{17, 24.01}, {16.95, 24}}, PlotStyle -> Thick, Ticks -> None, GridLines -> Automatic, GridLinesStyle -> Directive[Dashed, Gray], Axes -> False]
Export[ $ScriptCommandLine[[3]] , plot, ImageSize->{1200}] 
