#!/usr/local/bin/MathematicaScript -script


list1 = ToExpression[ $ScriptCommandLine[[2]] ]

plot = ListPlot[list1, Joined -> True, PlotMarkers ->{Automatic, Large}, PlotRange -> {{15, 25}, {15, 25}}, PlotStyle -> Thick] 

Export[ $ScriptCommandLine[[3]] , plot, ImageSize->{1200}] 
