:Begin:
:Function: 		pot
:Pattern:			fhcygtestinner[Coord_]
:Arguments:		{Coord}
:ArgumentTypes:	{RealList}
:ReturnType:		Manual
:End:
:Evaluate:		fhcygtest[coord_]:=Module[{answer=fhcygtestinner[coord]},
					Return[{answer[[1]],
							Table[answer[[i]],{i,2,4}],
							Table[answer[[i]],{i,5,7}],
							Table[answer[[i]],{i,8,10}]
						}];
					]
