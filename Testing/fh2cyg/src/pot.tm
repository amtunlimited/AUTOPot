:Begin:
:Function: 		pot
:Pattern:			fh2cyginner[Coord_]
:Arguments:		{Coord}
:ArgumentTypes:	{RealList}
:ReturnType:		Manual
:End:
:Evaluate:		fh2cyg[coord_]:=Module[{answer=fh2cyginner[coord]},
					Return[{answer[[1]],
							Table[answer[[i]],{i,2,4}],
							Table[answer[[i]],{i,5,7}],
							Table[answer[[i]],{i,8,10}]
						}];
					]
