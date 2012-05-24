:Begin:
:Function: 		pot
:Pattern:			WinTestinner[Coord_]
:Arguments:		{Coord}
:ArgumentTypes:	{RealList}
:ReturnType:		Manual
:End:
:Evaluate:		WinTest[coord_]:=Module[{answer=WinTestinner[coord]},
					Return[{answer[[1]],
							Table[answer[[i]],{i,2,4}],
							Table[answer[[i]],{i,5,7}],
							Table[answer[[i]],{i,8,10}]
						}];
					]
