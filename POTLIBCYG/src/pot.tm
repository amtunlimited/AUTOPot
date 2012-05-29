:Begin:
:Function: 		pot
:Pattern:			FUNCinner[Coord_]
:Arguments:		{Coord}
:ArgumentTypes:	{RealList}
:ReturnType:		Manual
:End:
:Evaluate:		FUNC[coord_]:=Module[{answer=FUNCinner[coord]},
					Return[{answer[[1]],
							Table[answer[[i]],{i,2,4}],
							Table[answer[[i]],{i,5,7}],
							Table[answer[[i]],{i,8,10}]
						}];
					]
