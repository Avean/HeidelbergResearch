using ..Dictionaries

F = DictSinCosFun["Neumann"];
E = DictSinCosEig["Neumann"];

using ..LaplaceDiscretisation

SC.N.Trunc
SC.N.Full


SC.P.Trunc
SC.P.Full

SC.N.Trunc * SC.N.Full'