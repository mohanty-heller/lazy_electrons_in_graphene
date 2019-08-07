tstart = AbsoluteTime[];
nproc = 8;
LaunchKernels[nproc];

row = 1;

ClearSystemCache[];
meminit = MemoryInUse[];

tint = 0;
tfinal = 30000;
\[Tau] = 1;
trange = Table[t, {t, tint, tfinal, \[Tau]}];

kx = 0;
ky = 0;
kvec = {kx, ky};

n = 8(*number of atoms in row, even*);
m = n/2(*number of rows, even*);
a0 = 1;
a = .142/.0529 a0(*0.142 nm is equilibrium bond length*);
t0 = 2.7/13.605698065;(*Rydberg energy*)
c = 3.37;
\[HBar] = 1;
\[Omega]0 = 0.0146;
M = 10946.8;
kconst = 2.333419;
kB = 8.61733*^-5/13.60569;
Temp = 300;

lattice1 = {n a Sqrt[3]/2, 0};
lattice2 = {m a Sqrt[3]/2, m 3/2 a};


xmat = ConstantArray[0, {m, n}];
ymat = ConstantArray[0, {m, n}];
Do[x = (i - 1) Sqrt[3]/2 a;
 Do[xmat[[i]][[j]] = x; x = x + a Sqrt[3]/2;, {j, 1, n}], {i, 1, m}]
Do[y = (i - 1) 3/2 a;
  Do[If[Mod[j, 2] == 1, ymat[[i]][[j]] = y, 
     ymat[[i]][[j]] = y + a/2];, {j, 1, n}], {i, 1, m}];
flatxmat = Flatten[xmat];
flatymat = Flatten[ymat];

(*Manipulate[ListPlot[Join[Thread[{flatxmat,flatymat}],Thread[{\
flatxmat,flatymat}+lattice1]],PlotRange\[Rule]{{-1,20},{-1,20}},\
AspectRatio\[Rule]1,PlotMarkers\[Rule]Graphics[{Red,Disk[]},ImageSize\
\[Rule]10]],{t,0,10000\[Pi]}]*)

V = 0;
Do[Do[index = (i - 1) n + j;
   neighbors = {};
   If[i == 1 && Mod[j, 2] == 1, 
    Print[(m - 1) n + j + 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[(m - 1) n + j + 1]] + 
          m a Sqrt[3]/2)^2 + (flatymat[[index]] - 
          flatymat[[(m - 1) n + j + 1]] + 3/2 m a)^2]];
    V = V + (\[Sqrt](((var[[2 index - 1]] + flatxmat[[index]] + 
                m a Sqrt[3]/2 - flatxmat[[(m - 1) n + j + 1]] - 
                var[[2 ((m - 1) n + j + 1) - 1]]))^2 + ((var[[
                 2 index]] + flatymat[[index]] + 3/2 m a - 
                flatymat[[(m - 1) n + j + 1]] - 
                var[[2 ((m - 1) n + j + 1)]]))^2) - a)^2;, {}];
   If[i == m && Mod[j, 2] == 0, 
    Print[j - 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[j - 1]] - 
          m a Sqrt[3]/2)^2 + (flatymat[[index]] - flatymat[[j - 1]] - 
          3/2 m a)^2]];
    V = V + (\[Sqrt](((var[[2 index - 1]] + flatxmat[[index]] - 
                m a Sqrt[3]/2 - flatxmat[[j - 1]] - 
                var[[2 (j - 1) - 1]]))^2 + ((var[[2 index]] + 
                flatymat[[index]] - 3/2 m a - flatymat[[j - 1]] - 
                var[[2 (j - 1)]]))^2) - a)^2;, {}];
   If[j == 1, 
    Print[i n, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[i n]] + 
          n a Sqrt[3]/2)^2 + (flatymat[[index]] - flatymat[[i n]])^2]];
    V = V + (\[Sqrt](((var[[2 index - 1]] + flatxmat[[index]] - 
                flatxmat[[i n]] + n a Sqrt[3]/2 - 
                var[[2 i n - 1]]))^2 + ((var[[2 index]] + 
                flatymat[[index]] - flatymat[[i n]] - 
                var[[2 i n]]))^2) - a)^2;, 
    AppendTo[neighbors, (i - 1) n + (j - 1)]];
   If[j == n, 
    Print[(i - 1) n + 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[(i - 1) n + 1]] - 
          n a Sqrt[3]/2)^2 + (flatymat[[index]] - 
          flatymat[[(i - 1) n + 1]])^2]];
    V = V + (\[Sqrt](((var[[2 index - 1]] + flatxmat[[index]] - 
                flatxmat[[(i - 1) n + 1]] - n a Sqrt[3]/2 - 
                var[[2 ((i - 1) n + 1) - 1]]))^2 + ((var[[2 index]] + 
                flatymat[[index]] - flatymat[[(i - 1) n + 1]] - 
                var[[2 ( (i - 1) n + 1)]]))^2) - a)^2;, 
    AppendTo[neighbors, (i - 1) n + (j + 1)]];
   If[Mod[j, 2] == 0, 
    If[i == m, {}, AppendTo[neighbors, (i) n + j - 1]], 
    If[i == 1, {}, AppendTo[neighbors, (i - 2) n + j + 1]]];
   Print[neighbors];
   Do[V = 
      V + (\[Sqrt](((var[[2 index - 1]] + flatxmat[[index]] - 
                 flatxmat[[neighbors[[k]]]] - 
                 var[[2 neighbors[[k]] - 1]]))^2 + ((var[[2 index]] + 
                 flatymat[[index]] - flatymat[[neighbors[[k]]]] - 
                 var[[2 neighbors[[k]]]]))^2) - a)^2;, {k, 1, 
     Length[neighbors]}];, {j, 1, n}], {i, 1, m}];
V = V/4 kconst(*1/2 from 1/2kx^2,1/2 from the summation convention*);
U = ConstantArray[0, {2 m n, 2 m n}];
Do[Do[Clear[temp];
   temp = D[D[V, var[[i]]], var[[j]]];
   var = Flatten[Thread[{flatxmat, flatymat}]];
   U[[i]][[j]] = temp;
   Clear[var];, {j, 1, 2 m n}], {i, 1, 2 m n}];
P = Eigenvectors[U];
eigval = Eigenvalues[U] // Chop;
eigval[[2 m n]] = 0;
eigval[[2 m n - 1]] = 0;
eigval[[2 m n - 2]] = 0;
(*Do[If[Abs[eigval[[i]]]<10^-9,eigval[[i]]=0,{}],{i,1,Length[eigval]}]\
;*)
Udiag = DiagonalMatrix[eigval];
var0 = ConstantArray[0, {1, 2 m n}][[
   1]](*RandomReal[{0,1.5/10a},2m n]*);
var0 = Sum[
   If[eigval[[i]] == 0, ConstantArray[0, {2 m n}], 
    P[[i]] Sqrt[(2 kB*Temp)/eigval[[i]]]], {i, 1, 2 m n}];

timediag = {};
Do[If[eigval[[i]] == 0, AppendTo[timediag, 0], 
  AppendTo[timediag, Cos[Sqrt[eigval[[i]]/M] t]]], {i, 1, 2 m n}]
coords[t_] = Transpose[P].DiagonalMatrix[timediag].P.var0 // Chop;
flatx[t_] = Table[coords[t][[2 i - 1]] + flatxmat[[i]], {i, 1, m n}];
flaty[t_] = Table[coords[t][[2 i]] + flatymat[[i]], {i, 1, m n}];

dflatx[T_] = D[flatx[T], T];
dflaty[T_] = D[flaty[T], T];

Print["Lattice generated."];


(* Manipulate[ListPlot[Join[Thread[{flatx[t],flaty[t]}]],PlotRange\
\[Rule]{{-1,20},{-1,20}},AspectRatio\[Rule]1,PlotMarkers\[Rule]\
Graphics[{Red,Disk[]},ImageSize\[Rule]10]],{t,0,10000\[Pi]}] *)

ClearSystemCache[];
H0 = ConstantArray[0, {m n, 1}];
S0 = ConstantArray[0, {m n, 1}];
B0 = ConstantArray[0, {m n, 1}];
Clear[index];
Do[ParallelDo[index = (i - 1) n + j;
   neighbors = {};
   \[Rho] = 
    ConstantArray[0, {IntegerPart[(tfinal - tint)/\[Tau]] + 1}];
   If[i == 1 && Mod[j, 2] == 1, 
    Print[(m - 1) n + j + 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[(m - 1) n + j + 1]] + 
          m a Sqrt[3]/2)^2 + (flatymat[[index]] - 
          flatymat[[(m - 1) n + j + 1]] + 3/2 m a)^2]];
    Print[ToString[index], ", ", ToString[(m - 1) n + j + 1]];
    hopping = {};
    xdiff = 
     Table[flatx[t][[index]] - flatx[t][[(m - 1) n + j + 1]] + 
       m a Sqrt[3]/2, {t, tint, tfinal, \[Tau]}];
    ydiff = 
     Table[flaty[t][[index]] - flaty[t][[(m - 1) n + j + 1]] + 
       3/2 m a, {t, tint, tfinal, \[Tau]}];
    dxdtneighbor = 
     Table[dflatx[t][[(m - 1) n + j + 1]], {t, tint, 
       tfinal, \[Tau]}];
    dydtneighbor = 
     Table[dflaty[t][[(m - 1) n + j + 1]], {t, tint, 
       tfinal, \[Tau]}];
    R = Sqrt[((#[[1]])^2) + ((#[[2]])^2)] & /@ 
      Thread[{xdiff, ydiff}];
    Print["R", R];
    ClearSystemCache[];
    Print[MemoryInUse[] - meminit];
    hopping = (-t0 Exp[-c ( # /a - 1)]) Exp[I kvec.lattice2] & /@ R;
    Print["hopping,", hopping];
    Export[
     "~/outputs/gamma/H" <> ToString[index] <> "," <> 
      ToString[(m - 1) n + j + 1] <> ".mat", hopping];
    Print[" Time elapsed: ", AbsoluteTime[] - tstart];, {}];
   If[i == m && Mod[j, 2] == 0, 
    Print[j - 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[j - 1]] - 
          m a Sqrt[3]/2)^2 + (flatymat[[index]] - flatymat[[j - 1]] - 
          3/2 m a)^2]];
    hopping = {};
    xdiff = 
     Table[flatx[t][[index]] - flatx[t][[j - 1]] - m a Sqrt[3]/2, {t, 
       tint, tfinal, \[Tau]}];
    ydiff = 
     Table[flaty[t][[index]] - flaty[t][[j - 1]] - 3/2 m a, {t, tint, 
       tfinal, \[Tau]}];
    dxdtneighbor = Table[dflatx[t][[j]], {t, tint, tfinal, \[Tau]}];
    dydtneighbor = Table[dflaty[t][[j]], {t, tint, tfinal, \[Tau]}];
    R = Sqrt[((#[[1]])^2) + ((#[[2]])^2)] & /@ 
      Thread[{xdiff, ydiff}];
    ClearSystemCache[];
    Print[MemoryInUse[] - meminit];
    Print["R", R];
    hopping = Exp[-I kvec.lattice2] (-t0 Exp[-c ( # /a - 1)]) & /@ R;
    Print["hopping,", hopping];
    Export[
     "~/outputs/gamma/H" <> ToString[index] <> "," <> ToString[j - 1] <> 
      ".mat", hopping];
    Print[" Time elapsed: ", AbsoluteTime[] - tstart];, {}];
   If[j == 1, 
    Print[i n, "---", 
       Sqrt[(flatxmat[[index]] - flatxmat[[i n]] + 
            n a Sqrt[3]/2)^2 + (flatymat[[index]] - 
            flatymat[[i n]])^2]]
      hopping = {};
    xdiff = 
     Table[flatx[t][[index]] - flatx[t][[i n]] + n a Sqrt[3]/2, {t, 
       tint, tfinal, \[Tau]}];
    ydiff = 
     Table[flaty[t][[index]] - flaty[t][[i n]], {t, tint, 
       tfinal, \[Tau]}];
    dxdtneighbor = Table[dflatx[t][[j]], {t, tint, tfinal, \[Tau]}];
    dydtneighbor = Table[dflaty[t][[j]], {t, tint, tfinal, \[Tau]}];
    R = Sqrt[((#[[1]])^2) + ((#[[2]])^2)] & /@ 
      Thread[{xdiff, ydiff}];
    ClearSystemCache[];
    Print[MemoryInUse[] - meminit];
    Print["R", R];
    hopping = Exp[I kvec.lattice1] (-t0 Exp[-c ( # /a - 1)]) & /@ R;
    Print["hopping,", hopping];
    Export[
     "~/outputs/gamma/H" <> ToString[index] <> "," <> ToString[i n] <> 
      ".mat", hopping];
    Print[" Time elapsed: ", AbsoluteTime[] - tstart];, 
    AppendTo[neighbors, (i - 1) n + (j - 1)]];
   If[j == n, 
    Print[(i - 1) n + 1, "---", 
     Sqrt[(flatxmat[[index]] - flatxmat[[(i - 1) n + 1]] - 
          n a Sqrt[3]/2)^2 + (flatymat[[index]] - 
          flatymat[[(i - 1) n + 1]])^2]];
    hopping = {};
    xdiff = 
     Table[flatx[t][[index]] - flatx[t][[(i - 1) n + 1]] - 
       n a Sqrt[3]/2, {t, tint, tfinal, \[Tau]}];
    ydiff = 
     Table[flaty[t][[index]] - flaty[t][[(i - 1) n + 1]], {t, tint, 
       tfinal, \[Tau]}];
    dxdtneighbor = Table[dflatx[t][[j]], {t, tint, tfinal, \[Tau]}];
    dydtneighbor = Table[dflaty[t][[j]], {t, tint, tfinal, \[Tau]}];
    R = Sqrt[((#[[1]])^2) + ((#[[2]])^2)] & /@ 
      Thread[{xdiff, ydiff}];
    ClearSystemCache[];
    Print[MemoryInUse[] - meminit];
    Print["R", R];
    hopping = Exp[-I kvec.lattice1] (-t0 Exp[-c ( # /a - 1)]) & /@ R;
    Print["hopping,", hopping];
    Export[
     "~/outputs/gamma/H" <> ToString[index] <> "," <> 
      ToString[(i - 1) n + 1] <> ".mat", hopping];
    Print[" Time elapsed: ", AbsoluteTime[] - tstart];, 
    AppendTo[neighbors, (i - 1) n + (j + 1)]];
   If[Mod[j, 2] == 0, 
    If[i == m, {}, AppendTo[neighbors, (i) n + j - 1]], 
    If[i == 1, {}, AppendTo[neighbors, (i - 2) n + j + 1]]];
   Do[Print[ToString[index], ", ", ToString[neighbors[[k]]]];
    hopping = {};
    xdiff = 
     Table[flatx[t][[index]] - flatx[t][[neighbors[[k]]]], {t, tint, 
       tfinal, \[Tau]}];
    ydiff = 
     Table[flaty[t][[index]] - flaty[t][[neighbors[[k]]]], {t, tint, 
       tfinal, \[Tau]}];
    dxdtneighbor = 
     Table[dflatx[t][[neighbors[[k]]]], {t, tint, tfinal, \[Tau]}];
    dydtneighbor = 
     Table[dflaty[t][[neighbors[[k]]]], {t, tint, tfinal, \[Tau]}];
    R = Sqrt[((#[[1]])^2) + ((#[[2]])^2)] & /@ 
      Thread[{xdiff, ydiff}];
    ClearSystemCache[];
    Print[MemoryInUse[] - meminit];
    Print["R", R];
    hopping = (-t0 Exp[-c ( # /a - 1)]) & /@ R;
    Print["hopping,", hopping];
    Export[
     "~/outputs/gamma/H" <> ToString[index] <> "," <> 
      ToString[neighbors[[k]]] <> ".mat", hopping];
    Print[" Time elapsed: ", AbsoluteTime[] - tstart];, {k, 1, 
     Length[neighbors]}];
   , {j, 1, n}], {i, row, row}];
Print["Hamiltonian generated."]
Quit[]