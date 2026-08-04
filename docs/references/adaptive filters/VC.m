(*
	Implementation of Schlicht's (1989)
	time Varying Coefficients (VC) estimation algorithm.

	Reference: Schlicht's paper can be downloaded from
		http://www.semverteilung.vwl.uni-muenchen.de

	Implemented by Johannes Ludsteck,
		Institute of Employment Research (IAB) Nuremberg, Germany

	The Package comes without any warranty...
*)

BeginPackage["VC`VC`",
  {"LinearAlgebra`MatrixManipulation`",
   "Statistics`NormalDistribution`",
   "Utilities`FilterOptions`"}];

Unprotect[{VCEstimate, CIPlot, Panel, CILevel, Columns, RegressorNames, VarianceRatios}];

(* Usage Information -------------------------------------------- *)

VCEstimate::usage = "VCEstimate[x,y] computes Schlicht's (1989) time varying
coefficients estimate of the linear model \n y[t] = x[t] b[t] + u[t]. 
with b[t] = b[t-1] + v[t]. y is the dependent 
variable, x the a matrix of regressors (each regressor occupies one column in x). 
VCEstimate returns the List {coeff, sdb, sdu, sdi} where coeff is the matrix of
coefficient time paths, sdb the matrix of the respective coefficient 
standard deviations, sdu the estimated standard deviation of u, 
and sdi the vector of the respective random walk residual standard deviations. \n\n
The estimation method is determined by the option EstimationMethod. 
If option Panel-> atimevar is given, data are interpreted as panel data. 
atimevar must be a sorted vector of timeindices (conformable with y and x). 
If option VarianceRatios-> val is given, val is used as starting values for the 
moment equations. If val is a scalar, it is expanded to a 
list {{0.5 val, 1.5 val}, {0.5 val, 1.5 val}, ...} of equal starting ranges for all coefficients. 
Alternatively, val may begiven as list or pairs
{{s1lo, s1hi}, {s2lo, s2hi}, ...} conformable to the number of regressors.
I.e. two starting values silo < sihi have to given for each coefficient.
Note that the constant is not included 
automatically as regressor by VCEstimate.";

CIPlot::usage ="CIPlot[coeff,standarddev]renders confidence interval plots for the
VCEstimate return values coeff and standarddev. Options:CILevel (default is 0.95),
and RegressorNames (default is Automatic).";

Options[VCEstimate] = {Panel -> False, VarianceRatios -> 1.0,
      Method -> {"NelderMead", "PostProcess" -> False}, MaxIterations -> 500};

Options[CIPlot] = {CILevel->0.95, RegressorNames->Automatic, ImageSize->200, Columns->2};

(* Function definitions ---- *)

Begin["`Private`"]

Tp[x_] := Transpose[x];

VCEstimate[x_, y_, opt___Rule] :=
      Module[{M, iS, P, Q, X, tX, tXy, a, T, n, sumopt,
       sdb, vu, vr, t, u, vi, Pa, ag, pan, rsol, rstart,
       critfun, cfunval, nb = SelectedNotebook[], argmin, neval = 0},

     sumopt = setOptions[{opt}];
     pan = Panel /. sumopt;
     {T, t, n} = setDimensions[y, x, pan];
     rstart = setVarRat[VarianceRatios /. sumopt, n];
     P = makeP[n, t];
     X = makeX[x, pan];
     tXX = (tX = Tp[X]).X;
     tXy = tX . y;
     r = Table[Unique[],{n}];

     critfun[r_/; VectorQ[r,NumericQ]]:= (
       iS = makeS[1.0/r, t];
       M = tXX + Tp[P].iS.P;
       a = LinearSolve[M, tXy, Method->Cholesky];
       Pa = P.a; u = y - X.a;
       Q = u.u + Pa.iS.Pa;
       Log[Det[M]] + (T-n) Log[Q] + (T-1) Tr[Log[r]]);

     argmin = r/. Last[NMinimize[{cfunval = critfun[r],
                Sequence @@ Thread[Power[10.0,-10.0] <= r]},
       MapThread[Prepend,{rstart,r}],
       FilterOptions[NMinimize, sumopt],
       StepMonitor :> monitor[neval++, nb, r, cfunval]]];

     critfun[argmin];

     vu = Q / (T-n);
     sdb = Partition[
       Sqrt[MapIndexed[Part[#1, First[#2]] &, vu*Inverse[M]]], n];

     {Tp[Partition[a, n]], Tp[sdb], Sqrt[vu], Sqrt[vu argmin]}];


(* Helper functions ------------------------------------------------------ *)

makeQ[t_] := SparseArray[
  {{i_, i_} -> -1, {i_, j_} /; j == i + 1 -> 1},
  {t - 1,t}];

makeP[n_, t_] := SparseArray[
  {{i_, i_} -> -1, {i_, j_} /;
  j == i + n -> 1},
  {(t - 1) n, t n}];

(* Kronecker Product *)
kronprod[a_?MatrixQ, b_?MatrixQ] := BlockMatrix[Map[#b &, a, {2}]];

makePi[n_, T_] := With[{Q = makeQ[T]},
  Table[SparseArray[
    kronprod[Q,{ReplacePart[Table[0, {n}], 1, i]}]],
  {i, n}]];

makeS[vars_List, t_] := With[{n = Length[vars]},
      SparseArray[{i_, i_} :> 
          vars[[1 + Mod[i - 1, n]]], {(t - 1) n, (t - 1) n}]];

makeX[x_, paninfo_] := If[paninfo === False,
  makeTimeSeriesX[x], makePanelX[x, paninfo]];

makeTimeSeriesX[x_] := Module[{t, n},
  {t, n} = Dimensions[x];
  SparseArray[{i_, j_} /; (i - 1) n < j <= i n :> x[[i, j - (i - 1) n]],
  {t, t n}]];

makePanelX[x_, paninfo_] := Module[{ind, per, i, j, t, n},
    {t, n} = Dimensions[x];
    per = (# - Min[#]) &[paninfo];
    ind = MapIndexed[#1 + n {0, 1} Part[per, #2[[1]]] &,
        Table[{i, j}, {i, t}, {j, n}], {2}];
    SparseArray[Thread[Flatten[ind, 1] -> Flatten[x]]]];

setDimensions[y_, x_, paninfo_] := Module[{T, t, n},
    checkDimensions[y, x, paninfo];
    {T, n} = Dimensions[x];
    t = If[paninfo === False, T, Length[Split[paninfo]]];
    {T, t, n}];

setVarRat[vq_, n_] := Module[{r},
      Which[NumberQ[vq], Table[{0.5 vq, 1.5 vq}, {n}],
        Length[vq] == n, vq,
        True, Message[VCEstimate::"varrat"]; Abort[]]];

checkDimensions[y_, x_, paninfo_] := Block[{a := Abort[]},
      If[VectorQ[y, NumericQ] && MatrixQ[x, NumericQ],
        If[Length[y] != Length[x], Message[VCEstimate::dimxy]; a],
        Message[VCEstimate::"typexy"]; a];
      If[paninfo =!= False,
        If[! VectorQ[paninfo, NumericQ], Message[VCEstimate::"pantype"]; a,
          If[Length[paninfo] != Length[x], Message[VCEstimate::"pandim"]; a]]]];

(*
setOptions[opt_] := Block[{comp},
    comp = Complement[
        Map[First, opt],
        Map[First, Join[Options[VCEstimate],Options[NMinimize]]]];
    If[comp =!= {}, Scan[Message[VCEstimate::"opt", #] &, comp]; Abort[]];
    {Panel, VarianceRatios} /. opt /. Options[VCEstimate]];
*)
    
setOptions[useropt_]:= Block[{isect, defopt, netopt, sumopt},
    defopt = Options[VCEstimate];
    isect = Intersection[Map[First, useropt], Map[First, defopt]];
    netopt = DeleteCases[defopt, Alternatives @@ Map[Rule[#, _] &, isect]];
    sumopt = Union[useropt, netopt]];

monitor[neval_, nb_, r_, crit_]:=
     NotebookWrite[nb,
       Cell["Iteration: " <> ToString[neval] <>
         " Objective function value: " <> ToString[NumberForm[-crit,10]] <> "\n" <>
         "Variance Ratios: " <> ToString[AccountingForm[r]], "Output"], All];

(* Routine for plotting coefficient time paths together with their
   confidence intervals ------------------------------------------------ *)

CIPlot[coeff_?MatrixQ, sd_?MatrixQ, opts___Rule] := Module[
    {labels, quantile, level, names, cols,
     size, i, dp, pj, d, plist, gapfill, n = Length[coeff]},

    {level, cols, size, names} =
      {CILevel, Columns, ImageSize, RegressorNames} /. {opts} /. Options[CIPlot];

    pj = PlotJoined -> True;
    dp = DisplayFunction :> Identity;
    d  = PlotStyle -> Dashing[{0.02, 0.02}];
    quantile = Quantile[NormalDistribution[0.,1.], 1 - (1-level)/2];
    labels = If[names =!= Automatic, names,
        Table["Coeff. " <> ToString[i], {i, Length[coeff]}]];

    plist = MapThread[
      Show[
        ListPlot[#1, PlotLabel -> #3, pj, dp],
        ListPlot[#1 - quantile #2, pj, dp, d],
        ListPlot[#1 + quantile #2, pj, dp, d]] &,
      {coeff, sd, labels}];

   cols = Min[cols, n];
   gapfill = Graphics[{RGBColor[1,1,1], Rectangle[{0,0},{1,1}]}];
   plist = Join[plist, Table[gapfill, {Mod[n, cols]}]];

   Show[GraphicsArray[Partition[plist, cols]],
        ImageSize -> size * cols]];

(* Error Messages ----------------------------------------------------------------- *)
VCEstimate::dimxy = "Dimensions of arguments y and x do not match.
Check dimensions and retry.";

VCEstimate::"typexy" = "x and/or y are not a numeric matrix/vector.
Check them and retry.";

VCEstimate::"pantype" = "time variable tvar in option Panel->tvar is not a numeric vector.
Check tvar and retry";

VCEstimate::"pandim" = "time variable tvar in option Panel->tvar is not conformable
with the number of observations. Check tvar andretry.";

VCEstimate::"varrat" = "Dimensions of starting values for variance ratios do not
match dimension of the regressor matrix. Check dimensionsand retry.";

VCEstimate::"opt" = "`1` is not an valid option for VCEstimate. Check and retry.";

VCEstimate::"conv" = "Convergence achieved after `1` iterations.";

VCEstimate::"noconv" = "Convergence not achieved after `1` iterations.
Returned results are not reliable. Check the input data and the model and retry.";

SetAttributes[
  {VCEstimate, CIPlot, Panel, CILevel, Columns, RegressorNames, VarianceRatios, EstimationMethod},
  {Protected, ReadProtected}];

End[];
EndPackage[];

