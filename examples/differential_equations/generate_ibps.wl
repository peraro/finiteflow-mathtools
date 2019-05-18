(* ::Package:: *)

<<LiteIBP`


<<"topology.wl"


(* GetSeeds generates the seed integrals for each sector, and for each type of identity *)
rank=4;
maxdots=1;
GetSeeds["IBP"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["LI"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["SR"][sector_]:=LIBPGenerateSeeds[sector,{0,0},0,rank];
GetSeeds["Map"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots+1},0,rank];


(* You can decide here how many parallel kernels you want to use to generate IBPs,
   but it does not work anymore in Mathematica 11.3 and later versions.  *)
If[$VersionNumber<11.3,
  nkernels = 4;
  LaunchKernels[4];
  ParallelEvaluate[<<"init.m"];
  ParallelNeeds["LiteIBP`"];
  DistributeDefinitions[GetSeeds];
]


(* generate identities and put them in the "ibps/" directory *)
LIBPFastGenIds[family,GetSeeds,"Directory"->"ibps"]


(* an overcomplete list of UT integrals, where preferred integrals should appear towards the end of the list *)
(* the list was provided by Simone Zoia, with the help of a program by Pascal Wasser. *)
utints = {ut[1]->s23 j[family,0,1,1,0,1,1,1,0,0]+s23 j[family,0,1,1,1,1,0,1,0,0]-2 s23 j[family,0,1,1,1,1,1,1,0,-1]-s23^2 j[family,0,1,1,1,1,1,1,0,0]-s23 j[family,1,0,1,1,1,0,1,0,0]-s12 s23 j[family,1,1,0,1,1,1,1,0,0]-s23 j[family,1,1,1,0,1,1,1,-1,0]-s23^2 j[family,1,1,1,0,1,1,1,0,0]+s13 s23 j[family,1,1,1,1,1,1,0,0,0]+s23 j[family,1,1,1,1,1,1,1,-1,-1]+(-s12 s23-s13 s23) j[family,1,1,1,1,1,1,1,-1,0]+s23^2 j[family,1,1,1,1,1,1,1,0,-1]-s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[2]->s12 j[family,0,1,1,0,1,1,1,0,0]-s23 j[family,0,1,1,1,1,0,1,0,0]+2 s23 j[family,0,1,1,1,1,1,1,0,-1]+s23^2 j[family,0,1,1,1,1,1,1,0,0]+s23 j[family,1,0,1,1,1,0,1,0,0]+s12 s23 j[family,1,1,0,1,1,1,1,0,0]+s23 j[family,1,1,1,0,1,1,1,-1,0]+s23^2 j[family,1,1,1,0,1,1,1,0,0]-s13 s23 j[family,1,1,1,1,1,1,0,0,0]-s23 j[family,1,1,1,1,1,1,1,-1,-1]+(s12 s23+s13 s23) j[family,1,1,1,1,1,1,1,-1,0]-s23^2 j[family,1,1,1,1,1,1,1,0,-1]+s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[3]->s13 j[family,0,1,1,0,1,1,1,0,0]+s23 j[family,0,1,1,1,1,0,1,0,0]-2 s23 j[family,0,1,1,1,1,1,1,0,-1]-s23^2 j[family,0,1,1,1,1,1,1,0,0]-s23 j[family,1,0,1,1,1,0,1,0,0]-s12 s23 j[family,1,1,0,1,1,1,1,0,0]+(-s12-s13-s23) j[family,1,1,1,0,1,1,1,-1,0]+s13 s23 j[family,1,1,1,1,1,1,0,0,0]+s23 j[family,1,1,1,1,1,1,1,-1,-1]+(-s12 s23-s13 s23) j[family,1,1,1,1,1,1,1,-1,0]-s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[4]->s12 j[family,0,1,1,1,0,1,1,0,0]-s23 j[family,0,1,1,1,1,1,1,0,-1]-s23^2 j[family,0,1,1,1,1,1,1,0,0]-1/2 s12 j[family,1,1,0,1,0,1,1,0,0]-s12 s23 j[family,1,1,0,1,1,1,1,0,0]+1/2 (-s12-s13) j[family,1,1,1,1,0,1,1,-1,0]-1/2 s12 s23 j[family,1,1,1,1,0,1,1,0,0]+(-s12 s23-s13 s23) j[family,1,1,1,1,1,1,1,-1,0]-s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[5]->s23 j[family,0,1,1,1,0,1,1,0,0]-s23 j[family,0,1,1,1,1,1,1,0,-1]-s23^2 j[family,0,1,1,1,1,1,1,0,0]-s12 s23 j[family,1,1,0,1,1,1,1,0,0]+(-s12 s23-s13 s23) j[family,1,1,1,1,1,1,1,-1,0]-s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[6]->s13 j[family,0,1,1,1,0,1,1,0,0]+s23 j[family,0,1,1,1,1,1,1,0,-1]+s23^2 j[family,0,1,1,1,1,1,1,0,0]+s12 s23 j[family,1,1,0,1,1,1,1,0,0]+(s12 s23+s13 s23) j[family,1,1,1,1,1,1,1,-1,0]+s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[7]->(s13 s23-s23^2) j[family,0,1,1,1,1,1,1,0,0]-s12 s23 j[family,1,1,0,1,1,1,1,0,0]+(-s12 s23-s13 s23) j[family,1,1,1,1,1,1,1,-1,0]-s12 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[8]->(s12 s23+s23^2) j[family,1,1,1,0,1,1,1,0,0]-s23^2 j[family,1,1,1,1,1,1,1,0,-1],ut[9]->s13 s23^2 j[family,1,1,1,1,1,1,1,0,0],ut[10]->(s12 s23+s23^2) j[family,0,1,1,1,1,1,1,0,0],ut[11]->s23 j[family,1,1,0,0,1,1,1,0,0]-s23 j[family,1,1,0,1,1,1,1,0,-1],ut[12]->(s12 s23+s13 s23) j[family,1,1,1,1,1,1,0,0,0],ut[13]->(s12 s23+s13 s23) j[family,1,1,1,1,1,0,1,0,0],ut[14]->s23^2 j[family,1,0,1,1,1,1,1,0,0],ut[15]->s13 s23 j[family,1,1,1,1,0,1,1,0,0],ut[16]->s13 s23 j[family,1,1,1,0,1,1,1,0,0],ut[17]->s13 s23 j[family,1,1,0,1,1,1,1,0,0],ut[18]->(s12+s13) j[family,1,1,0,1,1,0,1,0,0],ut[19]->(s12+s13) j[family,0,1,1,1,1,0,1,0,0],ut[20]->s13 j[family,1,1,0,1,0,1,1,0,0],ut[21]->s12 j[family,1,1,0,0,1,1,1,0,0],ut[22]->s23 j[family,1,0,1,1,0,1,1,0,0],ut[23]->s23 j[family,1,0,1,0,1,1,1,0,0],ut[24]->((-1+2 eps) (-2+3 eps) (-1+3 eps) j[family,0,0,1,1,0,0,1,0,0])/(2 eps^3 s23)};


(* write the additional identities with the definition of the UT integrals *)
Export["ibps/utints_def.m",utints];
Export["ibps/ids_family_utints.mx",(-#[[1]]+#[[2]])&/@utints,"MX"];
(* ...and the list of standard integrals appearing in them *)
Export["ibps/ints_family_utints.mx",LIBPIntegralsIn[utints],"MX"];


(* serialize in JSON format *)
LIBPSerializeFastIds[
  FileNames["ibps/ids_family_*.mx"],
  s23->1,
  {eps,s12,s13},
  "ExtraInts"->(First/@utints)]
