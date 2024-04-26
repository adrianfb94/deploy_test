/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2017 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include "cdo_int.h"
#include "error.h"
#include "modules.h"
#include "operator_help.h"
#include <cdi.h>

#include <dirent.h>
#include <dlfcn.h>
#include <regex>
#include <set>
#include <string>
// for std::sort()
#include <algorithm>


/* \cond */
void *Adisit(void *argument);
void *Afterburner(void *argument);
void *Arith(void *argument);
void *Arithc(void *argument);
void *Arithdays(void *argument);
void *Arithlat(void *argument);
void *Cat(void *argument);
void *CDItest(void *argument);
void *CDIread(void *argument);
void *CDIwrite(void *argument);
void *Change(void *argument);
void *Change_e5slm(void *argument);
void *Cloudlayer(void *argument);
void *CMOR(void *argument);
void *CMOR_lite(void *argument);
void *CMOR_table(void *argument);
void *Collgrid(void *argument);
void *Command(void *argument);
void *Comp(void *argument);
void *Compc(void *argument);
void *Complextorect(void *argument);
void *Cond(void *argument);
void *Cond2(void *argument);
void *Condc(void *argument);
void *Consecstat(void *argument);
void *Copy(void *argument);
void *Deltat(void *argument);
void *Deltime(void *argument);
void *Derivepar(void *argument);
void *Detrend(void *argument);
void *Diff(void *argument);
void *Distgrid(void *argument);
void *Duplicate(void *argument);
void *Echam5ini(void *argument);
void *Enlarge(void *argument);
void *Enlargegrid(void *argument);
void *Ensstat(void *argument);
void *Ensstat3(void *argument);
void *Ensval(void *argument);
void *Eofcoeff(void *argument);
void *Eofcoeff3d(void *argument);
void *EOFs(void *argument);
void *EOF3d(void *argument);
void *EstFreq(void *argument);
void *Expr(void *argument);
void *FC(void *argument);
void *Filedes(void *argument);
void *Fillmiss(void *argument);
void *Filter(void *argument);
void *Fldrms(void *argument);
void *Fldstat(void *argument);
void *Fldstat2(void *argument);
void *Fourier(void *argument);
void *Gengrid(void *argument);
void *Gradsdes(void *argument);
void *Gridboxstat(void *argument);
void *Gridcell(void *argument);
void *Gridsearch(void *argument);
void *Harmonic(void *argument);
void *Histogram(void *argument);
void *Importamsr(void *argument);
void *Importbinary(void *argument);
void *Importcmsaf(void *argument);
void *Importobs(void *argument);
void *Info(void *argument);
void *Input(void *argument);
void *Intgrid(void *argument);
void *Intgridtraj(void *argument);
void *Intlevel(void *argument);
void *Intlevel3d(void *argument);
void *Inttime(void *argument);
void *Intntime(void *argument);
void *Intyear(void *argument);
void *Invert(void *argument);
void *Invertlev(void *argument);
void *Isosurface(void *argument);
void *Log(void *argument);
void *MapReduce(void *argument);
void *Maskbox(void *argument);
void *Mastrfu(void *argument);
void *Math(void *argument);
void *Merge(void *argument);
void *Mergegrid(void *argument);
void *Mergetime(void *argument);
void *Merstat(void *argument);
void *Monarith(void *argument);
void *Mrotuv(void *argument);
void *Mrotuvb(void *argument);
void *Ninfo(void *argument);
void *Nmldump(void *argument);
void *Output(void *argument);
void *Outputgmt(void *argument);
void *Pack(void *argument);
void *Pardup(void *argument);
void *Pinfo(void *argument);
void *Pressure(void *argument);
void *Regres(void *argument);
void *Remap(void *argument);
void *Remapeta(void *argument);
void *Replace(void *argument);
void *Replacevalues(void *argument);
void *Rotuv(void *argument);
void *Rhopot(void *argument);
void *Runpctl(void *argument);
void *Runstat(void *argument);
void *Samplegridicon(void *argument);
void *Seascount(void *argument);
void *Seaspctl(void *argument);
void *Seasstat(void *argument);
void *Selbox(void *argument);
void *Selgridcell(void *argument);
void *Select(void *argument);
void *Selvar(void *argument);
void *Seloperator(void *argument);
void *Selrec(void *argument);
void *Seltime(void *argument);
void *Set(void *argument);
void *Setattribute(void *argument);
void *Setbox(void *argument);
void *Setgatt(void *argument);
void *Setgrid(void *argument);
void *Sethalo(void *argument);
void *Setmiss(void *argument);
void *Setpartab(void *argument);
void *Setrcaname(void *argument);
void *Settime(void *argument);
void *Setzaxis(void *argument);
void *Shiftxy(void *argument);
void *Showinfo(void *argument);
void *Showattribute(void *argument);
void *Sinfo(void *argument);
void *Smooth(void *argument);
void *Sort(void *argument);
void *Sorttimestamp(void *argument);
void *Specinfo(void *argument);
void *Spectral(void *argument);
void *Spectrum(void *argument);
void *Split(void *argument);
void *Splitrec(void *argument);
void *Splitsel(void *argument);
void *Splittime(void *argument);
void *Splityear(void *argument);
void *Subtrend(void *argument);
void *Tee(void *argument);
void *Template1(void *argument);
void *Template2(void *argument);
void *Test(void *argument);
void *Test2(void *argument);
void *Testdata(void *argument);
void *Tests(void *argument);
void *Timsort(void *argument);
void *Timcount(void *argument);
void *Timcumsum(void *argument);
void *Timpctl(void *argument);
void *Timselpctl(void *argument);
void *Timselstat(void *argument);
void *XTimstat(void *argument);
void *Timstat(void *argument);
void *Timstat2(void *argument);
void *Timstat3(void *argument);
void *Tinfo(void *argument);
void *Tocomplex(void *argument);
void *Transpose(void *argument);
void *Trend(void *argument);
void *Trms(void *argument);
void *Tstepcount(void *argument);
void *Vargen(void *argument);
void *Varrms(void *argument);
void *Vertintml(void *argument);
void *Vertintap(void *argument);
void *Vertstat(void *argument);
void *Vertcum(void *argument);
void *Vertwind(void *argument);
void *Verifygrid(void *argument);
void *Wind(void *argument);
void *Writegrid(void *argument);
void *Writerandom(void *argument);
void *YAR(void *argument);
void *Yearmonstat(void *argument);
void *Ydayarith(void *argument);
void *Ydaypctl(void *argument);
void *Ydaystat(void *argument);
void *Ydrunpctl(void *argument);
void *Ydrunstat(void *argument);
void *Yhourarith(void *argument);
void *Yhourstat(void *argument);
void *Ymonarith(void *argument);
void *Ymonpctl(void *argument);
void *Ymonstat(void *argument);
void *Yseaspctl(void *argument);
void *Yseasstat(void *argument);
void *Zonstat(void *argument);

void *EcaCfd(void *argument);
void *EcaCsu(void *argument);
void *EcaCwdi(void *argument);
void *EcaCwfi(void *argument);
void *EcaEtr(void *argument);
void *EcaFd(void *argument);
void *EcaGsl(void *argument);
void *EcaHd(void *argument);
void *EcaHwdi(void *argument);
void *EcaHwfi(void *argument);
void *EcaId(void *argument);
void *EcaSu(void *argument);
void *EcaTr(void *argument);
void *EcaTg10p(void *argument);
void *EcaTg90p(void *argument);
void *EcaTn10p(void *argument);
void *EcaTn90p(void *argument);
void *EcaTx10p(void *argument);
void *EcaTx90p(void *argument);

void *EcaCdd(void *argument);
void *EcaCwd(void *argument);
void *EcaRr1(void *argument);
void *EcaPd(void *argument);
void *EcaR75p(void *argument);
void *EcaR75ptot(void *argument);
void *EcaR90p(void *argument);
void *EcaR90ptot(void *argument);
void *EcaR95p(void *argument);
void *EcaR95ptot(void *argument);
void *EcaR99p(void *argument);
void *EcaR99ptot(void *argument);
void *EcaRx1day(void *argument);
void *EcaRx5day(void *argument);
void *EcaSdii(void *argument);

void *Fdns(void *argument);
void *Strwin(void *argument);
void *Strbre(void *argument);
void *Strgal(void *argument);
void *Hurr(void *argument);

// void *Hi(void *argument);
void *Wct(void *argument);

void *Magplot(void *argument);
void *Magvector(void *argument);
void *Maggraph(void *argument);

// HIRLAM_EXTENSIONS
void *Selmulti(void *argument);   // "selmulti", "delmulti"
void *WindTrans(void *argument);  // "uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"
void *Samplegrid(void *argument); // "samplegrid", "subgrid"


/* clang-format off */
#define  AdisitOperators        {"adisit", "adipot"}
#define  AfterburnerOperators   {"after"}
#define  ArithOperators         {"add",  "sub",  "mul",  "div", "min", "max", "atan2", "setmiss"}
#define  ArithcOperators        {"addc", "subc", "mulc", "divc", "mod"}
#define  ArithdaysOperators     {"muldpm", "divdpm", "muldpy", "divdpy", "muldoy"}
#define  ArithlatOperators      {"mulcoslat", "divcoslat"}
#define  CatOperators           {"cat"}
#define  CDItestOperators       {"ncopy"}
#define  CDIreadOperators       {"cdiread"}
#define  CDIwriteOperators      {"cdiwrite"}
#define  ChangeOperators        {"chcode", "chtabnum", "chparam", "chname", "chunit", "chlevel", "chlevelc", "chlevelv", "chltype"}
#define  Change_e5slmOperators  {"change_e5slm", "change_e5lsm", "change_e5mask"}
#define  CloudlayerOperators    {"cloudlayer"}
#define  CMOROperators          {"cmor"}
#define  CMORliteOperators      {"cmorlite"}
#define  CMORtableOperators     {"dump_cmor_table", "conv_cmor_table"}
#define  CollgridOperators      {"collgrid"}
#define  CommandOperators       {"command", "com", "cmd"}
#define  CompOperators          {"eq",  "ne",  "le",  "lt",  "ge",  "gt"}
#define  CompcOperators         {"eqc", "nec", "lec", "ltc", "gec", "gtc"}
#define  ComplextorectOperators {"complextorect"}
#define  CondOperators          {"ifthen",  "ifnotthen"}
#define  Cond2Operators         {"ifthenelse"}
#define  CondcOperators         {"ifthenc", "ifnotthenc"}
#define  ConsecstatOperators    {"consects", "consecsum"}
#define  CopyOperators          {"copy", "selall", "szip"}
#define  DeltatOperators        {"deltat"}
#define  DeltimeOperators       {"delday", "del29feb"}
#define  DeriveparOperators     {"gheight", "sealevelpressure"}
#define  DetrendOperators       {"detrend"}
#define  DiffOperators          {"diff", "diffp", "diffn", "diffc"}
#define  DistgridOperators      {"distgrid"}
#define  DuplicateOperators     {"duplicate"}
#define  Echam5iniOperators     {"import_e5ml", "import_e5res", "export_e5ml", "export_e5res"}
#define  EnlargeOperators       {"enlarge"}
#define  EnlargegridOperators   {"enlargegrid"}
#define  EnsstatOperators       {"ensrange", "ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensvar1", "ensstd", "ensstd1", "enspctl"}
#define  Ensstat3Operators      {"ensrkhistspace", "ensrkhisttime", "ensroc"}
#define  EnsvalOperators        {"enscrps", "ensbrs"}
#define  EofcoeffOperators      {"eofcoeff"}
#define  Eofcoeff3dOperators    {"eofcoeff3d"}
#define  EOFsOperators          {"eof", "eofspatial", "eoftime"}
#define  EOF3dOperators         {"eof3d","eof3dspatial","eof3dtime"}
#define  EstFreqOperators       {"estfreq"}
#define  ExprOperators          {"expr", "exprf", "aexpr", "aexprf"}
#define  FCOperators            {"fc2sp", "sp2fc", "fc2gp", "gp2fc"}
#define  FiledesOperators       {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vct2", "codetab", \
                                 "vlist", "partab", "partab2", "spartab"}
#define  FillmissOperators      {"fillmiss", "fillmiss2"}
#define  FilterOperators        {"bandpass", "highpass", "lowpass"}
#define  FldrmsOperators        {"fldrms"}
#define  FldstatOperators       {"fldrange", "fldmin", "fldmax", "fldsum", "fldmean", "fldavg", "fldstd", "fldstd1", "fldvar", "fldvar1", "fldpctl"}
#define  FldcorOperators        {"fldcor"}
#define  FldcovarOperators      {"fldcovar"}
#define  FourierOperators       {"fourier"}
#define  GengridOperators       {"gengrid"}
#define  GradsdesOperators      {"gradsdes", "dumpmap"}
#define  GridboxstatOperators   {"gridboxrange", "gridboxmin", "gridboxmax", "gridboxsum", "gridboxmean", "gridboxavg", "gridboxstd", "gridboxstd1", "gridboxvar", "gridboxvar1"}
#define  GridcellOperators      {"gridarea", "gridweights", "gridmask", "griddx", "griddy"}
#define  GridsearchOperators    {"testpointsearch", "testcellsearch"}
#define  HarmonicOperators      {"harmonic"}
#define  HistogramOperators     {"histcount", "histsum", "histmean", "histfreq"}
#define  ImportamsrOperators    {"import_amsr"}
#define  ImportbinaryOperators  {"import_binary"}
#define  ImportcmsafOperators   {"import_cmsaf"}
#define  ImportobsOperators     {"import_obs"}
#define  InfoOperators          {"info", "infop", "infon", "infoc", "xinfon", "map"}
#define  InputOperators         {"input", "inputsrv", "inputext"}
#define  IntgridOperators       {"intgridbil", "intpoint", "interpolate", "boxavg", "thinout"}
#define  IntgridtrajOperators   {"intgridtraj"}
#define  IntlevelOperators      {"intlevel", "intlevelx"}
#define  Intlevel3dOperators    {"intlevel3d", "intlevelx3d"}
#define  InttimeOperators       {"inttime"}
#define  IntntimeOperators      {"intntime"}
#define  IntyearOperators       {"intyear"}
#define  InvertOperators        {"invertlat", "invertlon", "invertlatdes", "invertlondes", "invertlatdata", "invertlondata"}
#define  InvertlevOperators     {"invertlev"}
#define  IsosurfaceOperators    {"isosurface"}
#define  LogOperators           {"dumplogs", "daylogs", "monlogs", "dumplogo", "snamelogo", "scalllogo", "smemlogo", "stimelogo", "sperclogo"}
#define  MapReduceOperators     {"reducegrid"}
#define  MaskboxOperators       {"masklonlatbox", "maskindexbox"}
#define  MaskregionOperators    {"maskregion"}
#define  MastrfuOperators       {"mastrfu"}
#define  MathOperators          {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", \
                                 "cos", "tan", "asin", "acos", "atan", "pow", "reci"}
#define  MergeOperators         {"merge"}
#define  MergegridOperators     {"mergegrid"}
#define  MergetimeOperators     {"mergetime"}
#define  MerstatOperators       {"merrange", "mermin", "mermax", "mersum", "mermean", "meravg", "merstd", "merstd1", "mervar", "mervar1", "merpctl"}
#define  MonarithOperators      {"monadd", "monsub", "monmul", "mondiv"}
#define  MrotuvOperators        {"mrotuv"}
#define  MrotuvbOperators       {"mrotuvb"}
#define  NinfoOperators         {"nyear", "nmon", "ndate", "ntime", "ncode", "npar", "nlevel", "ngridpoints", "ngrids"}
#define  NmldumpOperators       {"nmldump", "kvldump"}
#define  OutputOperators        {"output", "outputint", "outputsrv", "outputext", "outputf", "outputts", \
                                 "outputfld", "outputarr", "outputxyz"}
#define  OutputtabOperators     {"outputtab"}
#define  OutputgmtOperators     {"gmtxyz", "gmtcells", "outputcenter2", "outputcentercpt", \
                                 "outputboundscpt", "outputvector", "outputtri", "outputvrml"}
#define  PackOperators          {"pack"}
#define  PardupOperators        {"pardup", "parmul"}
#define  PinfoOperators         {"pinfo", "pinfov"}
#define  PressureOperators      {"pressure_fl", "pressure_hl", "deltap"}
#define  RegresOperators        {"regres"}
#define  RemapOperators         {"remap"}
#define  RemapbilOperators      {"remapbil", "genbil"}
#define  RemapbicOperators      {"remapbic", "genbic"}
#define  RemapnnOperators       {"remapnn", "gennn"}
#define  RemapdisOperators      {"remapdis", "gendis"}
#define  RemapyconOperators     {"remapycon", "genycon"}
#define  RemapconOperators      {"remapcon", "gencon"}
#define  Remapcon2Operators     {"remapcon2", "gencon2"}
#define  RemaplafOperators      {"remaplaf", "genlaf"}
#define    RemapgridOperators   {"remapsum"}
#define  RemapetaOperators      {"remapeta", "remapeta_s", "remapeta_z"}
#define  ReplaceOperators       {"replace"}
#define  ReplacevaluesOperators {"setvals", "setrtoc", "setrtoc2"}
#define  RhopotOperators        {"rhopot"}
#define  RotuvOperators         {"rotuvb"}
#define  RunpctlOperators       {"runpctl"}
#define  RunstatOperators       {"runrange", "runmin", "runmax", "runsum", "runmean", "runavg", "runstd", "runstd1", "runvar", "runvar1"}
#define  SamplegridiconOperators {"samplegridicon"}
#define  SeascountOperators     {"seascount"}
#define  SeaspctlOperators      {"seaspctl"}
#define  SeasstatOperators      {"seasrange", "seasmin", "seasmax", "seassum", "seasmean", "seasavg", "seasstd", "seasstd1", "seasvar", "seasvar1"}
#define  SelboxOperators        {"sellonlatbox", "selindexbox"}
#define  SelgridcellOperators   {"selgridcell", "delgridcell"}
#define  SelectOperators        {"select", "delete"}
#define  SelvarOperators        {"selparam", "selcode", "selname", "selstdname", "sellevel", "sellevidx", "selgrid", \
                                 "selzaxis", "selzaxisname", "seltabnum", "delparam", "delcode", "delname", "selltype"}
#define  SeloperatorOperators   {"seloperator"}
#define  SelrecOperators        {"selrec"}
#define  SeltimeOperators       {"seltimestep", "selyear", "selseason", "selmonth", "selday", "selhour", "seldate", \
                                 "seltime", "selsmon"}
#define  SetOperators           {"setcode", "setparam", "setname", "setunit", "setlevel", "setltype", "settabnum"}
#define  SetattributeOperators  {"setattribute"}
#define  SetboxOperators        {"setclonlatbox", "setcindexbox"}
#define  SetgattOperators       {"setgatt", "setgatts"}
#define  SetgridOperators       {"setgrid", "setgridtype", "setgridarea", "setgridmask", "unsetgridmask", "setgridnumber", "setgriduri", "usegridnumber"}
#define  SethaloOperators       {"sethalo", "tpnhalo"}
#define  SetmissOperators       {"setmissval", "setctomiss", "setmisstoc", "setrtomiss", "setvrange"}
#define  SetmisstonnOperators   {"setmisstonn", "setmisstodis"}
#define  SetcodetabOperators    {"setcodetab"}
#define  SetpartabOperators     {"setpartabc", "setpartabp", "setpartabn"}
#define  SetrcanameOperators    {"setrcaname"}
#define  SettimeOperators       {"setyear", "setmon", "setday", "setdate", "settime", "settunits", \
                                 "settaxis", "settbounds", "setreftime", "setcalendar", "shifttime"}
#define  SetzaxisOperators      {"setzaxis", "genlevelbounds"}
#define  ShiftxyOperators       {"shiftx", "shifty"}
#define  ShowinfoOperators      {"showyear", "showmon", "showdate", "showtime", "showtimestamp", "showcode", "showunit", \
                                 "showparam", "showname", "showstdname", "showlevel", "showltype", "showformat", "showgrid", "showatts", "showattsglob"}
#define  ShowattributeOperators {"showattribute", "showattsvar"}
#define  SinfoOperators         {"sinfo", "sinfop", "sinfon", "sinfoc", "seinfo", "seinfop", "seinfon", "seinfoc"}
#define  SmoothOperators        {"smooth", "smooth9"}
#define  SortOperators          {"sortcode", "sortparam", "sortname", "sortlevel"}
#define  SorttimestampOperators {"sorttimestamp", "sorttaxis"}
#define  SpecinfoOperators      {"specinfo"}
#define  SpectralOperators      {"gp2sp", "gp2spl", "sp2gp", "sp2gpl", "sp2sp", "spcut"}
#define  SpectrumOperators      {"spectrum"}
#define  SplitOperators         {"splitcode", "splitparam", "splitname", "splitlevel", "splitgrid", "splitzaxis", "splittabnum"}
#define  SplitrecOperators      {"splitrec"}
#define  SplitselOperators      {"splitsel"}
#define  SplittimeOperators     {"splithour", "splitday", "splitmon", "splitseas"}
#define  SplityearOperators     {"splityear", "splityearmon"}
#define  SubtrendOperators      {"subtrend"}
#define  TeeOperators           {"tee"}
#define  Template1Operators     {"template1"}
#define  Template2Operators     {"template2"}
#define  TestOperators          {"test"}
#define  Test2Operators         {"test2"}
#define  TestdataOperators      {"testdata"}
#define  TestsOperators         {"normal", "studentt", "chisquare", "beta", "fisher"}
#define  TimsortOperators       {"timsort"}
#define  TimcountOperators      {"timcount"}
#define    YearcountOperators   {"yearcount"}
#define    MoncountOperators    {"moncount"}
#define    DaycountOperators    {"daycount"}
#define    HourcountOperators   {"hourcount"}
#define  TimcumsumOperators     {"timcumsum"}
#define  TimpctlOperators       {"timpctl"}
#define    YearpctlOperators    {"yearpctl"}
#define    MonpctlOperators     {"monpctl"}
#define    DaypctlOperators     {"daypctl"}
#define    HourpctlOperators    {"hourpctl"}
#define  TimselpctlOperators    {"timselpctl"}
#define  TimselstatOperators    {"timselrange", "timselmin", "timselmax", "timselsum", "timselmean", "timselavg", "timselvar", "timselvar1", "timselstd", "timselstd1"}
#define  XTimstatOperators      {"xtimmin",  "xtimmax",  "xtimsum",  "xtimmean",  "xtimavg",  "xtimvar",  "xtimvar1",  "xtimstd",  "xtimstd1", \
                                 "xyearmin", "xyearmax", "xyearsum", "xyearmean", "xyearavg", "xyearvar", "xyearvar1", "xyearstd", "xyearstd1", \
                                 "xmonmin",  "xmonmax",  "xmonsum",  "xmonmean",  "xmonavg",  "xmonvar",  "xmonvar1",  "xmonstd",  "xmonstd1"}
#define  TimstatOperators       {"timrange",  "timmin",  "timmax",  "timsum",  "timmean",  "timavg",  "timvar",  "timvar1",  "timstd",  "timstd1"}
#define    YearstatOperators    {"yearrange", "yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearvar1", "yearstd", "yearstd1"}
#define    MonstatOperators     {"monrange",  "monmin",  "monmax",  "monsum",  "monmean",  "monavg",  "monvar",  "monvar1",  "monstd",  "monstd1"}
#define    DaystatOperators     {"dayrange",  "daymin",  "daymax",  "daysum",  "daymean",  "dayavg",  "dayvar",  "dayvar1",  "daystd",  "daystd1"}
#define    HourstatOperators    {"hourrange", "hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourvar1", "hourstd", "hourstd1"}
#define  TimcorOperators        {"timcor"}
#define  TimcovarOperators      {"timcovar"}
#define  Timstat3Operators      {"meandiff2test", "varquot2test"}
#define  TinfoOperators         {"tinfo"}
#define  TocomplexOperators     {"retocomplex", "imtocomplex"}
#define  TransposeOperators     {"transxy"}
#define  TrendOperators         {"trend"}
#define  TrmsOperators          {"trms"}
#define  TstepcountOperators    {"tstepcount"}
#define  VargenOperators        {"random", "const", "sincos", "coshill", "for", "topo", "temp", "mask", "stdatm"}
#define  VarrmsOperators        {"varrms"}
#define  VertintmlOperators     {"ml2pl", "ml2hl", "ml2plx", "ml2hlx", "ml2pl_lp", "ml2hl_lp", "ml2plx_lp", "ml2hlx_lp"}
#define  VertintapOperators     {"ap2pl", "ap2plx", "ap2pl_lp", "ap2plx_lp", "ap2hl", "ap2hlx"}
#define  VertstatOperators      {"vertrange", "vertmin", "vertmax", "vertsum", "vertint", "vertmean", "vertavg", "vertstd", "vertstd1", "vertvar", "vertvar1"}
#define  VertcumOperators       {"vertcum", "vertcumhl"}
#define  VertwindOperators      {"vertwind"}
#define  VerifygridOperators    {"verifygrid"}
#define  WindOperators          {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl", "dv2ps"}
#define  WritegridOperators     {"writegrid"}
#define  WriterandomOperators   {"writerandom"}
#define  YAROperators           {"yarbil", "yarnn", "yarcon"}
#define  YearmonstatOperators   {"yearmonmean", "yearmonavg"}
#define  YdayarithOperators     {"ydayadd", "ydaysub", "ydaymul", "ydaydiv"}
#define  YdaypctlOperators      {"ydaypctl"}
#define  YdaystatOperators      {"ydayrange", "ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydaystd", "ydaystd1", "ydayvar", "ydayvar1"}
#define  YdrunpctlOperators     {"ydrunpctl"}
#define  YdrunstatOperators     {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunstd", "ydrunstd1", "ydrunvar", "ydrunvar1"}
#define  YhourarithOperators    {"yhouradd", "yhoursub", "yhourmul", "yhourdiv"}
#define  YhourstatOperators     {"yhourrange", "yhourmin", "yhourmax", "yhoursum", "yhourmean", "yhouravg", "yhourstd", "yhourstd1", "yhourvar", "yhourvar1"}
#define  YmonarithOperators     {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
#define  YseasarithOperators    {"yseasadd", "yseassub", "yseasmul", "yseasdiv"}
#define  YmonpctlOperators      {"ymonpctl"}
#define  YmonstatOperators      {"ymonrange", "ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonstd", "ymonstd1", "ymonvar", "ymonvar1"}
#define  YseaspctlOperators     {"yseaspctl"}
#define  YseasstatOperators     {"yseasrange", "yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasstd", "yseasstd1", "yseasvar", "yseasvar1"}
#define  ZonstatOperators       {"zonrange", "zonmin", "zonmax", "zonsum", "zonmean", "zonavg", "zonstd", "zonstd1", "zonvar", "zonvar1", "zonpctl"}

#define  EcaCfdOperators        {"eca_cfd"}
#define  EcaCsuOperators        {"eca_csu"}
#define  EcaCwfiOperators       {"eca_cwfi"}
#define  EcaHwdiOperators       {"eca_hwdi"}
#define  EcaEtrOperators        {"eca_etr"}
#define  EcaFdOperators         {"eca_fd"}
#define  EcaGslOperators        {"eca_gsl"}
#define  EcaHdOperators         {"eca_hd"}
#define  EcaCwdiOperators       {"eca_cwdi"}
#define  EcaHwfiOperators       {"eca_hwfi"}
#define  EcaIdOperators         {"eca_id"}
#define  EcaSuOperators         {"eca_su"}
#define  EcaTrOperators         {"eca_tr"}
#define  EcaTg10pOperators      {"eca_tg10p"}
#define  EcaTg90pOperators      {"eca_tg90p"}
#define  EcaTn10pOperators      {"eca_tn10p"}
#define  EcaTn90pOperators      {"eca_tn90p"}
#define  EcaTx10pOperators      {"eca_tx10p"}
#define  EcaTx90pOperators      {"eca_tx90p"}

#define  EcaCddOperators        {"eca_cdd"}
#define  EcaCwdOperators        {"eca_cwd"}
#define  EcaRr1Operators        {"eca_rr1"}
/*
#define  EcaR10mmOperators      {"eca_r10mm"}
#define  EcaR20mmOperators      {"eca_r20mm"}
*/
#define  EcaPdOperators         {"eca_pd", "eca_r10mm", "eca_r20mm"}
#define  EcaR75pOperators       {"eca_r75p"}
#define  EcaR75ptotOperators    {"eca_r75ptot"}
#define  EcaR90pOperators       {"eca_r90p"}
#define  EcaR90ptotOperators    {"eca_r90ptot"}
#define  EcaR95pOperators       {"eca_r95p"}
#define  EcaR95ptotOperators    {"eca_r95ptot"}
#define  EcaR99pOperators       {"eca_r99p"}
#define  EcaR99ptotOperators    {"eca_r99ptot"}
#define  EcaRx1dayOperators     {"eca_rx1day"}
#define  EcaRx5dayOperators     {"eca_rx5day"}
#define  EcaSdiiOperators       {"eca_sdii"}

#define  FdnsOperators          {"fdns"}

#define  StrwinOperators        {"strwin"}
#define  StrbreOperators        {"strbre"}
#define  StrgalOperators        {"strgal"}
#define  HurrOperators          {"hurr"}

#define  HiOperators            {"hi"}
#define  WctOperators           {"wct"}

#define  MagplotOperators       {"contour", "shaded", "grfill"}
#define  MagvectorOperators     {"vector"}
#define  MaggraphOperators      {"graph"}

// HIRLAM_EXTENSIONS
#define  SelmultiOperators      {"selmulti", "delmulti", "changemulti"}
#define  WindTransOperators     {"uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"}
#define  SamplegridOperators    {"samplegrid", "subgrid"}


/* clang-format on */
/* \endcond */

/**
 * @param a pointer to a string/substring
 * @param b pointer to a string/substring
 * @param alen length of string a
 * @param blen lenght of string b
 * @retval true if a is similar to b
 * @retval false if a is not similar to b
 *
 * Recursive function for finding substrings of a operator name that match other
 * operators.
 */

static
bool similar(const char *a, const char *b, unsigned long alen, unsigned long blen) {
    if (alen > 2 && blen > 2 && strstr(b, a))
        return true;

    while (*a && *b && *a == *b) {
        a++;
        b++;
    }
    if (!*a && !*b)
        return true;
    /*
      printf("%d %d %s %s\n", alen, blen, a, b);
    */
    if (alen >= 2 && blen >= 1 && *a && similar(a + 1, b, alen - 2, blen - 1))
        return true;

    if (alen >= 1 && blen >= 2 && *b && similar(a, b + 1, alen - 1, blen - 2))
        return true;

    return false;
}

/**
 * @param original string tested for similarity to \p other
 * @param other string that \p original will be compared to
 * @retval true if original and other are similar
 * @retval false if not
 *
 * Wrapper function for #similar() to parse c++ strings to c strings
 */
static
bool similar(std::string original, std::string other) {
    return (similar(original.c_str(), other.c_str(), original.size(), other.size()));
}

/**
 * @param operatorName operator name
 * @retval true if #modules_map contains \p operatorName
 * @retval false if not
 */
static
bool operator_name_exists(std::string operatorName) {
    if (modules_map.find(operatorName) != modules_map.end()) {
        return true;
    }
    if (aliases.find(operatorName) != aliases.end()) {
        return true;
    }
    return false;
}

/**
 * @param moduleName module name
 * @retval true if #modules contains \a moduleName
 * @retval false if not
 */
static
bool module_map_contains(std::string moduleName) {
    if (modules.find(moduleName) != modules.end()) {
        return true;
    } else {
        Error("Module %s not found", moduleName.c_str());
    }
    return false;
}

/***
 * function for finding similar operator names for the given string
 * @param operatorName operator name to find similar operators for
 * @returns A string with all found names. The string is seqmented into lines
 * with a max lenght of 75 characters
 */
static
std::string find_similar(std::string operatorName) {
    std::string found_similar_operators = "";
    unsigned long lines = 1;
    unsigned long line_length = 105;
    if (operatorName != "") {
        // Searching for simlar operator names in operator to module map
        for (auto str : modules_map) {
            if (similar(string2lower(operatorName), str.first)) {
                if (found_similar_operators.size() + str.first.size() > lines * line_length) {
                    found_similar_operators += "\n";
                    lines++;
                }
                found_similar_operators += str.first;
                found_similar_operators += " ";
            }
        }
        // Searching for similar operator names in aliases to original map
        for (auto str : aliases) {
            if (similar(string2lower(operatorName), str.first)) {
                if (found_similar_operators.size() + str.first.size() > lines * line_length) {
                    found_similar_operators += "\n";
                    lines++;
                }
                found_similar_operators += str.first;
                found_similar_operators += " ";
            }
        }
    }
    return found_similar_operators;
}

/**
 * @param operatorName operator name.
 * @retval true if \p operatorName exists.
 * @retval false if \p operatorName is not in #modules_map
 *
 * Checks if given \p operatorName is in #modules_map. Else returns false.

 * Checks if \p operatorName is not a file.

 * If no matching operator is found checks for similar operators using
 find_similar().
 *
 *  \note If \p operatorName is a file name the program will exit.
 */
static
bool check_operator(std::string operatorName) {
    if (operator_name_exists(operatorName)) {
        return true;
    } else if (operatorName == "")
        Error("Operator name missing!");

    else {
        // Checking if the operatorname is an existing file name
        FILE *fp = fopen(operatorName.c_str(), "r");
        if (fp) {
            fclose(fp);
            fprintf(stderr, "Use commandline option -h for help.");
            Error("Operator missing, %s is a file on disk!", operatorName.c_str());
        }
        // Operator is no filename
        // Checking for similar operators
        fprintf(stderr, "Operator >%s< not found!\n", operatorName.c_str());
        fprintf(stderr, "Similar operators are:\n");
        std::string found_similar_operators = find_similar(operatorName);

        if (found_similar_operators.size() > 0) {
          std::cerr << found_similar_operators << std::endl;
        } else {
          fprintf(stderr, "(not found)\n");
        }
        exit(EXIT_FAILURE);
    }
    return false;
}

/***
 * Adds a module and its operators to cdo.
 * Adds the module to modules
 * Adds the operators of modules to modules_map
 * @param new_module newly constructed module
 * @note: if an error happens while adding the new module cdo will exit.
 */
static
void add_module(std::string module_name, modules_t new_module) {
    if (modules.find(module_name) == modules.end()) {
        modules[module_name] = new_module;
        for (std::string operatorName : new_module.operators) {
            // if the operator name is not already in the map or in the aliases
            if (!operator_name_exists(operatorName)) {
                modules_map[operatorName] = module_name;
            } else {
                Error("Tried to add operator but the operator name already exists");
            }
        }
    } else {
        Error("Module %s name already exists", module_name.c_str());
    }
}

/**
 * adds an key value pair to #modules_map with alias as key and originals name
 * as value
 * @param alias new alias to be added
 * @param original original operator name
 */
static
int add_alias(std::string alias, std::string original) {
    auto iter_original = modules_map.find(original);
    auto iter_alias = aliases.find(alias);

    if (iter_alias != aliases.end()) {
        Warning("alias %s could not be added: it already exists", alias.c_str());
        return -1;
    }

    if (iter_original == modules_map.end()) {
        Error("alias %s could not be added: operator %s does not exist", alias.c_str(),
              original.c_str());
        return -2;
    }
    if (modules_map.find(alias) != modules_map.end()) {
      Error("alias %s could not be added: alias name already exists as an operator", alias.c_str());
    }
    aliases[alias] = original;

    return 0;
}
/* clang-format off */
// stream in  -1 means: unlimited number of input streams
// stream out -1 means: usage of obase
/***
 * Initializes all hardcoded modules.
 */
void init_modules()
{
  
/*                             function        help function      operator names          mode number     num streams
                                                                                                  type       in out      */
  add_module("Adisit"        , {Adisit        , AdisitHelp        , AdisitOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Afterburner"   , {Afterburner   , AfterburnerHelp   , AfterburnerOperators   , 1 , CDI_REAL , -1 , 1  });
  add_module("Arith"         , {Arith         , ArithHelp         , ArithOperators         , 1 , CDI_REAL , 2  , 1  });
  add_module("Arithc"        , {Arithc        , ArithcHelp        , ArithcOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Arithdays"     , {Arithdays     , ArithdaysHelp     , ArithdaysOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Arithlat"      , {Arithlat      , {}                , ArithlatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Cat"           , {Cat           , CopyHelp          , CatOperators           , 1 , CDI_REAL , -1 , 1  });
  add_module("CDItest"       , {CDItest       , {}                , CDItestOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("CDIread"       , {CDIread       , {}                , CDIreadOperators       , 1 , CDI_REAL , 1  , 0  });
  add_module("CDIwrite"      , {CDIwrite      , {}                , CDIwriteOperators      , 1 , CDI_REAL , 0  , 1  });
  add_module("Change"        , {Change        , ChangeHelp        , ChangeOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Change_e5slm"  , {Change_e5slm  , {}                , Change_e5slmOperators  , 0 , CDI_REAL , 1  , 1  });
  add_module("Cloudlayer"    , {Cloudlayer    , {}                , CloudlayerOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("CMOR"          , {CMOR          , CMORHelp          , CMOROperators          , 1 , CDI_REAL , 1  , 0  });
  add_module("CMOR_lite"     , {CMOR_lite     , CMORliteHelp      , CMORliteOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("CMOR_table"    , {CMOR_table    , {}                , CMORtableOperators     , 1 , CDI_REAL , 0  , 0  });
  add_module("Collgrid"      , {Collgrid      , CollgridHelp      , CollgridOperators      , 1 , CDI_REAL , -1 , 1  });
  add_module("Command"       , {Command       , {}                , CommandOperators       , 0 , CDI_REAL , 1  , 0  });
  add_module("Comp"          , {Comp          , CompHelp          , CompOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Compc"         , {Compc         , CompcHelp         , CompcOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Complextorect" , {Complextorect , {}                , ComplextorectOperators , 1 , CDI_COMP , 1  , 2  });
  add_module("Cond"          , {Cond          , CondHelp          , CondOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Cond2"         , {Cond2         , Cond2Help         , Cond2Operators         , 1 , CDI_REAL , 3  , 1  });
  add_module("Condc"         , {Condc         , CondcHelp         , CondcOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Consecstat"    , {Consecstat    , ConsecstatHelp    , ConsecstatOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Copy"          , {Copy          , CopyHelp          , CopyOperators          , 1 , CDI_REAL , -1 , 1  });
  add_module("Deltat"        , {Deltat        , {}                , DeltatOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Deltime"       , {Deltime       , {}                , DeltimeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Derivepar"     , {Derivepar     , DeriveparHelp     , DeriveparOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Detrend"       , {Detrend       , DetrendHelp       , DetrendOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Diff"          , {Diff          , DiffHelp          , DiffOperators          , 1 , CDI_REAL , 2  , 0  });
  add_module("Distgrid"      , {Distgrid      , DistgridHelp      , DistgridOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Duplicate"     , {Duplicate     , DuplicateHelp     , DuplicateOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Echam5ini"     , {Echam5ini     , {}                , Echam5iniOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Enlarge"       , {Enlarge       , EnlargeHelp       , EnlargeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Enlargegrid"   , {Enlargegrid   , {}                , EnlargegridOperators   , 0 , CDI_REAL , 1  , 1  });
  add_module("Ensstat"       , {Ensstat       , EnsstatHelp       , EnsstatOperators       , 1 , CDI_REAL , -1 , 1  });
  add_module("Ensstat3"      , {Ensstat3      , Ensstat2Help      , Ensstat3Operators      , 1 , CDI_REAL , -1 , 1  });
  add_module("Ensval"        , {Ensval        , EnsvalHelp        , EnsvalOperators        , 1 , CDI_REAL , -1 , 1  });
  add_module("Eofcoeff"      , {Eofcoeff      , EofcoeffHelp      , EofcoeffOperators      , 1 , CDI_REAL , 2  , -1 });
  add_module("Eofcoeff3d"    , {Eofcoeff3d    , EofcoeffHelp      , Eofcoeff3dOperators    , 1 , CDI_REAL , 2  , -1 });
  add_module("EOFs"          , {EOFs          , EOFsHelp          , EOFsOperators          , 1 , CDI_REAL , 1  , 2  });
  add_module("EOF3d"         , {EOF3d         , EOFsHelp          , EOF3dOperators         , 1 , CDI_REAL , 1  , 2  });
  add_module("EstFreq"       , {EstFreq       , {}                , EstFreqOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Expr"          , {Expr          , ExprHelp          , ExprOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("FC"            , {FC            , {}                , FCOperators            , 1 , CDI_REAL , 1  , 1  });
  add_module("Filedes"       , {Filedes       , FiledesHelp       , FiledesOperators       , 1 , CDI_BOTH , 1  , 0  });
  add_module("Fillmiss"      , {Fillmiss      , {}                , FillmissOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Filter"        , {Filter        , FilterHelp        , FilterOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Fldrms"        , {Fldrms        , {}                , FldrmsOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Fldstat"       , {Fldstat       , FldstatHelp       , FldstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Fldstatcor"    , {Fldstat2      , FldcorHelp        , FldcorOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Fldstatvar"    , {Fldstat2      , FldcovarHelp      , FldcovarOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Fourier"       , {Fourier       , {}                , FourierOperators       , 1 , CDI_COMP , 1  , 1  });
  add_module("Gengrid"       , {Gengrid       , {}                , GengridOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("Gradsdes"      , {Gradsdes      , GradsdesHelp      , GradsdesOperators      , 1 , CDI_REAL , 1  , 0  });
  add_module("Gridboxstat"   , {Gridboxstat   , GridboxstatHelp   , GridboxstatOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Gridcell"      , {Gridcell      , GridcellHelp      , GridcellOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Gridsearch"    , {Gridsearch    , {}                , GridsearchOperators    , 0 , CDI_REAL , 0  , 0  });
  add_module("Harmonic"      , {Harmonic      , {}                , HarmonicOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Histogram"     , {Histogram     , HistogramHelp     , HistogramOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Importamsr"    , {Importamsr    , ImportamsrHelp    , ImportamsrOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Importbinary"  , {Importbinary  , ImportbinaryHelp  , ImportbinaryOperators  , 1 , CDI_REAL , 1  , 1  });
  add_module("Importcmsaf"   , {Importcmsaf   , ImportcmsafHelp   , ImportcmsafOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Importobs"     , {Importobs     , {}                , ImportobsOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Info"          , {Info          , InfoHelp          , InfoOperators          , 1 , CDI_BOTH , -1 , 0  });
  add_module("Input"         , {Input         , InputHelp         , InputOperators         , 1 , CDI_REAL , 0  , 1  });
  add_module("Intgrid"       , {Intgrid       , {}                , IntgridOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Intgridtraj"   , {Intgridtraj   , {}                , IntgridtrajOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Intlevel"      , {Intlevel      , IntlevelHelp      , IntlevelOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Intlevel3d"    , {Intlevel3d    , Intlevel3dHelp    , Intlevel3dOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Inttime"       , {Inttime       , InttimeHelp       , InttimeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Intntime"      , {Intntime      , InttimeHelp       , IntntimeOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Intyear"       , {Intyear       , IntyearHelp       , IntyearOperators       , 1 , CDI_REAL , 2  , -1 });
  add_module("Invert"        , {Invert        , InvertHelp        , InvertOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Invertlev"     , {Invertlev     , InvertlevHelp     , InvertlevOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Isosurface"    , {Isosurface    , {}                , IsosurfaceOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Log"           , {Log           , {}                , LogOperators           , 0 , CDI_REAL , 1  , 0  });
  add_module("MapReduce"     , {MapReduce     , MapReduceHelp     , MapReduceOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Maskbox"       , {Maskbox       , MaskboxHelp       , MaskboxOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Maskregion"    , {Maskbox       , MaskregionHelp    , MaskregionOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Mastrfu"       , {Mastrfu       , MastrfuHelp       , MastrfuOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Math"          , {Math          , MathHelp          , MathOperators          , 1 , CDI_BOTH , 1  , 1  });
  add_module("Merge"         , {Merge         , MergeHelp         , MergeOperators         , 1 , CDI_REAL , -1 , 1  });
  add_module("Mergetime"     , {Mergetime     , MergeHelp         , MergetimeOperators     , 1 , CDI_REAL , -1 , 1  });
  add_module("Mergegrid"     , {Mergegrid     , MergegridHelp     , MergegridOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Merstat"       , {Merstat       , MerstatHelp       , MerstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Monarith"      , {Monarith      , MonarithHelp      , MonarithOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Mrotuv"        , {Mrotuv        , {}                , MrotuvOperators        , 1 , CDI_REAL , 1  , 2  });
  add_module("Mrotuvb"       , {Mrotuvb       , {}                , MrotuvbOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("Ninfo"         , {Ninfo         , NinfoHelp         , NinfoOperators         , 1 , CDI_BOTH , 1  , 0  });
  add_module("Nmldump"       , {Nmldump       , {}                , NmldumpOperators       , 0 , CDI_REAL , 0  , 0  });
  add_module("Output"        , {Output        , OutputHelp        , OutputOperators        , 1 , CDI_REAL , -1 , 0  });
  add_module("Outputtab"     , {Output        , OutputtabHelp     , OutputtabOperators     , 1 , CDI_REAL , -1 , 0  });
  add_module("Outputgmt"     , {Outputgmt     , OutputgmtHelp     , OutputgmtOperators     , 1 , CDI_REAL , 1  , 0  });
  add_module("Pack"          , {Pack          , {}                , PackOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Pardup"        , {Pardup        , {}                , PardupOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Pinfo"         , {Pinfo         , {}                , PinfoOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Pressure"      , {Pressure      , {}                , PressureOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Regres"        , {Regres        , RegresHelp        , RegresOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Remap"         , {Remap         , RemapHelp         , RemapOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapbil"      , {Remap         , RemapbilHelp      , RemapbilOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapbic"      , {Remap         , RemapbicHelp      , RemapbicOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapnn"       , {Remap         , RemapnnHelp       , RemapnnOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapdis"      , {Remap         , RemapdisHelp      , RemapdisOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapycon"     , {Remap         , RemapyconHelp     , RemapyconOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapcon"      , {Remap         , RemapconHelp      , RemapconOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapcon2"     , {Remap         , Remapcon2Help     , Remapcon2Operators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remaplaf"      , {Remap         , RemaplafHelp      , RemaplafOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapgrid"     , {Remap         , {}                , RemapgridOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapeta"      , {Remapeta      , RemapetaHelp      , RemapetaOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Replace"       , {Replace       , ReplaceHelp       , ReplaceOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("Replacevalues" , {Replacevalues , ReplacevaluesHelp , ReplacevaluesOperators , 1 , CDI_REAL , 1  , 1  });
  add_module("Rhopot"        , {Rhopot        , RhopotHelp        , RhopotOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Rotuv"         , {Rotuv         , RotuvbHelp        , RotuvOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Runpctl"       , {Runpctl       , RunpctlHelp       , RunpctlOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Runstat"       , {Runstat       , RunstatHelp       , RunstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Samplegridicon", {Samplegridicon, {}                , SamplegridiconOperators, 1,  CDI_REAL,  1  , 2  });
  add_module("Seascount"     , {Seascount     , {}                , SeascountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Seaspctl"      , {Seaspctl      , SeaspctlHelp      , SeaspctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Seasstat"      , {Seasstat      , SeasstatHelp      , SeasstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Selbox"        , {Selbox        , SelboxHelp        , SelboxOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Selgridcell"   , {Selgridcell   , {}                , SelgridcellOperators   , 1 , CDI_BOTH , 1  , 1  });
  add_module("Select"        , {Select        , SelectHelp        , SelectOperators        , 1 , CDI_BOTH , -1 , 1  });
  add_module("Selvar"        , {Selvar        , SelvarHelp        , SelvarOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Selrec"        , {Selrec        , SelvarHelp        , SelrecOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Seloperator"   , {Seloperator   , {}                , SeloperatorOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Seltime"       , {Seltime       , SeltimeHelp       , SeltimeOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Set"           , {Set           , SetHelp           , SetOperators           , 1 , CDI_BOTH , 1  , 1  });
  add_module("Setattribute"  , {Setattribute  , SetattributeHelp  , SetattributeOperators  , 1 , CDI_REAL , 1  , 1  });
  add_module("Setbox"        , {Setbox        , SetboxHelp        , SetboxOperators        , 1 , CDI_REAL , 1  , 1  });
  //add_module("Setgatt"       , {Setgatt       , SetgattHelp       , SetgattOperators       , 1 , CDI_BOTH , 1  , 1  });  
  add_module("Setgrid"       , {Setgrid       , SetgridHelp       , SetgridOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Sethalo"       , {Sethalo       , SethaloHelp       , SethaloOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Setmiss"       , {Setmiss       , SetmissHelp       , SetmissOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Setmisstonn"   , {Fillmiss      , SetmissHelp       , SetmisstonnOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Setcodetab"    , {Setpartab     , SetHelp           , SetcodetabOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Setpartab"     , {Setpartab     , SetpartabHelp     , SetpartabOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Setrcaname"    , {Setrcaname    , {}                , SetrcanameOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Settime"       , {Settime       , SettimeHelp       , SettimeOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Setzaxis"      , {Setzaxis      , SetzaxisHelp      , SetzaxisOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Shiftxy"       , {Shiftxy       , {}                , ShiftxyOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Showinfo"      , {Showinfo      , ShowinfoHelp      , ShowinfoOperators      , 1 , CDI_BOTH , 1  , 0  });
  add_module("Showattribute" , {Showattribute , ShowattributeHelp , ShowattributeOperators , 1 , CDI_REAL , 1  , 0  });
  add_module("Sinfo"         , {Sinfo         , SinfoHelp         , SinfoOperators         , 1 , CDI_BOTH , -1 , 0  });
  add_module("Smooth"        , {Smooth        , SmoothHelp        , SmoothOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Sort"          , {Sort          , {}                , SortOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Sorttimestamp" , {Sorttimestamp , {}                , SorttimestampOperators , 1 , CDI_REAL , -1 , 1  });
  add_module("Specinfo"      , {Specinfo      , {}                , SpecinfoOperators      , 1 , CDI_REAL , 0  , 0  });
  add_module("Spectral"      , {Spectral      , SpectralHelp      , SpectralOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Spectrum"      , {Spectrum      , {}                , SpectrumOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Split"         , {Split         , SplitHelp         , SplitOperators         , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splitrec"      , {Splitrec      , SplitHelp         , SplitrecOperators      , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splitsel"      , {Splitsel      , SplitselHelp      , SplitselOperators      , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splittime"     , {Splittime     , SplittimeHelp     , SplittimeOperators     , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splityear"     , {Splityear     , SplittimeHelp     , SplityearOperators     , 1 , CDI_BOTH , 1  , -1 });
  add_module("Subtrend"      , {Subtrend      , SubtrendHelp      , SubtrendOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Tee"           , {Tee           , TeeHelp           , TeeOperators           , 1 , CDI_REAL , 2  , 1  });
  add_module("Template1"     , {Template1     , {}                , Template1Operators     , 0 , CDI_REAL , 1  , 1  });
  add_module("Template2"     , {Template2     , {}                , Template2Operators     , 0 , CDI_REAL , 1  , 1  });
  add_module("Test"          , {Test          , {}                , TestOperators          , 0 , CDI_REAL , 1  , 1  });
  add_module("Test2"         , {Test2         , {}                , Test2Operators         , 0 , CDI_REAL , 2  , 1  });
  add_module("Testdata"      , {Testdata      , {}                , TestdataOperators      , 0 , CDI_REAL , 1  , 1  });
  add_module("Tests"         , {Tests         , {}                , TestsOperators         , 0 , CDI_REAL , 1  , 1  });
  add_module("Timcount"      , {Timcount      , {}                , TimcountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Yearcount"     , {Timcount      , {}                , YearcountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Moncount"      , {Timcount      , {}                , MoncountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Daycount"      , {Timcount      , {}                , DaycountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Hourcount"     , {Timcount      , {}                , HourcountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timcumsum"     , {Timcumsum     , TimcumsumHelp     , TimcumsumOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timpctl"       , {Timpctl       , TimpctlHelp       , TimpctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Yearpctl"      , {Timpctl       , YearpctlHelp      , YearpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Monpctl"       , {Timpctl       , MonpctlHelp       , MonpctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Daypctl"       , {Timpctl       , DaypctlHelp       , DaypctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Hourpctl"      , {Timpctl       , HourpctlHelp      , HourpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Timselpctl"    , {Timselpctl    , TimselpctlHelp    , TimselpctlOperators    , 1 , CDI_REAL , 3  , 1  });
  add_module("Timsort"       , {Timsort       , TimsortHelp       , TimsortOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Timselstat"    , {Timselstat    , TimselstatHelp    , TimselstatOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("XTimstat"      , {XTimstat      , {}                , XTimstatOperators      , 0 , CDI_BOTH , 1  , 1  });
  add_module("Timstat"       , {Timstat       , TimstatHelp       , TimstatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Yearstat"      , {Timstat       , YearstatHelp      , YearstatOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Monstat"       , {Timstat       , MonstatHelp       , MonstatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Daystat"       , {Timstat       , DaystatHelp       , DaystatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Hourstat"      , {Timstat       , HourstatHelp      , HourstatOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timcor"        , {Timstat2      , TimcorHelp        , TimcorOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Timscorvar"    , {Timstat2      , TimcovarHelp      , TimcovarOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Timstat3"      , {Timstat3      , {}                , Timstat3Operators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Tinfo"         , {Tinfo         , {}                , TinfoOperators         , 1 , CDI_BOTH , 1  , 0  });
  add_module("Tocomplex"     , {Tocomplex     , {}                , TocomplexOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Transpose"     , {Transpose     , {}                , TransposeOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Trend"         , {Trend         , TrendHelp         , TrendOperators         , 1 , CDI_REAL , 1  , 2  });
  add_module("Trms"          , {Trms          , {}                , TrmsOperators          , 0 , CDI_REAL , 2  , 1  });
  add_module("Tstepcount"    , {Tstepcount    , {}                , TstepcountOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Vargen"        , {Vargen        , VargenHelp        , VargenOperators        , 1 , CDI_REAL , 0  , 1  });
  add_module("Varrms"        , {Varrms        , {}                , VarrmsOperators        , 0 , CDI_REAL , 2  , 1  });
  add_module("Vertintml"     , {Vertintml     , VertintmlHelp     , VertintmlOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertintap"     , {Vertintap     , VertintapHelp     , VertintapOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertstat"      , {Vertstat      , VertstatHelp      , VertstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertcum"       , {Vertcum       , {}                , VertcumOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertwind"      , {Vertwind      , {}                , VertwindOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Verifygrid"    , {Verifygrid    , {}                , VerifygridOperators    , 1 , CDI_REAL , 1  , 0  });
  add_module("Wind"          , {Wind          , WindHelp          , WindOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Writegrid"     , {Writegrid     , {}                , WritegridOperators     , 1 , CDI_REAL , 1  , 1  }); // no cdi output
  add_module("Writerandom"   , {Writerandom   , {}                , WriterandomOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("YAR"           , {YAR           , {}                , YAROperators           , 0 , CDI_REAL , 1  , 1  });
  add_module("Yearmonstat"   , {Yearmonstat   , YearmonstatHelp   , YearmonstatOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Ydayarith"     , {Ydayarith     , YdayarithHelp     , YdayarithOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Ydaypctl"      , {Ydaypctl      , YdaypctlHelp      , YdaypctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Ydaystat"      , {Ydaystat      , YdaystatHelp      , YdaystatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Ydrunpctl"     , {Ydrunpctl     , YdrunpctlHelp     , YdrunpctlOperators     , 1 , CDI_REAL , 3  , 1  });
  add_module("Ydrunstat"     , {Ydrunstat     , YdrunstatHelp     , YdrunstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Yhourarith"    , {Yhourarith    , YhourarithHelp    , YhourarithOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Yhourstat"     , {Yhourstat     , YhourstatHelp     , YhourstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Ymonarith"     , {Ymonarith     , YmonarithHelp     , YmonarithOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Yseasarith"    , {Ymonarith     , YseasarithHelp    , YseasarithOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Ymonpctl"      , {Ymonpctl      , YmonpctlHelp      , YmonpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Ymonstat"      , {Ymonstat      , YmonstatHelp      , YmonstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Yseaspctl"     , {Yseaspctl     , YseaspctlHelp     , YseaspctlOperators     , 1 , CDI_REAL , 3  , 1  });
  add_module("Yseasstat"     , {Yseasstat     , YseasstatHelp     , YseasstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Zonstat"       , {Zonstat       , ZonstatHelp       , ZonstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCfd"        , {EcaCfd        , EcaCfdHelp        , EcaCfdOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCsu"        , {EcaCsu        , EcaCsuHelp        , EcaCsuOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCwdi"       , {EcaCwdi       , EcaCwdiHelp       , EcaCwdiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaCwfi"       , {EcaCwfi       , EcaCwfiHelp       , EcaCwfiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaEtr"        , {EcaEtr        , EcaEtrHelp        , EcaEtrOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaFd"         , {EcaFd         , EcaFdHelp         , EcaFdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaGsl"        , {EcaGsl        , EcaGslHelp        , EcaGslOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaHd"         , {EcaHd         , EcaHdHelp         , EcaHdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaHwdi"       , {EcaHwdi       , EcaHwdiHelp       , EcaHwdiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaHwfi"       , {EcaHwfi       , EcaHwfiHelp       , EcaHwfiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaId"         , {EcaId         , EcaIdHelp         , EcaIdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaSu"         , {EcaSu         , EcaSuHelp         , EcaSuOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaTr"         , {EcaTr         , EcaTrHelp         , EcaTrOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaTg10p"      , {EcaTg10p      , EcaTg10pHelp      , EcaTg10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTg90p"      , {EcaTg90p      , EcaTg90pHelp      , EcaTg90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTn10p"      , {EcaTn10p      , EcaTn10pHelp      , EcaTn10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTn90p"      , {EcaTn90p      , EcaTn90pHelp      , EcaTn90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTx10p"      , {EcaTx10p      , EcaTx10pHelp      , EcaTx10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTx90p"      , {EcaTx90p      , EcaTx90pHelp      , EcaTx90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaCdd"        , {EcaCdd        , EcaCddHelp        , EcaCddOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCwd"        , {EcaCwd        , EcaCwdHelp        , EcaCwdOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaRr1"        , {EcaRr1        , EcaRr1Help        , EcaRr1Operators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaPd"         , {EcaPd         , EcaPdHelp         , EcaPdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaR75p"       , {EcaR75p       , EcaR75pHelp       , EcaR75pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR75ptot"    , {EcaR75ptot    , EcaR75ptotHelp    , EcaR75ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR90p"       , {EcaR90p       , EcaR90pHelp       , EcaR90pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR90ptot"    , {EcaR90ptot    , EcaR90ptotHelp    , EcaR90ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR95p"       , {EcaR95p       , EcaR95pHelp       , EcaR95pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR95ptot"    , {EcaR95ptot    , EcaR95ptotHelp    , EcaR95ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR99p"       , {EcaR99p       , EcaR99pHelp       , EcaR99pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR99ptot"    , {EcaR99ptot    , EcaR99ptotHelp    , EcaR99ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaRx1day"     , {EcaRx1day     , EcaRx1dayHelp     , EcaRx1dayOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaRx5day"     , {EcaRx5day     , EcaRx5dayHelp     , EcaRx5dayOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaSdii"       , {EcaSdii       , EcaSdiiHelp       , EcaSdiiOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Fdns"          , {Fdns          , FdnsHelp          , FdnsOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Strwin"        , {Strwin        , StrwinHelp        , StrwinOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Strbre"        , {Strbre        , StrbreHelp        , StrbreOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Strgal"        , {Strgal        , StrgalHelp        , StrgalOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Hurr"          , {Hurr          , HurrHelp          , HurrOperators          , 1 , CDI_REAL , 1  , 1  });
  // add_module("Hi"         , { Hi           , {}                , HiOperators            , 1 , CDI_REAL , 3  , 1   });
  add_module("Wct"           , {Wct           , WctHelp           , WctOperators           , 1 , CDI_REAL , 2  , 1  });
  add_module("Magplot"       , {Magplot       , MagplotHelp       , MagplotOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Magvector"     , {Magvector     , MagvectorHelp     , MagvectorOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Maggraph"      , {Maggraph      , MaggraphHelp      , MaggraphOperators      , 1 , CDI_REAL , -1 , 1  });
  // HIRLAM_EXTENSIONS
  add_module( "Samplegrid"   , { Samplegrid   , SamplegridHelp    , SamplegridOperators    , 1 , CDI_REAL , 1  , 1 });
  add_module( "Selmulti "    , { Selmulti     , SelmultiHelp      , SelmultiOperators      , 1 , CDI_REAL , 1  , 1 });
  add_module( "WindTrans"    , { WindTrans    , WindTransHelp     , WindTransOperators     , 1 , CDI_REAL , 1  , 1 });

  init_aliases();
}

/**
 * Initializes all hardcoded aliases
 */
void init_aliases()
{
  add_alias("afterburner"     , "after");
  add_alias("anomaly"         , "ymonsub");
  add_alias("deltap_fl"       , "deltap");
  add_alias("diffv"           , "diffn");
  add_alias("covar0"          , "timcovar");
  add_alias("covar0r"         , "fldcovar");
  add_alias("gather"          , "collgrid");
  add_alias("geopotheight"    , "gheight");
  add_alias("globavg"         , "fldavg");
  add_alias("import_grads"    , "import_binary");
  add_alias("infos"           , "sinfo");
  add_alias("infov"           , "infon");
  add_alias("intgrid"         , "intgridbil");
  add_alias("log"             , "ln");
  add_alias("lmean"           , "ymonmean");
  add_alias("lmmean"          , "ymonmean");
  add_alias("lmavg"           , "ymonavg");
  add_alias("lmstd"           , "ymonstd");
  add_alias("lsmean"          , "yseasmean");
  add_alias("chvar"           , "chname");
  add_alias("nvar"            , "npar");
  add_alias("outputkey"       , "outputtab");
  add_alias("vardes"          , "codetab");
  add_alias("pardes"          , "codetab");
  add_alias("selvar"          , "selname");
  add_alias("delvar"          , "delname");
  add_alias("selindex"        , "selgridcell");
  add_alias("remapcon1"       , "remaplaf");
  add_alias("remapdis1"       , "remapnn");
  add_alias("scatter"         , "distgrid");
  add_alias("showvar"         , "showname");
  add_alias("selgridname"     , "selgrid");
  add_alias("setvar"          , "setname");
  add_alias("setpartabv"      , "setpartabn");
  add_alias("setpartab"       , "setcodetab");
  add_alias("sinfov"          , "sinfon");
  add_alias("sortvar"         , "sortname");
  add_alias("splitvar"        , "splitname");
  add_alias("sort"            , "timsort");
  add_alias("eca_r1mm"        , "eca_rr1");
  add_alias("fpressure"       , "pressure_fl");
  add_alias("hpressure"       , "pressure_hl");
  add_alias("ensrkhist_space" , "ensrkhistspace");
  add_alias("ensrkhist_time"  , "ensrkhisttime");
  add_alias("gridverify"      , "verifygrid");
  add_alias("outputcenter"    , "gmtxyz");
  add_alias("outputbounds"    , "gmtcells");
  add_alias("selseas"         , "selseason");
  add_alias("selmon"          , "selmonth");
}
/* clang-format on */

/***
 * returns the module of given operator
 *
 * parameters:
 *      std::string operatorName -> name of the operator
 * return value:
 *      std::string -> name of the operators module
 */
static
std::string get_module_name_to(std::string operatorName) {
    // if not true the programm will exit in function check_operator
    if (check_operator(operatorName)) {
        if (modules_map.count(operatorName) > 0) {
            return modules_map[operatorName];
        } else if (aliases.count(operatorName) > 0) {
            return modules_map[aliases[operatorName]];
        }
    } else {
        // TODO: see if depricated since no operator can be added without a
        // module
        Error("Module for >%s< not installed", operatorName.c_str());
    }
    // Only to quell the warning. Function should never reach this.
    return "";
}

/**
 * @fn void *(*operatorModule(const char *operatorName))(void *)

 * returns the function of the module of \a operatorName
 * @param operatorName name of the operator
 * @returns the function of the #modules_t of \a operatorName
 */
void *(*operatorModule(const char *operatorName))(void *) {
    std::string operator_name = std::string(operatorName);
    return modules[get_module_name_to(operator_name)].func;
}

/***
 * returns help for given operator name
 * if there is no help a empty vector will be returned
 * @param operatorName operator name
 * @return vector of strings containing the help
 */
std::vector<std::string> operatorHelp(std::string operatorName) {
    std::string operator_name = std::string(operatorName);
    return modules[get_module_name_to(operator_name)].help;
}

/*** 
 * Returns the number of input streams a operator takes.
 * returns -1 for a unlimited number of input streams.
 * @param operatorName operator name
 * @retval short
 */
int operatorStreamInCnt(const char *operatorName) {
    std::string operator_name = std::string(operatorName);
    return modules[get_module_name_to(operator_name)].streamInCnt;
}

/***
 * Returns the number of output streams a given operator has.
 * returns -1 if obase is used
 * @param operatorName operator name
 * @return 1 for CDI_REAL, 2 for CDI_COMP (complex), 3 for CDI_REAL and CDI_COMP
 * @retval short
 */
int operatorStreamOutCnt(const char *operatorName) {
    std::string operator_name = std::string(operatorName);
    return modules[get_module_name_to(operator_name)].streamOutCnt;
}

/*** 
 * Returns the number type this operator works with
 * @param operatorName operator name
 * @reval short
 */
int operatorStreamNumber(const char *operatorName) {
    std::string operator_name = std::string(operatorName);
    return modules[get_module_name_to(operator_name)].number;
}

/***
 * Creates a sorted vector with all operator names and alisases excluding all modules that are
 * marked as internal
 * @return a sorted std::vector containing all operator names and aliases excluding all operators
 * which modules are marked as internal
 */
static
std::vector<std::string> get_sorted_operator_name_list() {

    std::vector<std::string> names;
    for (std::pair<std::string, std::string> operator_module_names_pair : modules_map) {
        if (modules[operator_module_names_pair.second].mode == 1) {
            names.push_back(operator_module_names_pair.first);
        }
    }
    // adding operators names from alias_map
    for (std::pair<std::string, std::string> alias : aliases) {
        names.push_back(alias.first);
    }
    std::sort(names.begin(), names.end());
    return names;
}

std::vector<std::string> get_no_output_operator_list()
{
 std::vector<std::string> names;
    for (std::pair<std::string, std::string> operator_module_names_pair : modules_map) {
        if (modules[operator_module_names_pair.second].mode == 1 
                && modules[operator_module_names_pair.second].streamOutCnt == 0) 
        {
            names.push_back(operator_module_names_pair.first);
        }
    }
    // adding operators names from alias_map
    std::string original;
    for (std::pair<std::string, std::string> alias : aliases) {
        original = alias.second;
        if(modules[modules_map[original]].mode == 1
                && modules[modules_map[original]].streamOutCnt == 0){
            names.push_back(alias.first);
        }
    }
    std::sort(names.begin(), names.end());
    return names;

}

void operatorPrintAll(void) {
    int number_of_chars = 0;
    std::string tab = "   ";
    int tab_width = tab.size();
    // using a set because it sorts the operators alphabetically on its own
    std::vector<std::string> sorted_operator_names = get_sorted_operator_name_list();

    std::cout << tab;
    for (auto operatorName : sorted_operator_names) {
        if (number_of_chars > 85) {
            number_of_chars = tab_width;
            std::cerr << std::endl << tab;
        }

        std::cerr << " " << operatorName;
        number_of_chars += 1 + operatorName.size();
    }

    std::cerr << std::endl;
}

#ifdef CUSTOM_MODULES
/***
  loads all custom modules in a specified folder
  @param folder_path custom modules folder
*/
#ifdef CUSTOM_MODULES
void load_custom_modules(std::string folder_path) {
    DIR *dir = opendir(folder_path.c_str());
    std::string file_path;
    std::regex library_regex("(.*\\.so)");
    if (dir != NULL) {
        struct dirent *ent = readdir(dir);
        while (ent != NULL) {
            if (std::regex_match(ent->d_name, library_regex)) {
                file_path = folder_path + "/" + ent->d_name;
                load_custom_module(file_path);
            }
            ent = readdir(dir);
        }
    } else {
        std::cerr << "Could not find " << folder_path << "for loading custom modules" << std::endl;
    }
}

/***
 * Loads a custom module from given path.
 * Modules must contain a (TODO: rename function) init_custom_module function
 * Program exits if a module could not be loaded.
 * @param module file path
 */
void load_custom_module(std::string file_path) {
    void (*init_custom_module)();
    void *lib_handle = dlopen(file_path.c_str(), RTLD_LAZY);
    custom_modules_lib_handles.push_back(lib_handle);
    if (!lib_handle) {
        std::cerr << "Cannot open library: " << dlerror() << std::endl;
        return;
    }

    dlerror();
    init_custom_module = (void (*)())dlsym(lib_handle, "init_custom_module");
    const char *dlsym_error = dlerror();

    if (dlsym_error) {
        std::cerr << "Cannot load symbol 'init_custom_module': " << dlsym_error << std::endl;
        dlclose(lib_handle);
        return;
    }

    init_custom_module();
    std::cout << "loaded custom module from '" << file_path << "'" << std::endl;
}
#endif
/***
  closes the handles for the loaded custum modules
*/
void close_library_handles() {
    for (unsigned long i = 0; i < custom_modules_lib_handles.size(); i++) {
        dlclose(custom_modules_lib_handles[i]);
    }
}
#endif
	
// helper function for setting the spacing in operatorPrintList
std::string get_spacing_for(std::string str) {
    int max = 16;
    std::string spacing = "";
    for (int i = str.size(); i <= max; i++) {
        spacing += " ";
    }
    return spacing;
}
/***
 * Prints all operator names and their short descriptions
 * Aliases are listed and point to their original operator name.
 * If the module is not documented the description is empty
 * If a module has only one operator the short module description is listed
 * If the operator is not documented the description is empty
 */
void operatorPrintList(bool print_no_output) {
    std::vector<std::string> output_list ;
    if(print_no_output)
    {
        output_list = get_no_output_operator_list();
    }
    else
    {
        output_list = get_sorted_operator_name_list();
    }
    std::vector<std::string> help;
    unsigned long list_length = output_list.size();
    unsigned long cur_help_idx;
    std::string line;
    std::string description;
    bool help_contains_name;

    for (unsigned long out_list_idx = 0; out_list_idx < list_length; out_list_idx++) {
        std::string current_op_name = output_list[out_list_idx];
        help = modules[get_module_name_to(current_op_name)].help;
        if (aliases.find(current_op_name) != aliases.end()) {

            output_list[out_list_idx] +=
                std::string(get_spacing_for(current_op_name) + "--> " + aliases[current_op_name]);
        }
        else if (!help.empty()) {
            description = "";

            unsigned long operator_section = 0;
            cur_help_idx = 0;
            //search for operator section
            while (operator_section == 0 && cur_help_idx < help.size() - 1) {
                line = help[++cur_help_idx];
                if (line.find("OPERATORS") != std::string::npos) {
                    operator_section = cur_help_idx;
                }
            }
            //if no operator section is found
            if (operator_section == 0) {
                cur_help_idx = 0;
                line = help[0];
                std::string name_section = help[0];
                help_contains_name = false;
                //search for the operator name in the description
                while (!line.empty()) {
                    line = help[++cur_help_idx];
                    if (line.find(current_op_name) != std::string::npos) {
                        help_contains_name = true;
                    }
                    name_section += line;
                }
                //if the name was found save description for later use
                if (help_contains_name) {
                    description = name_section.substr(name_section.find_first_of('-') + 2,
                                                      name_section.size());
                }
            } else {

                line = help[++operator_section];
                //search the operator section for current operator line
                while (line.find(current_op_name + " ") == std::string::npos && !line.empty() &&
                       operator_section < help.size()) {
                    line = help[++operator_section];
                    ;
                }
                //if operator line found save description for later use
                if (!line.empty() && operator_section < help.size()) {
                    auto op_name_start = line.find_first_not_of(" \t");

                    description = line.substr(
                        line.find_first_not_of(" \t", op_name_start + current_op_name.size()),
                        line.size());
                }
            }
            //add spaceing and saving output line to the output list
            output_list[out_list_idx] += get_spacing_for(current_op_name) + description;
        }
    }
    //print generated output list
    for (std::string str : output_list) {
        std::cout << str << std::endl;
    }
}

bool is_alias(char * operatorName)
{
    return (aliases.find(std::string(operatorName)) != aliases.end());
}

char* get_original(char * operatorName)
{
    char* original = NULL;
    if(is_alias(operatorName)){
        std::string opName = aliases[std::string(operatorName)];
        original = (char*)realloc(operatorName, opName.size());
        strcpy(original, opName.c_str());
    }
    else{
        Error("%s is not an alias", operatorName);
    }
    return original; 
}


modules_t &getModule(const std::string &operator_name)
{
    return modules[get_module_name_to(operator_name)];
}
