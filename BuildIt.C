void BuildIt(TString filstr){
 
  //  gSystem->Load("/eic/u/zhengl/eicsmear-dev/BeAGLE-install/lib/libeicsmear");
  //  gSystem->AddIncludePath(" -I/eic/u/zhengl/eicsmear-dev/BeAGLE-install/include/");
  gSystem->Load("/afs/rhic.bnl.gov/eic/lib/libeicsmear.so");
  TString fullstr = "outForPythiaMode/" + filstr;
  BuildTree(fullstr);
}
