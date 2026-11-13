# R reads the frozen bytes directly; this avoids this Windows build's C++ path canonicalization bug.
read_frozen_model<-function(path){
 stopifnot(file.exists(path));xgboost::xgb.load.raw(readBin(path,'raw',n=file.size(path)))
}
