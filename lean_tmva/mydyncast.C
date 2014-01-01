#include "TMVA/MethodBase.h"

using namespace TMVA;

MethodBase* dyncast(IMethod* met) {
    return dynamic_cast<MethodBase*>(met);
}
