/*
 * nsht_legmat_api.c
 *
 * Code generation for function 'nsht_legmat_api'
 *
 * C source code generated on: Thu Aug 14 15:05:14 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nsht_legmat.h"
#include "nsht_legmat_api.h"
#include "nsht_legmat_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo b_emlrtRTEI = { 1, 1, "nsht_legmat_api", "" };

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static real_T c_emlrt_marshallIn(const mxArray *L, const char_T *identifier);
static real_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static void e_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void emlrt_marshallIn(const mxArray *thetas, const char_T *identifier,
  emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static real_T f_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static void info_helper(ResolvedFunctionInfo info[37]);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  e_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T c_emlrt_marshallIn(const mxArray *L, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(emlrtAlias(L), &thisId);
  emlrtDestroyArray(&L);
  return y;
}

static real_T d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv3[2];
  boolean_T bv0[2];
  int32_T i4;
  static const boolean_T bv1[2] = { FALSE, TRUE };

  int32_T iv4[2];
  for (i4 = 0; i4 < 2; i4++) {
    iv3[i4] = 1 + 204799 * i4;
    bv0[i4] = bv1[i4];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv3, bv0, iv4);
  ret->size[0] = iv4[0];
  ret->size[1] = iv4[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void emlrt_marshallIn(const mxArray *thetas, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(emlrtAlias(thetas), &thisId, y);
  emlrtDestroyArray(&thetas);
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv2[2] = { 0, 0 };

  const mxArray *m2;
  y = NULL;
  m2 = mxCreateNumericArray(2, (int32_T *)&iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u->data);
  mxSetDimensions((mxArray *)m2, u->size, 2);
  emlrtAssign(&y, m2);
  return y;
}

static real_T f_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void info_helper(ResolvedFunctionInfo info[37])
{
  info[0].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[0].name = "error";
  info[0].dominantType = "char";
  info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/lang/error.m";
  info[0].fileTimeLo = 1319697566U;
  info[0].fileTimeHi = 0U;
  info[0].mFileTimeLo = 0U;
  info[0].mFileTimeHi = 0U;
  info[1].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[1].name = "length";
  info[1].dominantType = "double";
  info[1].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  info[1].fileTimeLo = 1303117406U;
  info[1].fileTimeHi = 0U;
  info[1].mFileTimeLo = 0U;
  info[1].mFileTimeHi = 0U;
  info[2].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!norm_factor_first";
  info[2].name = "mtimes";
  info[2].dominantType = "double";
  info[2].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[2].fileTimeLo = 1289483692U;
  info[2].fileTimeHi = 0U;
  info[2].mFileTimeLo = 0U;
  info[2].mFileTimeHi = 0U;
  info[3].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!norm_factor_first";
  info[3].name = "mrdivide";
  info[3].dominantType = "double";
  info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[3].fileTimeLo = 1357915548U;
  info[3].fileTimeHi = 0U;
  info[3].mFileTimeLo = 1319697566U;
  info[3].mFileTimeHi = 0U;
  info[4].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[4].name = "rdivide";
  info[4].dominantType = "double";
  info[4].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[4].fileTimeLo = 1346481588U;
  info[4].fileTimeHi = 0U;
  info[4].mFileTimeLo = 0U;
  info[4].mFileTimeHi = 0U;
  info[5].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[5].name = "eml_scalexp_compatible";
  info[5].dominantType = "double";
  info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  info[5].fileTimeLo = 1286786396U;
  info[5].fileTimeHi = 0U;
  info[5].mFileTimeLo = 0U;
  info[5].mFileTimeHi = 0U;
  info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[6].name = "eml_div";
  info[6].dominantType = "double";
  info[6].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  info[6].fileTimeLo = 1313319010U;
  info[6].fileTimeHi = 0U;
  info[6].mFileTimeLo = 0U;
  info[6].mFileTimeHi = 0U;
  info[7].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!norm_factor_first";
  info[7].name = "sqrt";
  info[7].dominantType = "double";
  info[7].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  info[7].fileTimeLo = 1343801586U;
  info[7].fileTimeHi = 0U;
  info[7].mFileTimeLo = 0U;
  info[7].mFileTimeHi = 0U;
  info[8].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  info[8].name = "eml_error";
  info[8].dominantType = "char";
  info[8].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  info[8].fileTimeLo = 1343801558U;
  info[8].fileTimeHi = 0U;
  info[8].mFileTimeLo = 0U;
  info[8].mFileTimeHi = 0U;
  info[9].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  info[9].name = "eml_scalar_sqrt";
  info[9].dominantType = "double";
  info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  info[9].fileTimeLo = 1286786338U;
  info[9].fileTimeHi = 0U;
  info[9].mFileTimeLo = 0U;
  info[9].mFileTimeHi = 0U;
  info[10].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!norm_factor_first";
  info[10].name = "abs";
  info[10].dominantType = "double";
  info[10].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[10].fileTimeLo = 1343801566U;
  info[10].fileTimeHi = 0U;
  info[10].mFileTimeLo = 0U;
  info[10].mFileTimeHi = 0U;
  info[11].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[11].name = "eml_scalar_abs";
  info[11].dominantType = "double";
  info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  info[11].fileTimeLo = 1286786312U;
  info[11].fileTimeHi = 0U;
  info[11].mFileTimeLo = 0U;
  info[11].mFileTimeHi = 0U;
  info[12].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!norm_factor_first";
  info[12].name = "mpower";
  info[12].dominantType = "double";
  info[12].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  info[12].fileTimeLo = 1286786442U;
  info[12].fileTimeHi = 0U;
  info[12].mFileTimeLo = 0U;
  info[12].mFileTimeHi = 0U;
  info[13].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  info[13].name = "power";
  info[13].dominantType = "double";
  info[13].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  info[13].fileTimeLo = 1348163130U;
  info[13].fileTimeHi = 0U;
  info[13].mFileTimeLo = 0U;
  info[13].mFileTimeHi = 0U;
  info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  info[14].name = "eml_scalar_eg";
  info[14].dominantType = "double";
  info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[14].fileTimeLo = 1286786396U;
  info[14].fileTimeHi = 0U;
  info[14].mFileTimeLo = 0U;
  info[14].mFileTimeHi = 0U;
  info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  info[15].name = "eml_scalexp_alloc";
  info[15].dominantType = "double";
  info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  info[15].fileTimeLo = 1352388860U;
  info[15].fileTimeHi = 0U;
  info[15].mFileTimeLo = 0U;
  info[15].mFileTimeHi = 0U;
  info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  info[16].name = "floor";
  info[16].dominantType = "double";
  info[16].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  info[16].fileTimeLo = 1343801580U;
  info[16].fileTimeHi = 0U;
  info[16].mFileTimeLo = 0U;
  info[16].mFileTimeHi = 0U;
  info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  info[17].name = "eml_scalar_floor";
  info[17].dominantType = "double";
  info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  info[17].fileTimeLo = 1286786326U;
  info[17].fileTimeHi = 0U;
  info[17].mFileTimeLo = 0U;
  info[17].mFileTimeHi = 0U;
  info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  info[18].name = "eml_scalar_eg";
  info[18].dominantType = "double";
  info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[18].fileTimeLo = 1286786396U;
  info[18].fileTimeHi = 0U;
  info[18].mFileTimeLo = 0U;
  info[18].mFileTimeHi = 0U;
  info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  info[19].name = "mtimes";
  info[19].dominantType = "double";
  info[19].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[19].fileTimeLo = 1289483692U;
  info[19].fileTimeHi = 0U;
  info[19].mFileTimeLo = 0U;
  info[19].mFileTimeHi = 0U;
  info[20].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!sine_factor";
  info[20].name = "sin";
  info[20].dominantType = "double";
  info[20].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  info[20].fileTimeLo = 1343801586U;
  info[20].fileTimeHi = 0U;
  info[20].mFileTimeLo = 0U;
  info[20].mFileTimeHi = 0U;
  info[21].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  info[21].name = "eml_scalar_sin";
  info[21].dominantType = "double";
  info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  info[21].fileTimeLo = 1286786336U;
  info[21].fileTimeHi = 0U;
  info[21].mFileTimeLo = 0U;
  info[21].mFileTimeHi = 0U;
  info[22].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!sine_factor";
  info[22].name = "mtimes";
  info[22].dominantType = "double";
  info[22].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[22].fileTimeLo = 1289483692U;
  info[22].fileTimeHi = 0U;
  info[22].mFileTimeLo = 0U;
  info[22].mFileTimeHi = 0U;
  info[23].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!sine_factor";
  info[23].name = "abs";
  info[23].dominantType = "double";
  info[23].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[23].fileTimeLo = 1343801566U;
  info[23].fileTimeHi = 0U;
  info[23].mFileTimeLo = 0U;
  info[23].mFileTimeHi = 0U;
  info[24].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!sine_factor";
  info[24].name = "mpower";
  info[24].dominantType = "double";
  info[24].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  info[24].fileTimeLo = 1286786442U;
  info[24].fileTimeHi = 0U;
  info[24].mFileTimeLo = 0U;
  info[24].mFileTimeHi = 0U;
  info[25].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[25].name = "mtimes";
  info[25].dominantType = "double";
  info[25].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[25].fileTimeLo = 1289483692U;
  info[25].fileTimeHi = 0U;
  info[25].mFileTimeLo = 0U;
  info[25].mFileTimeHi = 0U;
  info[26].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[26].name = "mpower";
  info[26].dominantType = "double";
  info[26].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  info[26].fileTimeLo = 1286786442U;
  info[26].fileTimeHi = 0U;
  info[26].mFileTimeLo = 0U;
  info[26].mFileTimeHi = 0U;
  info[27].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  info[27].name = "eml_error";
  info[27].dominantType = "char";
  info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  info[27].fileTimeLo = 1343801558U;
  info[27].fileTimeHi = 0U;
  info[27].mFileTimeLo = 0U;
  info[27].mFileTimeHi = 0U;
  info[28].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[28].name = "mrdivide";
  info[28].dominantType = "double";
  info[28].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[28].fileTimeLo = 1357915548U;
  info[28].fileTimeHi = 0U;
  info[28].mFileTimeLo = 1319697566U;
  info[28].mFileTimeHi = 0U;
  info[29].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[29].name = "sqrt";
  info[29].dominantType = "double";
  info[29].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  info[29].fileTimeLo = 1343801586U;
  info[29].fileTimeHi = 0U;
  info[29].mFileTimeLo = 0U;
  info[29].mFileTimeHi = 0U;
  info[30].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m";
  info[30].name = "abs";
  info[30].dominantType = "double";
  info[30].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[30].fileTimeLo = 1343801566U;
  info[30].fileTimeHi = 0U;
  info[30].mFileTimeLo = 0U;
  info[30].mFileTimeHi = 0U;
  info[31].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!Kostelec_recursion_scaled";
  info[31].name = "mtimes";
  info[31].dominantType = "double";
  info[31].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[31].fileTimeLo = 1289483692U;
  info[31].fileTimeHi = 0U;
  info[31].mFileTimeLo = 0U;
  info[31].mFileTimeHi = 0U;
  info[32].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!Kostelec_recursion_scaled";
  info[32].name = "mrdivide";
  info[32].dominantType = "double";
  info[32].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[32].fileTimeLo = 1357915548U;
  info[32].fileTimeHi = 0U;
  info[32].mFileTimeLo = 1319697566U;
  info[32].mFileTimeHi = 0U;
  info[33].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!Kostelec_recursion_scaled";
  info[33].name = "sqrt";
  info[33].dominantType = "double";
  info[33].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  info[33].fileTimeLo = 1343801586U;
  info[33].fileTimeHi = 0U;
  info[33].mFileTimeLo = 0U;
  info[33].mFileTimeHi = 0U;
  info[34].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!Kostelec_recursion_scaled";
  info[34].name = "mpower";
  info[34].dominantType = "double";
  info[34].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  info[34].fileTimeLo = 1286786442U;
  info[34].fileTimeHi = 0U;
  info[34].mFileTimeLo = 0U;
  info[34].mFileTimeHi = 0U;
  info[35].context =
    "[E]C:/Users/u5472257/Documents/GitHub/dMRI_Sampling_Scheme/nsht_legmat.m!Kostelec_recursion_scaled";
  info[35].name = "cos";
  info[35].dominantType = "double";
  info[35].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  info[35].fileTimeLo = 1343801572U;
  info[35].fileTimeHi = 0U;
  info[35].mFileTimeLo = 0U;
  info[35].mFileTimeHi = 0U;
  info[36].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  info[36].name = "eml_scalar_cos";
  info[36].dominantType = "double";
  info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  info[36].fileTimeLo = 1286786322U;
  info[36].fileTimeHi = 0U;
  info[36].mFileTimeLo = 0U;
  info[36].mFileTimeHi = 0U;
}

const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  ResolvedFunctionInfo info[37];
  ResolvedFunctionInfo u[37];
  int32_T i;
  const mxArray *y;
  int32_T iv1[1];
  ResolvedFunctionInfo *r0;
  const char * b_u;
  const mxArray *b_y;
  const mxArray *m1;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  uint32_T c_u;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  nameCaptureInfo = NULL;
  info_helper(info);
  for (i = 0; i < 37; i++) {
    u[i] = info[i];
  }

  y = NULL;
  iv1[0] = 37;
  emlrtAssign(&y, mxCreateStructArray(1, iv1, 0, NULL));
  for (i = 0; i < 37; i++) {
    r0 = &u[i];
    b_u = r0->context;
    b_y = NULL;
    m1 = mxCreateString(b_u);
    emlrtAssign(&b_y, m1);
    emlrtAddField(y, b_y, "context", i);
    b_u = r0->name;
    c_y = NULL;
    m1 = mxCreateString(b_u);
    emlrtAssign(&c_y, m1);
    emlrtAddField(y, c_y, "name", i);
    b_u = r0->dominantType;
    d_y = NULL;
    m1 = mxCreateString(b_u);
    emlrtAssign(&d_y, m1);
    emlrtAddField(y, d_y, "dominantType", i);
    b_u = r0->resolved;
    e_y = NULL;
    m1 = mxCreateString(b_u);
    emlrtAssign(&e_y, m1);
    emlrtAddField(y, e_y, "resolved", i);
    c_u = r0->fileTimeLo;
    f_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m1) = c_u;
    emlrtAssign(&f_y, m1);
    emlrtAddField(y, f_y, "fileTimeLo", i);
    c_u = r0->fileTimeHi;
    g_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m1) = c_u;
    emlrtAssign(&g_y, m1);
    emlrtAddField(y, g_y, "fileTimeHi", i);
    c_u = r0->mFileTimeLo;
    h_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m1) = c_u;
    emlrtAssign(&h_y, m1);
    emlrtAddField(y, h_y, "mFileTimeLo", i);
    c_u = r0->mFileTimeHi;
    i_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m1) = c_u;
    emlrtAssign(&i_y, m1);
    emlrtAddField(y, i_y, "mFileTimeHi", i);
  }

  emlrtAssign(&nameCaptureInfo, y);
  emlrtNameCapturePostProcessR2012a(emlrtAlias(nameCaptureInfo));
  return nameCaptureInfo;
}

void nsht_legmat_api(const mxArray * const prhs[3], const mxArray *plhs[2])
{
  emxArray_real_T *thetas;
  emxArray_real_T *P;
  emxArray_real_T *Sc;
  real_T L;
  real_T m;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&thetas, 2, &b_emlrtRTEI, TRUE);
  emxInit_real_T(&P, 2, &b_emlrtRTEI, TRUE);
  emxInit_real_T(&Sc, 2, &b_emlrtRTEI, TRUE);

  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAlias(prhs[0]), "thetas", thetas);
  L = c_emlrt_marshallIn(emlrtAliasP(prhs[1]), "L");
  m = c_emlrt_marshallIn(emlrtAliasP(prhs[2]), "m");

  /* Invoke the target function */
  nsht_legmat(thetas, L, m, P, Sc);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(P);
  plhs[1] = emlrt_marshallOut(Sc);
  Sc->canFreeData = FALSE;
  emxFree_real_T(&Sc);
  P->canFreeData = FALSE;
  emxFree_real_T(&P);
  thetas->canFreeData = FALSE;
  emxFree_real_T(&thetas);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (nsht_legmat_api.c) */
