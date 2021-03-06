/******************************************
Copyright (c) 2018, Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#ifndef ALL_PREDICTORS_H
#define ALL_PREDICTORS_H

#include "predict_func_type.h"
#include "clustering.h"

#include <vector>
using std::vector;

namespace CMSat {

void fill_pred_funcs();

//return value must be indexed by cluster
const vector<keep_func_type>& get_short_pred_keep_funcs(size_t conf);

//return value must be indexed by cluster
const vector<keep_func_type>& get_long_pred_keep_funcs(size_t conf);

bool short_pred_func_exists(size_t conf);
bool long_pred_func_exists(size_t conf);

}

#endif //ALL_PREDICTORS_H
