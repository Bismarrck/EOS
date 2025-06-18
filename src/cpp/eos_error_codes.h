//
// Created by Bismarrck on 2025/6/18.
//

#ifndef EOS_ERROR_CODES_H
#define EOS_ERROR_CODES_H

// Error codes
constexpr int EOS_SUCCESS = 0;
constexpr int EOS_ERROR_FILE_NOT_FOUND = 1;
constexpr int EOS_ERROR_FILE_PARSE = 2;
constexpr int EOS_ERROR_INVALID_DIMENSIONS = 3;
constexpr int EOS_ERROR_TFD_LOAD_A = 11;
constexpr int EOS_ERROR_TFD_LOAD_B = 12;
constexpr int EOS_ERROR_TFD_NOT_INIT = 13;
constexpr int EOS_ERROR_UNKNOWN_EOS_ID = 21;
constexpr int EOS_ERROR_INVALID_EOS_TYPE = 22;
constexpr int EOS_ERROR_PARAMS_NOT_LOADED = 23;
constexpr int EOS_ERROR_ANALYTIC_DISPATCH = 24;

#endif //EOS_ERROR_CODES_H
