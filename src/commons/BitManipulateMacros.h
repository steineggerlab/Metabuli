//
// Created by matchy233 on 7/3/20.
//

#ifndef ADCLASSIFIER2_BITMANIPULATEMACROS_H
#define ADCLASSIFIER2_BITMANIPULATEMACROS_H

#define SET_END_FLAG(num)           (0x8000U | (num))
#define IS_LAST_15_BITS(num)        (0x8000U & (num))
#define GET_15_BITS(num)            (0x7fffU & (num))
#define GET_3_BITS(num)             (0X7U    & (num))
#define GET_2_BITS(num) (0X3U & (num))
#define DECODE_15_BITS(diff, num)   ((diff) | (GET_15_BITS(num)))

#endif //ADCLASSIFIER2_BITMANIPULATEMACROS_H
