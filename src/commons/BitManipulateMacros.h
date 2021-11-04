//
// Created by matchy233 on 7/3/20.
//

#ifndef ADCLASSIFIER2_BITMANIPULATEMACROS_H
#define ADCLASSIFIER2_BITMANIPULATEMACROS_H

#define SET_END_FLAG(num)           (0x8000U | (num))
#define IS_LAST_15_BITS(num)        (0x8000U & (num))
#define GET_15_BITS(num)            (0x7fffU & (num))
#define GET_3_BITS(num)             (0X7U    & (num))
#define GET_2_BITS_1(num) (signed char)(0x3U & (num))
#define GET_2_BITS_2(num) (signed char)(0xCu & (num))
#define GET_2_BITS_3(num) (signed char)(0x30u & (num))
#define GET_2_BITS_4(num) (signed char)(0XC0u & (num))
#define GET_2_BITS_5(num) (signed char)(0X300U & (num))
#define GET_2_BITS_6(num) (signed char)(0XC00U & (num))
#define GET_2_BITS_7(num) (signed char)(0X3000U & (num))
#define GET_2_BITS_8(num) (signed char)(0XC000U & (num))

#define DECODE_15_BITS(diff, num)   ((diff) | (GET_15_BITS(num)))

#endif //ADCLASSIFIER2_BITMANIPULATEMACROS_H
