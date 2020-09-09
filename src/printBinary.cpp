//
// Created by KJB on 05/08/2020.
//
#include "printBinary.h"

void print_binary64(int n, uint64_t x){
    if(n==0){
        return;
    }
    print_binary64(n - 1, x >> 1);
    putc((x & 0x01)? '1':'0', stdout);

}

void print_binary16(int n, uint16_t x){
    if(n==0){
        return;
    }
    print_binary16(n - 1, x >> 1);
    putc((x & 0x01)? '1':'0', stdout);
}