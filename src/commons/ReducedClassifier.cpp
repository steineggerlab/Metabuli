//
// Created by 김재범 on 2022/06/28.
//

#include "ReducedClassifier.h"

ReducedClassifier::ReducedClassifier(LocalParameters & par)
: Classifier(par){
    MARKER = 0Xffffffff;
    MARKER = ~MARKER;
    bitsForCodon = 4;
}