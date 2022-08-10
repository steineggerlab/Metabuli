//
// Created by 김재범 on 2022/06/28.
//

#include "ReducedClassifier.h"

ReducedClassifier::ReducedClassifier(LocalParameters & par)
: Classifier(par){
    setMarker(0Xffffffff);
    setNumOfBitsForCodon(4);
}