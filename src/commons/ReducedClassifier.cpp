//
// Created by 김재범 on 2022/06/28.
//

#include "ReducedClassifier.h"

ReducedClassifier::ReducedClassifier(LocalParameters & par, const vector<TaxID> & taxIdList)
: Classifier(par, taxIdList){
    setMarker(0Xffffffff);
    setNumOfBitsForCodon(4);
}