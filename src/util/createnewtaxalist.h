#ifndef METABULI_CREATE_NEW_TAXA_LIST_H
#define METABULI_CREATE_NEW_TAXA_LIST_H

#include <string>
#include "TaxonomyWrapper.h"

int createnewtaxalist(TaxonomyWrapper * oldTaxonomy,
                      TaxonomyWrapper * newTaxonomy,
                      const std::string & newNodesFileName,
                      const std::string & newNamesFileName, 
                      const std::string & accession2taxidFileName, 
                      std::vector<NewTaxon> & newTaxaList,
                      std::map<std::string, TaxID> & newAccessions,
                      std::vector<std::string> & unmappedAccessions);

#endif //METABULI_CREATE_NEW_TAXA_LIST_H