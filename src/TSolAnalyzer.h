#ifndef __TSOLANALYZER_H
#define __TSOLANALYZER_H

#include "THaAnalyzer.h"

class TSolAnalyzer : public THaAnalyzer {
    public:
	TSolAnalyzer() {;}
	virtual ~TSolAnalyzer() {;}

    private:
	Int_t f;

    public:
	ClassDef(TSolAnalyzer,0)
};

#endif//__TSOLANALYZER_H
