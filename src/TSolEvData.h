#ifndef __TSOLEVDATA_H
#define __TSOLEVDATA_H

#include "THaEvData.h"

class TSolEVIOFile;

class TSolEvData : public THaEvData {
    public:
	TSolEvData();
	~TSolEvData();

	void Clear();

	int LoadEvent( const int*, THaCrateMap *);
	int LoadEvent( TSolEVIOFile *, THaCrateMap *);

    private:

    public:
	ClassDef(TSolEvData,0)
};

#endif//__TSOLEVDATA_H
