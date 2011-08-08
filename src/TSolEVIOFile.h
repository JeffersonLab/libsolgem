#ifndef __TSOLEVIOFILE_H
#define __TSOLEVIOFILE_H

#include "TObject.h"

#include "evioFileChannel.hxx"
#include "evioUtil.hxx"

class TSolEVIOFile : public TObject {

    public:
	TSolEVIOFile();
	TSolEVIOFile( const char *name );
	virtual ~TSolEVIOFile();

	void SetFilename( const char *name );
	Int_t Open();

	Int_t ReadNextEvent();

    private:
	char fFilename[255];
	evio::evioFileChannel *fChan;

    public:
	ClassDef(TSolEVIOFile, 0)
};

#endif//__TSOLEVIOFILE_H
