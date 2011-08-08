#include "TSolEVIOFile.h"

TSolEVIOFile::TSolEVIOFile(){
    fFilename[0] = '\0';
    fChan = NULL;

    return;
}

TSolEVIOFile::TSolEVIOFile(const char *f){
    SetFilename(f);
    return;
}

TSolEVIOFile::~TSolEVIOFile(){ return; }

void TSolEVIOFile::SetFilename( const char *f ){
   strcpy( fFilename, f );
   return;
}

Int_t TSolEVIOFile::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }

    fChan = new evio::evioFileChannel(fFilename, "r", 0 );

    if( !fChan ){ return 0; }

    // This is returns nothing
    fChan->open();

    return 1;
}

Int_t TSolEVIOFile::ReadNextEvent(){
    if( !fChan->read() ){
	return 0;
    }

    evio::evioDOMTree EDT(fChan);

    // Read in header info
    evio::evioDOMNodeListP evntNl = EDT.getNodeList(evio::tagNumEquals(1, 1));

}
