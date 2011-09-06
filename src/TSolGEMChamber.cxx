#include "TClonesArray.h"
#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
{
    printf("I'm a GEM chamber named %s\n", name );

    fPlanes = new TClonesArray("TSolGEMPlane", 2);

    return;
}
