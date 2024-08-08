#ifndef GENERATE_INCLUDE
#define GENERATE_INCLUDE

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For potential constants.
#include "../include/potentials/laradjiRevalee.h"
#include "../include/potentials/laradjiSpangler.h"
#include "../include/potentials/laradjiPoursoroush.h"

//For certain lipid models
#include "../models/brushModels.h"

//For certain nanoparticle models
#include "../models/nanoModels.h"

//For certain lipid models
#include "../models/lipidModels.h"

//For Yu's Janus particle
#include "../include/algorithms/tesselaSphere.h"
#include "../models/janus.h"

#define CUTOFF 2.0
#define RMIN 1.0

#define HEAD 2
#define TAIL 3
#define NANOPARTICLE 4
#define HEAD2 6
#define TAIL2 7


#endif
