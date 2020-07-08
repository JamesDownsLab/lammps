//
// Created by james on 07/07/2020.
//

#include "fix_randomforce.h"
#include "error.h"
#include <cstring>
#include "force.h"
#include "random_mars.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include <cmath>
#include "domain.h"
#include "region.h"
#include <iostream>
using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE, CONSTANT, EQUAL, ATOM};

FixRandomForce::FixRandomForce(LAMMPS *lmp, int narg, char **arg) :
        Fix(lmp, narg, arg),
        fvalue(NULL), interval(NULL), seed(NULL), random(NULL), forceConstant(NULL)
{
    if (narg < 3) error->all(FLERR, "Illegal fix random force command");

//    if (strstr(arg[3], "v_") == arg[3]) {
//        int n = strlen(&arg[3][2]) + 1;
//        fstr = new char[n];
//        strcpy(fstr, &arg[3][2]);
//    } else {
//        fstart = force->numeric(FLERR, arg[3]);
//        fstyle = CONSTANT;
//    }

    fstart = force->numeric(FLERR, arg[3]);
    fend = force->numeric(FLERR, arg[4]);

    interval = force->inumeric(FLERR, arg[5]);
    seed = force->inumeric(FLERR, arg[6]);

    if (interval <= 0) error->all(FLERR, "Fix random force interval must be > 0");
    if (seed <= 0) error->all(FLERR, "Fix random force seed must be > 0");

    // initialize Marsaglia RNG with processor-unique seed
    random = new RanMars(lmp, seed + comm->me);

    iregion = -1;
}

FixRandomForce::~FixRandomForce()
{
    delete random;
    delete [] fstr;
}

int FixRandomForce::setmask() {
    int mask = 0;
    mask |= POST_FORCE;
    return mask;
}

void FixRandomForce::init()
{
    dt = update->dt;
    forceConstant = sqrt(1/(dt*interval));
}

void FixRandomForce::post_force(int vflag)
{
    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    fvalue = sqrt(fstart + delta * (fend-fstart));

    if (update->ntimestep % interval) return;

    Region *region = NULL;
    if (iregion >= 0) {
        region = domain->regions[iregion];
        region->prematch();
    }

    double unwrap[3];
    for (int i = 0; i < nlocal; i++){
        if (mask[i] & groupbit) {
            if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
            domain->unmap(x[i], image[i], unwrap);
            f[i][0] += forceConstant*fvalue*random->gaussian();
            f[i][1] += forceConstant*fvalue*random->gaussian();
        }
    }

}