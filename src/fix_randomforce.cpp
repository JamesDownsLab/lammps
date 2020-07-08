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
#include "group.h"
#include <iostream>
using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE, CONSTANT, EQUAL, ATOM};

FixRandomForce::FixRandomForce(LAMMPS *lmp, int narg, char **arg) :
        Fix(lmp, narg, arg),
        fvalue(NULL), interval(NULL), seed(NULL), random(NULL), forceConstant(NULL)
{
    if (narg < 3) error->all(FLERR, "Illegal fix random force command");

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

    bigint count = group->count(igroup);

    if (update->ntimestep % interval) return;

    double fx_all = 0;
    double fy_all = 0;
    double unwrap[3];
    for (int i = 0; i < nlocal; i++){
        if (mask[i] & groupbit) {
            domain->unmap(x[i], image[i], unwrap);
            double fx = forceConstant*fvalue*random->gaussian();
            double fy = forceConstant*fvalue*random->gaussian();
            f[i][0] += fx;
            f[i][1] += fx;
            fx_all += fx;
            fy_all += fy;
        }
    }
    // Set total force to zero
    for (int i = 0; i<nlocal; i++){
        if (mask[i] & groupbit) {
            f[i][0] -= fx_all/count;
            f[i][1] -= fy_all/count;
        }
    }

}