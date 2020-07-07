
#ifdef FIX_CLASS

FixStyle(randomforce,FixRandomForce)

#else

#ifndef LAMMPS_FIX_RANDOMFORCE_H
#define LAMMPS_FIX_RANDOMFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRandomForce : public Fix {
public:
    FixRandomForce(class LAMMPS *, int, char **);

    virtual ~FixRandomForce();

    int setmask();

    void init();

    virtual void post_force(int);


private:
    double fvalue;
    int interval, seed, iregion;
    char *fstr;
    int fstyle;
    double dt;

    class RanMars *random;
};

}


#endif //LAMMPS_FIX_RANDOMFORCE_H
#endif
